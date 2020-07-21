#' @title Realod source
#' @description Function to reaload libraries and source utils.R
#'
#' @return null
reload_source <- function() {
  if (!require("dplyr")) install.packages("dplyr")
  library("dplyr")
  if (!require("tidyr")) install.packages("tidyr")
  library("tidyr")
  if (!require("ggplot2")) install.packages("ggplot2")
  library("ggplot2")
  if (!require("purrr")) install.packages("purrr")
  library("purrr")
  if (!require("readr")) install.packages("readr")
  library("readr")
  if (!require("lubridate")) install.packages("lubridate")
  library("lubridate")
  if (!require("RColorBrewer")) install.packages("RColorBrewer")
  library("RColorBrewer")
  if (!require("knitr")) install.packages("knitr")
  library("knitr")
  
  source("source/utils.R")
}

#' @title Run Stan analysis with random effects
#' @description Runs the Stan seroprevalence model with household random effects
#'
#' @param model_script Stan model script
#' @param dat list with data to run the model
#' @param analysis which analysis to run among different seropositivity thresholds (see details)
#' @param coef_eqn formula in character format expressing the probability of seropositivity on the logit scale
#' @param pos_control number of predicted positives in the control dataset
#' @param neg_control number of predicted negatives in the control dataset
#' @param control_tp number of true positives in the control dataset
#' @param control_tn number of true negatives in the control dataset
#' @param pop_age_cats age categories
#' @param n_cores number of cores to use for parallel computation
#' @param sex_ref reference category for sex
#' @param age_ref reference category for age
#' @param wk_ref reference category for serosurvey week
#' @param redo redo fit or load pre-computed posteriors if available
#' 
#' @details The analysis parameter needs to be one of:
#' - ei: main analysis with manufacturer threshold
#' - gva: custom Geneva threshold
#' - ifa: alternative thershold
#' 
#' @return a list with parameter posteriors and results
run_analysis_stan_re <- function(model_script,
                                 dat,
                                 analysis = "ei",
                                 coef_eqn,
                                 pos_control,
                                 neg_control,
                                 control_tp,
                                 control_fp,
                                 pop_age_cats,
                                 n_cores = detectCores() - 2,
                                 sex_ref = 0,
                                 age_ref = "[20,50)",
                                 wk_ref = "2",
                                 redo = F,
                                 chains,
                                 iter,
                                 warmup,
                                 control, 
                                 ...) {
  ## Prepare data for analysis
  ana_suffix <- case_when(analysis == "ei" ~ "",
                          analysis == "gva" ~ "_gva",
                          analysis == "ifa" ~ "_ifa",
                          T ~ as.character(NA))
  if (is.na(ana_suffix))
    stop("Unkonwn analysis parameter, must be one of ei, gva, ifa")
  
  # Set analysis data
  ana_dat <- as_tibble(dat)
  
  for (var in c("pos", "neg", "ind")) {
    ana_dat[[var]] <- as.numeric(dat[[paste0(var, ana_suffix)]])
  }
  
  # Set model matrix
  X <- model.matrix(as.formula(paste("~", coef_eqn)), data = ana_dat)
  
  # name of the stan output file 
  stan_out_file <- paste0("results/stan_fit_", analysis,".rds")
  
  if (!file.exists(stan_out_file) | redo) {
    
    re_model <- stan_model(model_script)
    
    # Unique household ids from 1 to N
    u_hh_ids <- unique(ana_dat$household_id)
    ana_dat$u_household_id <- map_dbl(ana_dat$household_id, ~which(u_hh_ids == .))
    
    # Run stan model
    stan_est <- sampling(re_model,
                         data = list(
                           N_survey = nrow(ana_dat),
                           H = length(u_hh_ids),
                           hh = ana_dat$u_household_id,
                           p_vars = ncol(X),
                           X = X,
                           survey_pos = ana_dat$pos,
                           N_pos_control = pos_control,
                           control_tp = control_tp,
                           N_neg_control = neg_control,
                           control_fp = control_fp
                         ),
                         chains = chains,
                         iter = iter,
                         warmup = warmup,
                         control = control)
    
    saveRDS(stan_est, stan_out_file)
  } else {
    cat("Loading pre-computed posteriors from ", stan_out_file, "\n")
    stan_est <- readRDS(stan_out_file)
  }
  
  ## extract log-likelihoods
  stan_ll <- loo::extract_log_lik(stan_est) %>% loo::loo()
  
  ## extract parameters
  beta <- extract(stan_est, pars = "beta")[[1]]
  sigma <- extract(stan_est, pars = "sigma_h")[[1]]
  
  pop_cat_mat <- pop_age_cats %>%
    model.matrix(as.formula(paste("~", coef_eqn)), data = .)
  
  ## compute estimates by age category
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  
  pop_cat_p <- foreach(i = 1:nrow(pop_cat_mat),
                       .combine = rbind, 
                       .inorder = F, 
                       .packages = c("tidyverse", "foreach")) %dopar% 
    { 
      foreach(j = 1:nrow(beta), 
              .combine = rbind, 
              .inorder = T) %do% 
        {
          
          # Compute probability integrating across in household random effects
          prob <- integrate(function(x) {
            plogis(qnorm(
              x, beta[j, , drop = F] %*% t(pop_cat_mat[i, , drop = F]),
              sigma[j]
            ))
          }, 0, 1)[[1]]
          
          tibble(
            age_cat = pop_age_cats$age_cat[i],
            Sex = pop_age_cats$Sex[i],
            week = pop_age_cats$week[i],
            pop = pop_age_cats$pop[i],
            seropos = prob
          ) %>%
            mutate(sim = j)
        }
    }
  parallel::stopCluster(cl)
  
  ## overall estimate
  overall_re <- pop_cat_p %>%
    filter(week != 1) %>%
    mutate(
      var = "Overall (excluding week 1)",
      val = ""
    ) %>%
    group_by(sim, var, val) %>%
    summarize(p = weighted.mean(seropos, pop)) %>%
    ungroup()
  
  # weekyl estimates
  wk_re <- pop_cat_p %>%
    mutate(var = "Week") %>%
    rename(val = week) %>%
    group_by(sim, var, val) %>%
    summarize(p = weighted.mean(seropos, pop)) %>%
    ungroup()
  
  ## find age specific probabilities in order to make relative risks
  age_re <- pop_cat_p %>%
    filter(Sex == sex_ref, week == wk_ref) %>%
    mutate(var = "Age") %>%
    rename(val = age_cat) %>%
    group_by(sim, var, val) %>%
    summarize(p = weighted.mean(seropos, pop)) %>%
    ungroup()
  
  # sex-specific probabilities
  sex_re <- pop_cat_p %>%
    filter(age_cat == age_ref, week == wk_ref) %>%
    mutate(var = "Sex") %>%
    rename(val = Sex) %>%
    group_by(sim, var, val) %>%
    summarize(p = weighted.mean(seropos, pop)) %>%
    ungroup()
  
  # restuls list
  res <- list(
    beta = extract(stan_est, pars = "beta")[[1]],
    model_mtx = X,
    sigma_h = extract(stan_est, pars = "sigma_h")[[1]],
    sens = extract(stan_est, pars = "sens")[[1]],
    spec = extract(stan_est, pars = "spec")[[1]],
    obs = nrow(ana_dat),
    pos = sum(ana_dat$pos),
    neg = sum(ana_dat$neg),
    ind = sum(ana_dat$ind),
    stan_ll = stan_ll,
    pop_cat_p = pop_cat_p,
    subset_est = bind_rows(
      overall_re,
      wk_re,
      sex_re,
      age_re
    )
  )
  
  return(res)
}

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
#' - ifa: alternative threshold
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
                                 age_ref = "[20,40)",
                                 wk_ref = "2",
                                 rha_name_ref = "WRHA",
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
  if (!dir.exists('results')) {dir.create('results')}
  stan_out_file <- paste0("results/stan_fit_", analysis,".rds")

  if (!file.exists(stan_out_file) | redo) {

    set_cmdstan_path("C:/cmdstan/")
    re_model <- cmdstan_model(model_script)

    # Unique household ids from 1 to N
    u_hh_ids <- unique(ana_dat$household_id)
    ana_dat$u_household_id <- map_dbl(ana_dat$household_id, ~which(u_hh_ids == .))

    # Run stan model
    cmdstan_est <- re_model$sample(data = list(
                           N_survey = nrow(ana_dat),
                           H = length(u_hh_ids),
                           n_time = max(sero_dat$time_point),
                           R = max(sero_dat$age_cat1),
                           hh = ana_dat$u_household_id,
                           p_vars = ncol(X),
                           X = X,
                           time_point = sero_dat$time_point,
                           survey_pos = ana_dat$pos,
                           N_pos_control = pos_control,
                           control_tp = control_tp,
                           N_neg_control = neg_control,
                           control_fp = control_fp,
                           rha = sero_dat$rha,
                           age = sero_dat$age_cat1,
                           diagnostic_file = "diag.csv"
                         ),
                         chains = chains,
                         parallel_chains = chains,
                         iter_sampling = iter - warmup,
                         iter_warmup = warmup,
                         adapt_delta = control[[1]],
                         max_treedepth = control[[2]],
                         refresh = 100)

    stan_est <- rstan::read_stan_csv(cmdstan_est$output_files())
    saveRDS(stan_est, stan_out_file)
  } else {
    cat("Loading pre-computed posteriors from ", stan_out_file, "\n")
    stan_est <- readRDS(stan_out_file)
  }
  ## extract log-likelihoods
  stan_ll <- stan_est %>% loo::loo()

  ## extract parameters
  beta <- extract(stan_est, pars = "beta")[[1]]
  sigma <- extract(stan_est, pars = "sigma_h")[[1]]
  sigmas_timepoints <- extract(stan_est, pars = "sigma_t")$sigma_t
  means_t <- extract(stan_est, pars = "means_t")$means_t
  z_timepoints <- extract(stan_est, pars = "eta_t")$eta_t
  n_timepoints <- dim(z_timepoints)[2]

  post_mat <- cbind(beta, sigma, sigmas_timepoints)

  pop_cat_mat <- pop_age_cats %>%
    model.matrix(as.formula(paste("~", coef_eqn)), data = .)

  
  N_rows <- nrow(beta)
  N_preds <- ncol(beta)

  ## compute estimates by age category
  #Parallel for MacOS
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  #Parallel for Linux
  #doParallel::registerDoParallel(n_cores)

    pop_cat_p <- foreach(i = 1:nrow(pop_cat_mat),
                       .combine = bind_rows,
                       .inorder = T,
                       .packages = c("tidyverse", "foreach")) %dopar%
    {
      
      output <- apply(as.array(1:N_rows), 1, function(posterior){
        print(posterior)
        integrate(function(x) {
            plogis(qnorm(
              x, post_mat[posterior, 1:N_preds, drop = F] %*% t(pop_cat_mat[i, 1:N_preds, drop = F]) +
               means_t[posterior, pop_age_cats$week[i]] +
               post_mat[posterior, N_preds + 2, drop = F] * z_timepoints[posterior, pop_age_cats$week[i], pop_age_cats$age[i]],
              post_mat[posterior, N_preds + 1]
            ))
          }, 0, 1,rel.tol=1e-5)[[1]]
      })

      out <- tibble(
            age_cat = pop_age_cats$age_cat[i],
            Sex = pop_age_cats$Sex[i],
            week = pop_age_cats$week[i],
            pop = pop_age_cats$pop[i],
            rha = pop_age_cats$rha[i],
            rha_name = pop_age_cats$rha_name[i],
            seropos = output
      ) %>%
      mutate(sim = 1:n())
      validate_tibble(out)
      return(out)
    }
  # parallel::stopCluster(cl)
  stopImplicitCluster()

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
    mutate(var = "Week",
            week = as.character(week)) %>%
    rename(val = week) %>%
    group_by(sim, var, val) %>%
    summarize(p = weighted.mean(seropos, pop)) %>%
    ungroup()

  # weekly estimates by age
  wk_age_re <- pop_cat_p %>%
    mutate(var = "Week_Age",
           val = paste0(week, "_", age_cat1)) %>%
    group_by(sim, var, val) %>%
    summarize(p = weighted.mean(seropos, pop)) %>%
    ungroup()

  ## find age specific probabilities in order to make relative risks
  age_re <- pop_cat_p %>%
    filter(Sex == sex_ref, week == wk_ref, rha_name == rha_name_ref) %>%
    mutate(var = "Age") %>%
    rename(val = age_cat) %>%
    group_by(sim, var, val) %>%
    summarize(p = weighted.mean(seropos, pop)) %>%
    ungroup()

  # sex-specific probabilities
  sex_re <- pop_cat_p %>%
    filter(age_cat == age_ref, week == wk_ref, rha_name == rha_name_ref) %>%
    mutate(var = "Sex") %>%
    rename(val = Sex) %>%
    group_by(sim, var, val) %>%
    summarize(p = weighted.mean(seropos, pop)) %>%
    ungroup()

  # RHA-specific probabilities
  rha_re <- pop_cat_p %>%
    filter(age_cat == age_ref, week == wk_ref) %>%
    mutate(var = "RHA") %>%
    rename(val = rha_name) %>%
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
      wk_age_re,
      sex_re,
      age_re,
      rha_re
    )
  )

  return(res)
}

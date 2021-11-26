//
// This Stan program defines a model for adjusting a predicted
// seroincidence by the sensitivity and specificity of the diagnostic using data from lab validation study.
// used for primary analyses of serocov-pop study in Geneva, Switzerland

// We have input data from a validation set (which determines the diagnostic performance) and the survey (from which we'd like to estimate seropos).
data {
    int<lower=1> N_survey; //numbber of participants in the survey
    int<lower=1> H; //number of households in the survey
    int<lower = 1> R; //number of age categories
    int<lower=0> survey_pos[N_survey]; //observation positive or negative
    int<lower=1, upper=H> hh[N_survey]; //household of observation
    int<lower=1> p_vars; //number of variables to adjust for
    int<lower = 1> n_time; //number of time periods
    matrix[N_survey, p_vars] X; //covariate model matrix (age, sex in these analyses)
    int<lower = 1, upper = n_time> time_point[N_survey]; //model matrix of the time periods in the analysis
    int<lower=0> N_pos_control; //number of positive controls in the validation data
    int<lower=0,upper=N_pos_control> control_tp; // number of true positive tests in the validation data
    int<lower=0> N_neg_control; // number of negative controls in the validation data
    int<lower=0,upper=N_neg_control> control_fp;// number of false positives by the diagnostic test in the validation study
    int<lower = 1, upper = R> age[N_survey]; //Age category a specific observation
}

parameters {
    real<lower=0, upper=1> spec; // specificity of the diagnostic test.
    real<lower=0, upper=1> sens; // sensitivity of the diagnostic test.
    vector[p_vars] beta; // fixed regression coefficients
    real<lower=0> sigma_h; // variability of household random effect
    vector[H] eta_h; // standard normals for the household random effect
    vector<lower = 0>[n_time] sigma_t; //variability of each time period across RHAs
    matrix[n_time, R] eta_t; //standard normals for the RHA x time effects
    vector[n_time] means_t; // mean for each time point, around which the RHAs are centered
    // real<lower=0> sigma_i;
    // vector[N_survey] eta_i;
}

transformed parameters {
  vector<lower=0, upper=1>[N_survey] p; // probability of seropositivity for an observation
{
  matrix[R, n_time] time_coeffs = (rep_matrix(sigma_t, R) .* eta_t)';
  vector[N_survey] time_effects;
  for(i in 1:N_survey)
    time_effects[i] = means_t[time_point[i]] + time_coeffs[age[i], time_point[i]];
  p = inv_logit(X * beta + time_effects + sigma_h*eta_h[hh]);
}
}

//  We observe 'survey_pos' cases as a bernoulli distribution based on the
//  survey with observations coming as the sum of the true
//  positive rate p*sens and false negative rate (1-p)*(1-spec).
model {
    target+= bernoulli_lpmf(survey_pos | p*sens+(1-p)*(1-spec));
    target+= binomial_lpmf(control_tp | N_pos_control, sens);
    target+= binomial_lpmf(control_fp | N_neg_control, 1-spec);
    eta_h ~ std_normal();
    for(i in 1:n_time)
      to_vector(eta_t[i, ]) ~ std_normal();
    sigma_t ~ std_normal();
    means_t ~ normal(0, 2);

    beta ~ std_normal(); // priors for coefficients
}

generated quantities {
  vector[N_survey] log_lik;

    for(i in 1:N_survey){
      log_lik[i] = bernoulli_lpmf(survey_pos[i] | p[i]*sens+(1-p[i])*(1-spec));
    }

}

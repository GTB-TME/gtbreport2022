/* non H, separate new/ret, global random effects */
functions{
#include function_lognormaladd.stan
}
data{
  real h_beta_mu;
  real h_beta_sigma;
  real h_alpha_sigma;
  real h_delta_sigma;
  real h_phi_sigma;
#include data_common.stan
}
#include transformeddata_common.stan

parameters{
#include parameters_common.stan
  real<lower=0> beta_sigma;
  real<lower=0> alpha_sigma;
}

transformed parameters {
#include transformedparameters_common.stan
}

model {
  /* ===== Prior likelihood ===== */
  beta_sigma ~ cauchy(0,h_beta_sigma);
  alpha_sigma  ~ cauchy(0,h_alpha_sigma);
  to_vector(beta) ~ normal(h_beta_mu, beta_sigma);
  to_vector(alpha) ~ normal(0, alpha_sigma);
  to_vector(phi) ~ normal(0,h_phi_sigma);
  to_vector(delta) ~ normal(0,h_delta_sigma);

  /* ==== Data likelihood ==== */
  if(J==2){                     /* J=2 */
#include datalikelihood_N_noH.stan
#include datalikelihood_R_noH.stan
  } else {                      /* J=4 */
#include datalikelihood_LN.stan
#include datalikelihood_LR.stan
#include datalikelihood_YN.stan
#include datalikelihood_YR.stan
  }
}

generated quantities {
  /* the log-likelihood vector needed by loo */
  vector[NDP] log_lik;
  int k;
  k = 1;
  if(J==2){                     /* J=2 */
#include GQ_loo_loglike_noH.stan
  } else {                      /* J=4 */
#include GQ_loo_loglike_H.stan
  }
}

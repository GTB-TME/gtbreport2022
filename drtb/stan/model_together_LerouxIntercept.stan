/* non H, separate new/ret, global random effects */
functions{
#include function_lognormaladd.stan
#include function_sparse_Leroux.stan
}
data{
  real h_beta_mu;
  real h_beta_sigma;
  real h_alpha_sigma;
  real h_delta_sigma;
#include data_common.stan
}
#include transformeddata_common.stan

parameters{
  cholesky_factor_corr[2*(J-1)] L_Omega; /* count correlations */
  vector<lower=0>[2*(J-1)] omega;        /* RVs for Sigma */
  matrix[N,2*(J-1)] phi;                  /* spatial REs */
  matrix[N, 2*(J-1)] delta;                /* slope REs */
  matrix[P, 2*(J-1)] beta;                /* regression coefs */
  matrix[P, 2*(J-1)] alpha;                /* regression coefs */
  real<lower = 0, upper = 1> alphsp;   /* spatial parms */
}

transformed parameters {
#include transformedparameters_common.stan
}

model {
  L_Omega ~ lkj_corr_cholesky(2);
  omega ~ normal(0, 2.5);
  phi ~ sparse_leroux(alphsp,E,Dvec,lambda,
                      diag_pre_multiply(omega, L_Omega),
                      N,2*(J-1)+1,Nedges);
  /* ===== Prior likelihood ===== */
  to_vector(beta) ~ normal(h_beta_mu, h_beta_sigma);
  to_vector(alpha) ~ normal(0, h_alpha_sigma);
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

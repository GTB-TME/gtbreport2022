/* non H, separate new/ret, global random effects */

functions{
#include function_lognormaladd.stan
#include function_sparse_Leroux.stan
}
data{
  int ark;                      /* AR order>0 */
  real h_beta_mu;
  real h_beta_sigma;
  real h_alpha_sigma;
#include data_common.stan
}
#include transformeddata_common.stan

parameters{
  matrix[P, 2*(J-1)] beta;                /* regression coefs */
  matrix[P, 2*(J-1)] alpha;                /* regression coefs */
  matrix[N, 2*(J-1)] dpsi[T];                /* incremental REs */
  cholesky_factor_corr[2*(J-1)] L_Omega; /* count correlations */
  vector<lower=0>[2*(J-1)] omega;        /* RVs for Sigma */
  real<lower = 0, upper = 1> alphsp;   /* spatial parms */
  vector<lower = 0, upper = 1>[ark] alphti;   /* time AR parm */
}

transformed parameters {
#include transformedparameters_ARk.stan
}

model {

  /* ===== RE  priors ===== */
  L_Omega ~ lkj_corr_cholesky(2);
  omega ~ normal(0, 2.5);
  for( j in 1:T)
    dpsi[j] ~ sparse_leroux(alphsp,E,Dvec,lambda,
                            diag_pre_multiply(omega, L_Omega),
                            N,2*(J-1)+1,Nedges);


  /* ===== regression coef priors ===== */
  to_vector(beta) ~ normal(h_beta_mu, h_beta_sigma);
  to_vector(alpha) ~ normal(0, h_alpha_sigma);

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

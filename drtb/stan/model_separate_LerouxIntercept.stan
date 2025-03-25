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
  cholesky_factor_corr[(J-1)] L_Omegan; /* count correlations */
  vector<lower=0>[(J-1)] omegan;        /* RVs for Sigma */
  cholesky_factor_corr[(J-1)] L_Omegar; /* count correlations */
  vector<lower=0>[(J-1)] omegar;        /* RVs for Sigma */
  matrix[N,2*(J-1)] phi;                  /* spatial REs */
  matrix[N, 2*(J-1)] delta;                /* slope REs */
  matrix[P, 2*(J-1)] beta;                /* regression coefs */
  matrix[P, 2*(J-1)] alpha;                /* regression coefs */
  real<lower = 0, upper = 1> alphspn;   /* spatial parms */
  real<lower = 0, upper = 1> alphspr;   /* spatial parms */
}

transformed parameters {
#include transformedparameters_common.stan
}

model {
  L_Omegan ~ lkj_corr_cholesky(2);
  omegan ~ normal(0, 2.5);
  L_Omegar ~ lkj_corr_cholesky(2);
  omegar ~ normal(0, 2.5);
  phi[,1:(J-1)] ~ sparse_leroux(alphspn,E,Dvec,lambda,
                                diag_pre_multiply(omegan, L_Omegan),
                                N,(J-1)+1,Nedges);
  phi[,J:(2*(J-1))] ~ sparse_leroux(alphspr,E,Dvec,lambda,
                                    diag_pre_multiply(omegar, L_Omegar),
                                    N,(J-1)+1,Nedges);
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

/* common data description */

/* --- shared sizes --- */
int<lower=0> T;               /* number of sites */
int<lower=0> N;               /* number of sites */
int<lower=0> J;               /* number of cats */
int<lower=0> P;               /* number of parms, inc intercept */
int<lower=0> NDP;             /* number of RR data points */
matrix[N,P] X;                /* covariates */
/* random effects */
int<lower=0> Nedges;          /* number of edges */
matrix<lower = 0, upper = 1>[N, N] W; /* spatial matrix */
int<lower=0> NTedges;          /* number of edges */
matrix<lower = 0, upper = 1>[T, T] D; /* temporal matrix */


/* --- sizes of each case --- */
int<lower=0> N_L_new_rr;
int<lower=0> N_L_ret_rr;
int<lower=0> N_L_new_hr;
int<lower=0> N_L_ret_hr;
int<lower=0> N_L_new_rm;
int<lower=0> N_L_ret_rm;
int<lower=0> N_L_new_x;
int<lower=0> N_L_ret_x;
int<lower=0> N_L_new_dh0;
int<lower=0> N_L_new_dh1;
int<lower=0> N_L_new_dh2;
int<lower=0> N_L_ret_dh0;
int<lower=0> N_L_ret_dh1;
int<lower=0> N_L_ret_dh2;
int<lower=0> N_Y_new_rr;
int<lower=0> N_Y_ret_rr;
int<lower=0> N_Y_new_x;
int<lower=0> N_Y_ret_x;
int<lower=0> N_Y_new_dh0;
int<lower=0> N_Y_ret_dh0;
int<lower=0> N_Y_ret_dh1r;
int<lower=0> N_Y_new_dh1h;
int<lower=0> N_Y_ret_dh1h;
int<lower=0> N_Y_new_dh2;
int<lower=0> N_Y_ret_dh2;

/* --- surveillance data --- */
int<lower=0> L_new_rr_id[N_L_new_rr];
int<lower=0> L_new_rr_N[N_L_new_rr,2];
int L_new_rr_t[N_L_new_rr];

int<lower=0> L_ret_rr_id[N_L_ret_rr];
int<lower=0> L_ret_rr_N[N_L_ret_rr,2];
int L_ret_rr_t[N_L_ret_rr];

int<lower=0> L_new_hr_id[N_L_new_hr];
int<lower=0> L_new_hr_N[N_L_new_hr,2];
int L_new_hr_t[N_L_new_hr];

int<lower=0> L_ret_hr_id[N_L_ret_hr];
int<lower=0> L_ret_hr_N[N_L_ret_hr,2];
int L_ret_hr_t[N_L_ret_hr];

int<lower=0> L_new_rm_id[N_L_new_rm];
int<lower=0> L_new_rm_N[N_L_new_rm,2];
int L_new_rm_t[N_L_new_rm];

int<lower=0> L_ret_rm_id[N_L_ret_rm];
int<lower=0> L_ret_rm_N[N_L_ret_rm,2];
int L_ret_rm_t[N_L_ret_rm];

int<lower=0> L_new_x_id[N_L_new_x];
int<lower=0> L_new_x_N[N_L_new_x,2];
int L_new_x_t[N_L_new_x];

int<lower=0> L_ret_x_id[N_L_ret_x];
int<lower=0> L_ret_x_N[N_L_ret_x,2];
int L_ret_x_t[N_L_ret_x];

int<lower=0> L_new_dh0_id[N_L_new_dh0];
int<lower=0> L_new_dh0_N[N_L_new_dh0,4];
int L_new_dh0_t[N_L_new_dh0];

int<lower=0> L_new_dh1_id[N_L_new_dh1];
int<lower=0> L_new_dh1_N[N_L_new_dh1,3];
int L_new_dh1_t[N_L_new_dh1];

int<lower=0> L_new_dh2_id[N_L_new_dh2];
int<lower=0> L_new_dh2_N[N_L_new_dh2,2];
int L_new_dh2_t[N_L_new_dh2];

int<lower=0> L_ret_dh0_id[N_L_ret_dh0];
int<lower=0> L_ret_dh0_N[N_L_ret_dh0,4];
int L_ret_dh0_t[N_L_ret_dh0];

int<lower=0> L_ret_dh1_id[N_L_ret_dh1];
int<lower=0> L_ret_dh1_N[N_L_ret_dh1,3];
int L_ret_dh1_t[N_L_ret_dh1];

int<lower=0> L_ret_dh2_id[N_L_ret_dh2];
int<lower=0> L_ret_dh2_N[N_L_ret_dh2,2];
int L_ret_dh2_t[N_L_ret_dh2];

/* --- survey data --- */

int<lower=0> Y_new_rr_id[N_Y_new_rr];
real<upper=0> Y_new_rr_M[N_Y_new_rr];
real<lower=0> Y_new_rr_S[N_Y_new_rr];
int Y_new_rr_t[N_Y_new_rr];

int<lower=0> Y_ret_rr_id[N_Y_ret_rr];
real<upper=0> Y_ret_rr_M[N_Y_ret_rr];
real<lower=0> Y_ret_rr_S[N_Y_ret_rr];
int Y_ret_rr_t[N_Y_ret_rr];

int<lower=0> Y_new_x_id[N_Y_new_x];
real<upper=0> Y_new_x_M[N_Y_new_x];
real<lower=0> Y_new_x_S[N_Y_new_x];
int Y_new_x_t[N_Y_new_x];

int<lower=0> Y_ret_x_id[N_Y_ret_x];
real<upper=0> Y_ret_x_M[N_Y_ret_x];
real<lower=0> Y_ret_x_S[N_Y_ret_x];
int Y_ret_x_t[N_Y_ret_x];

int<lower=0> Y_new_dh0_id[N_Y_new_dh0];
real<upper=0> Y_new_dh0_M[N_Y_new_dh0,3];
real<lower=0> Y_new_dh0_S[N_Y_new_dh0,3];
int Y_new_dh0_t[N_Y_new_dh0];

int<lower=0> Y_ret_dh0_id[N_Y_ret_dh0];
real<upper=0> Y_ret_dh0_M[N_Y_ret_dh0,3];
real<lower=0> Y_ret_dh0_S[N_Y_ret_dh0,3];
int Y_ret_dh0_t[N_Y_ret_dh0];

int<lower=0> Y_ret_dh1r_id[N_Y_ret_dh1r];
real<upper=0> Y_ret_dh1r_M[N_Y_ret_dh1r,2];
real<lower=0> Y_ret_dh1r_S[N_Y_ret_dh1r,2];
int Y_ret_dh1r_t[N_Y_ret_dh1r];

int<lower=0> Y_new_dh1h_id[N_Y_new_dh1h];
real<upper=0> Y_new_dh1h_M[N_Y_new_dh1h,2];
real<lower=0> Y_new_dh1h_S[N_Y_new_dh1h,2];
int Y_new_dh1h_t[N_Y_new_dh1h];

int<lower=0> Y_ret_dh1h_id[N_Y_ret_dh1h];
real<upper=0> Y_ret_dh1h_M[N_Y_ret_dh1h,2];
real<lower=0> Y_ret_dh1h_S[N_Y_ret_dh1h,2];
int Y_ret_dh1h_t[N_Y_ret_dh1h];

int<lower=0> Y_new_dh2_id[N_Y_new_dh2];
real<upper=0> Y_new_dh2_M[N_Y_new_dh2];
real<lower=0> Y_new_dh2_S[N_Y_new_dh2];
int Y_new_dh2_t[N_Y_new_dh2];

int<lower=0> Y_ret_dh2_id[N_Y_ret_dh2];
real<upper=0> Y_ret_dh2_M[N_Y_ret_dh2];
real<lower=0> Y_ret_dh2_S[N_Y_ret_dh2];
int Y_ret_dh2_t[N_Y_ret_dh2];

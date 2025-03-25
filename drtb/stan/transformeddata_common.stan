transformed data {
  /* matrices for aggregating where data is missing */
  matrix [4,1] m4k = [[0],[0],[0],[1]];    /* keep 4 */
  matrix [4,2] m4 = [[1,0],[1,0],[1,0],[0,1]];/* k 4,aggregate rest */
  matrix [4,3] m24 = [[1,0,0],[0,1,0],[1,0,0],[0,0,1]]; /* k=24,a13 */
  matrix [4,2] m24k = [[0,0],[1,0],[0,0],[0,1]]; /* keep 24 */
  matrix [4,2] m34k = [[0,0],[0,0],[1,0],[0,1]]; /* keep 34 */

  matrix [4,2] m1234 = [[1,0],[1,0],[0,1],[0,1]]; /* two aggregates: 1&2, 3&4*/
  matrix [4,2] m1324 = [[1,0],[0,1],[1,0],[0,1]]; /* two aggregates: 1&3, 2&4 */
  matrix [4,2] m1423 = [[1,0],[0,1],[0,1],[1,0]]; /* two aggregates: 1&4, 2&3 */
  /* spatial stuff */
  real small = 1e-9;
  vector[N] zeros;
  vector[N] ones;
  vector[N] W_rowsums;
  row_vector[N] Dvec;
  int E[Nedges, 2];   // adjacency pairs
  vector[N] lambda;       // eigenvalues
  matrix[N,N] Wstar;
  real tz[T];       // times
  vector[T] tzeros;  /* times */
  /* temporal */
  vector[T] D_rowsums;
  row_vector[T] TDvec;
  int TE[NTedges, 2];   // adjacency pairs
  vector[T] Tlambda;       // eigenvalues
  matrix[T,T] Dstar;
  /* regional keys */
  int regkey[N];
  int regkeym[2*(J-1)*N];


  /* reduced forms of data for no-H modelling only */
  /* new */
  int<lower=0> L_new_dh0_N_noH[N_L_new_dh0,2];
  real<upper=0> Y_new_dh0_M_noH[N_Y_new_dh0];
  real<lower=0> Y_new_dh0_S_noH[N_Y_new_dh0];
  real<upper=0> Y_new_dh1h_M_noH[N_Y_new_dh1h];
  real<lower=0> Y_new_dh1h_S_noH[N_Y_new_dh1h];
  /* ret */
  int<lower=0> L_ret_dh0_N_noH[N_L_ret_dh0,2];
  real<upper=0> Y_ret_dh0_M_noH[N_Y_ret_dh0];
  real<lower=0> Y_ret_dh0_S_noH[N_Y_ret_dh0];
  real<upper=0> Y_ret_dh1h_M_noH[N_Y_ret_dh1h];
  real<lower=0> Y_ret_dh1h_S_noH[N_Y_ret_dh1h];
  /* new */
  for(j in 1:N_L_new_dh0){
    L_new_dh0_N_noH[j,1] = L_new_dh0_N[j,1]+L_new_dh0_N[j,2];
    L_new_dh0_N_noH[j,2] = L_new_dh0_N[j,3]+L_new_dh0_N[j,4];
  }
  for(j in 1:N_Y_new_dh0){
    Y_new_dh0_M_noH[j] = LNaddmu(Y_new_dh0_M[j,3],Y_new_dh0_S[j,3],
                                 Y_new_dh0_M[j,2],Y_new_dh0_S[j,2]);
    Y_new_dh0_S_noH[j] = LNaddsigma(Y_new_dh0_M[j,3],Y_new_dh0_S[j,3],
                                   Y_new_dh0_M[j,2],Y_new_dh0_S[j,2]);
  }
  for(j in 1:N_Y_new_dh1h){
    Y_new_dh1h_M_noH[j]=LNaddmu(Y_new_dh1h_M[j,1],Y_new_dh1h_S[j,1],
                                Y_new_dh1h_M[j,2],Y_new_dh1h_S[j,2]);
   Y_new_dh1h_S_noH[j]=LNaddsigma(Y_new_dh1h_M[j,1],Y_new_dh1h_S[j,1],
                                 Y_new_dh1h_M[j,2],Y_new_dh1h_S[j,2]);
  }
  /* ret */
  for(j in 1:N_L_ret_dh0){
    L_ret_dh0_N_noH[j,1] = L_ret_dh0_N[j,1]+L_ret_dh0_N[j,2];
    L_ret_dh0_N_noH[j,2] = L_ret_dh0_N[j,3]+L_ret_dh0_N[j,4];
  }
  for(j in 1:N_Y_ret_dh0){
    Y_ret_dh0_M_noH[j] = LNaddmu(Y_ret_dh0_M[j,3],Y_ret_dh0_S[j,3],
                                 Y_ret_dh0_M[j,2],Y_ret_dh0_S[j,2]);
    Y_ret_dh0_S_noH[j]= LNaddsigma(Y_ret_dh0_M[j,3],Y_ret_dh0_S[j,3],
                                   Y_ret_dh0_M[j,2],Y_ret_dh0_S[j,2]);
  }
  for(j in 1:N_Y_ret_dh1h){
    Y_ret_dh1h_M_noH[j]=LNaddmu(Y_ret_dh1h_M[j,1],Y_ret_dh1h_S[j,1],
                                Y_ret_dh1h_M[j,2],Y_ret_dh1h_S[j,2]);
   Y_ret_dh1h_S_noH[j]=LNaddsigma(Y_ret_dh1h_M[j,1],Y_ret_dh1h_S[j,1],
                                 Y_ret_dh1h_M[j,2],Y_ret_dh1h_S[j,2]);
  }

  /* regional key */
  for(i in 1:N){
    regkey[i] = 1;            /* AFR is default w/current coding */
    for(j in 2:P){
      if(X[i,j]==1){
        regkey[i] = j;        /* which region (if not AFR) */
      }
    }
  }
  for(j in 1:(2*(J-1))){
    for(i in 1:N){
      regkeym[(j-1)*N+i] = regkey[i];
    }
  }

  /* times */
  for(i in 1:T){
      tz[i] = 1.0*i;
  }
  tzeros = rep_vector(0, T);

  /* spatial stuff */
  for (i in 1:N) {
    W_rowsums[i] = sum(W[i, ]);
  }
  Dvec = W_rowsums';
  zeros = rep_vector(0, N);
  { // generate sparse representation for W
    int counter;
    counter = 1;
    // loop over upper triangular part of W to identify neighbor pairs
    for (i in 1:(N - 1)) {
      for (j in (i + 1):N) {
        if (W[i, j] == 1) {
          E[counter, 1] = i;
          E[counter, 2] = j;
          counter = counter + 1;
        }
      }
    }
  }
  Wstar = -W;
  for(n in 1:N)
    Wstar[n,n] = Wstar[n,n] + Dvec[n] + small;
  lambda = eigenvalues_sym(Wstar);


  /* temporal adjacency stuff */
  for (i in 1:T) {
    D_rowsums[i] = sum(D[i, ]);
  }
  TDvec = D_rowsums';
  { // generate sparse representation for D
    int counter;
    counter = 1;
    // loop over upper triangular part of D to identify neighbor pairs
    for (i in 1:(T - 1)) {
      for (j in (i + 1):T) {
        if (D[i, j] == 1) {
          TE[counter, 1] = i;
          TE[counter, 2] = j;
          counter = counter + 1;
        }
      }
    }
  }
  Dstar = -D;
  for(n in 1:T)
    Dstar[n,n] = Dstar[n,n] + TDvec[n] + small;
  Tlambda = eigenvalues_sym(Dstar);
}

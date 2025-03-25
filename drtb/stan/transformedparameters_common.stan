/* common transformed parameters block */
matrix[N,J] ntheta[T];            /* transformed probs: new cases */
matrix[N,J] nprops[T];
matrix[N,J] rtheta[T];            /* transformed probs: ret cases */
matrix[N,J] rprops[T];

/* NOTE new version with intercept/slope */
matrix[N,J-1] nincpt;
matrix[N,J-1] rincpt;
matrix[N,J-1] nslope;
matrix[N,J-1] rslope;

/* intercept/slope */
nincpt = X * beta[,1:(J-1)] + phi[,1:(J-1)];
nslope = X * alpha[,1:(J-1)] + delta[,1:(J-1)];
rincpt = X * beta[,J:(2*(J-1))] + phi[,J:(2*(J-1))];
rslope = X * alpha[,J:(2*(J-1))] + delta[,J:(2*(J-1))];

/* time */
for(j in 1:T){
  ntheta[j] = append_col(zeros, nincpt + (j-1) * nslope);
  rtheta[j] = append_col(zeros, rincpt + (j-1) * rslope);
}

for(j in 1:T){
  for (k in 1:N){
    nprops[j][k] = softmax(ntheta[j][k]')';
    rprops[j][k] = softmax(rtheta[j][k]')';
  }
}



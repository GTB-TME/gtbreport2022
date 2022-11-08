/**
 * Return the log probability of a proper conditional autoregressive (CAR) prior
 * with a sparse representation for the adjacency matrix. Adapted for Leroux priors
 * from Max Joseph's work.
 *
 * @param phi matrix containing the parameters with a CAR prior
 * @param alpha Dependence (usually spatial) parameter for the CAR prior (real)
 * @param E Sparse representation of adjacency matrix (int array of node pairs)
 * @param N Length of phi (int)
 * @param J columns in data
 * @param Nedges Number of adjacent pairs (int)
 * @param D Number of neighbors for each location (vector)
 * @param lambda Eigenvalues of W^* (vector)
 * @param L the cholesky for Sigma
 *
 * @return Log probability density of CAR prior up to additive constant
 */
real sparse_leroux_lpdf(matrix phi,  real alpha,
                        int[,] E, row_vector D,
                        vector lambda, matrix L,
                        int N, int J, int Nedges) {
  matrix[J-1,N] psi;
  vector[N] ldet_terms;
  vector[J-1] quad_terms;

  psi = mdivide_left_tri_low(L,phi'); /* psi = L \ phi' */
  /* phi' * ((1-alpha)*eye + alpha*Wstar) * phi; */
  for(k in 1:(J-1))
    quad_terms[k] = (1-alpha)*dot_self(psi[k]) +
      alpha*dot_product( (D .* psi[k]), psi[k] ) -
      2*alpha*dot_product( psi[k,E[,1]], psi[k,E[,2]] );

  /* (1-a)1 + aW* = (1-a)(1+a W* / (1-a)) */
  /* deta(K (x) Sigma) = (det K)^{J-1} (det Sigma)^N */
  for (i in 1:N) ldet_terms[i] = log1p(alpha * lambda[i] / (1-alpha));

  return 0.5 * ( +N * (J-1) * log(1-alpha)
                 +sum(ldet_terms)*(J-1)
                 -sum(quad_terms)) - N*sum(log(diagonal(L)));
}


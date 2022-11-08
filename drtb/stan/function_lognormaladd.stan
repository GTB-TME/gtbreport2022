/* LN parms for sum of two LN distributions */
real LNaddmu(real mu1,  real sigma1, real mu2, real sigma2){
  real E1; real E2; real E; real T2;
  E1 = exp(mu1 + sigma1^2/2);
  E2 = exp(mu2 + sigma2^2/2);
  E = E1 + E2;
  T2 = E1^2 * exp(sigma1^2) + E2^2 * exp(sigma2^2) + 2 * E1* E2;
  return log(E^2/sqrt(T2));
}

real LNaddsigma(real mu1,  real sigma1, real mu2, real sigma2){
  real E1; real E2; real E; real T2;
  E1 = exp(mu1 + sigma1^2/2);
  E2 = exp(mu2 + sigma2^2/2);
  E = E1 + E2;
  T2 = E1^2 * exp(sigma1^2) + E2^2 * exp(sigma2^2) + 2 * E1* E2;
  return sqrt(log(T2/E^2));
}

## add CIs
## helper functions
getCI1 <- function(x) binom.test(x[1],x[2],p=.025)$conf.int #k,N
getCI <- function(k,N) t(apply(cbind(k,N),1,getCI1))
getCI(c(5,5,5),c(10,10,10))
## function to add binomial CIs
MLH <- function(k,N) {
  if(length(k)!=length(N)) stop('k and N have different lengths!')
  k <- as.integer(k); N <- as.integer(N)
  mid <- lo <- hi <- rep(NA,length(k))
  who <- which(!is.na(k) & !is.na(N) & N>0)
  if(any(k[who]<0)){ stop('k<0!')}
  if(any(k[who]>N[who])){ stop('k>N!')}
  HL <- getCI(k[who],N[who])
  mid[who] <- k[who]/N[who]
  lo[who] <- HL[,1]; hi[who] <- HL[,2]
  list(mid*1e2,lo*1e2,hi*1e2)
}
MLH(c(5,5,5),c(10,10,10))


## PG vectorized lohi: see R/fun.R for lohi
vlohi <- Vectorize(lohi, c('ev', 'sd'))

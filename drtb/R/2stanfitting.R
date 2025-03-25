## this can be run as R --slave --vanilla --args <2stanfitting.R modelname N X
## 
## if N and X are absent, it will produce full outputs for modelname, including PSIS LOO estimates
## if N is present and X is absent, explicit omission of country N from the omitlist will be inferred
## if N>0 is present and X is present, omission of most recent year (dy=0) from country N will be run
## if N<=0 and X is present, all of the last dy=N years' data from all countries will be omitted
## if N<=0 and X is absent, an error will be thrown

args <- commandArgs(trailingOnly = TRUE)
## NOTE debugging
## args <- c('together_LerouxIntercept_2','-2','T')

variant <- as.character(args[1])
cat('variant = ',variant,'\n')
## variant <- 'separate_global_2'
## variant <- 'together_AR1_2'
## variant <- 'together_AR1_4'
nargs <- length(args)


## arg logic
tallstring <- '' #used if omitting all countries recent data
omitTstr <- ''
omittingT <- FALSE
omitting <- FALSE
if(nargs>1){
  omitting <- TRUE
  omitcno <- as.integer(args[2])
  cat('Will be omitting country number ',omitcno,'\n(<=0 zero means -n years for all countries)\n')
  if(omitcno<=0 & nargs<=2) stop('Needs additional flag to do rollback!')
  if(nargs>2){
    omittingT <- TRUE
    omitTstr <- 'T'
    cat('...only omitting most recent data points\n')
    if(omitcno<=0){
      omitcno <- -omitcno
      tallstring <- 'TA' #used as flag that doing this case
      cat('...omitting most recent years up to dy=',omitcno,' data points dropped for *all* countries!\n')
    }
  } else {
    omittingT <- FALSE
  }
} else {
  omitcno <- 0
}

print(omitcno)
print(c(omitting,omittingT))
print(c(omitTstr,tallstring))

## libraries
library(here)
library(glue)
library(rstan)
library(posterior)
library(loo)
library(data.table)
library(ggplot2)
gh <- function(x) glue(here(x))
rot45 <- theme(axis.text.x = element_text(angle = 45, hjust = 1))

## data
load(gh('drtb/data/D4S.Rdata'))
load(gh('drtb/data/isoidx5.Rdata'))
load(here('drtb/data/regkey.Rdata'))
load(here('drtb/data/RPD.Rdata'))
ckey <- data.table(iso3=isoidx,cid=1:length(isoidx))
RPD <- merge(RPD[year>=2000],regkey[,.(iso3,g_whoregion)],by='iso3')
## make temporal adjacency matrix 
DM <- matrix(0,nrow=D4S$T,ncol=D4S$T)
addones <- function(X,condition=function(x,y) abs(x-y)==1){
  for(i in 1:nrow(X)){
    for(j in 1:ncol(X)){
      if(condition(i,j))
        X[i,j] <- 1
    }
  }
  X
}
confn <- function(x,y) abs(x-y)<=3
DM <- addones(DM,confn)
diag(DM) <- 0
D4S$D <- DM
D4S$NTedges <- sum(DM)/2
## omit iso3 list
c.omit.v <- scan(gh('drtb/data/c.omit.vA.txt'),what = 'string') #NOTE change list here
c.omit.vT <- scan(gh('drtb/data/c.omit.vAt.txt'),what = 'string') #NOTE change list here
if(omittingT){
  c.omit.v <- c.omit.vT
}
if(omitting){ #get iso3 code to omit for this job
  if(tallstring==''){
    omitcn <- c.omit.v[omitcno]
    cat('omitting country ',omitcn,'\n')
  } else {
    omitcn <- omitcno #n years for output file names
  }
} else {
  omitcn <- ''
}


## load plotting functions
source(gh('drtb/R/utils/plotsetc.R'))
source(gh('drtb/R/utils/countryomitter.R'))

## prior scoping
## ilogit(2)
## ilogit(-2)
## softmax <- function(x) exp(x)/sum(exp(x))
## softmax(c(0,-2*(1.5+1)))

## === VARIANTS ===

getJ <- function(x) as.integer(substr(x,start=nchar(x),stop=nchar(x)))
getModelName <- function(x) substr(x,start=1,stop=nchar(x)-2)

## modelname <- getModelName('together_AR1_2')
## getJ('together_AR1_2') #set J
(modelname <- getModelName(variant))
smodelname <- modelname #stan modelname

## using ARk to specify k
if(grepl('_AR',modelname)){
  ## ark <- as.integer(substr(modelname,
  ##                          start = nchar(modelname),
  ##                          stop=nchar(modelname)))
  ark <- 1
  smodelname <- gsub(ark,'k',smodelname)
  D4S$ark <- ark
}


## common data:
D4S$J <- getJ(variant) #set J
D4S$h_beta_mu <- 1.5
D4S$h_beta_sigma <- 1
D4S$h_alpha_sigma <- 1/20

## different data:
wphi <- c("separate_global","separate_globalH",
          "separate_regional")
wdelt <- c(wphi,"separate_LerouxIntercept","together_LerouxIntercept")

if(modelname %in% wphi){
  D4S$h_phi_sigma <- 1
}

if(modelname %in% wdelt){
  D4S$h_delta_sigma <- 1/20
}

if(grepl('ARk',modelname)){
  D4S$ark <- 1
}

if(omitting){ #get iso3 code to omit for this job
  cat('modifying data \n')
  if(tallstring==''){ #LOOCV either country or most recent
    D4S <- countryomitter(omitcn,D4S,omittingT)
  } else { #omit most recent n years from everyone
    cat('..country loop ...\n')
    for(cn in ckey[,unique(iso3)]){ #country loop
      for(i in 0:omitcno){          #loop up to most recent time
        D4S <- countryomitter(cn,D4S,omittingT,i,verbose = FALSE) #omit i-th most recent year from country cn
      }
    }
  }
}
## stop()

## === compile model
fn <- gh('drtb/stan/model_{smodelname}.stan')
sm <- stan_model(file=fn)


## === run model:
niter <- 1e3
nchains <- 2
cores <- 1


tt1 <- system.time({
  samps <- sampling(sm,data = D4S,
                    iter = niter, chains = nchains,cores=cores,
                    verbose = FALSE)
})
tt1/60 #100


## save(samps,file=gh('stout/samps{variant}.Rdata'))
print(get_elapsed_time(samps)/60)
## load(file=gh('stout/samps{variant}.Rdata'))

fn <- gh('drtb/plots/t.{variant}{omitcn}{omitTstr}{tallstring}.txt')
cat(tt1[3]/60,file=fn)

K <- getRRprops(samps,D4S$J)
if(omitting){
  if(!omittingT){
    ## save means in most recent to evaluate wobble
    W <- K[dy==0,.(value=mean(value)),by=.(iso3,patients,dst)]
    save(W,file=gh('drtb/outdata/W.{variant}.{omitcn}.Rdata'))
    K <- K[iso3==omitcn]
  }
  fn <- gh('drtb/outdata/K{omitTstr}{tallstring}.{variant}.{omitcn}.Rdata')
  fn <- gsub('\\.\\.','\\.',fn) #simplest ends up with ".."; removed for backwards consistency
  save(K,file=fn)
  if(tallstring=='TA'){
    KO <- makeplotdata(K)
    save(KO,file=gh('drtb/outdata/KOTTA.{variant}.{omitcn}.Rdata'))
  }
} else { #full inference

  KO <- makeplotdata(K)
  save(KO,file=gh('drtb/outdata/KO.{variant}.Rdata'))

  fn <- gh('drtb/plots/r.{variant}')
  makeregionplots(KO,fn)

  ## === save summaries
  drws <- as_draws_df(samps)
  drws <- subset_draws(drws, c("beta", "alpha"))
  sdr <- summarize_draws(drws,mean,sd,rhat,ess_bulk)
  fn <- gh('drtb/outdata/sdr.{variant}.csv')
  fwrite(sdr,file=fn)

  nmz <- names(samps)
  nmz <- nmz[!grepl("\\[",nmz)] #only top level things
  smy <- summary(samps,pars=nmz)
  save(smy,file=gh('drtb/outdata/smy.{variant}.Rdata'))

  loo.out <- loo(samps, save_psis = FALSE)
  save(loo.out,file=gh('drtb/outdata/loo.out.{variant}.Rdata'))

  ll <- extract_log_lik(samps)
  waic.out <- waic(ll)
  save(waic.out,file=gh('drtb/outdata/waic.out.{variant}.Rdata'))

  lls <- rowSums(ll)
  dll <- list(data.likelihood.mean = mean(lls),
              data.likelihood.sd=sd(lls))
  save(dll,file=gh('drtb/outdata/dll.{variant}.Rdata'))
}

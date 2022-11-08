#' ---
#' title: Analysis of prevalence survey data
#' author: Pete Dodd
#' date: 3 July, 2022
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' 
#' #Pre-amble
#' (Last updated: `r Sys.Date()`)
#' 
#' This processes the prevalence survey data  to extract some patterns that can be used for age disaggregation
#'
#' N.B. This file should render to html with `rmarkdown::render('01prev_survey.R',output_dir='../html')` or from the command line with `R -q -e "rmarkdown::render(\"01prev_survey.R\",output_dir=\"../html\")"`
#'
#' #Initial inspection of prevalence survey data
#'
#' Relevant libraries for loading and manipulating data:
#'
rm(list=ls())
library(here)
source(here('disaggregation/R/utilities.R'))  #most dependencies in here

#' Data inputs loaded here
load(here('data/prev.rda'))          #prevalence survey data
load(here("data/pop.rda"))           #population data
load(here("data/tb.rda"))            #notification data
load(here("data/est.rda"))           #used for incidence method

#' Read in the updated prevalence survey data data
names(prev)

#' Some initial plots, first of all only looking at all sexes and skipping some of the young age categories and age categories that duplicate older adults:
ggplot(data=prev[sex=='a' &
                 age.group %ni% c('10plus',
                                  '15plus',
                                  '5plus',
                                  '50plus',
                                  '3monthplus'),],
       aes(x=age.group,y=prev,col=case.type)) +
  geom_errorbar(aes(ymin=prev.lo,ymax=prev.hi),width=0) + 
  geom_point() +
  ylab('Prevalence per 100,000')+
  facet_wrap(~iso3,scales = 'free') + 
  rot45

#' Note that:
#'
#' * there are data on both bacteriologically confirmed and smear-positive cases, sometimes together, sometimes not
#' * some countries use very different age-groups in reporting
#'
#'
#' ## Ad hoc changes to data
#'
#' We will make some changes to the data to deal with these problems
#'
#' To deal with the non-standard age categories, where countries use broader age categories, we will assume the same measurement in each subgroup that aggregates to a larger group. This ensures that aggregation would yield the same point estimate regardless of demography.
#'
#' This approach will apply to GMB, RWA:

## GMB
tmp <- prev[iso3=='GMB' & case.type=='b']
tmp <- tmp[age.group!='15plus',]
tmp2 <- copy(tmp)
tmp$age.group <- c('15_24','35_44','55_64')
tmp2$age.group <- c('25_34','45_54','65plus')
prev <- prev[!(iso3=='GMB' & sex=='a'),]    #drop these records
prev <- rbind(prev,tmp,tmp2)                  #join back on modified

## RWA
tmp <- prev[iso3=='RWA' & case.type=='b']
tmp <- tmp[age.group!='15plus',]
tmp2 <- copy(tmp)
tmp$age.group <- c('15_24','35_44','55_64')
tmp2$age.group <- c('25_34','45_54','65plus')
prev <- prev[!(iso3=='RWA' & sex=='a')]    #drop these records
prev <- rbind(prev,tmp,tmp2)                  #join back on modified

#' Further, PHL and KHM have some additional age categories from previous surveys, which will simply drop:
## KMH: we can keep only the most recent survey
prev <- prev[!(iso3=='KHM' & !year==2011)]
## PHL: again keep only the most recent survey
prev <- prev[!(iso3=='PHL' & year!=2016)]

#' For TZA, their only survey was done with smear. Relable this as bacteriologically confirmed for the age results. Then drop the age results that are smear-only (which will also drop the older Chinese surveys)
prev <- prev[!(iso3=='TZA' & case.type=='b')] #drop the overall B+
prev <- prev[!(iso3=='TZA' & age.group=='15plus' & sex=='a')] #drop the overall smr+
prev[iso3=='TZA' & sex=='a',case.type:='b']

#' Check new surveys
(newiso <- prev[year>=2018,unique(iso3)])
prev[iso3=='LSO']                       #OK
prev <- prev[!(iso3=='MMR' & year!=2018)] #drop older survey
prev[iso3=='MMR']
prev[iso3=='MOZ']                         #OK
prev[iso3=='NAM']                         #OK
prev[iso3=='NPL']                         #OK
prev[iso3=='SWZ']                         #OK
prev <- prev[!(iso3=='VNM' & year!=2018)] #drop older survey
prev[iso3=='VNM']                         #OK
prev[iso3=='ZAF']                         #OK


#'
#' ## Other preparation
#' 
#' Focus on age patterns and get rid of funny groupings and smr-only:
#' 
PR <- prev[sex=='a' &
        age.group %ni% c('10plus','15plus','5plus','50plus','3monthplus','all') &
        case.type=='b']
PR$case.type <- NULL                    #since no longer true

#' Now introduce a RR for being prevalence with respect to the youngest age category as a reference point
PR[,iso3.year:=paste(iso3,year)]        #introduce factor for each survey
PR[,RR:=.SD[,prev]/.SD[age.group=='15_24',prev],by=iso3.year] #calculate RR
PR[,prevse:=(prev.hi-prev.lo)/3.92]                           #prevalence SE

#' Want to include errors from reference category and others for $Q=X/Y$ as $var(Q)/\bar{Q}^2=(var(X)/\bar{X}^2+\var(Y)/\bar{Y}^2)^2$
PR[,RRse:=.SD[age.group=='15_24',prevse/prev],by=iso3.year] #reference group errors
PR[age.group=='15_24',RRse:=0]          #no error in reference group
PR[,RRse:=RR*sqrt(prevse^2/prev^2 + RRse^2)] #combine errors
PR[,RRlo:=pmax(RR - 1.96*RRse,0)]
PR[,RRhi:=RR + 1.96*RRse]

#' Check the prevalence data again
GP <- ggplot(data=PR,
             aes(x=age.group,y=prev)) +
  geom_errorbar(aes(ymin=prev.lo,ymax=prev.hi),width=0) + 
  geom_point() + 
  facet_wrap(~iso3.year,scales = 'free') +
  rot45
GP

ggsave(GP,file=here('disaggregation/output/prevsplits/PSinput.pdf'),
       w=26,h=15)

#' # Patterns with respect to age (first pass)
#'
#' Inspect the pattern of RRs
dodge <- position_dodge(.1)
ggplot(data=PR,
       aes(x=age.group,y=RR,col=iso3.year,group=iso3.year,ymin=RRlo,ymax=RRhi)) +
  geom_hline(yintercept=1,alpha=.2) +
  geom_errorbar(width=0,position=dodge) + 
  geom_point(position=dodge) + geom_line()+
  facet_wrap(~iso3.year)+
  ylab('Prevalence risk ratio') +
  rot45

#' Merge in WHO regions
regkey <- unique(pop[,.(iso3,g.whoregion)])
PR <- merge(PR,regkey,by='iso3',all.x=TRUE)
PR$age.group <- factor(PR$age.group)    #this wasn't a factor before


#' In terms of simple predictors of this pattern, there is some degree of cohesiveness by region:`g_whoregion`:
ggplot(data=PR,
       aes(x=age.group,y=RR,col=iso3.year,group=iso3.year)) +
  geom_errorbar(aes(ymin=RRlo,ymax=RRhi),width=0,position=dodge) + 
  geom_point(position=dodge) + geom_line()+
  facet_wrap(~g.whoregion)+
  geom_hline(yintercept=1,col='grey') +
  ylab('Prevalence risk ratio') + 
  rot45

#' Introduce the log of these RRs
PR[,logRR:=log(RR)]

#'
#' # Hierarchical approach to age patterns
#'
#' There are several shortcomings of the above approach:
#' 
#' * doesn't use measurement uncertainty
#' * ignores country specific data and only uses pooled
#' * no way of predicting for region without prevalence survey
#' * doesn't include prevalence surveys with non-standard age patterns
#'
#' Try to fix the first two by using a Normal Inverse Wishart (NIW) prior for unknown means and covariances in each region, together with a normal error model.
#'
#' This is described in more detail in this peer-reviewed artice:
#' 
#' https://pubmed.ncbi.nlm.nih.gov/33624797/
#'
#' The normal situation is to use an NIW as a prior when the data $y$ is generated from a MVN with unknown mean and variance:
#'
#' $$\mu, \Sigma \sim NIW(\mu_0,\lambda,\Psi,\nu) = NIW(\xi)$$
#' $$y_i \sim MVN(\mu,\Sigma)$$
#'
#' Then we have
#'
#' $$\mu, \Sigma | Y=\{y_i\}_{i=1}^n \sim NIW(\mu_n,\lambda_n,\Psi_n,\nu_n)$$
#'
#' where the update rules for the parameters $\xi$ are available on [wikipedia](https://en.wikipedia.org/wiki/Normal-inverse-Wishart_distribution), whose notation I use (expect for $\xi$ for all parameters collectively).
#'
#' Our situation is not this simple because we don't observe true prevalence and therefore true RRs, and so only have measurements for $z_i$ with known variances $E_i$. However we can use a Gibbs sampler that alternates:
#'
#' * Sample $Y|Z \& \mu, \Sigma$, which is normal (see below)
#' * Sample $\mu,\Sigma | Y$, which is NIW (by above)
#'
#' By completing the square, the Normal distribution of the first step above is that in a country, $i$, $y_i$ is MVN with
#'
#' * mean $(z_iE_i^{-1}+\mu\Sigma^{-1})S_i$
#' * variance $S_i$
#'
#' where $S_i^{-1} = E_i^{-1} + \Sigma^{-1}$. Below, we write some code to implement such a sampler.
#'
#'
#'
#' Normal inverse Wishart sampler (see link above):
getaNIW <- function(mu,lam,Psi,nu) {
    Sigh <- riwish(nu, Psi)
    muout <- mvrnorm(n=1,mu=mu,Sigma=Sigh) #single draw
    list(mu=muout,Sig=Sigh)
}
## getaNIW(mu0,1,Psi,nu)                   #test

#' We're going to work with  $Y$, which is $(n \times p)$ (underlying true RRs) where n=number of obervations and p=dimension. We want a function that can sample $Y$  given a noisy observation of $z$, with known measurement errors, and $\mu$ & $\Sigma$. For us $p=5$. A sampling function for single row/observation:
getaY <- function(z,z.se,mu,invSig){       #z a row of Z
    p <- length(z)
    Einv <- diag(1/z.se^2)              #observation precisions in country
    invS <- Einv + invSig
    S <- solve(invS)
    offset <- (z %*% Einv + mu %*% invSig) %*% S
    mvrnorm(n=1,mu=offset,Sigma=S)
}
## getaY(c(1,2),c(.1,.1),mu0,Psi)          #test

#' A function that generates a whole set of Y given a `cbind`-ed input of log(RR) and SEs 
getY <- function(Zboth,mu,Sig){
    invSig <- solve(Sig)
    ## Zboth <- cbind(Z,Z.se)
    p <- ncol(Zboth)/2
    ans <- apply(Zboth,
                 1,
                 FUN = function(x) getaY(x[1:p],x[(p+1):(2*p)],mu,invSig) )
    t(ans)
}
## getY( cbind( diag(rep(1,2)), matrix(.1,2,2) ), mu=mu0, Sig=Psi) #test
## getY(cbind(Z,Z.se),rep(0,5),diag(rep(1,5)))                     #test
## getaY(Z[3,],Z.se[3,],rep(0,5),diag(rep(1,5)))                   #test

#' The main function to do Gibbs sampling:
getGibbs <- function(Z,Z.se,nchain=1e2,
                     init=list(nu=5,lam=1,mu=rep(0,5),Psi=diag(rep(1,5)))){
    nu <- init$nu; lam <- init$lam; mu <- init$mu; Psi <- init$Psi;
    Zboth <- cbind(Z,Z.se)
    n <- nrow(Z); p <- ncol(Z)
    ## containers for keeping data
    mustore <- ystore <- pstore <- list()
    ## this bit of conjugate update is always the same so not below
    nu.n <- nu + n
    lam.n <- lam + n
    ## sample initial mu, Sig from prior:
    pr <- getaNIW(mu=mu,lam=lam,Psi=Psi,nu=nu)
    mu <- pr$mu ; Sig <- pr$Sig
    for(i in 1:nchain){
        if(!i%%1e2)cat('i=',i,'\n')
        ## sample Y | Z, mu, Sig
        Y <- getY(Zboth,mu=mu,Sig=Sig)
        ## sample mu Sig | Y
        ## conjugate update
        ybar <- colMeans(Y)
        mu.n <- (lam*mu + n*ybar)/(lam + n)
        ymybar <- Y - matrix(ybar,nrow=n,ncol=p,byrow = TRUE)
        S <- t(ymybar) %*% ymybar       #sum over data - pxp
        ybarmmuS <- (ybar-mu) %*% t(ybar-mu) # pxp
        Psi.n <- Psi + S + ybarmmuS*lam*n/lam.n
        pr <- getaNIW(mu=mu.n,lam=lam.n,Psi=Psi.n,nu=nu.n)
        mu <- pr$mu ; Sig <- pr$Sig
        ## store
        mustore[[i]] <- mu; ystore[[i]] <- Y; pstore[[i]] <- mvrnorm(1,mu,Sig)
        ## Sigstore[[i]] <- Sig
    }
    list(muz=mustore,Yz=ystore,ps=pstore)
}

save(PR,file=here('disaggregation/output/prevsplits/PR.Rdata')) #save out data to be used
load(file=here('disaggregation/output/prevsplits/PR.Rdata'))

#'
#' ## Application to AFR
#' 
#' Now try this out on some real data. Focus on AFR region:

## a function to get data from a region:
getZdata <- function(reg){
    Z.df <- dcast(data=PR[g.whoregion==reg,],formula=iso3.year ~ age.group,value.var='RR')
    Z.se.df <- dcast(data=PR[g.whoregion==reg,],formula=iso3.year ~ age.group,value.var='RRse')
    Z.df$`15_24` <- Z.se.df$`15_24` <- NULL                    #drop ref cat
    Z <- as.matrix(Z.df[,-1])
    Z.se <- as.matrix(Z.se.df[,-1])
    Z.se <- 0.5 * log((Z+Z.se)/(Z-Z.se))
    Z <- log(Z)
    list(Z=Z,Z.se=Z.se,Z.df=Z.df)
}

## make data in required form for Africa
ZZ <- getZdata('AFR')
Z <- ZZ$Z; Z.se <- ZZ$Z.se; Z.df <- ZZ$Z.df

#' Apply the Gibbs sampler to this:
psi0 <- matrix(1,5,5)+diag(rep(.5,5))
out <- getGibbs(Z=Z,Z.se=Z.se,nchain=1e3,
                init=list(nu=10,lam=1,mu=rep(0,5),Psi=psi0))

#' Introduce a container for storing all of these results
#'
agedisag <- list()

#' Assess convergence of group means:
MU <- do.call('rbind',out$mu)
colnames(MU) <- names(Z.df)[-1]
lattice::xyplot(mcmc(MU),layout=c(2,3))

#' and of individual country means.
YZ <- do.call('rbind',lapply(out$Yz,function(x)x[1,])) #first country
colnames(YZ) <- names(Z.df)[-1]
lattice::xyplot(mcmc(YZ),layout=c(2,3))

#' and of posterior predictions.
PZ <- do.call('rbind',out$ps)
colnames(PZ) <- names(Z.df)[-1]
lattice::xyplot(mcmc(PZ),layout=c(2,3))
lattice::splom(PZ)

#' store for later:
agedisag[['AFR']] <- list(isoz = unlist(lapply(strsplit(unlist(Z.df[,1]),split=" "),function(x)x[1])),Ylist = out$Yz,Pmat=PZ)


#' To compare against the data, here are a couple of functions to reformat and plot:
makeplotdata <- function(LoV,isoy=NA,reg='AFR',n=NULL){
    df <- do.call('rbind',LoV)
    colnames(df) <- colnames(Z)
    if(!is.null(n)) if(n<nrow(df)) df <- df[sample(nrow(df),n),]
    df <- as.data.frame(df)
    df[['15_24']] <- 0
    df$n <- 1:nrow(df)
    dfm <- reshape2::melt(df,id.vars='n')
    names(dfm)[2:3] <- c('age.group','logRR')
    dfm$age.group <- factor(dfm$age.group,
                            levels=c('15_24',"25_34","35_44",
                                     "45_54","55_64","65plus"),
                            ordered=TRUE)
    if(!is.null(nrow(isoy))) isoy <- unlist(isoy)[1]
    dfm$iso3.year <- isoy
    dfm$g.whoregion <- reg
    dfm$wt <- 1
    dfm
}

makeregionplot <- function(xtradata,reg='AFR',alph=.1){
 PP <-  ggplot(data=PR[g.whoregion==reg],
               aes(x=age.group,y=logRR,col=iso3.year,fill=iso3.year,
                   group=g.whoregion)) + 
     facet_wrap(~g.whoregion)+
     ylab('log prevalence risk ratio') +
     theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
     stat_summary(data=xtradata,geom="ribbon",alpha=.2,col=NA,
                  fun.min = function(x) quantile(x, 0.025),
                  fun.max = function(x) quantile(x, 0.975)) +
     stat_summary(data=xtradata,geom="line", fun=median) +
     geom_point()#position=dodge)
 PP
}


#' Have a look at outputs for a specfic country with a prevalence survey:
cn <- 1                                 #specifc country
topym <- makeplotdata(lapply(out$Yz,function(x)x[cn,]),isoy=Z.df[cn,1],reg='AFR')
makeregionplot(topym,'AFR',alph=.1)

#' and then a different country with a prevalence survey:
cn <- 11                                 #specifc country
topym <- makeplotdata(lapply(out$Yz,function(x)x[cn,]),isoy=Z.df[cn,1],reg='AFR')
makeregionplot(topym,'AFR')

#' A country in the region without a prevalence survey would use the predictive patterns, which have more variability:
PZdfm <- makeplotdata(out$ps,reg='AFR')
makeregionplot(PZdfm,'AFR')
save(PZdfm,file=here('disaggregation/output/prevsplits/PZdfm.Rdata'))

## lattice::splom(PZ)


#'
#' ## Corresponding checks in SEA
#'
#'
#'

ZZ <- getZdata('SEA')
Z <- ZZ$Z; Z.se <- ZZ$Z.se; Z.df <- ZZ$Z.df

#' Apply the Gibbs sampler to this:
out <- getGibbs(Z=Z,Z.se=Z.se,nchain=1e3,
                init=list(nu=10,lam=1,mu=rep(0,5),Psi=psi0))

#' Assess convergence of group means:
MU <- do.call('rbind',out$mu)
colnames(MU) <- names(Z.df)[-1]
lattice::xyplot(mcmc(MU),layout=c(2,3))

#' and of individual country means.
YZ <- do.call('rbind',lapply(out$Yz,function(x)x[1,])) #first country
colnames(YZ) <- names(Z.df)[-1]
lattice::xyplot(mcmc(YZ),layout=c(2,3))

#' and of posterior predictions.
PZ <- do.call('rbind',out$ps)
colnames(PZ) <- names(Z.df)[-1]
lattice::xyplot(mcmc(PZ),layout=c(2,3))

#' store for later:
agedisag[['SEA']] <- list(isoz = unlist(lapply(strsplit(unlist(Z.df[,1]),split=" "),function(x)x[1])),Ylist = out$Yz,Pmat=PZ)

#' Have a look at outputs for a specfic country with a prevalence survey:
cn <- 1                                 #specifc country
topym <- makeplotdata(lapply(out$Yz,function(x)x[cn,]),isoy=Z.df[cn,1],reg='SEA')
makeregionplot(topym,'SEA',alph=.1)

#' and then a different country with a prevalence survey:
cn <- 4                                 #specifc country
topym <- makeplotdata(lapply(out$Yz,function(x)x[cn,]),isoy=Z.df[cn,1],reg='SEA',2e2)
makeregionplot(topym,'SEA')

#' A country in the region without a prevalence survey would use the predictive patterns, which have more variability:
PZdfm <- makeplotdata(out$ps,reg='SEA')
makeregionplot(PZdfm,'SEA',alph=5e-2)


#'
#' ## Corresponding checks in WPR
#' 
#' 
ZZ <- getZdata('WPR')
Z <- ZZ$Z; Z.se <- ZZ$Z.se; Z.df <- ZZ$Z.df

#' Apply the Gibbs sampler to this:
out <- getGibbs(Z=Z,Z.se=Z.se,nchain=1e3,
                init=list(nu=10,lam=1,mu=rep(0,5),Psi=psi0))

#' Assess convergence of group means:
MU <- do.call('rbind',out$mu)
colnames(MU) <- names(Z.df)[-1]
lattice::xyplot(mcmc(MU),layout=c(2,3))

#' and of individual country means.
YZ <- do.call('rbind',lapply(out$Yz,function(x)x[1,])) #first country
colnames(YZ) <- names(Z.df)[-1]
lattice::xyplot(mcmc(YZ),layout=c(2,3))

#' and of posterior predictions.
PZ <- do.call('rbind',out$ps)
colnames(PZ) <- names(Z.df)[-1]
lattice::xyplot(mcmc(PZ),layout=c(2,3))

#' store for later:
agedisag[['WPR']] <- list(isoz = unlist(lapply(strsplit(unlist(Z.df[,1]),split=" "),function(x)x[1])),Ylist = out$Yz,Pmat=PZ)

#' Have a look at outputs for a specfic country with a prevalence survey:
cn <- 1                                 #specifc country
topym <- makeplotdata(lapply(out$Yz,function(x)x[cn,]),isoy=Z.df[cn,1],reg='WPR')
makeregionplot(topym,'WPR',alph=.1)

#' and then a different country with a prevalence survey:
cn <- 3                                 #specifc country
topym <- makeplotdata(lapply(out$Yz,function(x)x[cn,]),isoy=Z.df[cn,1],reg='WPR',2e2)
makeregionplot(topym,'WPR')

#' A country in the region without a prevalence survey would use the predictive patterns, which have more variability:
PZdfm <- makeplotdata(out$ps,reg='WPR')
makeregionplot(PZdfm,'WPR')

#'
#' ## Corresponding checks in EMR
#' 
#' 
ZZ <- getZdata('EMR')
Z <- ZZ$Z; Z.se <- ZZ$Z.se; Z.df <- ZZ$Z.df

#' Apply the Gibbs sampler to this (using more informative prior -- see nu):
out <- getGibbs(Z=Z,Z.se=Z.se,nchain=1e3,
                init=list(nu=10,lam=1,mu=rep(0,5),Psi=psi0))

#' Assess convergence of group means:
MU <- do.call('rbind',out$mu)
colnames(MU) <- names(Z.df)[-1]
lattice::xyplot(mcmc(MU),layout=c(2,3))

#' and of individual country means.
YZ <- do.call('rbind',lapply(out$Yz,function(x)x[1,])) #first country
colnames(YZ) <- names(Z.df)[-1]
lattice::xyplot(mcmc(YZ),layout=c(2,3))

#' and of posterior predictions.
PZ <- do.call('rbind',out$ps)
colnames(PZ) <- names(Z.df)[-1]
lattice::xyplot(mcmc(PZ),layout=c(2,3))

#' store for later:
agedisag[['EMR']] <- list(isoz = unlist(lapply(strsplit(unlist(Z.df[,1]),split=" "),function(x)x[1])),Ylist = out$Yz,Pmat=PZ)


#' Have a look at outputs for a specfic country with a prevalence survey:
cn <- 1                                 #specifc country
topym <- makeplotdata(lapply(out$Yz,function(x)x[cn,]),isoy=Z.df[cn,1],reg='EMR')
makeregionplot(topym,'EMR',alph=.1)

#' and then a different country with a prevalence survey:
cn <- 2                                 #specifc country
topym <- makeplotdata(lapply(out$Yz,function(x)x[cn,]),isoy=Z.df[cn,1],reg='EMR',2e2)
makeregionplot(topym,'EMR')

#' A country in the region without a prevalence survey would use the predictive patterns, which have more variability:
PZdfm <- makeplotdata(out$ps,reg='EMR')
makeregionplot(PZdfm,'EMR')

#'
#' # Patterns with respect to sex
#'
ggplot(data=prev[sex%in%c('f','m'),],
       aes(x=sex,y=prev,col=case.type,pch=age.group)) +
    geom_errorbar(aes(ymin=prev.lo,ymax=prev.hi),width=0) +
    geom_point() +
    facet_wrap(~factor(paste0(iso3,'.',year)),scales = 'free')


#' These data give a number of direct estimates of MF ratios:
#'
##make MF ratios
MFD <- prev[sex%in%c('f','m'),]
MFD[,iso3.year:=paste(iso3,year)]
MFD <- MFD[iso3.year %ni% c('CHN 1990','CHN 2000','BGD 2008',
                            'IDN 2004','KOR 1990','MMR 1994'),] #keep most recent
MFD[,table(case.type,iso3)]

MFD[iso3=='TZA',case.type:='b'] #keep this as only available
MFD <- MFD[case.type!='s',]                        #ditch rest
## drop korea as different age group and no precision
MFD <- MFD[age.group=='15plus',]
MFD[,MF:=prev/.SD[sex=='f',prev],by=iso3.year]
MFD[,MFse:=(prev.hi-prev.lo)/(3.92*prev)] #fractional SEs
MFD[,MFse:=MF*sqrt(sum(MFse^2)),by=iso3.year] #sum and multiply
MFD <- MFD[sex=='m',]                      #drop female
MFD[,sig2:=log(1+MFse^2/MF^2)]                #for LN
MFD[,mu:=log(MF)-sig2/2]                      #for LN
setkey(MFD,'iso3')
kable(MFD[,.(iso3.year,prev,MF,MFse,mu,sig2)])
#' In this last step, we have calculated parameters for modelling the MF RRs as following a log-normal distribution.

#'
#'
#' In addition, Katherine Horton's [systematic review](http://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1002119) includes regional estimates:

KH <- data.table(g.whoregion=c('SEA','WPR','AFR','EMR','AMR','GLO'),
                 MF=c(3.37,1.93,1.73,1.49,0.74,2.21),
                 MFlo=c(2.83,1.37,1.54,1.2,0.42,1.92),
                 MFhi=c(4.01,2.72,1.95,1.84,1.33,2.54))
KH[,sig2:=log(1+(MFhi-MFlo)^2/(3.92*MF)^2)]
KH[,mu:=log(MF)-sig2/2]                 #LN parms
setkey(KH,'g.whoregion')
kable(KH)

#' However, she also reports an age dependence of the MF ratio (Figure 5). Let the age and sex ratios (with respect to the youngest females) be $\rho_{a,s}$. For a particular age group $MF_a=\rho_{a,M}/\rho_{a,F}$. For simplicity, let $\rho_{a,F}=\rho_a$, and let $\rho_0=1$. The RRs we model above depend on the demography:
#'
#' $$RR_a = \frac{(MF_aN_{a,m}+N_{a,f})\rho_{a}/N_a}{(MF_0N_{0,m}+N_{0,f})/N_0}$$
#' $$MF = \frac{\sum_aMF_aN_{a,m}\rho_{a}N_m}{\sum_aN_{a,f}\rho_{a}N_f}$$
#'
#'
#' At first pass, assume a multiplicative effect, i.e. $MF_a=MF_0$. Then given an $RR_a$ and $MF$ for a country, need to solve these equations for $\rho_a$ and $MF_0$:
#'
#' $$RR_a = \frac{(MF_0.N_{a,m}+N_{a,f})\rho_{a}/N_a}{(MF_0N_{0,m}+N_{0,f})/N_0}$$
#' $$MF = MF_0\frac{\sum_aN_{a,m}\rho_{a}/N_m}{\sum_aN_{a,f}\rho_{a}/N_f}$$
#'
#' Can approximately solve this by iteration. Having obtained $MF_0$ and $\rho_a$, the age/sex disaggregation of (adult) prevalence is:
#'
#'
#' $$P_{a,s}=\frac{(MF_0)^{\delta_{s,m}}\rho_a N_{a,s}}{\sum_{b,t}(MF_0)^{\delta_{t,m}}\rho_b N_{b,t}}$$
#' 
#'
#' # Obtaining age/sex disaggregations for prevalence
#'
#' In this section we develop a set of functions to generate age/sex for a given country.
#'
#' This will need a few extra data files:

agz <- c('1524','2534','3544','4554','5564','65') #some age ranges for use below

#' Check incidence sources
est[,unique(source.inc)]

#' First sampling a set of probabilities for a given country (in adults)
SampleRRs <- function(n=1,iso3,verbose=TRUE){
  if(verbose) cat('Country=',iso3,'\n')
  reg <- as.character(est[iso3,g.whoregion][1]) #look up region
  if(verbose) cat('Region=',reg,'\n')
  srcinc <- est[iso3]
  srcinc <- srcinc[year==estyr,source.inc]
  print(srcinc)
  if(is.na(srcinc) |
     !grepl('Standard adjustment',srcinc) & reg %ni% c('EUR','AMR')){
    ##o/w use notifications
    ## first, sample some MFs
    if(iso3 %in% MFD$iso3){             #does the country have a prev survey?
      if(verbose) cat('Prevalence survey used for MF ratio. \n')
      mfs <- MFD[iso3,rlnorm(n=n,meanlog=mu,sdlog=sqrt(sig2))] #sample MFs
    } else {              #use Horton regional data
      if(verbose) cat('SR/MA used for MF ratio. \n')
      mfs <- KH[reg,rlnorm(n=n,meanlog=mu,sdlog=sqrt(sig2))]
    }
    ## sample RRa's
    if( reg %in% c('AFR','SEA','WPR','EMR') ){ #hierarchical model
      if(verbose) cat('Hierarchical model used for age RRs. \n')
      if(iso3 %in% agedisag[[reg]]$isoz){     #RR from prevalence survey
        if(verbose) cat('(Prevalence survey detected for age RRs.) \n')
        k <- which(agedisag[[reg]]$isoz==iso3)
        RRs <- lapply(agedisag[[reg]]$Ylist,function(x)x[k,]) #rows for this iso3
        RRs <- do.call('rbind',RRs)
      } else {                        #RR from regional prediction
        RRs <- agedisag[[reg]]$Pmat
      }
      RRs <- RRs[sample(1:nrow(RRs),n,replace=TRUE),] #sample n
      RRs <- exp(RRs)                 #transform
      if(n>1){
        RRs <- cbind(rep(1,n),RRs)      #add reference
      } else RRs <- c(1,RRs)
    } else {
      RRs <- matrix(nrow=n,ncol=6)
      cat('Oops! No method currently defined!\n')
    }
    list(MFs=mfs,RRs=RRs)
  } else {       # ---------- use notifications!
    if(verbose) cat('Using notification data for age/sex disaggregation.\n')
    notes <- matrix(nrow=2,ncol=6)
    yr <- estyr
    if(is.na(sum(tb[iso3==iso3&year==yr,paste0('newrel.m',agz),with=FALSE]))){
      yr <- (estyr-1)
      cat('Check notification data! Using last year due to NAs!\n')
    }
    here <- tb[iso3]
    tmp <- unlist(here[year==yr,paste0('newrel.m',agz),with=FALSE])
    notes[1,] <- tmp
    tmp <- unlist(here[year==yr,paste0('newrel.f',agz),with=FALSE])
    notes[2,] <- tmp
    notes <- notes/sum(notes)
    colnames(notes) <- c('15_24','25_34','35_44','45_54','55_64','65plus')
    rownames(notes) <- c('M','F')
    list(P=notes)
  }
}

## checks
test <- SampleRRs(1,'PAK')              #EMR, PS
test <- SampleRRs(1,'AFG')              #EMR, no PS
test <- SampleRRs(1,'ZWE')              #AFR, PS
test <- SampleRRs(1,'ZAF')              #AFR, no PS
test <- SampleRRs(1,'GBR')              #notifications
test <- SampleRRs(1,'ABW')              #std adj
test <- SampleRRs(1,'BOL')              #AMR problem before as was looking for PS


#' In the case where the disaggregation is not from notifications, this function accounts for demography and adapts the underlying MF and age RRs to match the observed marginals, assuming no interaction (see above). 
#'
correctRR <- function(cn,MF,RR,verbose=FALSE){
  Demo <- matrix(nrow=2,ncol=6)       #demography
  Demo[1,] <- unlist(pop[iso3==cn & year==estyr,paste0('e.pop.m',agz),with=FALSE])
  Demo[2,] <- unlist(pop[iso3==cn & year==estyr,paste0('e.pop.f',agz),with=FALSE])
  atot <- colSums(Demo)               #totals by age
  mf <- MF                            #underlying MF 1st guess
  rho <- RR                           #underlying RR 1st guess
  for(i in 1:10){                     #seems to work
    RRa <- (mf*Demo[1,] + Demo[2,])*rho/atot 
    RRa <- RRa/RRa[1] # based on application to demog
    MFc <- mf * sum(rho*Demo[1,]) / sum(rho*Demo[2,])
    mf <- MF * sum(rho*Demo[2,]) / sum(rho*Demo[1,]) #rescale
    rho <- RR*(atot/atot[1])*(mf*Demo[1,1]+Demo[2,1])/(mf*Demo[1,]+Demo[2,])
  }
  if(verbose){ print(MFc); print(RRa);}
  list(mf=mf,rho=rho)
}

## and a wrapped version for operating on multiple inputs:
correctRRs <- function(cn,MFs,RRs){
    if(length(MFs)==1) return(list(correctRR(cn,MFs,RRs)))
    lapply(1:length(MFs),FUN=function(i) correctRR(cn,MFs[i],RRs[i,]) )
}


test <- SampleRRs(1,'PAK')              #EMR, PS
test$MFs
test$RRs
correctRR('PAK',test$MFs[1],test$RRs,verbose = TRUE) 

test <- SampleRRs(10,'PAK')              #EMR, PS
correctRR('PAK',test$MFs[1],test$RRs[1,],verbose = TRUE)
correctRRs('PAK',test$MFs[1],test$RRs[1,])
test3 <- correctRRs('PAK',test$MFs,test$RRs)

#' This function converts these prevalence disaggregations to incidence disaggregations
#'
## for a specific country on a single sample
getaProp <- function(mf,rho,Demo){
  Demo[1,] <- mf*rho*Demo[1,]
  Demo[2,] <- rho*Demo[2,]
  Demo <- Demo/sum(Demo)
  colnames(Demo) <- c('15_24','25_34','35_44','45_54','55_64','65plus')
  rownames(Demo) <- c('M','F')
  Demo
}

## for a specific country on multiple samples
getProps <- function(n=1,cn,verbose=TRUE,shape=FALSE){
  smps <- SampleRRs(n=n,iso3=cn,verbose=verbose)     #observed RRs
  if('P' %in% names(smps)){
    ans <- smps$P
    if(shape){
      ans <- as.data.table(ans)
      ans[,rep:=1]; ans[,sex:=c('M','F')]
      ans <- melt(ans,id=c('rep','sex'))
    }
  } else {
    L <- correctRRs(cn=cn,MFs = smps$MFs,RRs=smps$RRs) #list underlying RRs for work
    Demo <- matrix(nrow=2,ncol=6)       #demography
    Demo[1,] <- unlist(pop[iso3==cn & year==estyr,
                           paste0('e.pop.m',agz),with=FALSE])
    Demo[2,] <- unlist(pop[iso3==cn & year==estyr,
                           paste0('e.pop.f',agz),with=FALSE])
    ans <- lapply(L,function(x) getaProp(mf=x$mf,rho=x$rho,Demo=Demo))
    if(shape){   #convert from list to data.table
      ans <- do.call('rbind',ans)
      rownames(ans) <- 1:nrow(ans)    #to avoid warning
      ans <- cbind(rep=rep(1:(nrow(ans)/2),each=2),
                   sex=rep(1:2,(nrow(ans)/2)),ans)
      ans <- as.data.frame(ans)
      ans <- reshape2::melt(ans,id=c('rep','sex'))
      ans$sex <- factor(c('M','F')[ans$sex])
      ans <- as.data.table(ans)
    }
  }
  ans
}

## test
getProps(n=50,'PAK',verbose=TRUE,shape=TRUE)
getProps(n=5,'GBR',verbose=TRUE,shape=TRUE)
getProps(n=5,'ABW',verbose=TRUE,shape=TRUE)

#'
#' ## Applying this method to countries
#'
#' 
#' Loop through countries in simple way, collect the data and compute mean and SEs:
#'
prevsplits <- prevsplitsL <- prevsrc <- list() #containers
for(cn in unique(est$iso3)){
  if(nrow(est[iso3==cn & year==estyr,])>0){ #skip things without estimates
      txt <- capture.output( Pdt <- getProps(n=1e2,
                                             cn=cn,
                                             verbose = TRUE,
                                             shape=TRUE) )
    print(prevsrc[[cn]] <- txt[1:5])
    Pdt[,iso3:=cn]
    prevsplitsL[[cn]] <- Pdt
    Pdt <- Pdt[,.(prop=mean(value),prop.se=sd(value)),
               by=.(sex,age.group=variable,iso3)]
    prevsplits[[cn]] <- Pdt
  }
}

prevsplits <- do.call('rbind',prevsplits)
prevsplitsL <- do.call('rbind',prevsplitsL)
prevsplits[is.na(prop),unique(iso3)]

cat(prevsplits[is.na(prop),unique(iso3)],
    file=here('disaggregation/output/prevsplits/zproblems.txt'))

#' save out results:
save(file=here('disaggregation/output/prevsplits/prevsplits.Rdata'),
     prevsplits)
save(file=here('disaggregation/output/prevsplits/prevsrc.Rdata'),prevsrc)
fn <- here('disaggregation/output/prevsplits/datal')
if(!file.exists(fn)) dir.create(fn)
save(file=here(fn,'prevsplitsL.Rdata'),prevsplitsL)

#' A function for outputting a given country's age/sex split
#'
countryplot <- function(cn){
  tmp <- prevsplits[iso3==cn,]
  tmp[is.na(prop.se),prop.se:=0]
  tmp[is.na(prop),prop:=0]
  ans <- ggplot(data=tmp,aes(x=age.group,y=prop,fill=sex)) +
    coord_flip() +
    geom_bar(data=tmp[sex=='M',],stat='identity',aes(y=prop))+
    geom_bar(data=tmp[sex=='F',],stat='identity',aes(y=-prop)) +
    geom_errorbar(data=tmp[sex=='M',],aes(ymin=prop-1.96*prop.se,ymax=prop+1.96*prop.se),width=0)+
    geom_errorbar(data=tmp[sex=='F',],aes(ymin=-prop-1.96*prop.se,ymax=-prop+1.96*prop.se),width=0) + ggtitle(paste0('Adult TB prevalence split: ',cn))
  ans
}

## test
countryplot('PAK')
countryplot('GBR')

#' Various country plots are now saved in directory prevsplitplots
#' loop & save out (SLOW)
fn <- here('disaggregation/output/prevsplits/plots')
if(!file.exists(fn)) dir.create(fn)
for(cn in unique(prevsplits$iso3)){
  if(!is.na(prevsplits[iso3==cn,sum(prop)])){
    plt <- countryplot(cn)
    fn <- glue(here("disaggregation/output/prevsplits/plots")) +
        "/" + cn + ".pdf"
    suppressMessages(ggsave(filename = fn,plot=plt))
  }
}





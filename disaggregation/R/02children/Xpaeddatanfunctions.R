## first go at putting together some estimates
rm(list=ls())
## packages---------------------------------------
library(here)
library(ggplot2);theme_set(theme_bw())  #change theme in release version
library(ggrepel)
library(scales)
library(grid)
library(data.table)
library(lhs)

## some requisite data-------------------------------------------
load(here('disaggregation/R/02children/data/LAT.Rdata'))  #latitude
load(here('data/pop.rda'))    #demography
load(here('data/unaids.rda'))     #HIV prev

## update here
estyr <- 2021                              # update
load(here('disaggregation/R/02children/data/eprev.Rdata')) #TB prevalence old
load(here('disaggregation/R/02children/data/estold.rda'))#TB est old
estold <- copy(est)
load(here('data/est.rda'))



## make BCG data if not there
BCGfn <- here('disaggregation/R/02children/data/BCG.Rdata')
if(!file.exists(BCGfn)){
  BCG <- fread(here('disaggregation/R/02children/data',
                    'wuenic2021rev_web-update.csv'),skip=1)
  BCG <- BCG[,.(iso3=V2,bcg=as.numeric(V5))]
  save(BCG,file=BCGfn)
} else{
  load(BCGfn)
}


##  update to old prevalence estimates
rat <- est[year==estyr][,.(iso3,inc)]
estold <- estold[year==2016][,.(iso3,oldinc=inc)]
rat <- merge(rat,estold)
rat[,rat:=inc/oldinc]
rat[!is.finite(rat),rat:=1]
rat[,summary(rat)]
eprev <- merge(eprev,rat[,.(iso3,rat)])
eprev[,c('prev','prev.sd','prev.lo','prev.hi'):= .(prev*rat,prev.sd*rat,prev.lo*rat,prev.hi*rat)] #scale
eprev[,rat:=NULL]
eprev[,year:=estyr]                      #to avoid errors

## TB prevalence
TBr <- eprev

## new demography -- NB slightly less granular than previous
Dr <- pop[year==estyr,.(country,
                       `<1`=(e.pop.f04+e.pop.m04)/5,
                       `1`=(e.pop.f04+e.pop.m04)/5,
                       `2-4`=3*(e.pop.f04+e.pop.m04)/5,
                       `5-9`=(e.pop.f514+e.pop.m514)/2,
                       `10-14`=(e.pop.f514+e.pop.m514)/2,
                       iso3)]

## key for regions
WHOkey <- unique(pop[,.(country,iso2,iso3,g_whoregion=g.whoregion)])

## bcg
bcg <- merge(BCG[,.(iso3,coverage=bcg)],WHOkey,by='iso3',all.y=TRUE)
bcgm <- bcg[,mean(coverage,na.rm=TRUE)]
bcg[is.na(coverage),coverage:=bcgm]     #NA to mean


## HIV prevalence
H <- unaids[year==estyr,.(iso3,hivu15=hiv014,hivlo=hiv014.lo,
                         hivhi=hiv014.hi,pop014=pop014)]

## ## --- modifications of data ---
setkey(Dr,'iso3')
setkey(TBr,'iso3')
setkey(H,'iso3')
setkey(bcg,'iso3')


## gamma: var=shape/rate^2, mean=shape/rate (sd=sqrt(shape)/rate)
## shape=(mean/sd)^2, rate = mean/sd^2
TBr[,shpz:=(prev/((prev.hi-prev.lo)/3.92))^2]
TBr[,rtz:=prev/((prev.hi-prev.lo)/3.92)^2]

## correct
miss <- which(is.na(TBr$shpz))
TBr[miss,]                              #0 prev
miss <- which(is.na(TBr$rtz))
TBr[miss,]
TBr[is.na(shpz),shpz:=0]
TBr[is.na(rtz),rtz:=1]                  #prevalence O(1)

## HIV
## H <- merge(H,TBr[,.(iso3)],all.y=TRUE)
H[,shpzh:=(hivu15/((hivhi-hivlo)/3.92))^2]
H[,rtzh:=hivu15/((hivhi-hivlo)/3.92)^2]
H[!is.finite(rtzh),rtzh:=Inf]
H[!is.finite(shpzh),shpzh:=1]
## H

dim(Dr);dim(bcg);dim(H); dim(TBr)
## ---- do some merging to ensure aligned -----
## LAT
Dr <- merge(Dr,LAT[,c('iso3','LAT')],all.x=TRUE,by='iso3')
summary(Dr$LAT)
Dr[is.na(LAT),.(iso3,country)]          #missing LAT, will need to revisit these
Dr <- Dr[!is.na(LAT)]                  #180 countries
## bcg
Dr <- merge(Dr,bcg[,.(iso3,coverage)],all.x=TRUE,by='iso3') #merge BCG
tmp <- Dr[is.na(coverage),.(iso3,country)]
if(nrow(tmp)>0){     #probably no big contributors; use mean
  print(tmp)
  mbcg <- mean(Dr$coverage,na.rm=TRUE); Dr[is.na(coverage),coverage:=mbcg]
}
## HIV
Dr <- merge(Dr,H,all.x=TRUE,by='iso3')
Dr[is.na(hivu15),hivu15:=0]             #zero where missing
## TB
Dr <- merge(Dr,TBr,all.x=TRUE,by='iso3')

## extra NA handling after merge
Dr[!is.finite(rtzh),rtzh:=1e6]
Dr[!is.finite(shpzh),shpzh:=1]
Dr[!is.finite(hivlo),hivlo:=0]
Dr[!is.finite(hivhi),hivhi:=0]
Dr[!is.finite(pop014),
   pop014:= as.integer(`<1` + `1` + `2-4` + `5-9` + `10-14`)]

## NOTE make sure lined up!

## making parameters--------------------------------------------
## s = E(1-E)/V - 1; a = sE, b=s(1-E); V = (u-l)/4
getAB <- function(E,L,U){
    V <- (U-L)/4
    sz <- E*(1-E)/V - 1
    a <- sz*E
    b <- sz*(1-E)
    return(list(a=a,b=b))
}

smpm <- 0.477# taken as previous value

(NH <- Dr[,sum(hivu15>0)])                   #number of HIV countries

## progression to disease
## prob TB
Ez <- c(0.5,0.25,0.05,0.02,0.15)
uz <- 1.25*Ez                           #+/- 25%
lz <- 0.75*Ez
drs1 <- getAB(Ez,lz,uz)$a
drs2 <- getAB(Ez,lz,uz)$b

## progression to DTB (proportion of TB)
Ezd <- c(0.3,0.35,0.1,0.125,0.0167)
uzd <- c(0.4,0.5,0.125,0.15625,0.021)
lzd <- c(0.2,0.2,0.075,0.09375,0.0125)
drs1d <- getAB(Ezd,lzd,uzd)$a
drs2d <- getAB(Ezd,lzd,uzd)$b

## protection against DTB
vds1 <- 1.25
vds2 <- 2.5

## fraction of protection for PTB
vs1 <- 4
vs2 <- 1

## beta risk
bmlg <- 1.677891 #log(4)
bsd <- 0.3714703 #1

## HIV IRR -- 
## 7.9 (95%CI: 4.5-13.7)
## hmlg <- log(7.9)
## hsd <- abs(log(4.5)-log(13.7))/3.92
## changed to new value -- see mortality piece
hmlg <- 4.36
hsd <- 1.14


## metaparameters for generating the replicate parameters
metaparm <- list(drs1 = drs1,drs2 = drs2,
                 drs1d = drs1d,drs2d = drs2d, #progression
                 vds1=vds1,vds2=vds2,vs1=vs1,vs2=vs2, #BCG
                 bmlg=bmlg,bsd=bsd,                   #beta
                 hmlg=hmlg,hsd=hsd,                   #HIV IRR
                 latprop = 0.41,  #propn variation due to latitude
                 relsmr = 0.23   #relative inf of smr-
                 )



## note: this relies on same configuration among inputs
## country, <1,1,2-4,5-9,10-14
## bcg in percent: country, coverage
## change to expect DT input
getIncidence <- function(Dr,parms,ARI,HIV,method='lat'){
    bcg <- Dr[,coverage]
    LAT <- Dr[,LAT]
    PTBu <- PTBv <- DTBu <-
        DTBv <- as.data.frame(Dr[,.(iso3,`<1`,`1`,`2-4`,`5-9`,`10-14`)])
    L <- rep(1,dim(Dr)[1])                #latitude
    if(method=='lat'){L <- L - parms$latprop * (1-LAT/90)} #proportion of protection left @ lat
    pp <- 1-parms$vBCGp                   #PTB protection by age
    dp <- 1-parms$vBCGd                   #DTB protection by age
    for(i in 1:5){
        PTBv[,i+1] <- PTBv[,i+1] * parms$prPTB[i] * (1-pp[i]*L) * (bcg/100)
        PTBu[,i+1] <- PTBu[,i+1] * parms$prPTB[i] * (1-bcg/100)
        DTBv[,i+1] <- DTBv[,i+1] * parms$prDTB[i] * (1-dp[i]*L) * (bcg/100)
        DTBu[,i+1] <- DTBu[,i+1] * parms$prDTB[i] * (1-bcg/100) 
    }
    DTBu[,2:6] <- ARI * (DTBu[,2:6] + DTBv[,2:6]) * 1e0 * (1-HIV+HIV*parms$hivIRR)
    PTBu[,2:6] <- ARI * (PTBu[,2:6] + PTBv[,2:6]) * 1e0 * (1-HIV+HIV*parms$hivIRR)
    return(list(PTB=PTBu,DTB=DTBu))
}

## for getting the LTBI prevalence and exposures
getLnE <- function(Dr,ARI,method){
  atrisk <- Dr
  if(method=='prev'){                   #convert to prevalence
    ariz <- matrix(c( 1-exp(-.5*ARI),
                     1-exp(-1.5*ARI),
                     1-exp(-4.0*ARI),
                     1-exp(-7.5*ARI),
                     1-exp(-12.5*ARI)),
                   nrow = dim(Dr)[1],ncol=5,byrow=FALSE)
  } else {                              #actually exposure
    ariz <- matrix(rep(ARI,5),ncol=5,nrow=dim(Dr)[1],byrow=FALSE)
  }
  ## print(ariz)
  atrisk[,2:6] <- atrisk[,2:6] * ariz * 1e3
  return(atrisk)
}

## indtest <- getLnE(Dr[Dr$country=='India',], 0.5*1e-2, method='prev')
## sum(indtest[,-1]*1e-6)                  #13M

## some natural history parameters
## age groups 0-1,1-2,2-5,5-10,10-15
## example - lots of bits over-written
parm <- list(prPTB=c(0.35,
                     0.15,
                     0.05,
                     0.02,
                     0.15), #risk of progression to PTB
             prDTB=c(0.15,
                     0.035,
                     0.005,
                     0.005,
                     0.005), #risk of progression TBM/miliary
             vBCGd=rep(0.2,5),  #HR for diseminated w BCG
             vBCGp=rep(0.3,5),   #HR for pulmonary w BCG
             hivIRR=20,
             cohabrate=0.3,     #beta hh
             latprop=.41       #latitude variability
             )

## ## ## check 1000s or not
## test <- getIncidence(Dr,parm, .01,Dr[,hivu15])



## varying gradient from household sizes
NC <- dim(Dr)[1]



## --------------
## version using LHS etc
## --------------

## ..............

LHN <- 22+NC+NH
LHsample <- randomLHS(100,LHN)                #latin hypercube sample

## no pars = 5.progt + 5.progD + 5.progProtD + 5.progProtF + beta  + hiv = 22
## TB prevalence = NC
## NH?
## convert quantiles into samples from the distributions
makeParameters2 <- function(LH,mp){
    nr <- dim(LH)[1]
    ## progression
    suppressWarnings({                  #numerical from qbeta
        p1z <- qbeta(LH[,1:5],
                     shape1=matrix(mp$drs1,ncol=5,nrow=nr,byrow=TRUE),
                     shape2=matrix(mp$drs2,ncol=5,nrow=nr,byrow=TRUE))
        p2z <- qbeta(LH[,1:5 + 5],
                     shape1=matrix(mp$drs1d,ncol=5,nrow=nr,byrow=TRUE),
                     shape2=matrix(mp$drs2d,ncol=5,nrow=nr,byrow=TRUE))
        ## BCG
        vdz <- qbeta(LH[,1:5 + 10],
                     shape1=matrix(mp$vds1,ncol=5,nrow=nr,byrow=TRUE),
                     shape2=matrix(mp$vds2,ncol=5,nrow=nr,byrow=TRUE))
        vpz <- qbeta(LH[,1:5 + 15],
                     shape1=matrix(mp$vs1,ncol=5,nrow=nr,byrow=TRUE),
                     shape2=matrix(mp$vs2,ncol=5,nrow=nr,byrow=TRUE))
    })
    ## vpz <- 1-(1-vpz) * vdz
    vpz <- 1-vpz * (1-vdz) #HRp = 1-frac*(1-HRd)
    ## HIV IRR
    irrz <- qlnorm(LH[,21],meanlog = mp$hmlg,sd=mp$hsd) #HIV
    ## betas
    ##different by country/random effects? could be argued for lots
    bz <- qlnorm(LH[,22],meanlog=mp$bmlg,sd=mp$bsd)
    ## HIV and TB prevalence data
    tb <- qgamma(LH[,23:(22+NC)],
                 shape = matrix(Dr$shpz,ncol=NC,nrow=nr,byrow=TRUE),
                 rate = matrix(Dr$rtz,ncol=NC,nrow=nr,byrow=TRUE))
    hiv <- qgamma(LH[,(23+NC):(22+NC+NH)],
                  shape = matrix(Dr[hivu15>0,shpzh],
                                 ncol=NH,nrow=nr,byrow=TRUE),
                  rate = matrix(Dr[hivu15>0,rtzh],
                                ncol=NH,nrow=nr,byrow=TRUE))
    ## safety cap at 1
    hiv[is.na(hiv)] <- 0; hiv[hiv>1] <- .99;tb[tb>1e5] <- 1e5
    return(list(p1=p1z,p2=p2z,vd=vdz,vp=vpz,irr=irrz,b=bz,
                tb=tb,h=hiv))
}

## test
## LHparms <- makeParameters2(LHsample,metaparm) #turn the raw LHS into parameters

## this loops through the list of parameter samples and returns PTB/DTB
## this data.table version is slightly slower but scales better than the dataframe one
loopAns2 <- function(LP,mp,method=c('comm','lat')){
    n <- length(LP$irr)
    k <- dim(Dr)[1]
    nr <- n*k
    rpl <- rep(1:n,each=nrow(Dr))
    ans <- data.frame(country=rep(Dr$country,n),
                      a=0.0,b=0.0,s=0.0,d=0.0,f=0.0,
                      replicate=rpl, iso3=rep(Dr$iso3,n),
                      ARI=method[1],LAT=method[2])
    names(ans)[2:6] <- c("<1","1","2-4","5-9","10-14")
    ans <- as.data.table(ans)
    ans2 <- copy(ans)
    ans4 <- data.table(country=rep(Dr$country,n),ariz=0.0,
                       replicate=rpl,ARI=method[1],
                       LAT=method[2],
                       iso3=rep(Dr$iso3,n)) #recording ARIs
    for(i in 1:n){
        if(!i%%100)print(paste(method[1],method[2],':',i))
        parms <- list(prPTB=LP$p1[i,]*(1-LP$p2[i,]),
                      prDTB=LP$p1[i,]*LP$p2[i,],vBCGd=LP$vd[i,],
                      vBCGp=LP$vp[i,],hivIRR=LP$irr[i],latprop = 0.41)
        hiv <- Dr$hivu15
        hiv[hiv>0] <- LP$h[i,] # random sample for isos with HIV
        ## make ARI
        ari <- LP$tb[i,] / 1e5
        vg <- 1
        if( method[1]=='comm'){
            ari <- ari * LP$b[i] * smpm
        } else {
            stop('deprecated method!')
        }
        tmp <- getIncidence(Dr,parms,ari,hiv,
                            method=method[2]) #change here
        ans[1:k +(i-1)*k,(2:6):=tmp$PTB[,2:6]]
        ans2[1:k +(i-1)*k,(2:6):=tmp$DTB[,2:6]]
        ans4[1:k +(i-1)*k,ariz:=ari]
      }
    return(list(PTB=ans,DTB=ans2,ARI=ans4))
}


## function for aggregating and returned data
make2cat <- function(L){#aggregate results
    ans <- as.data.frame(L$PTB)
    tmp <- as.data.frame(L$DTB)
    ans[,2:6] <- ans[,2:6] + tmp[,2:6]  #PTB+DTB
    ans[,2:6] <- ans[,2:6] + tmp[,2:6]  #PTB+DTB
    ## melted version
    n <- length(ans$country)
    ans2 <-  data.frame(country=rep(ans$country,2),
                        iso3=rep(ans$iso3,2),value=0,
                        age = c(rep('0-4',n),rep('5-14',n)),
                        replicate=ans$replicate,ARI=ans$ARI,LAT=ans$LAT)
    ans2$value <- c(rowSums(ans[,2:3]),rowSums(ans[,4:6]))
    ## returning
    return(ans2)
}


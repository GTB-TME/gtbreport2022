## ======== child mortality estimates prior  ===========
## run in this directory
rm(list=ls())
## libraries
library(lhs)
library(data.table)
library(scales)
library(ggplot2)
library(here)

## ------- load data from preparation.R -----------
load(here('disaggregation/R/02children/mortinput/MBB.Rdata'))
load(here('data/tb.rda'))
load(here('disaggregation/output/incsplits/data/incsplit.Rdata'))
WHOkey <- unique(tb[,.(iso3,g_whoregion=g.whoregion)])

## -------functions--------
oddit <- function(x) x/(1-x)
ioddit <- function(x) x/(1+x)
logit <- function(x) log(oddit(x))
ilogit <- function(x) ioddit(exp(x))



##function for extracting lognormal parz from mid point and variance
getLNparms <- function(mid,var){
    mu <- log(mid)                      #mid as median
    x <- (1+sqrt(1+4*var/mid^2))/2
    list(mu=mu,sig=sqrt(log(x)))
}


## take x,y U[0,1]
getHAORs <- function(x,y){
    x <- qnorm(x)
    y <- qnorm(y)                       #normals
    z <- matrix(c(x,y),ncol=2,nrow=length(x))
    z <- z %*% hivartOR$U      #make these correlated normals
    z <- z + matrix(hivartOR$mn,ncol=2,nrow=nrow(z),byrow=TRUE) #proper return
    colnames(z) <- c('HIV','ART')
    exp(z)
}



IQR <- function(x) round(quantile(x,c(.25,.75)),digits=2)
pIQR <- function(x) paste0('[',paste(IQR(x),collapse=', '),']')
rmn <- function(x) round(mean(x),digits = 2)

## ------- parameter work and verification ------------

## --- distribution for CFR with tx
## splitting by age
## Y: 1·9% (0·5%, 7.1%)
## O: 0·8% (0·3%, 2.1%)
ontxY <- getLNparms(1.9*1e-2,(1e-2*(7.1-.5)/3.92)^2)
ontxO <- getLNparms(.8*1e-2,(1e-2*(2.1-.3)/3.92)^2)

## HIV/ART ORs
## revised version from frequentist approach
## mean =  2.6375681 -0.5683867
## V = [[0.2325509 -0.2325509],[-0.2325509  0.6367345]]

hivartOR <- list(mn = c(2.6375681, -0.5683867),
                 sg = matrix(c(0.2325509, -0.2325509,-0.2325509,  0.6367345),2,2))
hivartOR$U <- chol(hivartOR$sg)

## ## ## inspect if sensible
## png('test/test_CFR_ontx_Y.png')
## curve(dlnorm(x,ontxY$mu,ontxY$sig),from=0.,to=.1,n=1e3,col=2,xlab='CFR on tx (0-4)',ylab='') 
## abline(v=1.9*1e-2);abline(v=.5*1e-2,lty=2);abline(v=7.1*1e-2,lty=2);
## dev.off()

## png('test/test_CFR_ontx_O.png')
## curve(dlnorm(x,ontxO$mu,ontxO$sig),from=0.,to=.05,n=1e3,col=2,xlab='CFR on tx (5-14)',ylab='')
## abline(v=.8*1e-2);abline(v=.3*1e-2,lty=2);abline(v=2.1*1e-2,lty=2);
## dev.off()

## test
Z <- getHAORs(runif(1e4),runif(1e4))
Yz0 <- qlnorm(runif(1e4),ontxY$mu,ontxY$sig)
Oz0 <- qlnorm(runif(1e4),ontxO$mu,ontxO$sig)
head(Yz0)
head(ioddit(Z[,'HIV']*oddit(Yz0)))

## HIV/ART

Yz <- ioddit(Z[,'HIV']*oddit(Yz0))

## PP <- qplot(x=Yz) +
##     ggtitle(paste0('on tx, 0-4,HIV+/ART-\nmean=',rmn(Yz),', IQR=',pIQR(Yz))) +
##     geom_vline(xintercept=14.3*1e-2,col=2) +
##     geom_vline(xintercept=7.4*1e-2,col=2,lty=2) +
##     geom_vline(xintercept=24.1*1e-2,col=2,lty=2) +
##     annotate('text',x=0.25,y=1500,label='HIV+CFR in pre-ART US data\n(all ages)',col=2) +
##     xlab('CFR')

## ggsave(PP,file='test/test_CFR_ontx_YH.png')

## Oz <- ioddit(Z[,'HIV']*oddit(Oz0))

## PP <- qplot(x=Oz) +
##     ggtitle(paste0('on tx, 5-14,HIV+/ART-\nmean=',rmn(Oz),', IQR=',pIQR(Oz))) +
##     geom_vline(xintercept=14.3*1e-2,col=2) +
##     geom_vline(xintercept=7.4*1e-2,col=2,lty=2) +
##     geom_vline(xintercept=24.1*1e-2,col=2,lty=2) +
##     annotate('text',x=0.25,y=1500,label='HIV+CFR in pre-ART US data\n(all ages)',col=2) + xlab('CFR')

## ggsave(PP,file='test/test_CFR_ontx_OH.png')

## Yz <- ioddit(Z[,'ART']*Z[,'HIV']*oddit(Yz0))

## PP <- qplot(x=Yz) +
##     ggtitle(paste0('on tx, 0-4,HIV+/ART+\nmean=',rmn(Yz),', IQR=',pIQR(Yz))) +
##     geom_vline(xintercept=3.4*1e-2,col=2) +
##     geom_vline(xintercept=.7*1e-2,col=2,lty=2) +
##     geom_vline(xintercept=9.6*1e-2,col=2,lty=2) +
##     annotate('text',x=0.25,y=1500,label='HIV+CFR in post-ART US data\n(all ages)',col=2) + xlab('CFR')


## ggsave(PP,file='test/test_CFR_ontx_YHA.png')

## Oz <- ioddit(Z[,'ART']*Z[,'HIV']*oddit(Oz0))

## PP <- qplot(x=Oz) +
##     ggtitle(paste0('on tx, 5-14,HIV+/ART+\nmean=',rmn(Oz),', IQR=',pIQR(Oz))) +
##     geom_vline(xintercept=3.4*1e-2,col=2) +
##     geom_vline(xintercept=.7*1e-2,col=2,lty=2) +
##     geom_vline(xintercept=9.6*1e-2,col=2,lty=2) +
##     annotate('text',x=0.25,y=1500,label='HIV+CFR in pre-ART US data\n(all ages)',col=2) + xlab('CFR')

## ggsave(PP,file='test/test_CFR_ontx_OHA.png')


## --- distribution for CFR w/o tx
## 0-4: 43.6%; 95% CI: 36.8%, 50.6%
## 5-14: 14.9%; 95% CI: 11.5%, 19.1%). 

notxY <- getLNparms(43.6*1e-2,(1e-2*(50.6-36.8)/3.92)^2)
notxO <- getLNparms(14.9*1e-2,(1e-2*(19.1-11.5)/3.92)^2)

## png('test/test_CFR_notxY.png')
## curve(dlnorm(x,notxY$mu,notxY$sig),from=.3,to=.6,n=1e3,col=2,xlab='CFR no tx (0-4)',ylab='') 
## abline(v=43.6*1e-2);abline(v=50.6*1e-2,lty=2);abline(v=36.8*1e-2,lty=2);
## dev.off()

## png('test/test_CFR_notxO.png')
## curve(dlnorm(x,notxO$mu,notxO$sig),from=0.,to=.3,n=1e3,col=2,xlab='CFR no tx (5-15)',ylab='') 
## abline(v=14.9*1e-2);abline(v=11.5*1e-2,lty=2);abline(v=19.1*1e-2,lty=2);
## dev.off()

## NB--- this is BETA not lnormal
## from roulette elicitation -- stored in data/ez.txt
## see data/elicit.R for methods todo

## ez <- scan('data/ez.txt')
ez <- c(77.13050, 11.10817, 15.18683, 12.87500, 19.59083,  6.89700, 10.43383, 11.08417)

## HO
notxHO <- list(a=ez[5],b=ez[6])
quantile(rbeta(n=1e4,notxHO$a,notxHO$b),probs=c(.025,.5,.975))

## png('test/test_CFR_notx_HO.png')
## curve(dbeta(x,notxHO$a,notxHO$b),from=0,to=1,n=1e3,col=2,
##       xlab='CFR no tx (HIV+/ART-/5-14)') 
## dev.off()

## HY
notxHY <- list(a=ez[1],b=ez[2])
quantile(rbeta(n=1e4,notxHY$a,notxHY$b),probs=c(.025,.5,.975))

## png('test/test_CFR_notx_HY.png')
## curve(dbeta(x,notxHY$a,notxHY$b),from=0,to=1,n=1e3,col=2,
##       xlab='CFR no tx (HIV+/ART-/0-4)') 
## dev.off()


## HAY
notxHAY <- list(a=ez[3],b=ez[4])
quantile(rbeta(n=1e4,notxHAY$a,notxHAY$b),probs=c(.025,.5,.975))

## png('test/test_CFR_notx_HAY.png')
## curve(dbeta(x,notxHAY$a,notxHAY$b),from=0,to=1,n=1e3,col=2,
##       xlab='CFR no tx (HIV+/ART+/0-4)') 
## dev.off()

## HAO
## convert gamma to beta
notxHAO <- list(a=ez[7],b=ez[8])
quantile(rbeta(n=1e4,notxHAO$a,notxHAO$b),probs=c(.025,.5,.975))

## png('test/test_CFR_notx_HAO.png')
## curve(dbeta(x,notxHAO$a,notxHAO$b),from=0,to=1,n=1e3,col=2,
##       xlab='CFR no tx (HIV+/ART+/5-14)') 
## dev.off()


## --- HIV/ART distributions
## ART HR: 0.30 [0.21 - 0.39]
artp <- getLNparms(0.3,((0.39-0.21)/3.92)^2)

## from UNAIDS meta-analysis. 
hivp <- list(mu=4.364, sig=1.14)        #REMA

## png('test/test_IRR_HIV.png')
## curve(dlnorm(x,hivp$mu,hivp$sig),from=1,to=150,n=1e3,col=2,xlab='IRR for UN HIV')
## dev.off()

## png('test/test_HR_ART.png')
## curve(dlnorm(x,artp$mu,artp$sig),from=.1,to=.6,n=1e3,col=2,xlab='HR for ART') 
## abline(v=0.3);abline(v=0.21,lty=2);abline(v=0.39,lty=2);
## dev.off()


## to cope with potential zeros
QG <- function(r,shape,scale){
  z <- !is.finite(scale) | !is.finite(shape) | !shape>0 | !scale>0
  ans <- r
  ans[z] <- 0
  ans[!z] <-  qgamma( r[!z],shape=shape[!z],scale=scale[!z])
  ans
}

## quantile functions for distributions
BURD <- function(x){
    ans <- QG(x,shape=MBB$k,scale=MBB$S) #NB for length(x)=nrow(MBB) !
    ans[BB$k==0] <- 0                  #safety for QG
    ans
}

## not on treatment
CFRnotxHY <- function(x) qbeta(x,notxHY$a,notxHY$b) #quantile function
CFRnotxHO <- function(x) qbeta(x,notxHO$a,notxHO$b) #quantile function
CFRnotxHAY <- function(x) qbeta(x,notxHAY$a,notxHAY$b) #quantile function
CFRnotxHAO <- function(x) qbeta(x,notxHAO$a,notxHAO$b) #quantile function
CFRnotxY <- function(x) qlnorm(x,notxY$mu,notxY$sig) #quantile function
CFRnotxO <- function(x) qlnorm(x,notxO$mu,notxO$sig) #quantile function

## on treatment
CFRontxY <- function(x) qlnorm(x,ontxY$mu,ontxY$sig) #quantile function
CFRontxO <- function(x) qlnorm(x,ontxO$mu,ontxO$sig) #quantile function

## see getHAORs above to convert these two into relevant for HIV/ART states on tx

IRRhiv <- function(x) qlnorm(x,hivp$mu,hivp$sig) #quantile function
HRart <- function(x) qlnorm(x,artp$mu,artp$sig) #quantile function

BURDH <- function(x) {
    ans <- QG(x,shape=MBB$hk,scale=MBB$hS) #NB for length(x)=nrow(BB) !
    ans[is.na(ans)] <- 0                     #safety
    ans[ans>1] <- .999                       #safety
    ans
}

BURDA <- function(x) {
    ## warns when ak==Inf
    ans <- QG(x,shape=MBB$ak,scale=MBB$aS) #NB for length(x)=nrow(BB) !
    ans[is.na(ans)] <- 0
    ans[ans>1] <- .999
    ans
}

test <- BURDA(runif(nrow(MBB)))

## vectorized beta sampling
vrbeta <- function(AB){
    apply(X=AB,MARGIN=1,FUN=function(x)rbeta(1,x[1],x[2]))
}

## vrbeta(cbind(rep(2,10),rep(5,10)))

## ## check
## hist(CFRnotxHAO(runif(1e4)))
## ## hist(CFRtx(runif(1e4)))

## ## hist(CFRtxH(runif(1e4)))
## ## hist(CFRnotxH(runif(1e4)))
## hist(CFRnotxY(runif(1e4)))
## hist(CFRnotxO(runif(1e4)))

## hist(IRRhiv(runif(1e4)))
## hist(HRart(runif(1e4)))

## ## some NAs in notification data
## cat(as.character(BB$country[is.na(BB$notifs)]),file='test/clean_NAnotifs.txt')
## sum(is.na(BB$notifs))                   #0
## cat(as.character(BB$country[is.na(BB$notifsY)]),file='test/clean_NAnotifsY.txt')
## sum(is.na(BB$notifsY))                   #0
## cat(as.character(BB$country[is.na(BB$notifsO)]),file='test/clean_NAnotifsO.txt')
## sum(is.na(BB$notifsO))                   #0

## no change...
MBB$notifs[is.na(MBB$notifs)] <- 0
MBB$notifsY[is.na(MBB$notifsY)] <- 0
MBB$notifsO[is.na(MBB$notifsO)] <- 0



## -------------- ,,, --------

## gwho <- unique(as.character(RESGs$g_whoregion))
txstat <- c('On tx (HIV-)','Off tx (HIV-)',
            'On tx (HIV+)','Off tx (HIV+)')


## ## ---------------- same but with age now --------------
## change sensitivity analyses here!!
useHIVneg <- FALSE                      #use the elicited data or be conservative
sensNote <- FALSE                       #sensitivity for private sectors?
## stuff for note=tx sensitivity analysis - from top 10
senscns <- c('IND','CHN','NGA','IDN','COD','BGD','PAK','MOZ','TZA','VNM')
senspriv <- c(46,
              23,
              10,
              56,
              0,
              7,
              32,
              0,
              0,
              9)*1e-2   #proportions in private sector
senswho <- rep(NA,length(senscns))
for(i in 1:length(senscns))
    senswho[i] <- which(as.character(MBB$iso3)==senscns[i])
F <- rep(1,nrow(MBB))                    #mult factor for notifications
if(sensNote)
    F[senswho] <- 1/(1-senspriv)            #mult inflation factor

## n replicates
nrep <- 1e3*10                                   #NOTE
RESA <- matrix(NA,nrow = nrow(MBB)*nrep,ncol=12) #for results
RESA[,1] <- 1:nrow(MBB)

## need runif for: burden, 2xCFR = 3
L <- randomLHS(n=nrow(RESA),k=15)         #latin hypercube sample
buryo <- incsplit[MBB[,as.character(iso3)]][age_group=='0_4',.(inc=sum(inc)),by=iso3]
burol <- incsplit[MBB[,as.character(iso3)]][age_group=='5_14',.(inc=sum(inc)),by=iso3]

## L is a latin hypercube sample on U[0,1]^15
## MBB contains the country data - calx done world-block-wise 
for(i in 1:nrep){
    if(!i%%1e3) print(i)
    now <- (i-1)*nrow(MBB) + 1:nrow(MBB)  # replicate block for all countries
    ## get parameters
    ## burden <- BURD(L[now,1])            # burden
    cfrz1y <- CFRontxY(L[now,2])            # CFR on tx 0-5
    cfrz1o <- CFRontxO(L[now,3])            # CFR on tx 5-15
    cfrz0y <- CFRnotxY(L[now,4])        # CFR off tx 0-5
    cfrz0o <- CFRnotxO(L[now,5])       # CFR off tx 5-15
    proph <- BURDH(L[now,6])            # burden HIV
    irrz <- IRRhiv(L[now,7])            # IRR HIV
    propa <- BURDA(L[now,8])            # prev ART
    hrz <- HRart(L[now,9])              # HR ART
    ## new HIV on tx via ORs
    Z <- getHAORs(L[now,10],L[now,11])    #ORs for CFRs HIV/ART (2 cols w/ correlation)
    cfrzh1Y <- ioddit(Z[,'HIV']*oddit(cfrz1y))        # CFR HIV on tx, ART-, Y
    cfrzh1O <- ioddit(Z[,'HIV']*oddit(cfrz1o))       # CFR HIV on tx, ART-, O
    cfrzh1aY <- ioddit(Z[,'HIV']*Z[,'ART']*oddit(cfrz1y))        # CFR HIV on tx, ART+, Y
    cfrzh1aO <- ioddit(Z[,'HIV']*Z[,'ART']*oddit(cfrz1o))       # CFR HIV on tx, ART+, O
    ## new HIV off tx elicited
    if(!useHIVneg){
        cfrzh0Y <- CFRnotxHY(L[now,12])        # CFR HIV no tx, ART-, Y
        cfrzh0O <- CFRnotxHO(L[now,13])        # CFR HIV no tx, ART-, O
        cfrzh0aY <- CFRnotxHAY(L[now,14])        # CFR HIV no tx, ART+, Y
        cfrzh0aO <- CFRnotxHAO(L[now,15])        # CFR HIV no tx, ART+, O
    } else{
        cfrzh0Y <- CFRnotxY(L[now,12])        # CFR HIV no tx, ART-, Y
        cfrzh0O <- CFRnotxO(L[now,13])        # CFR HIV no tx, ART-, O
        cfrzh0aY <- CFRnotxY(L[now,14])        # CFR HIV no tx, ART+, Y
        cfrzh0aO <- CFRnotxO(L[now,15])        # CFR HIV no tx, ART+, O
    }
    ## calculations
    ## disaggregation: inc = inc0 * (nhiv + IRR*hiv*(nart + HR*art)) - 
    tots <- (1-proph) + irrz*proph*(1-propa + hrz*propa)
    tbhiv <- 1-(1-proph)/tots                    #proportion of TB in HIV
    tbart <- hrz*propa/(1-propa+hrz*propa)       #proportion of HIV-TB in ART
    ## young <- vrbeta(MBB[,c('aa','ab')])  #proportion young
    ## deaths
    ## --young
    DeathsOnTxY <- F*MBB$notifsY * cfrz1y * (1-tbhiv) #F is notification inflation =1 unless priv
    DeathsOffTxY <- pmax(buryo$inc-F*MBB$notifsY,0) * cfrz0y * (1-tbhiv)
    DeathsOnTxHY <- F*MBB$notifsY * (cfrzh1aY * tbart + cfrzh1Y * (1-tbart)) * tbhiv 
    DeathsOffTxHY <- pmax(buryo$inc-F*MBB$notifsY,0) * tbhiv * (cfrzh0aY * tbart + cfrzh0Y * (1-tbart)) #averages CFR
    ## --old
    DeathsOnTxO <- F*MBB$notifsO * cfrz1o * (1-tbhiv)
    DeathsOffTxO <- pmax(burol$inc-F*MBB$notifsO,0) * cfrz0o * (1-tbhiv)
    DeathsOnTxHO <- F*MBB$notifsO * (cfrzh1aO * tbart + cfrzh1O * (1-tbart)) * tbhiv
    DeathsOffTxHO <- pmax(burol$inc-F*MBB$notifsO,0) * tbhiv * (cfrzh0aO * tbart + cfrzh0O * (1-tbart)) #averages CFR
    ## record results
    RESA[now,2] <- DeathsOnTxY + DeathsOnTxO
    RESA[now,3] <- DeathsOffTxY + DeathsOffTxO
    RESA[now,4] <- DeathsOnTxHY + DeathsOnTxHO
    RESA[now,5] <- DeathsOffTxHY + DeathsOffTxHO
    RESA[now,6] <- DeathsOnTxY 
    RESA[now,7] <- DeathsOffTxY
    RESA[now,8] <- DeathsOnTxHY
    RESA[now,9] <- DeathsOffTxHY
    RESA[now,10] <- proph
    RESA[now,11] <- tbhiv
    RESA[now,12] <- (DeathsOnTxHY + DeathsOnTxHO + DeathsOffTxHY + DeathsOffTxHO) / (DeathsOnTxY + DeathsOnTxO + DeathsOffTxY + DeathsOffTxO + DeathsOnTxHY + DeathsOnTxHO + DeathsOffTxHY + DeathsOffTxHO + 1e-6)
}

## neaten
colnames(RESA) <- c('id','DeathsOnTx','DeathsOffTx',
                    'DeathsOnTxH','DeathsOffTxH',
                    'DeathsOnTxY','DeathsOffTxY',
                    'DeathsOnTxHY','DeathsOffTxHY',
                    'HIV','HIVinTB','HIVinDeaths')
RESA <- as.data.frame(RESA)
RESA$iso3 <- MBB$iso3
RESA$country <- MBB$country
RESA$g_whoregion<- MBB$g.whoregion
RESA$replicate <- rep(1:nrep,each=nrow(MBB))


## convert for post-processing
RESA <- as.data.table(RESA)
RESA

summary(RESA)                           #fix

##  deaths by country
RESAc <- RESA[,list(deaths=mean(DeathsOnTx+DeathsOffTx+DeathsOnTxH+DeathsOffTxH)),by=list(iso3,g_whoregion)]



RESATOT <- RESA[,list(total = sum(DeathsOnTx,na.rm=TRUE)+
                          sum(DeathsOffTx,na.rm=TRUE)+
                          sum(DeathsOnTxH,na.rm=TRUE)+
                          sum(DeathsOffTxH,na.rm=TRUE)),
              by=list(replicate)]       #global totals

RESATOTS <- RESATOT[,list(mid=median(total),
                          lo=quantile(total,probs =.025),
                          hi=quantile(total,probs=.975))] #summarize totals

## inspect
if(!sensNote & !useHIVneg){
    png(here('disaggregation/R/02children/plots/global_hist.png'))
    hist(RESATOT$total,xlim=c(0e5,3.5e5),breaks=20)
    abline(v=RESATOTS$mid,col=2,lwd=3);abline(v=RESATOTS$lo,lty=2);abline(v=RESATOTS$hi,lty=2)
    dev.off()
}

## regional totals
RESAG <- RESA[,list(OnTx=sum(DeathsOnTx),OffTx=sum(DeathsOffTx),
                  OnTxH=sum(DeathsOnTxH),OffTxH=sum(DeathsOffTxH)),
            by=list(replicate,g_whoregion)]

RESAGs <- RESAG[,list(OnTx=quantile(OnTx,.5),
                      OffTx=quantile(OffTx,.5),
                      OnTxH=quantile(OnTxH,.5),
                      OffTxH=quantile(OffTxH,.5),
                      OnTx0=quantile(OnTx,.025),
                      OffTx0=quantile(OffTx,.025),
                      OnTxH0=quantile(OnTxH,.025),
                      OffTxH0=quantile(OffTxH,.025),
                      OnTx1=quantile(OnTx,.975),
                      OffTx1=quantile(OffTx,.975),
                      OnTxH1=quantile(OnTxH,.975),
                      OffTxH1=quantile(OffTxH,.975)
                    ),
            by=list(g_whoregion)]

## global tables
RESAT <- RESA[,list(OnTx=sum(DeathsOnTx),OffTx=sum(DeathsOffTx),
                  OnTxH=sum(DeathsOnTxH),OffTxH=sum(DeathsOffTxH)),
            by=list(replicate)]

RESATs <- RESAT[,list(OnTx=quantile(OnTx,.5),
                      OffTx=quantile(OffTx,.5),
                      OnTxH=quantile(OnTxH,.5),
                      OffTxH=quantile(OffTxH,.5),
                      OnTx0=quantile(OnTx,.025),
                      OffTx0=quantile(OffTx,.025),
                      OnTxH0=quantile(OnTxH,.025),
                      OffTxH0=quantile(OffTxH,.025),
                      OnTx1=quantile(OnTx,.975),
                      OffTx1=quantile(OffTx,.975),
                      OnTxH1=quantile(OnTxH,.975),
                      OffTxH1=quantile(OffTxH,.975))]

RESAGnc <- RESA[,list(total=sum(DeathsOnTx+DeathsOffTx+
                                DeathsOnTxH+DeathsOffTxH)),
            by=list(replicate,g_whoregion)]

RESAGncs <- RESAGnc[,list(total=median(total),
                          total0=quantile(total,.025),
                        total1=quantile(total,.975)),
                  by=list(g_whoregion)]

RESAGnca <- RESA[,list(total=sum(DeathsOnTx+DeathsOffTx+
                                 DeathsOnTxH+DeathsOffTxH)),
            by=list(replicate)]


gwho <- unique(as.character(RESAGs$g_whoregion))
gwho <- sort(gwho)
txstat <- c('On tx (HIV-)','Off tx (HIV-)','On tx (HIV+)',
            'Off tx (HIV+)')


RESAGsm <- melt(RESAGs[,1:5,with=FALSE],id='g_whoregion')

## reshaping for tables
RESAGsl <- melt(RESAGs[,c(1,6:9),with=FALSE],id='g_whoregion')
RESAGsh <- melt(RESAGs[,c(1,10:13),with=FALSE],id='g_whoregion')

LO <- dcast(RESAGsl,g_whoregion ~ variable) 
MID <- dcast(RESAGsm,g_whoregion ~ variable)
HI <- dcast(RESAGsh,g_whoregion ~ variable)
gwho <- (as.character(MID$g_whoregion))
LO <- LO[,-1]
MID <- MID[,-1]
HI <- HI[,-1]
MID <- cbind(MID,total=RESAGncs[order(g_whoregion),total])
LO <- cbind(LO,total=RESAGncs[order(g_whoregion),total0])
HI <- cbind(HI,total=RESAGncs[order(g_whoregion),total1])
gm <- c(unlist(RESATs[,1:4,with=FALSE]),total=median(RESAGnca$total))
gl <- c(unlist(RESATs[,5:8,with=FALSE]),
        total=quantile(RESAGnca$total,.025))
gh <- c(unlist(RESATs[,9:12,with=FALSE]),
        total=quantile(RESAGnca$total,.975))

MID <- as.data.table(rbind(as.matrix(MID),gm))
LO <- as.data.table(rbind(as.matrix(LO),gl))
HI <- as.data.table(rbind(as.matrix(HI),gh))

## youngs ----------
RESATy <- RESA[,list(OnTx=sum(DeathsOnTxY),OffTx=sum(DeathsOffTxY),
                  OnTxH=sum(DeathsOnTxHY),OffTxH=sum(DeathsOffTxHY)),
            by=list(replicate)]


RESAGy <- RESA[,list(OnTx=sum(DeathsOnTxY),OffTx=sum(DeathsOffTxY),
                  OnTxH=sum(DeathsOnTxHY),OffTxH=sum(DeathsOffTxHY)),
            by=list(replicate,g_whoregion)]

RESAGsy <- RESAGy[,list(OnTx=quantile(OnTx,.5),
                        OffTx=quantile(OffTx,.5),
                        OnTxH=quantile(OnTxH,.5),
                        OffTxH=quantile(OffTxH,.5),
                        OnTx0=quantile(OnTx,.025),
                        OffTx0=quantile(OffTx,.025),
                        OnTxH0=quantile(OnTxH,.025),
                        OffTxH0=quantile(OffTxH,.025),
                        OnTx1=quantile(OnTx,.975),
                        OffTx1=quantile(OffTx,.975),
                        OnTxH1=quantile(OnTxH,.975),
                        OffTxH1=quantile(OffTxH,.975)
                    ),
            by=list(g_whoregion)]


RESATsy <- RESATy[,list(OnTx=quantile(OnTx,.5),OffTx=quantile(OffTx,.5),
                   OnTxH=quantile(OnTxH,.5),OffTxH=quantile(OffTxH,.5),
                   OnTx0=quantile(OnTx,.025),OffTx0=quantile(OffTx,.025),
                   OnTxH0=quantile(OnTxH,.025),
                   OffTxH0=quantile(OffTxH,.025),
                   OnTx1=quantile(OnTx,.975),OffTx1=quantile(OffTx,.975),
                   OnTxH1=quantile(OnTxH,.975),
                   OffTxH1=quantile(OffTxH,.975))]

RESAGncy <- RESA[,list(total=sum(DeathsOnTxY+DeathsOffTxY+
                                 DeathsOnTxHY+DeathsOffTxHY)),
            by=list(replicate,g_whoregion)]

RESAGncsy <- RESAGncy[,list(total=median(total),
                            total0=quantile(total,.025),
                        total1=quantile(total,.975)),
                  by=list(g_whoregion)]

RESAGncay <- RESA[,list(total=sum(DeathsOnTxY+DeathsOffTxY+
                                  DeathsOnTxHY+DeathsOffTxHY)),
            by=list(replicate)]

gwho <- unique(as.character(RESAGsy$g_whoregion))
txstat <- c('On tx (HIV-)','Off tx (HIV-)',
            'On tx (HIV+)','Off tx (HIV+)')


RESAGsmy <- melt(RESAGsy[,1:5,with=FALSE],id='g_whoregion')


## reshaping for tables
RESAGsly <- melt(RESAGsy[,c(1,6:9),with=FALSE],id='g_whoregion')
RESAGshy <- melt(RESAGsy[,c(1,10:13),with=FALSE],id='g_whoregion')

LO <- dcast(RESAGsly,g_whoregion ~ variable) 
MID <- dcast(RESAGsmy,g_whoregion ~ variable)
HI <- dcast(RESAGshy,g_whoregion ~ variable)
gwho <- as.character(MID$g_whoregion)
LO <- LO[,-1]
MID <- MID[,-1]
HI <- HI[,-1]

MID <- cbind(MID,RESAGncsy[order(g_whoregion),total])
LO <- cbind(LO,RESAGncsy[order(g_whoregion),total0])
HI <- cbind(HI,RESAGncsy[order(g_whoregion),total1])
gm <- c(unlist(RESATsy[,1:4,with=FALSE]),median(RESAGncay$total))
gl <- c(unlist(RESATsy[,5:8,with=FALSE]),quantile(RESAGncay$total,.025))
gh <- c(unlist(RESATsy[,9:12,with=FALSE]),quantile(RESAGncay$total,.975))

MID <- as.data.table(rbind(as.matrix(MID),gm))
LO <- as.data.table(rbind(as.matrix(LO),gl))
HI <- as.data.table(rbind(as.matrix(HI),gh))

## compare young/old
RESATb <- RESA[,list(OnTx=sum(DeathsOnTx),OffTx=sum(DeathsOffTx),
                     OnTxH=sum(DeathsOnTxH),OffTxH=sum(DeathsOffTxH),
                     OnTxY=sum(DeathsOnTxY),OffTxY=sum(DeathsOffTxY),
                  OnTxHY=sum(DeathsOnTxHY),OffTxHY=sum(DeathsOffTxHY)),
            by=list(replicate)]

## mosaic
RESATB <- RESATb[,list(OnTxY = median(OnTxY),
                       OnTxO = median(OnTx-OnTxY),
                       OffTxY = median(OffTxY),
                       OffTxO = median(OffTx-OffTxY))]

RESATB <- data.frame(Treatment=rep(c('On','Off'),each=2),
                     Age=c('0-4','5-14'),Number=unlist(RESATB[1,]))

## for mosaic
RESATBM <- RESATb[,list(OnTxY = median(OnTxY),
                        OnTxO = median(OnTx-OnTxY),
                        OffTxY = median(OffTxY),
                        OffTxO = median(OffTx-OffTxY),
                        OnTxHY = median(OnTxHY),
                        OnTxHO = median(OnTxH-OnTxHY),
                        OffTxHY = median(OffTxHY),
                        OffTxHO = median(OffTxH-OffTxHY)
                       )
                  ]


## ------- Mosaic data ---------

mosdat <- data.frame(cat = names(RESATBM),
                     value = unname(unlist(RESATBM)),
                     age = rep(c(' age <5',' age 5 to <15'),4),
                     tbtx = rep(c('Treated for TB','Treated for TB',
                                  'Untreated for TB',
                                  'Untreated for TB'),2),
                     hiv = rep(c('HIV -ve','HIV +ve'),each=4))

pcuntx <- 1e2*sum(mosdat$value[mosdat$tbtx=='Untreated for TB']) / sum(mosdat$value)
fn <- paste0('pcuntx',useHIVneg,'_priv',sensNote,'.txt')
fn <- here('disaggregation/R/02children/plots',fn)
cat(pcuntx,file=fn)

pchiv <- 1e2*sum(mosdat$value[mosdat$hiv=='HIV +ve']) / sum(mosdat$value)
fn <- paste0('pchiv',useHIVneg,'_priv',sensNote,'.txt')
fn <- here('disaggregation/R/02children/plots',fn)
cat(pchiv,file=fn)


## ------- country data ---------

mid <- function(x) median(x)
ub <- function(x) quantile(x,.975)
lb <- function(x) quantile(x,.025)
uq <- function(x) quantile(x,.75)
lq <- function(x) quantile(x,.25)
midr <- function(x) round(mid(x))
ubr <- function(x)  round(ub(x))
lbr <- function(x)  round(lb(x))

RESAc <- RESA[,list(age=mid(DeathsOnTx+DeathsOffTx+DeathsOnTxH+
                            DeathsOffTxH),
                    aget=uq(DeathsOnTx+DeathsOffTx+DeathsOnTxH+
                            DeathsOffTxH),
                    ageb=lq(DeathsOnTx+DeathsOffTx+DeathsOnTxH+
                            DeathsOffTxH)),
              by=iso3]
RESAc <- RESAc[order(age,decreasing=TRUE),]
RESAc$iso3 <- factor(RESAc$iso3,levels=RESAc$iso3,ordered = TRUE)
tmp <- RESAc[1:20,]
tmp <- merge(tmp,WHOkey,by='iso3',all.y=FALSE)
tmp$country <- factor(tmp$country,levels=rev(tmp$country),ordered=TRUE)
100*sum(RESAc[1:10,age])/sum(RESAc[,age])

pcgs <- 100*(RESAc[1:10,age])/sum(RESAc[,age])
fn <- paste0('pcgs',useHIVneg,'_priv',sensNote,'.txt')
fn <- here('disaggregation/R/02children/plots',fn)
cat(pcgs,file=fn)

fn <- paste0('top10countries',useHIVneg,'_priv',sensNote,'.csv')
fn <- here('disaggregation/R/02children/plots',fn)
write.csv(file=fn,cbind(RESAc[1:10,],pcgs=round(pcgs),cpcgs=round(cumsum(pcgs))))


RESAC <- RESA[,.(deaths = midr(DeathsOnTx+DeathsOffTx+
                               DeathsOnTxH+DeathsOffTxH),
                 deathslo = lbr(DeathsOnTx+DeathsOffTx+
                                DeathsOnTxH+DeathsOffTxH),
                 deathshi = ubr(DeathsOnTx+DeathsOffTx+
                                DeathsOnTxH+DeathsOffTxH),
                 deathsY = midr(DeathsOnTxY+DeathsOffTxY+
                                DeathsOnTxHY+DeathsOffTxHY),
                 deathsYlo = lbr(DeathsOnTxY+DeathsOffTxY+
                                 DeathsOnTxHY+DeathsOffTxHY),
                 deathsYhi = ubr(DeathsOnTxY+DeathsOffTxY+
                                 DeathsOnTxHY+DeathsOffTxHY),
                 deathsH = midr(DeathsOnTxH+DeathsOffTxH),
                 deathsHlo = lbr(DeathsOnTxH+DeathsOffTxH),
                 deathsHhi = ubr(DeathsOnTxH+DeathsOffTxH),
                 ## young no HIV+/-
                 deathsHY = midr(DeathsOnTxHY+DeathsOffTxHY),
                 deathsHYlo = lbr(DeathsOnTxHY+DeathsOffTxHY),
                 deathsHYhi = ubr(DeathsOnTxHY+DeathsOffTxHY),
                 deathsNHY = midr(DeathsOnTxY+DeathsOffTxY),
                 deathsNHYlo = lbr(DeathsOnTxY+DeathsOffTxY),
                 deathsNHYhi = ubr(DeathsOnTxY+DeathsOffTxY),
                 ## old no HIV+/-
                 deathsHO = midr(DeathsOnTxH+DeathsOffTxH-
                                 DeathsOnTxHY+DeathsOffTxHY),
                 deathsHOlo = lbr(DeathsOnTxH+DeathsOffTxH-
                                  DeathsOnTxHY-DeathsOffTxHY),
                 deathsHOhi = ubr(DeathsOnTxH+DeathsOffTxH-
                                  DeathsOnTxHY-DeathsOffTxHY),
                 deathsNHO = midr(DeathsOnTx+DeathsOffTx-
                                  DeathsOnTxY+DeathsOffTxY),
                 deathsNHOlo = lbr(DeathsOnTx+DeathsOffTx-
                                   DeathsOnTxY-DeathsOffTxY),
                 deathsNHOhi = ubr(DeathsOnTx+DeathsOffTx-
                                   DeathsOnTxY-DeathsOffTxY),
                 ## previous stuff
                 DeathsOnTx = midr(DeathsOnTx),
                 DeathsOnTxlo = lbr(DeathsOnTx),
                 DeathsOnTxhi = ubr(DeathsOnTx),
                 DeathsOnTxY = midr(DeathsOnTxY),
                 DeathsOnTxYlo = lbr(DeathsOnTxY),
                 DeathsOnTxYhi = ubr(DeathsOnTxY),
                 DeathsOffTx = midr(DeathsOffTx),
                 DeathsOffTxlo = lbr(DeathsOffTx),
                 DeathsOffTxhi = ubr(DeathsOffTx),
                 DeathsOffTxY = midr(DeathsOffTxY),
                 DeathsOffTxYlo = lbr(DeathsOffTxY),
                 DeathsOffTxYhi = ubr(DeathsOffTxY),
                 HIV = midr(1e3*HIV),HIVinTB = midr(1e3*HIVinTB),
                 HIVinDeaths = midr(1e3*HIVinDeaths),
                 HIVlo = lbr(1e3*HIV),HIVinTBlo = lbr(1e3*HIVinTB),
                 HIVinDeathslo = lbr(1e3*HIVinDeaths),
                 HIVhi = ubr(1e3*HIV),HIVinTBhi = ubr(1e3*HIVinTB),
                 HIVinDeathshi = ubr(1e3*HIVinDeaths)
                 ),
              by=list(iso3,g_whoregion)]


RESAC <- merge(RESAC,WHOkey,by='iso3')

KM <- RESAC
save(KM,file=here('disaggregation/R/02children/data/KM.Rdata'))

## NB this is per MILLE not per CENT
if(!sensNote & !useHIVneg){
  fn <- paste0('countries',useHIVneg,'_priv',sensNote,'.csv')
  fn <- here('disaggregation/R/02children/plots',fn)
  write.csv(file=fn,RESAC)
  RESACH <- RESAC[HIV>0,list(iso3,HIV,HIVinTB,HIVinDeaths)]
  fn <- paste0('countries_HIVpm',useHIVneg,'_priv',sensNote,'.csv')
  fn <- here('disaggregation/R/02children/plots',fn)
  write.csv(file=fn,RESACH)
}

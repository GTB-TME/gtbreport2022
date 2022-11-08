library(here)
library(glue)
library(loo)
library(data.table)
library(ggplot2)
library(ggrepel)
library(ggpubr)
gh <- function(x) glue(here(x))
rot45 <- theme(axis.text.x = element_text(angle = 45, hjust = 1))
source(gh('drtb/R/utils/plotsetc.R'))


## ================== STATS ================
vrs <- dir(path='drtb/outdata',pattern='loo')
vrs <- gsub('loo\\.out\\.','',vrs)
vrs <- gsub('\\.Rdata','',vrs)
vrs <- vrs[!grepl('noH',vrs)]

## LOO stats
L <- list()
for(vnt in vrs){
  fn <- gh('drtb/outdata/loo.out.{vnt}.Rdata')
  load(fn)
  L[[vnt]] <- loo.out
}
length(L)



(LR <- loo_compare(L))
LRD <- as.data.table(LR)
LRD[,model:=rownames(LR)]


## WAIC stats
W <- list()
for(vnt in vrs){
  fn <- gh('drtb/outdata/waic.out.{vnt}.Rdata')
  load(fn)
  W[[vnt]] <- waic.out
}
length(W)

loo_compare(W) #similar


## log-like
LL <- list()
for(vnt in vrs){
  fn <- gh('drtb/outdata/dll.{vnt}.Rdata')
  load(fn)
  LL[[vnt]] <- dll
}
length(LL)
LL[[1]]

LD <- rbindlist(LL)
LD[,model:=names(LL)]

## join
LRD <- merge(LRD,LD,by='model')
LRD <- LRD[rev(order(elpd_diff))]

nmz <- c('model','elpd_diff','se_diff','p_loo','elpd_loo','data.likelihood.mean')
rst <- setdiff(names(LRD),nmz)
nmz <- c(nmz,rst)

setcolorder(LRD,nmz)

LRD

fwrite(LRD,file=gh('drtb/outdata/compare.csv'))


LRD[,.(model,elpd_diff,elpd_loo,data.likelihood.mean)]


## =================== PLOTS ========================

load(here('drtb/data/regkey.Rdata'))
load(here('drtb/data/RPD.Rdata'))
load(here('drtb/data/HBC.Rdata'))

RPD <- merge(RPD[year>=2000],regkey,by='iso3')
RPD <- merge(RPD,HBC,by='iso3')
hbc <- HBC[g.hbmdr==TRUE,unique(iso3)] #list of HBC30 DR

## RPD <- RPD[!(iso3=='SOM' & year==2020)] #NOTE drop

## most recent data
MR <- RPD[g.hbmdr==TRUE,
          .SD[year==max(year),.(year,RR.mid,RR.lo,RR.hi)],
          by=.(iso3,g_whoregion,patients)]
## element_rect(fill="red")

## by checking compare.csv
vrtz <- c('together_LerouxInterceptLerouxSlope_2',
          'separate_LerouxInterceptLerouxSlope_2',
          'together_LerouxInterceptLerouxSlope_4')
svo <- TRUE #saveout
GPeg <- LP <- list()
## loop over variants
for(vnt in vrtz){

  ## load
  load(gh('drtb/outdata/KO.{vnt}.Rdata'))
  KO <- merge(KO,regkey)

  ## all time
  GP <- ggplot(RPD[iso3 %in% hbc],
               aes(year,RR.mid,
                   ymin=RR.lo,ymax=RR.hi,
                   shape=type,col=patients,lty=type))+
    geom_point()+geom_pointrange()+
    geom_line(data=KO[iso3 %in% hbc])+
    geom_ribbon(data=KO[iso3 %in% hbc],aes(fill=patients),
                alpha=0.2,col=NA)+
    scale_y_sqrt()+
    facet_wrap(~country,scales='free',nrow=6,ncol=5)+
    rot45+theme_light()+
    theme(strip.background = element_blank())+
    theme(strip.text = element_text(colour = 'black'))+
    xlab('Year') + ylab('Rifampicin resistance (%), square root scale')+
    geom_vline(xintercept = 2015,col='darkblue',lty=2)
  ## GP

  ## ggsave(GP,file=gh('neatplots/HBC30allt.{vnt}.pdf'),h=18,w=15)
  if(svo)
    ggsave(GP,file=gh('drtb/plots/HBC30allt.{vnt}.png'),h=18,w=15)

  eg <- 'KGZ'
  GPeg[[vnt]] <- ggplot(RPD[iso3 == eg],
               aes(year,RR.mid,
                   ymin=RR.lo,ymax=RR.hi,
                   shape=type,col=patients,lty=type))+
    geom_point()+geom_pointrange()+
    geom_line(data=KO[iso3 == eg])+
    geom_ribbon(data=KO[iso3 == eg],aes(fill=patients),
                alpha=0.2,col=NA)+
    scale_y_sqrt()+
    facet_wrap(~country,scales='free',nrow=6,ncol=5)+
    rot45+theme_light()+
    theme(strip.background = element_blank())+
    theme(strip.text = element_text(colour = 'black'))+
    xlab('Year') + ylab('Rifampicin resistance (%), square root scale')+
    geom_vline(xintercept = 2015,col='darkblue',lty=2)



  ## restrict by time
  GP <- ggplot(RPD[iso3 %in% hbc & year >=2015],
               aes(year,RR.mid,
                   ymin=RR.lo,ymax=RR.hi,
                   shape=type,col=patients,lty=type))+
    geom_point()+geom_pointrange()+
    geom_line(data=KO[iso3 %in% hbc  & year >=2015])+
    geom_ribbon(data=KO[iso3 %in% hbc  & year >=2015],aes(fill=patients),
                alpha=0.2,col=NA)+
    scale_y_sqrt()+
    facet_wrap(~country,scales='free',nrow=6,ncol=5)+
    rot45+theme_light()+
    theme(strip.background = element_blank())+
    theme(strip.text = element_text(colour = 'black'))+
    xlab('Year') + ylab('Rifampicin resistance (%), square root scale')+
    geom_vline(xintercept = 2015,col='darkblue',lty=2)
  ## GP

  ## ggsave(GP,file=gh('neatplots/HBC30all15.{vnt}.pdf'),h=18,w=15)
  if(svo)
    ggsave(GP,file=gh('drtb/plots/HBC30all15.{vnt}.png'),h=18,w=15)

  ## most recent & difference
  MRB <- merge(MR,KO[,.(iso3,year,patients,
                        e.RR.mid=RR.mid,
                        e.RR.lo=RR.lo,e.RR.hi=RR.hi)],
               by=c('iso3','year','patients'),
               all.x=TRUE,all.y=FALSE)
  MRB[patients=='ret',patients:='retreatment']

  GP <-
  ggplot(MRB,aes(x=RR.mid,xmin=RR.lo,xmax=RR.hi,label=iso3,
                 y=e.RR.mid,ymin=e.RR.lo,ymax=e.RR.hi))+
    geom_errorbar()+geom_errorbarh()+geom_point()+
    facet_wrap(~patients)+
    geom_abline(intercept = 0,slope=1,col=2)+
    geom_text_repel()+
    theme_light()+
    theme(strip.background = element_blank())+
    theme(strip.text = element_text(colour = 'black'))+
    xlab('Most recent empirical measure of rifampicin resistance (%, square root)')+
    ylab('Corresponding model estimate of rifampicin resistance (%, square root)')
  GP <- GP + scale_y_sqrt() + scale_x_sqrt()
  LP[[vnt]] <- GP

  f <- 0.7
  if(svo)
    ggsave(GP,file=gh('drtb/plots/Compare.{vnt}.png'),h=10*f,w=20*f)

}
## element_rect(fill="red")


GA <- ggarrange(plotlist = LP,nrow=3,ncol=1,labels = LETTERS[1:3])
ggsave(GA,file=gh('drtb/plots/Compare.BOTH.png'),h=10*f*3,w=20*f)


GA <- ggarrange(plotlist = GPeg,nrow=1,ncol=3,labels = LETTERS[1:3],
                common.legend = TRUE,legend='top')

ggsave(GA,file=gh('drtb/plots/Compare.RE.png'),h=5,w=12)




## =================== CrossValidation ========================
## NOTE not including the CV data in main repo as rather large
## load data
load(gh('drtb/data/D4S.Rdata'))
load(here('drtb/data/regkey.Rdata'))
load(here('drtb/data/RPD.Rdata'))
## data prep
## ckey <- data.table(iso3=isoidx,cid=1:length(isoidx))
RPD <- merge(RPD[year>=2000],regkey[,.(iso3,g_whoregion)])
RPD[,max(year),by=type]
maxy <- max(RPD$year)

c.omit.v <- scan(gh('drtb/data/c.omit.vA.txt'),what = 'string')#NOTE change list here

RPD <- RPD[iso3 %in% c.omit.v,
           .(iso3,g_whoregion,patients,dy=maxy-year,RR.mid,RR.lo,RR.hi)]

## log normal
getms <- function(M,L,H){
  M[M==0] <- 1 #safety
  s <- (H-L)/3.92
  list(mu=log(M/(sqrt(M^2+s^2))),
       sig1=sqrt(log(1+s^2/M^2)))
}

RPD[,c('mu','sig'):=getms(RR.mid,RR.lo,RR.hi)] #add log-normal parms
RPD[!is.finite(mu)]

rootlst <- c('together_LerouxInterceptLerouxSlope_2',
             'separate_LerouxInterceptLerouxSlope_2',
             'together_LerouxIntercept_4',
             'separate_global_4')

meann <- function(x) mean(x[is.finite(x)])

SGA <- SRA <- SCA <- list()
for(i in 1:length(rootlst)){
  cat('----- starting work on ',i,'...\n')
  ## loop read-in data
  L <- list()
  for(cn in c.omit.v){
    fn <- gh('drtb/outdata/K.{rootlst[i]}.{cn}.Rdata')
    if(file.exists(fn)){
      load(fn)
      L[[cn]] <- K
    }
  }
  K <- rbindlist(L)
  K <- K[dst=='RR']
  K <- merge(K,RPD,by=c('iso3','patients','dy'),
        all.x = FALSE,all.y=FALSE)
  ## K <- merge(K,RPD,by=c('iso3','patients','dy'),
  ##            all.x = TRUE,all.y=FALSE)
  ## K <- K[!is.na(RR.mid)] #not all years have data
  K[,pc:=1e2*value]      #work with %ages
  K[,AE:=abs(RR.mid-pc)] #absolute error
  K[,APE:=1e2*abs(RR.mid-pc)/pc] #absolute % error
  K[,inout:=ifelse(pc<RR.hi,1,0) * ifelse(pc>RR.lo,1,0)]
  K[,neglogp:=-dlnorm(pc,mu,sig,log=TRUE)]
  cls <- c('coverage','MAE','MAPE','H','U')
  SC <- K[,.(
    coverage=1e2*mean(inout,na.rm=TRUE), #coverage
    MAE=abs(mean(pc)-RR.mid),    #MAE
    MAPE=1e2*abs(1-RR.mid/mean(pc)),    #MAPE
    H=mean(neglogp,na.rm=TRUE),  #H
    U=IQR(pc)  #uncertainty
    ),by=.(iso3,patients,g_whoregion)] #mean over times/reps in each iso
  SR <- SC[,lapply(.SD,meann),.SDcols=cls,by=.(g_whoregion,patients)]
  SG <- SC[,lapply(.SD,meann),.SDcols=cls,by=patients]
  SG[,model:=rootlst[i]]
  SR[,model:=rootlst[i]]
  SC[,model:=rootlst[i]]
  G <- K[,.(H=mean(neglogp)),by=.(id,patients)]
  G <- G[,.(H.sd=sd(H)),by=patients]
  SG <- merge(SG,G,by='patients')
  ## gather
  SGA[[i]] <- SG
  SRA[[i]] <- SR
  SCA[[i]] <- SC
  cat('----- done work on ',i,'...\n')
}
SGA <- rbindlist(SGA)
SRA <- rbindlist(SRA)
SCA <- rbindlist(SCA)

save(SGA,file=gh('drtb/outdata/SGA.Rdata'))
save(SRA,file=gh('drtb/outdata/SRA.Rdata'))
save(SCA,file=gh('drtb/outdata/SCA.Rdata'))

(tn <- SGA[patients=='new',
           .(model,H,H.sd,MAPE,MAE,U,C=coverage)][order(H)])
(tr <- SGA[patients=='ret',
           .(model,H,H.sd,MAPE,MAE,U,C=coverage)][order(H)])
fwrite(tn,file=gh('drtb/outdata/tn.csv'))
fwrite(tr,file=gh('drtb/outdata/tr.csv'))


(tnr <- SRA[patients=='new',
            .(g_whoregion,model,
              H,MAPE,MAE,U,C=coverage)][order(g_whoregion,H)])
(trr <- SRA[patients=='ret',
            .(g_whoregion,model,
              H,MAPE,MAE,U,C=coverage)][order(g_whoregion,H)])
fwrite(tnr,file=gh('drtb/outdata/tnr.csv'))
fwrite(trr,file=gh('drtb/outdata/trr.csv'))


tn[model %in% c('together_LerouxInterceptLerouxSlope_2','separate_LerouxInterceptLerouxSlope_2')]
tr[model %in% c('together_LerouxInterceptLerouxSlope_2','separate_LerouxInterceptLerouxSlope_2')]
tnr[model %in% c('together_LerouxInterceptLerouxSlope_2','separate_LerouxInterceptLerouxSlope_2')]
trr[model %in% c('together_LerouxInterceptLerouxSlope_2','separate_LerouxInterceptLerouxSlope_2')]


## separate: lowest new MAE, ret H MAPE

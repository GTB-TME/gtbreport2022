## this is for reporting data things and re-organizing to work with Stan
## libraries/utilities used in both
library(here)
library(glue)
library(data.table)
gh <- function(x) glue(here(x))

## used only in plots data section
library(ggplot2)
library(ggrepel)
rot45 <- theme(axis.text.x = element_text(angle = 45, hjust = 1))

## helper functions for binomial CIs
source(gh('drtb/R/utils/CIhelpers.R'))

## load processed DRS data
load(gh('drtb/data/LB.Rdata'))
load(gh('drtb/data/YB.Rdata'))

## from r_rlt:
##   rr (R-res regardless of H)
## 
## from dst_rlt
##   dst_rlt_hr (H-res regardless of R),
##   dst_rlt_rr (R-res regardless of H),
##   dr_h_nr (H-res not R-res),
##   dr_r_nh (R-res not H-res),
##   mdr

## ===================================================
## --- plots data ---
## for RR
## surveillance
LB[case %in% c('r','r.h','r.m'),c('RR.k','RR.N'):=.(rr,r_rlt)]
LB[case %in% c('dh0','dh2'),c('RR.k','RR.N'):=.(dr_r_nh + mdr,dst_rlt)]
LB[case %in% c('dh0.x','dh2.x'),
   c('RR.k','RR.N'):=.(dr_r_nh + mdr+xpert_dr_r,dst_rlt+xpert)]

LR <- LB[,.(id,patients,type='surveillance',iso3,year,
            case,RR.k,RR.N)]
LR <- LR[!is.na(RR.k) & !is.na(RR.N)]

LR[,c('RR.mid','RR.lo','RR.hi'):=MLH(RR.k,RR.N)]


## surveys
YB[case %in% c('r','r.hnr','r.hr'),
   c('RR.mid','RR.lo','RR.hi'):=.(rr_pcnt,rr_pcnt_lo,rr_pcnt_hi)]


YB[case %in% c('dh0','dh1.h','dh0.x'),
   c('se1','se2'):=.((dr_r_nh_pcnt_hi - dr_r_nh_pcnt_lo)/3.92,
   (mdr_pcnt_hi - mdr_pcnt_lo)/3.92)]
YB[,se:=sqrt(se1^2+se2^2)]
YB[,md:=dr_r_nh_pcnt+mdr_pcnt]

YB[case %in% c('dh0','dh1.h'),
   c('RR.mid','RR.lo','RR.hi'):=.(md,
                                  pmax(0,md-1.96*se),
                                  pmin(100,md+1.96*se))]

## use Xpert as it is
YB[case %in% c('x','r.hr.x','dh2.x'),
   c('RR.mid','RR.lo','RR.hi'):=.(xpert_dr_r_pcnt,
                                  xpert_dr_r_pcnt_lo,
                                  xpert_dr_r_pcnt_hi)]

## var-weighted average of xpert nad rest
YB[case %in% c('dh0.x'),se1:=se] #from above
YB[case %in% c('dh0.x'),
   se2:=(xpert_dr_r_pcnt_hi-xpert_dr_r_pcnt_lo)/3.92]
YB[case %in% c('dh0.x'),se:=sqrt(se1^2+se2^2)]
YB[case %in% c('dh0.x'),md:=(se1/se)^2*md + (se2/se)^2*xpert_dr_r_pcnt]
YB[case %in% c('dh0.x'),
   c('RR.mid','RR.lo','RR.hi'):=.(md,
                                  pmax(0,md-1.96*se),
                                  pmin(100,md+1.96*se))]
YR <- YB[,.(id,patients,type='survey',iso3,year,
            case,RR.mid,RR.lo,RR.hi)]
YR <- YR[!is.na(RR.mid)]
YR[,c('RR.k','RR.N'):=NA]

## join
RPD <- rbind(LR,YR,use.names = TRUE)

save(RPD,file=here('drtb/data/RPD.Rdata'))

## quick look
load(file=here('drtb/data/RPD.Rdata'))
load(file=here('drtb/data/R.Rdata'))

RPD <- merge(RPD,unique(R[,.(iso3,g_whoregion)]),all.x=TRUE,by='iso3')

for(reg in unique(RPD[,g_whoregion])){
  print(reg)
  GP <- ggplot(RPD[g_whoregion==reg],aes(year,RR.mid,
                       ymin=RR.lo,ymax=RR.hi,
                       shape=type,col=patients,lty=type))+
    geom_line()+geom_point()+geom_pointrange()+
    scale_y_sqrt()+
    facet_wrap(~iso3,scales='free')+
    rot45
  ggsave(GP,file=gh('drtb/plots/r_{reg}.pdf'),w=12,h=15)
}



## ==================== STAN DATA ===============================
## We assume that each point measured (for category $j$ at site $i$) is log-normally distributed
## 
## $$\pi_{i,j} \sim \mathrm{LogNorm}(\mu_{i,j},\sigma_{i,j}) $$
## 
## where $\mu_{i,j}$ and $\sigma_{i,j}$ are determined by
##
## $$\mu_{i,j} = \log\left(\frac{p_{i,j}}{\sqrt{1+v_{i,j}/p_{i,j}^2}}\right)$$
## $$\sigma^2_{i,j} = \log\left(1+v_{i,j}/p_{i,j}^2\right)$$
## $$v_{i,j} = \left(\frac{\pi^H_{i,j}-\pi^L_{i,j}}{2\times 1.96}\right)^2$$

## --- stan data
getmusig <- function(lo,mid,hi){
  mid <- mid/1e2; lo <- lo/1e2; hi <- hi/1e2 #convert pcnt to number
  mid[mid==0] <- 1e-6                   #safety for zeros
  v <- ((hi-lo)/3.92)^2
  v[v==0] <- 1e-6
  mu <- log(mid/sqrt(1+v/mid^2))
  sig2 <- log(1+v/mid^2)
  list(mu=mu,sig=sqrt(sig2))
}


## load 5 nearest country structure (generated above)
load(here('drtb/data/W5.Rdata'))
load(here('drtb/data/isoidx5.Rdata'))
load(here('data/cty.rda'))
dim(cty); length(isoidx); dim(W)
regkey <- cty[iso3 %in% isoidx]

for(i in 1:length(isoidx)) regkey[iso3==isoidx[i],id:=i]
regkey[is.na(id)] #none
## regkey <- regkey[!is.na(id)]

dim(regkey)
regkey[,range(id)]
## (missed <- setdiff(1:length(isoidx),regkey$id))
## isoidx[missed]
## regkey[id==missed]

## former soviet union countries:
## armenia,azer,bel,est,geo,kaz,kyr,lat,lith,mold,rus,taj,turk,ukr,uzb
## regkey[grepl('Uz',country)]
FSUlist <- c('ARM','AZE','BLR','EST','GEO','KAZ','KGZ','LVA','LTU',
             'MDA','RUS','TJK','TKM','UKR','UZB')
length(FSUlist)
regkey[,region:=g.whoregion]
regkey[iso3 %in% FSUlist,region:='FSU']
regkey[iso3=='RUS'] #check
regmm <- model.matrix(data=regkey,~region)
dim(regmm)
colnames(regmm)

## need spatial data to decide order
## may need redoing; see DRHMM 1st
## see also how region is coded
load(gh('drtb/data/LB.Rdata'))
load(gh('drtb/data/YB.Rdata'))

## restrict if needed
LB <- LB[year>=2000]
YB <- YB[year>=2000]

## maximum time
maxy <- max(LB$year,YB$year)
LB[,dy:=as.integer(maxy-year)]
YB[,dy:=as.integer(maxy-year)]

## merge in id from regkey
LB <- merge(LB,regkey[,.(iso3,cid=id)],by='iso3',all.x=TRUE)
YB <- merge(YB,regkey[,.(iso3,cid=id)],by='iso3',all.x=TRUE)

## === surveillance data
## need to think about new/ret - keep searate
LD <- list()

## --- [rr,r_rlt] Xrr = [whoregion,dt]

## new RR
css <- c('r','r.h','r.m') #cases
(tmp <- LB[patients=='new' & case %in% css,.(cid,dy,rr,r_rlt)])
tmp <- na.omit(tmp)
## include in stan data
LD[['L.new.rr.id']] <- as.array(tmp[,cid])
LD[['L.new.rr.N']] <- cbind(tmp[,r_rlt-rr],tmp[,rr])
LD[['L.new.rr.t']] <- as.array(tmp$dy)
## ret RR
(tmp <- LB[patients=='ret' & case %in% css,.(cid,dy,rr,r_rlt)])
tmp <- na.omit(tmp)
## include in stan data
LD[['L.ret.rr.id']] <- as.array(tmp[,cid])
LD[['L.ret.rr.N']] <- cbind(tmp[,r_rlt-rr],tmp[,rr])
LD[['L.ret.rr.t']] <- as.array(tmp$dy)

## --- RH

## new r.h
css <- c('r.h') #cases
(tmp <- LB[patients=='new' & case %in% css,.(cid,dy,dst_rlt_hr,dst_rlt)])
tmp <- tmp[!is.na(dst_rlt)]
## tmp[is.na(dst_rlt_hr),dst_rlt_hr:=0] #not needed
## include in stan data
LD[['L.new.hr.id']] <- as.array(tmp[,cid])
LD[['L.new.hr.N']] <- cbind(tmp[,dst_rlt-dst_rlt_hr],tmp[,dst_rlt_hr])
LD[['L.new.hr.t']] <- as.array(tmp$dy)
## new r.h
(tmp <- LB[patients=='ret' & case %in% css,.(cid,dy,dst_rlt_hr,dst_rlt)])
tmp <- tmp[!is.na(dst_rlt)]
## tmp[is.na(dst_rlt_hr),dst_rlt_hr:=0] #not needed
## include in stan data
LD[['L.ret.hr.id']] <- as.array(tmp[,cid])
LD[['L.ret.hr.N']] <- cbind(tmp[,dst_rlt-dst_rlt_hr],tmp[,dst_rlt_hr])
LD[['L.ret.hr.t']] <- as.array(tmp$dy)


## --- RM

## new r.m
css <- c('r.m') #cases
(tmp <- LB[patients=='new' & case %in% css,.(cid,dy,mdr,dst_rlt)])
## include in stan data
LD[['L.new.rm.id']] <- as.array(tmp[,cid])
LD[['L.new.rm.N']] <- cbind(tmp[,dst_rlt-mdr],tmp[,mdr])
LD[['L.new.rm.t']] <- as.array(tmp$dy)
## ret r.m
(tmp <- LB[patients=='ret' & case %in% css,.(cid,dy,mdr,dst_rlt)])
## include in stan data
LD[['L.ret.rm.id']] <- as.array(tmp[,cid])
LD[['L.ret.rm.k']] <- tmp[,mdr]
LD[['L.ret.rm.N']] <- cbind(tmp[,dst_rlt-mdr],tmp[,mdr])
LD[['L.ret.rm.t']] <- as.array(tmp$dy)


## --- [xpert_dr_r,xpert]

## new x
css <- unique(grep('x',LB$case,value=TRUE))
(tmp <- LB[patients=='new' & case %in% css,.(cid,dy,xpert_dr_r,xpert)])
tmp <- tmp[xpert>0]     #NOTE exclude cases with 0 denom
## include in stan data
LD[['L.new.x.id']] <- as.array(tmp[,cid])
LD[['L.new.x.N']] <- cbind(tmp[,xpert-xpert_dr_r],tmp[,xpert_dr_r])
LD[['L.new.x.t']] <- as.array(tmp$dy)
## ret x
(tmp <- LB[patients=='ret' & case %in% css,.(cid,dy,xpert_dr_r,xpert)])
tmp <- tmp[xpert>0]     #NOTE exclude cases with 0 denom
X <- regmm[tmp$cid,] #regional information
## include in stan data
LD[['L.ret.x.id']] <- as.array(tmp[,cid])
LD[['L.ret.x.N']] <- cbind(tmp[,xpert-xpert_dr_r],tmp[,xpert_dr_r])
LD[['L.ret.x.t']] <- as.array(tmp$dy)

## --- DH cases
## NOTE convention to keep order for preserved buckets
## dst_rlt: [hnr,rnh,mdr], [dst_rlt_hr]
## dh0, dh1,dh2

## -- new
css <- unique(grep('dh',LB$case,value=TRUE))
(TMP <- LB[patients=='new' & case %in% css,
           .(cid,dy,dr_h_nr,dr_r_nh,mdr,dst_rlt,case)])

## dh0
css <- c('dh0','dh0.x')
tmp <- TMP[case %in% css]
## include in stan data
LD[['L.new.dh0.id']] <- as.array(tmp[,cid])
LD[['L.new.dh0.N']] <- as.matrix(tmp[,.(dst_rlt-dr_h_nr-dr_r_nh-mdr,
                                        dr_h_nr,dr_r_nh,mdr)])
LD[['L.new.dh0.t']] <- as.array(tmp$dy)
## dh1 (RNH missing)
css <- c('dh1','dh1.x')
tmp <- TMP[case %in% css]
## include in stan data
LD[['L.new.dh1.id']] <- as.array(tmp[,cid])
LD[['L.new.dh1.N']] <- as.matrix(tmp[,.(dst_rlt-dr_h_nr-mdr,
                                        dr_h_nr,mdr)])
LD[['L.new.dh1.t']] <- as.array(tmp$dy)
## dh2
css <- c('dh2','dh2.x')
tmp <- TMP[case %in% css]
## include in stan data
LD[['L.new.dh2.id']] <- as.array(tmp[,cid])
LD[['L.new.dh2.N']] <- as.matrix(tmp[,.(dst_rlt-mdr,mdr)])
LD[['L.new.dh2.t']] <- as.array(tmp$dy)

## -- ret
css <- unique(grep('dh',LB$case,value=TRUE))
(TMP <- LB[patients=='ret' & case %in% css,
           .(cid,dy,dr_h_nr,dr_r_nh,mdr,dst_rlt,case)])

## dh0
css <- c('dh0','dh0.x')
tmp <- TMP[case %in% css]
## include in stan data
LD[['L.ret.dh0.id']] <- as.array(tmp[,cid])
LD[['L.ret.dh0.N']] <- as.matrix(tmp[,.(dst_rlt-dr_h_nr-dr_r_nh-mdr,
                                        dr_h_nr,dr_r_nh,mdr)])
LD[['L.ret.dh0.t']] <- as.array(tmp$dy)
## dh1
css <- c('dh1','dh1.x')
tmp <- TMP[case %in% css]
## include in stan data
LD[['L.ret.dh1.id']] <- as.array(tmp[,cid])
LD[['L.ret.dh1.N']] <- as.matrix(tmp[,.(dst_rlt-dr_h_nr-mdr,
                                        dr_h_nr,mdr)])
LD[['L.ret.dh1.t']] <- as.array(tmp$dy)
## dh2
css <- c('dh2','dh2.x')
tmp <- TMP[case %in% css]
## include in stan data
LD[['L.ret.dh2.id']] <- as.array(tmp[,cid])
LD[['L.ret.dh2.N']] <- as.matrix(tmp[,.(dst_rlt-mdr,mdr)])
LD[['L.ret.dh2.t']] <- as.array(tmp$dy)

## === survey data
YD <- list()

## [M,S]x:
## --- [rr]
## new
css <- c('r','r.hr','r.hnr','r.hr.x')
(tmp <- YB[patients=='new' & case %in% css,
           .(cid,dy,case,rr_pcnt,rr_pcnt_lo,rr_pcnt_hi)])
MS <- getmusig(tmp$rr_pcnt_lo,tmp$rr_pcnt,tmp$rr_pcnt_hi)
YD[['Y.new.rr.id']] <- as.array(tmp[,cid])
YD[['Y.new.rr.M']] <- MS$mu
YD[['Y.new.rr.S']] <- MS$sig
YD[['Y.new.rr.t']] <- as.array(tmp$dy)

## ret
css <- c('r','r.hr','r.hnr','r.hr.x')
(tmp <- YB[patients=='ret' & case %in% css,
           .(cid,dy,case,rr_pcnt,rr_pcnt_lo,rr_pcnt_hi)])
MS <- getmusig(tmp$rr_pcnt_lo,tmp$rr_pcnt,tmp$rr_pcnt_hi)
YD[['Y.ret.rr.id']] <- as.array(tmp[,cid])
YD[['Y.ret.rr.M']] <- MS$mu
YD[['Y.ret.rr.S']] <- MS$sig
YD[['Y.ret.rr.t']] <- as.array(tmp$dy)

## --- [xpert_dr_pcnt]
## new
css <- unique(grep('x',YB$case,value=TRUE))
(tmp <- YB[patients=='new' & case %in% css,
           .(cid,dy,case,xpert_dr_r_pcnt,
             xpert_dr_r_pcnt_lo,xpert_dr_r_pcnt_hi)])
MS <- getmusig(tmp$xpert_dr_r_pcnt_lo,tmp$xpert_dr_r_pcnt,
               tmp$xpert_dr_r_pcnt_hi)
YD[['Y.new.x.id']] <- as.array(tmp[,cid])
YD[['Y.new.x.M']] <- MS$mu
YD[['Y.new.x.S']] <- MS$sig
YD[['Y.new.x.t']] <- as.array(tmp$dy)

## ret
(tmp <- YB[patients=='ret' & case %in% css,
           .(cid,dy,case,xpert_dr_r_pcnt,
             xpert_dr_r_pcnt_lo,xpert_dr_r_pcnt_hi)])
MS <- getmusig(tmp$xpert_dr_r_pcnt_lo,tmp$xpert_dr_r_pcnt,
               tmp$xpert_dr_r_pcnt_hi)
YD[['Y.ret.x.id']] <- as.array(tmp[,cid])
YD[['Y.ret.x.M']] <- MS$mu
YD[['Y.ret.x.S']] <- MS$sig
YD[['Y.ret.x.t']] <- as.array(tmp$dy)

## --- [hnr,rnh,mdr]

## -- dh0
## new
css <- unique(grep('dh0',YB$case,value=TRUE))
(tmp <- YB[patients=='new' & case %in% css,
           .(cid,dy,case,
             dr_h_nr_pcnt,dr_r_nh_pcnt,mdr_pcnt,
             dr_h_nr_pcnt_lo,dr_r_nh_pcnt_lo,mdr_pcnt_lo,
             dr_h_nr_pcnt_hi,dr_r_nh_pcnt_hi,mdr_pcnt_hi)])
tmp <- na.omit(tmp)
MS.hnr <- getmusig(tmp$dr_h_nr_pcnt_lo,tmp$dr_h_nr_pcnt,
                   tmp$dr_h_nr_pcnt_hi)
MS.rnh <- getmusig(tmp$dr_r_nh_pcnt_lo,tmp$dr_r_nh_pcnt,
                   tmp$dr_r_nh_pcnt_hi)
MS.mdr <- getmusig(tmp$mdr_pcnt_lo,tmp$mdr_pcnt,tmp$mdr_pcnt_hi)
YD[['Y.new.dh0.id']] <- as.array(tmp[,cid])
YD[['Y.new.dh0.M']] <- cbind(MS.hnr$mu,MS.rnh$mu,MS.mdr$mu)
YD[['Y.new.dh0.S']] <- cbind(MS.hnr$sig,MS.rnh$sig,MS.mdr$sig)
YD[['Y.new.dh0.t']] <- as.array(tmp$dy)
## ret
(tmp <- YB[patients=='ret' & case %in% css,
           .(cid,dy,case,
             dr_h_nr_pcnt,dr_r_nh_pcnt,mdr_pcnt,
             dr_h_nr_pcnt_lo,dr_r_nh_pcnt_lo,mdr_pcnt_lo,
             dr_h_nr_pcnt_hi,dr_r_nh_pcnt_hi,mdr_pcnt_hi)])
tmp <- na.omit(tmp)
MS.hnr <- getmusig(tmp$dr_h_nr_pcnt_lo,tmp$dr_h_nr_pcnt,
                   tmp$dr_h_nr_pcnt_hi)
MS.rnh <- getmusig(tmp$dr_r_nh_pcnt_lo,tmp$dr_r_nh_pcnt,
                   tmp$dr_r_nh_pcnt_hi)
MS.mdr <- getmusig(tmp$mdr_pcnt_lo,tmp$mdr_pcnt,tmp$mdr_pcnt_hi)
YD[['Y.ret.dh0.id']] <- as.array(tmp[,cid])
YD[['Y.ret.dh0.M']] <- cbind(MS.hnr$mu,MS.rnh$mu,MS.mdr$mu)
YD[['Y.ret.dh0.S']] <- cbind(MS.hnr$sig,MS.rnh$sig,MS.mdr$sig)
YD[['Y.ret.dh0.t']] <- as.array(tmp$dy)

## -- dh1r
## new
css <- unique(grep('dh1.r',YB$case,value=TRUE))
(tmp <- YB[patients=='new' & case %in% css,
           .(cid,dy,case,
             dr_h_nr_pcnt,dr_r_nh_pcnt,mdr_pcnt,
             dr_h_nr_pcnt_lo,dr_r_nh_pcnt_lo,mdr_pcnt_lo,
             dr_h_nr_pcnt_hi,dr_r_nh_pcnt_hi,mdr_pcnt_hi)])
## NOTE none!

## ret
(tmp <- YB[patients=='ret' & case %in% css,
           .(cid,dy,case,
             dr_h_nr_pcnt,mdr_pcnt,
             dr_h_nr_pcnt_lo,mdr_pcnt_lo,
             dr_h_nr_pcnt_hi,mdr_pcnt_hi)])
MS.hnr <- getmusig(tmp$dr_h_nr_pcnt_lo,tmp$dr_h_nr_pcnt,
                   tmp$dr_h_nr_pcnt_hi)
MS.mdr <- getmusig(tmp$mdr_pcnt_lo,tmp$mdr_pcnt,tmp$mdr_pcnt_hi)
YD[['Y.ret.dh1r.id']] <- as.array(tmp[,cid])
YD[['Y.ret.dh1r.M']] <- cbind(MS.hnr$mu,MS.mdr$mu)
YD[['Y.ret.dh1r.S']] <- cbind(MS.hnr$sig,MS.mdr$sig)
YD[['Y.ret.dh1r.t']] <- as.array(tmp$dy)

## -- dh1h
## new
css <- unique(grep('dh1.h',YB$case,value=TRUE))
(tmp <- YB[patients=='new' & case %in% css,
           .(cid,dy,case,
             dr_h_nr_pcnt,dr_r_nh_pcnt,mdr_pcnt,
             dr_h_nr_pcnt_lo,dr_r_nh_pcnt_lo,mdr_pcnt_lo,
             dr_h_nr_pcnt_hi,dr_r_nh_pcnt_hi,mdr_pcnt_hi)])
MS.rnh <- getmusig(tmp$dr_r_nh_pcnt_lo,tmp$dr_r_nh_pcnt,
                   tmp$dr_r_nh_pcnt_hi)
MS.mdr <- getmusig(tmp$mdr_pcnt_lo,tmp$mdr_pcnt,tmp$mdr_pcnt_hi)
YD[['Y.new.dh1h.id']] <- as.array(tmp[,cid])
YD[['Y.new.dh1h.M']] <- cbind(MS.rnh$mu,MS.mdr$mu)
YD[['Y.new.dh1h.S']] <- cbind(MS.rnh$sig,MS.mdr$sig)
YD[['Y.new.dh1h.t']] <- as.array(tmp$dy)

## ret
(tmp <- YB[patients=='ret' & case %in% css,
           .(cid,dy,case,
             dr_h_nr_pcnt,dr_r_nh_pcnt,mdr_pcnt,
             dr_h_nr_pcnt_lo,dr_r_nh_pcnt_lo,mdr_pcnt_lo,
             dr_h_nr_pcnt_hi,dr_r_nh_pcnt_hi,mdr_pcnt_hi)])
MS.rnh <- getmusig(tmp$dr_r_nh_pcnt_lo,tmp$dr_r_nh_pcnt,
                   tmp$dr_r_nh_pcnt_hi)
MS.mdr <- getmusig(tmp$mdr_pcnt_lo,tmp$mdr_pcnt,tmp$mdr_pcnt_hi)
YD[['Y.ret.dh1h.id']] <- as.array(tmp[,cid])
YD[['Y.ret.dh1h.M']] <- cbind(MS.rnh$mu,MS.mdr$mu)
YD[['Y.ret.dh1h.S']] <- cbind(MS.rnh$sig,MS.mdr$sig)
YD[['Y.ret.dh1h.t']] <- as.array(tmp$dy)

## -- dh2
## new
css <- unique(grep('dh2',YB$case,value=TRUE))
(tmp <- YB[patients=='new' & case %in% css,
           .(cid,dy,case,
             dr_h_nr_pcnt,dr_r_nh_pcnt,mdr_pcnt,
             dr_h_nr_pcnt_lo,dr_r_nh_pcnt_lo,mdr_pcnt_lo,
             dr_h_nr_pcnt_hi,dr_r_nh_pcnt_hi,mdr_pcnt_hi)])
MS.mdr <- getmusig(tmp$mdr_pcnt_lo,tmp$mdr_pcnt,tmp$mdr_pcnt_hi)
YD[['Y.new.dh2.id']] <- as.array(tmp[,cid])
YD[['Y.new.dh2.M']] <- MS.mdr$mu
YD[['Y.new.dh2.S']] <- MS.mdr$sig
YD[['Y.new.dh2.t']] <- as.array(tmp$dy)

## ret
(tmp <- YB[patients=='ret' & case %in% css,
           .(cid,dy,case,
             dr_h_nr_pcnt,dr_r_nh_pcnt,mdr_pcnt,
             dr_h_nr_pcnt_lo,dr_r_nh_pcnt_lo,mdr_pcnt_lo,
             dr_h_nr_pcnt_hi,dr_r_nh_pcnt_hi,mdr_pcnt_hi)])
MS.mdr <- getmusig(tmp$mdr_pcnt_lo,tmp$mdr_pcnt,tmp$mdr_pcnt_hi)
YD[['Y.ret.dh2.id']] <- as.array(tmp[,cid])
YD[['Y.ret.dh2.M']] <- MS.mdr$mu
YD[['Y.ret.dh2.S']] <- MS.mdr$sig
YD[['Y.ret.dh2.t']] <- as.array(tmp$dy)

## === join & save out
D4S <- c(LD,YD)
names(D4S) <- gsub("\\.","_",names(D4S))

## size of cases
idnmz <- grep('id',names(D4S),value=TRUE)
for(nm in idnmz){ #size of each case
  nm2 <- gsub('_id','',nm)
  nm2 <- paste0('N_',nm2)
  D4S[[nm2]] <- length(D4S[[nm]])
}

## add N, J, P, T
X <- regmm
D4S$X <- X
D4S$P <- ncol(X)
D4S$N <- nrow(X)
D4S$W <- W
D4S$Nedges <- sum(W)/2
tz <- sort(unique(c(YB$dy,LB$dy)))
D4S$T <- length(tz)

## this is to be altered depending on the model type
D4S$J <- 4 #all types


## NOTE issues corrected here:
summary(D4S$L_ret_dh0_N)
bad <- which(D4S$L_ret_dh0_N[,1]<0)
D4S$L_ret_dh0_N[bad,]
D4S$L_ret_dh0_N[bad,1] <- 0
isoidx[bad] #ISR - small numbers TODO check/query

## add number of RR data-points
NDP <- with(data=D4S,{
  ## new
  N_L_new_rr +
    N_L_new_x +
    N_L_new_dh0 +
    N_Y_new_rr +
    N_Y_new_x +
    N_Y_new_dh0 +
    N_Y_new_dh1h +
    ## ret
    N_L_ret_rr +
    N_L_ret_x +
    N_L_ret_dh0 +
    N_Y_ret_rr +
    N_Y_ret_x +
    N_Y_ret_dh0 +
    N_Y_ret_dh1h
})

D4S$NDP <- NDP

save(D4S,file=gh('drtb/data/D4S.Rdata'))


## also output table
names(D4S)
nmz <- grep("N_",names(D4S),value=TRUE)
rnmz <- gsub("N_","",nmz)
inmz <- paste0(rnmz,"_id")

cc <- cy <- rep(0,length(nmz))
for(i in 1:length(nmz)){
  cy[i] <- D4S[[nmz[i]]]
  cc[i] <- length(unique(D4S[[inmz[i]]]))
}

TD <- data.table(rcase=nmz,cc=cc,cy=cy)
TD[,c('ditch','data.type','patient.group','case'):=tstrsplit(rcase,"_")]
TD[,ditch:=NULL]
TD[,patient.group:=ifelse(patient.group=='new','new','retreatment')]
TD[,data.type:=ifelse(data.type=='L','surveillance','survey')]
TD <- TD[order(patient.group,data.type,case)]
TD <- TD[,.(data.type,patient.group,case,country.count=cc,country.years=cy)]

fwrite(TD,file=gh('drtb/data/TD.csv'))

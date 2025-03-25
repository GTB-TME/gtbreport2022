#' ---
#' title: Generating age/sex incidence disaggregation
#' author: Pete Dodd
#' date: 18 July, 2022
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---

#' #Pre-amble
#' (Last updated: `r Sys.Date()`)
#'
#' This file is for generating age/sex disaggregation of incidence.
#'
#' This is therefore to be run after the analysis of prevalence data and after the child TB data has been generated
#'
#' N.B. This file should render to html with `rmarkdown::render('03inc_split.R',output_dir='../html')` or from the command line with `R -q -e "rmarkdown::render(\"03inc_split.R\",output_dir=\"../html\")"`
#'
#'
#' #Read in necessary data
#'
#' Relevant libraries for loading and manipulating data:
#'
rm(list=ls())
library(here)                           # file locations
source(here('disaggregation/R/utilities.R'))

#' Various pre-requisite data
load(here('data/est.rda'))         # current burden estimates
load(here('data/tb.rda'))          # includes notifications
load(here('data/pop.rda'))          # population for denom
load(here('disaggregation/R/02children/data/K.Rdata')) # Child prior
load(here('disaggregation/output/prevsplits/prevsplits.Rdata')) # from prevalence split work
names(prevsplits)[2] <- 'age_group'

#' Flag for whether cairo device is functioning for outward plotting:
gypt <- TRUE
plotting <- TRUE                       #whether to plot ot not
set.seed(2345)                          #PRNG seed

## a country hash
hbcsh <- merge(unique(est[,.(iso3,g.hbc)]),
               unique(tb[,.(iso3,name=country)]))
hbcsh <- hbcsh[g.hbc==TRUE]; hbcsh[,g.hbc:=NULL]
## shorter names for graphs
hbcsh[iso3=='COD',name:="DR Congo"]
hbcsh[iso3=='PRK',name:="DPR Korea"]
hbcsh[iso3=='TZA',name:="UR Tanzania"]


#' Reading BB for kid splits prior
fn <- here('disaggregation/R/02children/data/BB.Rdata')
load(fn)

#' Having a look at using the model for kids prior. First needs some NAs filled
K <- merge(K,unique(tb[,.(iso3,g.whoregion)]),
           by='iso3',all.y = TRUE,all.x = TRUE)
K[,rav1:=mean(propu15b,na.rm = TRUE),by=g.whoregion]
K[,rav2:=mean(propu15b.sd,na.rm = TRUE),by=g.whoregion]
K[is.na(propu15b),propu15b:=rav1];
K[is.na(propu15b.sd),propu15b.sd:=rav2]
K[propu15b==0,propu15b.sd:=rav2]; K[propu15b==0,propu15b:=rav1]; 
K[,propu15:=propu15b]; K[,propu15.sd:=propu15b.sd]; #progression model!


#' First, for countries with `NA` in prevsplit (from method fails, especially in AMR)
#' create a global median of proportions, normalised:
names(prevsplits)[names(prevsplits)=='age.group'] <- 'age_group'
(missed <- prevsplits[is.na(prop),unique(iso3)])
glop <- prevsplits[,.(propm=median(prop,na.rm=TRUE), #global median
                      propm.se=median(prop.se,na.rm=TRUE)),
                   by=.(sex,age_group)]
tot <- glop[,sum(propm)]; glop[,propm:=propm/tot] #renormalize
prevsplits <- merge(prevsplits,glop,by=c('sex','age_group'),
                     all.x=TRUE,all.y=TRUE,allow.cartesian = TRUE)
prevsplits[iso3 %in% missed,
           c('prop','prop.se'):=list(propm,propm.se)] #fill out with global version
(missed <- prevsplits[is.na(prop.se),unique(iso3)])
prevsplits[iso3 %in% missed,prop.se:=propm.se] #fill out with global version
prevsplits[,c('propm','propm.se'):=NULL]

prevsplits$age_group <- factor(prevsplits$age_group,
                               levels=agz2,ordered=TRUE)
prevsplits$sex <- factor(prevsplits$sex,levels=c('M','F'),ordered=TRUE)
setkeyv(prevsplits,c('iso3','age_group','sex'))

#'
#' # Child sex ratios
#' 
#' Looking at sex ratios in the child notification data
#' Make data:
clz <- c('iso3','g.whoregion',
         outer(c('newrel.m','newrel.f'),kdzboth,paste0))
tmpw <- tb[year==estyr,..clz]
tmpw[,N:=newrel.m04 + newrel.f04 + newrel.m514 + newrel.f514 ]
tmpw[,table(is.na(N),is.na(newrel.m59))] #no rescue from new age gps
tmpw <- tmpw[!is.na(N)]
tmpw <- tmpw[N>0]


#' Use metafor and conduct some meta-analyses for younger and older age groups
#' 
#' {r fig.height = 10, fig.width = 5}
## ======== 0-4 ===== proportion M
dat <- tmpw[,.(id=iso3,xi=newrel.m04,
               ni=newrel.m04 + newrel.f04,g.whoregion)]
dat <- dat[dat$ni>0,]
dat <- dat[dat$ni>dat$xi,]
dat <- escalc(measure="PR", xi=xi, ni=ni, data=dat)
citmp <- t(apply(dat[,c('xi','ni')],1,function(x) binom.test(x[1], x[2])$conf.int))
dat$ci.lb <- citmp[,1]
dat$ci.ub <- citmp[,2]
res <- rma.glmm(measure="PLO", xi=xi, ni=ni, data=dat)
res.rma <- rma.glmm(measure="PLO", mods=~g.whoregion,
                    xi=xi, ni=ni, data=dat) #regional
## ## plot
## with(dat, forest(yi, ci.lb=ci.lb, ci.ub=ci.ub, refline=predict(res, transf=transf.ilogit)$pred,slab = id))
## addpoly(res, row=-1, transf=transf.ilogit)
## title('Proportion male 0-4')
## abline(h=0)
(YS <- predict(res, transf=transf.ilogit, digits=3))
## by region
YSR <- predict(res.rma,transf=transf.ilogit, digits=3)
YSR <- data.table(mid=YSR$pred,hi=YSR$cr.ub,lo=YSR$cr.lb,
                  g.whoregion=dat$g.whoregion)
(YSR <- YSR[,.(mid=mid[1],hi=hi[1],lo=lo[1]),by=g.whoregion])

fwrite(YSR,file=here('disaggregation/R/02children/data/YSR.csv'))


#' {r fig.height = 10, fig.width = 5}
## ======== 5-14 ===== proportion M
dat <- tmpw[,.(id=iso3,xi=newrel.m514,ni=newrel.m514 + newrel.f514,
               g.whoregion)]
dat <- dat[dat$ni>0,]
dat <- dat[dat$ni>dat$xi,]
dat <- escalc(measure="PR", xi=xi, ni=ni, data=dat)
citmp <- t(apply(dat[,c('xi','ni')],1,function(x) binom.test(x[1], x[2])$conf.int))
dat$ci.lb <- citmp[,1]
dat$ci.ub <- citmp[,2]
res <- rma.glmm(measure="PLO", xi=xi, ni=ni, data=dat)
res.rma <- rma.glmm(measure="PLO", mods=~g.whoregion,
                    xi=xi, ni=ni, data=dat) #regional
## ## plot
## with(dat, forest(yi, ci.lb=ci.lb, ci.ub=ci.ub, refline=predict(res, transf=transf.ilogit)$pred,slab = id))
## addpoly(res, row=-1, transf=transf.ilogit)
## title('Proportion male 5-14')
## abline(h=0)
(OS <- predict(res, transf=transf.ilogit, digits=3))
## by region
OSR <- predict(res.rma,transf=transf.ilogit, digits=3)
## OSR <- predict(res,transf=transf.ilogit, digits=3)
OSR <- data.table(mid=OSR$pred,hi=OSR$cr.ub,lo=OSR$cr.lb,
                  g.whoregion=dat$g.whoregion)
(OSR <- OSR[,.(mid=mid[1],hi=hi[1],lo=lo[1]),by=g.whoregion])

fwrite(OSR,file=here('disaggregation/R/02children/data/OSR.csv'))


#'
#' # Ad hoc changes to remove NAs where possible
#'
#' 

#'  ## First look at 0-14 groupings and NA as zero problems
#'
#' Boys:
see <- c('iso3','newrel.m014','newrel.m04',
         'newrel.m1524','newrel.m2534')
see2 <- c('iso3','newrel.m014','newrel.m04',
          'newrel.m514','newrel.m1524','newrel.m2534')
tb[year==estyr & !is.na(newrel.m014) & is.na(newrel.m04),..see]

#' countries with oddities
tb[year==estyr-1 & !is.na(newrel.m014) & is.na(newrel.m04),..see] #using 014
tb[year==estyr-1 & is.na(newrel.m014) & is.na(newrel.m04),..see]  #all NA


#' Record out:
(knotem21 <- tb[year==estyr & !is.na(newrel.m014) & is.na(newrel.m04),
                unique(iso3)])
(knotem20 <- tb[year==estyr-1 & !is.na(newrel.m014) & is.na(newrel.m04),
                unique(iso3)])

cat(knotem21,file=here('disaggregation/output/incsplits/data/knotem21.txt'))
cat(knotem20,file=here('disaggregation/output/incsplits/data/knotem20.txt'))

#' For countries where the 014 count is zero, fill out the separate categories
tb[year==estyr & !is.na(newrel.m014) &
   is.na(newrel.m04) & newrel.m014==0,..see]
tb[year==estyr & !is.na(newrel.m014) &
   is.na(newrel.m04) & newrel.m014==0,
   c('newrel.m04','newrel.m514'):=.(0,0)]

#' For countries where there is a total count, assign to categories by regional weight
tmp <- tb[year==estyr & !is.na(newrel.m04) & !is.na(newrel.m514),
          .(g.whoregion,newrel.m04,newrel.m514)]

ggplot(data=tmp[newrel.m04>0,
                .(frac=newrel.m04/(newrel.m04+newrel.m514)),
                by=g.whoregion],aes(x=frac)) +
  geom_histogram()+
  facet_wrap(~g.whoregion) +
    geom_vline(data=tmp[newrel.m04>0,
                        .(frac=mean(newrel.m04 /
                                    (newrel.m04+newrel.m514))),
                        by=g.whoregion],
               aes(xintercept=frac),col=2) +
    xlab('Fraction of boy notifications <5')

tmp <- tmp[newrel.m04>0,.(fracm=mean(newrel.m04 /
                                     (newrel.m04+newrel.m514))),
           by=g.whoregion]

fwrite(tmp,file=here('disaggregation/output/incsplits/data/knotem.csv'))

tb[year==estyr & !is.na(newrel.m014) & is.na(newrel.m04),..see]
work <- tb[year==estyr & !is.na(newrel.m014) & is.na(newrel.m04),iso3]

tb <- merge(tb,tmp,by='g.whoregion')

tb[year==estyr & !is.na(newrel.m014) & is.na(newrel.m04),
   newrel.m514 := floor((1-fracm) * newrel.m014)]
tb[year==estyr & !is.na(newrel.m014) & is.na(newrel.m04),
   newrel.m04 := newrel.m014-newrel.m514]

#' Check
tb[year==estyr & !is.na(newrel.m014) & is.na(newrel.m04),..see]
kable(tb[year==estyr & iso3 %in% work,..see2])


#' Girls:
see <- c('iso3','newrel.f014','newrel.f04',
         'newrel.f1524','newrel.f2534')
see2 <- c('iso3','newrel.f014','newrel.f04',
          'newrel.f514','newrel.f1524','newrel.f2534')
tb[year==estyr & !is.na(newrel.f014) & is.na(newrel.f04),..see]

#' Oddities
tb[year==estyr-1 & !is.na(newrel.f014) & is.na(newrel.f04),..see] #using 014
tb[year==estyr-1 & is.na(newrel.f014) & is.na(newrel.f04),..see]  #all NA

#' Record out:
(knotef21 <- tb[year==estyr & !is.na(newrel.f014) & is.na(newrel.f04),unique(iso3)])
(knotef20 <- tb[year==estyr-1 & !is.na(newrel.f014) & is.na(newrel.f04),unique(iso3)])

cat(knotef21,file=here('disaggregation/output/incsplits/data/knotef21.txt'))
cat(knotef20,file=here('disaggregation/output/incsplits/data/knotef20.txt'))

#' For countries where the 014 count is zero, fill out the separate categories
tb[year==estyr & !is.na(newrel.f014) & is.na(newrel.f04) &
   newrel.f014==0,..see]
tb[year==estyr & !is.na(newrel.f014) & is.na(newrel.f04) &
   newrel.f014==0,
   c('newrel.f04','newrel.f514'):=.(0,0)]

#' For countries where there is a total count, assign to categories by regional weight
tmp <- tb[year==estyr & !is.na(newrel.f04) & !is.na(newrel.f514),
          .(g.whoregion,newrel.f04,newrel.f514)]

ggplot(data=tmp[newrel.f04>0,.(frac=newrel.f04 /
                                   (newrel.f04+newrel.f514)),
                by=g.whoregion],
       aes(x=frac)) +
  geom_histogram()+
  facet_wrap(~g.whoregion) +
    geom_vline(data=tmp[newrel.f04>0,
                        .(frac=mean(newrel.f04 /
                                    (newrel.f04+newrel.f514))),
                        by=g.whoregion],
               aes(xintercept=frac),col=2) +
    xlab('Fraction girl notifications <5')

tmp <- tmp[newrel.f04>0,.(fracf=mean(newrel.f04 /
                                     (newrel.f04+newrel.f514))),
           by=g.whoregion]

fwrite(tmp,file=here('disaggregation/output/incsplits/data/knotef.csv'))

tb[year==estyr & !is.na(newrel.f014) & is.na(newrel.f04),..see]
work <- tb[year==estyr & !is.na(newrel.f014) & is.na(newrel.f04),iso3]

tb <- merge(tb,tmp,by='g.whoregion')

tb[year==estyr & !is.na(newrel.f014) & is.na(newrel.f04),
   newrel.f514 := floor((1-fracf) * newrel.f014)]
tb[year==estyr & !is.na(newrel.f014) & is.na(newrel.f04),
   newrel.f04 := newrel.f014-newrel.f514]

#' Check
tb[year==estyr & !is.na(newrel.f014) & is.na(newrel.f04),..see]
kable(tb[year==estyr & iso3 %in% work,..see2])


#'
#' # Other NAs that are fixable
#' 
#' Inspect:
am <- paste0('newrel.m',agz); af <- paste0('newrel.f',agz);
tbnow <- tb[year==estyr]
gotna <- rowSums(tbnow[,..am]) + rowSums(tbnow[,..af])
gotnas <- rowSums(is.na(tbnow[,..am])) + rowSums(is.na(tbnow[,..af]))
tots <- rowSums(tbnow[,..am],na.rm=TRUE) +
    rowSums(tbnow[,..af],na.rm=TRUE) #totals
tbnow[,tots:=tots]
toshowm <- c('iso3','tots',am); toshowf <- c('iso3','tots',af)
tbnow[is.na(gotna) & gotnas<16,..toshowm] 
tbnow[is.na(gotna) & gotnas<16,..toshowf]

#' Assume that if total is <100, NA is actually zero. Ie for
(naz <- tbnow[is.na(gotna) & gotnas<16 & tots<100,iso3])

cat(naz,file=here('disaggregation/output/incsplits/data/naz.txt'))

## loop through
ind <- c(am,af)
for(j in ind){
  set(tb, i = which(tb$year==estyr & tb$iso3 %in% naz & is.na(tb[[j]])),
      j = j, value = 0)
}

#' Check
#' 
toshowm <- toshowm[-2]
toshowf <- toshowf[-2]
tb[year==estyr & iso3 %in% naz,..toshowm]
tb[year==estyr & iso3 %in% naz,..toshowf]

#' Keep this version of `tb.Rdata` in case needed for CFR mortality
save(tb,file=here('disaggregation/output/incsplits/data/tbimp.Rdata'))

#' 
#' # Simplification of flow
#'
#' ## Pass1: countries using standard adjustments
#'
#' Make a list of countries to be done
isotodo <- est[,unique(iso3)]
length(isotodo)

#' Use standard adjustment for standard adjustment or CR with good CDR?
#'
#' CDRs of countries
cdrtab <- est[,.(iso3,year,inc.num=inc*pop)]
##newinc includes unk h
cdrtab <- merge(cdrtab,tb[,.(iso3,year,c.newunk,c.newinc)],
                by=c('iso3','year'))
cdrtab <- cdrtab[!is.na(c.newunk)]
cdrtab <- cdrtab[,.SD[year==max(year)],by=iso3]
dim(cdrtab)
cdrtab[,cdr:=c.newinc/inc.num*1e5]
cdrtab[,qplot(cdr)] + geom_vline(xintercept = 0.85,col=2)
cdrtab[!is.finite(cdr)]                      #
cdrtab[!is.finite(cdr),cdr:=1]
cdrtab[cdr>1]
setkey(cdrtab,iso3)

#' Lists of countries todo
est[year==estyr,unique(source.inc)]
est[year==estyr,table(source.inc)]
est[year==estyr & source.inc=="Model",iso3] #COVID disruption, directly
est[year==estyr & source.inc=="Regional model",iso3] #COVID disruption, more crudely
est[year==estyr & source.inc=="Current trends",iso3] #

isosf <- est[grepl("Standard adjustment",source.inc),
             as.character(unique(iso3))]
cdrtab[isosf]


#' Delete any existing plots in relevant folder!
#' 
fz <- list.files(path=here('disaggregation/output/incsplits/plots'),
                 pattern="*.pdf",
                 recursive = TRUE, full.names = TRUE)
if(length(fz)>0) do.call(file.remove,list(fz)) #clean up

#' _Pass 1_ will be countries that are either:
#' 
#'  * SF adjustment & $CDR>0.85$
#'  * fewer than 1000 notifications in total
#' 
isopass1 <- c(cdrtab[isosf][cdr>0.85,iso3], #SF + CDR>0.85
              cdrtab[c.newinc<1e3,iso3]     #few notifications
              )
isopass1 <- unique(isopass1)
isopass1

cat(isopass1,file=here('disaggregation/output/incsplits/data/pass1.txt'))

#' Create directory(s) if missing
drnm <- here('disaggregation/output/incsplits/plots')
if(!file.exists(drnm)) dir.create(drnm)
drnm <- here('disaggregation/output/incsplits/plots/pass1')
if(!file.exists(drnm)) dir.create(drnm)

#' Work up incsplits for pass 1 countries
H <- W <- 7
ISL <- hi30 <- pltlist <- noteissue <- naissue <- list()
## starting loop
for(cn in isopass1){
  tmp <- getNotes(cn)                   #most recent notifications
  if(tmp[,sum(is.na(newrel))] || tmp[,sum(newrel,na.rm=TRUE)]<50){
    tmp[is.na(newrel),newrel:=0] #correcting GRD, eg
    noteissue[[cn]] <- cn
  }
  cdr <- cdrtab[cn,cdr]
  tot.newrel <- tmp[,sum(newrel)]
  if(tot.newrel==0){
    tot.newrel <- 1
    cdr <- 1
  }
  tmp[,prop:=newrel/tot.newrel]
  tmp[,inc:=newrel/cdr]
  tmp[,c('iso3','prop.se','src','pass','problem'):=
         .(cn,NA,'Notifications x factor',1,FALSE)]
  tmp <- merge(tmp,AA[,.(age_group,age)],by='age_group')
  ISL[[cn]] <- tmp
  plt <- niplot(tmp)
  pltlist[[cn]] <- plt
  if(est[iso3==cn,g.hbc[1]]) hi30[[cn]] <- plt #add to hibc30
  if(plotting){
    fn <- drnm + glue('/{cn}1.pdf')
    if(gypt){
      ggsave(filename=fn,plot=plt,device = cairo_pdf,h=H,w=W)
    } else
      ggsave(filename=fn,plot=plt,h=H,w=W)
  }
}



#' Following countries with fewer than 50 notifications have some issues:
(noteissue <- unlist(noteissue))
cat(noteissue,file=here('disaggregation/output/incsplits/data/noteissue.txt'))

#' Still left todo
isotodo <- setdiff(isotodo,isopass1) #
length(isotodo)

#' combine & check
incsplit <- rbindlist(ISL)
incsplit[!is.finite(inc)] #DMA has NAs
incsplit[!is.finite(inc),c('inc','prop.se'):=0.0] #remove

incsplit[,sum(inc)]*1e-6
incsplit[age_group %in% kds,sum(inc)]/incsplit[,sum(inc)]*1e2
incsplit[,length(unique(iso3))]

#'
#' ## Pass 2: sampling against prior
#'
#' 
#' ## Function for country proportions
#'
#' Basically, this function uses the adult prevalence splits worked in the prevalence analysis and applies these after the child incidence has been split off. However, sometimes the notifications exceed incidence in some age group. This function will model the output of the prevalence analysis as a Dirichlet distribution, and sample from this in order to arrive at some samples that do not fail this test (i.e. rejection sampleing from this as a prior).
#'
#' Now the main function. Setting `asis=TRUE` will by-pass any attempt to do rejection sampling and return the (potentially faulty) incidence split. This will be returned in any case if it is not possible to fix, either because there are insufficient incident adult cases to exceed notifications (across all adult ages), or because there were no samples not rejected:

priorsample <- function(cn,
                        nsim=1e6,
                        silent=FALSE,
                        tol=1e-2,
                        adultadjust=FALSE
                        ){
  tmp <- getNotes(cn,silent=silent)                 #Notifications
  ages3 <- ifelse(nrow(tmp)==6,TRUE,FALSE)          #no adult split in notes
  if(tmp$year[1] < max(tb$year)){# if not most recent year then rescale
    EL <- est[cn][,.(iso3,year,inc,pop)]
    sf <- EL[year==max(tb$year),inc*pop]/EL[year==tmp$year[1],inc*pop]
    if(!is.finite(sf)) sf <- 1
    tmp$newrel <- tmp$newrel * sf
  }
  ## fetch adult child notifications, total incidence
  INC <- est[cn][year==estyr,inc*pop/1e5] #total incidence TODO check OK
  atmp <- tmp[age_group %ni% kds,]    #adult notifications
  tmp2 <- K[cn]
  ## --- make splits ----
  if(nsim>1){ ## replace with samples
    ## child splits
    lnp <- getLNparms(tmp2[,propu15],tmp2[,propu15.sd]^2)
    propu15 <- rlnorm(n=nsim,meanlog=lnp[1,1],sdlog = sqrt(lnp[2,1])) #adult/child
    if( cn %ni% BB$iso3 ){
      txt <- paste0(cn,' not found in BB.Rdata! You probably needed to delete this first.\n')
      stop(txt)
    }
    aa <- BB[iso3==cn,aa]; ab <- BB[iso3==cn,ab]
    suppressWarnings({propu5 <- rbeta(nsim,aa,ab)})#child age split
    tt <- getLNparms(YS$pred,(YS$ci.lb-YS$ci.ub)^2/3.92^2)
    pu5m <- rlnorm(nsim,tt[1,1],sqrt(tt[2,1])) #young child sex split
    tt <- getLNparms(OS$pred,(OS$ci.lb-OS$ci.ub)^2/3.92^2)
    po5m <- rlnorm(nsim,tt[1,1],sqrt(tt[2,1])) #old child sex split
    ## adult splits
    if(!adultadjust){
      lnp <- getLNparms(prevsplits[cn][,prop],prevsplits[cn][,prop.se]^2)
      adultprops <- apply(t(lnp),1,
                          function(x)rlnorm(n=nsim,meanlog=x[1],sdlog = sqrt(x[2])))
    } else {
      adultprops <- matrix(tmp$newrel[-c(1:4)],nrow=nsim,ncol=12,byrow=TRUE)
    }
    adultprops <- adultprops / rowSums(adultprops) #normalize
    ## make matrix nsim x age/sex of proportions
    PZ <- cbind(
      propu15 * propu5 * pu5m,   #M<5
      propu15 * propu5 * (1-pu5m), #F<5
      propu15 * (1-propu5) * po5m, #M514
      propu15 * (1-propu5) * (1-po5m), #F514
      (1-propu15) * adultprops
                )
    PZ <- PZ/rowSums(PZ)                #safety
  } else {stop('nsim must be >1!');}
  ## --- compare to notifications ---
  if(!ages3){
    NZ <- matrix(tmp$newrel,nrow=nsim,ncol=16,byrow=TRUE)
    NZ <- INC*PZ - NZ                     #hoping >0
  } else {
    NZ <- matrix(tmp$newrel,nrow=nsim,ncol=6,byrow=TRUE)
    PZS <- cbind( PZ[,1:4] ,
                 rowSums(PZ[,seq(5,15,by=2)]),
                 rowSums(PZ[,seq(5,15,by=2)+1]))
    NZ <- INC*PZS - NZ                     #hoping >0
  }
  keep <- apply(NZ,1,function(x) all(x>0,na.rm = TRUE)) #NA for cases like MOZ
  accept <- TRUE
  if(all(!keep)){                       #all failed
    if(!silent)cat('None accepted, using tolerance!\n')
    accept <- FALSE
    NZ[is.na(NZ)] <- 0
    NZ[NZ>0] <- 0                       #only penalize undershoot
    penalty <- rowSums(NZ^2)
    threshold <- quantile(penalty,tol) #best 1%
    keep <- penalty < threshold
    if(!silent)cat('Using best ',sum(keep),' samples\n')
  } else{ if(!silent)cat('Accepted ',sum(keep),' samples\n');}
  if(sum(keep)>0){                      #safety
    if(sum(keep)>1){
      PZ <- PZ[keep,]
      PZ <- colMeans(PZ)
    } else PZ <- PZ[keep,]
    PZ <- PZ/sum(PZ)
    if(ages3){
      tmp <- ni3cat(tmp)
      tmp[,note:='1 adult newrel acat!']
    } else {
      tmp[,note:='']
    }
    tmp[,prop:=PZ]
    PZ <- INC * PZ
    tmp[,inc:=PZ]
    tmp[,prop.se:=NA]
    if(accept) tmp[,src:='Rejection sampling from prior']
    if(!accept) tmp[,src:='Best samples from prior']
  } else {                              #error!
    tmp[,c('prop','inc','prop.se','src'):=NA]
  }
  tmp
}

#' test!
tmp <- priorsample('NGA') #
niplot(tmp)

#' problematic countries looked at
tmp <- priorsample('MOZ')
niplot(tmp)

#' look at 3 age plot
tmp[age_group%ni% kds,age_group:='15plus']
tmp3 <- tmp[,.(newrel=sum(newrel),inc=sum(inc)),
            by=.(age_group,sex,iso3,year)]
niplot(tmp3)

#' another example with only one adult age category
tmp <- priorsample('GMB')
niplot(tmp)

#' look at 3 age
tmp[age_group%ni% kds,age_group:='15plus']
tmp3 <- tmp[,.(newrel=sum(newrel),inc=sum(inc)),
            by=.(age_group,sex,iso3,year)]
niplot(tmp3)

#' Create directory if missing
drnm <- here('disaggregation/output/incsplits/plots/pass2')
if(!file.exists(drnm)) dir.create(drnm)


#' Looping through countries
#'
#' 
isotodo <- isotodo[!isotodo %in% c('SCG')] #premerge
for(cn in isotodo){
  cno <- priorsample(cn)
  cno$pass <- 2
  ISL[[cn]] <- copy(cno)
  if(cno$note[1]!='') cno[age_group %ni% kds,newrel:=NA] #1 adult age
  plt <- niplot(cno)
  pltlist[[cn]] <- plt
  if(est[iso3==cn,g.hbc[1]]) hi30[[cn]] <- plt #add to hibc30
  if(plotting){
    fn <- drnm + glue('/{cn}2.pdf')
    if(gypt){
      ggsave(filename=fn,plot=plt,device = cairo_pdf,w=W,h=H)
    } else
      ggsave(filename=fn,plot=plt,w=W,h=H)
  }
}

#' validity issues:
#' BGD, CHN, IND?, MAR, MMR?, PNG, PRK, TGO?, ZMB

#' 3 category case:
#' MOZ

#' Check UGA 
## UGA, 
(tmp <- getNotes('UGA',silent=FALSE))                 #Notifications
#' Check this and consider
mn <- paste0('newrel.m',agz)
fn <- paste0('newrel.f',agz)
bn <- c(fn,mn)
(v <- tb[iso3=='UGA' & year==2021][,..bn])
sum(v,na.rm=TRUE)
tb[iso3=='UGA' & year==2021,c.newinc]
#' No error


#' Save out the 1 adult age cat graphs
tmp <- copy(ISL[['MOZ']])
mozf <- niplot(tmp)
tmp[age_group%ni% kds,age_group:='15plus']
tmp3 <- tmp[,.(newrel=sum(newrel),inc=sum(inc)),
            by=.(age_group,sex,iso3,year)]

moz <- niplot(tmp3)
ggsave(moz,device = cairo_pdf,w=W,h=H,
       filename=here('disaggregation/output/incsplits/plots/pass2/MOZ2b.pdf'))

#' # Pass 3: Ad hoc reversion to adjustments
#'
#'  Oddities among the pass2 countries deserving special treatment:
#'
#'
#' 

## ## last year:
## isopass3 <- c(## 'BGR',#now in pass 1
##               'CHN','MAR','PNG','PRK','SYR','THA',
##               ## from last year (edgecases @ pass 2 this yr):
##               ## would also include UGA but missing newrel precludes
##               'MMR','IND','IDN','BGD')


## BGD, CHN, IND?, MAR,MMR? , PNG, PRK, TGO?, ZMB
isopass3 <- c('BGD',
              'CHN',
              'IND', #?
              'MAR',
              'MMR', #?
              'PAN',
              'PNG',
              'PRK',
              'TGO', #?
              'ZMB'
              )

#' Create directory if missing
drnm <- here('disaggregation/output/incsplits/plots/pass3')
if(!file.exists(drnm)) dir.create(drnm)

#' Sample following the adult pattern of notifications
for(cn in isopass3){
  cno <- priorsample(cn,adultadjust = TRUE)
  cno$pass <- 3
  ISL[[cn]] <- cno
  plt <- niplot(cno)
  pltlist[[cn]] <- plt
  if(est[iso3==cn,g.hbc[1]]) hi30[[cn]] <- plt #add to hibc30
  if(plotting){
    fn <- drnm + glue('/{cn}3.pdf')
    if(gypt){
      ggsave(filename=fn,plot=plt,device = cairo_pdf,w=W,h=H)
    } else
      ggsave(filename=fn,plot=plt,w=W,h=H)
  }
}

#' Create directory if missing
drnm <- here('disaggregation/output/incsplits/plots/pass2/0pass2fails')
if(!file.exists(drnm)) dir.create(drnm)


## move them into the failed directory
for( cn in isopass3){
  ofn <- here('disaggregation/output/incsplits/plots/pass2') + glue('/{cn}2.pdf')
  nfn <- here('disaggregation/output/incsplits/plots/pass2/0pass2fails') + glue('/{cn}2.pdf')
  file.rename(from=ofn,to=nfn)
}

#'
#'
#' # Finishing
#'
#' 

#' Create directory if missing
drnm <- here('disaggregation/output/incsplits/plots/data')
if(!file.exists(drnm)) dir.create(drnm)


#' Checks
length(hi30)
length(ISL)
length(pltlist)
save(file=here('disaggregation/output/incsplits/plots/data/pltlist.Rdata'),pltlist)
save(file=here('disaggregation/output/incsplits/plots/data/hi30.Rdata'),hi30)

## aggregations
incsplit <- rbindlist(ISL,fill=TRUE)
incsplit[!is.finite(inc)] #DMA has NAs
incsplit[!is.finite(inc),c('inc','prop.se'):=0.0] #remove

incsplit[,sum(inc,na.rm=TRUE)]*1e-6
incsplit[age_group %in% kds,sum(inc,na.rm=TRUE)]/incsplit[,sum(inc,na.rm=TRUE)]*1e2
incsplit[,length(unique(iso3))]

#' checks
## proportions
incsplit[,testsum:=sum(prop),by=iso3]
(tf <- unique(incsplit[abs(testsum-1)>1e-2,iso3]))
incsplit[iso3 %in% tf,.(iso3,newrel,prop,src)]
incsplit[iso3 %in% tf,all(newrel==0)]   #only failing as zero
setkey(incsplit,iso3)

## incidence
tmp <- incsplit[,.(testinc=sum(inc)),by=iso3]
tmp2 <- merge(est[year==estyr,.(iso3,inc)],
              pop[year==estyr,.(iso3,e.pop.num)],by='iso3')
tmp2[,inc.num:=inc*e.pop.num/1e5]
tmp <- merge(tmp,tmp2[,.(iso3,inc.num)],
             by='iso3',all.x = TRUE,all.y = FALSE)
incsplit <- merge(incsplit,tmp2[,.(iso3,inc.num)],
                  by='iso3',all.x = TRUE,all.y = FALSE)
(tf <- unique(tmp[abs(testinc/inc.num-1)>1e-2,iso3]))
incsplit[iso3 %in% tf,inc:=prop*inc.num,by=iso3]
incsplit[,inc.num:=NULL]
incsplit[,age:=NULL]


#' Save out data:
save(file=here('disaggregation/output/incsplits/data/incsplit.Rdata'),
     incsplit)



#' resume
#'
#' Load
load(file=here('disaggregation/output/incsplits/data/incsplit.Rdata'))
load(file=here('disaggregation/output/incsplits/plots/data/hi30.Rdata'))

#'
#' Plot of HBC 30 countries
## hbc <- as.character(est[g.hbc==TRUE,unique(iso3)])
hbc <- as.character(hbcsh[order(as.character(name)),iso3])
hbcn <- as.character(hbcsh[order(as.character(name)),name])
hbc <- c(t(matrix(hbc,byrow = TRUE,ncol=5)))         #re-order for plot
hbcn <- c(t(matrix(hbcn,byrow = TRUE,ncol=5)))         #re-order for plot

## pltlst <- list()                        #change order
## for(i in seq_along(hbc)) pltlst[[i]] <- hi30[[hbc[i]]]

#' Record which countries are not using most recent year
#' 

mryr <- list()
for(cn in tb[,unique(iso3)]){
  tmp <- getNotes(cn)[,.(iso3,year)]
  mryr[[cn]] <- tmp
}
mryr <- rbindlist(mryr)
mryr <- unique(mryr)
mryr[year==1980,year:=NA]
mryr[year!=estyr]
fwrite(mryr,file=here('disaggregation/output/incsplits/data/mryr.csv'))

#' find out if any HBC are old notification data
(hbcna <- mryr[iso3 %in% hbc & year<estyr,iso3])

#' and correct their plots
if(length(hbcna)>0){
  cat(hbcna,file=here('disaggregation/output/incsplits/data/hbcna.txt'))
  for(cn in hbcna){
    tmp <- incsplit[cn]
    tmp[,newrel:=NA]                        #remove notification data
    hi30[[cn]] <- niplot(tmp)
  }
}


## 30 HBC
nplst <- list()
for(i in 1:30){
  plt <- hi30[[hbc[i]]]
  plt <- plt + theme_classic()
  plt <- plt + theme(legend.position="none") + xlab('') + ylab('')
  if(!(i%%5==1)) plt <- plt + theme(axis.text.y=element_blank())
  ## if(i>6) plt <- plt + theme(axis.text.y=element_blank())
  plt <- plt + ggtitle(hbcn[i]) + theme(axis.text.x=element_text(size=7))
  nplst[[i]] <- plt
}


hi30[['MOZ']]
nplst

#' Save out plot
#' 
fn <- here('disaggregation/output/incsplits/aplots/HI30.pdf')
plt <- ggarrange(plotlist = nplst,ncol=5,nrow=6)
## ggexport(filename=fn,plot=plt,height=15,width=15*.75,device = cairo_pdf)
if(gypt){
  ggsave(filename=fn,plot=plt,height=15,width=15*.75,device = cairo_pdf)
} else
    ggsave(filename=fn,plot=plt,height=15,width=15*.75)

#' Width changes

wz <- rep(1,5);wz[1] <- 1.15; wz <- wz/1.15
fn <- here('disaggregation/output/incsplits/aplots/HI30b.pdf')
plt <- grid.arrange(nplst[[1]],nplst[[2]],nplst[[3]],
                    nplst[[4]],nplst[[5]],
                    nplst[[6]],nplst[[7]],nplst[[8]],nplst[[9]],
                    nplst[[10]],
                    nplst[[11]],nplst[[12]],nplst[[13]],nplst[[14]],
                    nplst[[15]],
                    nplst[[16]],nplst[[17]],nplst[[18]],nplst[[19]],
                    nplst[[20]],
                    nplst[[21]],nplst[[22]],nplst[[23]],nplst[[24]],
                    nplst[[25]],
                    nplst[[26]],nplst[[27]],nplst[[28]],nplst[[29]],
                    nplst[[30]],
                    ncol=5,widths=wz)
if(gypt){
  ggsave(filename=fn,plot=plt,,height=15,width=15*.75,device = cairo_pdf)
} else
  ggsave(filename=fn,plot=plt,,height=15,width=15*.75)


#' Regional and global splits

regsplt <- merge(incsplit,est[year==estyr,.(iso3,g.whoregion)],
                 by='iso3',
                 all.y = FALSE)
regsplt <- regsplt[!is.na(inc),
                   .(inc=sum(inc),newrel=sum(newrel,na.rm=TRUE)),
                   by=.(sex,age_group,g.whoregion)]
regsplt$age <- gsub('_','-',regsplt$age_group)
regsplt$age <- factor(regsplt$age,levels=agz3,ordered=TRUE)
regsplt$g.whoregion <- factor(regsplt$g.whoregion,
                              levels=sort(unique(regsplt$g.whoregion)))
regsplt$sex <- factor(regsplt$sex,ordered=FALSE)
fn <- here('disaggregation/output/incsplits/data/regsplt.Rdata')
save(regsplt,file=fn)    #save
load(file=fn)    #load



whoz <- sort(unique(regsplt$g.whoregion))
whozt <- c('Africa','The Americas',
           'Eastern Mediterranean','Europe','South-East Asia',
           'Western Pacific')
## sw <- c(1,4,2,5,3,6); whozt <- whozt[sw]; whoz <- whoz[sw]

nplst <- list()
for(i in seq_along(whoz)){
  reg <- regsplt[g.whoregion==whoz[i]]
  nplst[[i]] <- 
    ggplot(data=reg,aes(x=age,y=newrel,fill=sex)) +
    coord_flip() +
    geom_bar(data=reg[sex=='M'],stat='identity',aes(y=newrel))+
    geom_bar(data=reg[sex=='F'],stat='identity',aes(y=-newrel))+
    scale_y_continuous(labels = absspace) +
    geom_bar(data=reg[sex=='M'],stat='identity',
             aes(x=age,y=inc),fill='transparent',col=1)+
    geom_bar(data=reg[sex=='F'],stat='identity',
             aes(x=age,y=-inc),fill='transparent',col=1) +
    ggtitle(whozt[i]) + theme_classic()+
    theme(legend.position="none") + xlab('') + ylab('') +
    scale_x_discrete(labels=agz4)
} 


#' Save out plot for regions
fn <- here('disaggregation/output/incsplits/aplots/regional.pdf')
plt <- ggarrange(plotlist = nplst,ncol=3,nrow=2)
if(gypt){
  ggsave(plt,file=fn,height=10,width=15,device = cairo_pdf)
} else {
  ggsave(plt,file=fn,height=10,width=15)
}


#' Global calculation and graph
glob <- regsplt[!is.na(inc),
                .(inc=sum(inc),newrel=sum(newrel,na.rm=TRUE)),
                by=.(sex,age)]
glob$age <- factor(glob$age,levels=agz3,ordered=TRUE)
fn <- here('disaggregation/output/incsplits/data/glob.Rdata')
save(glob,file=fn)    #save
load(file=fn)    #load


plt <- ggplot(data=glob,aes(x=age,y=newrel,fill=sex)) +
  coord_flip() +
  geom_bar(data=glob[sex=='M'],stat='identity',aes(y=newrel))+
  geom_bar(data=glob[sex=='F'],stat='identity',aes(y=-newrel))+
  ylab('New & relapse cases') + xlab('Age group')+
  scale_y_continuous(labels = absspace) +
  geom_bar(data=glob[sex=='M'],stat='identity',
           aes(x=age,y=inc),fill='transparent',col=1)+
  geom_bar(data=glob[sex=='F'],stat='identity',
           aes(x=age,y=-inc),fill='transparent',col=1) +
  theme_classic() + 
  theme(legend.position="none") + xlab('') + ylab('') +
  scale_x_discrete(labels=agz4)

fn <- here('disaggregation/output/incsplits/aplots/global.pdf')
if(gypt){
  ggsave(plt,file=fn,height=7,width=7,device = cairo_pdf)
} else {
  ggsave(plt,file=fn,height=7,width=7)
}

#' Also save data relevant to plots elsewhere:
#'

## make denominator data in correct form
denom <- pop[year==estyr]
keep <- c('iso3',
          grep('e.pop.m',names(pop),value=TRUE),
          grep('e.pop.f',names(pop),value=TRUE))
keep <- keep[!grepl('plus',keep)]
keep <- keep[!grepl('014',keep)]
denom <- denom[,..keep]
denom <- melt(denom,id='iso3')
denom[,variable:=gsub("e\\.pop\\.","",variable)]
denom[,sex:=ifelse(grepl('m',variable),'M','F')]

denom[,age:=gsub('m|f','',variable)]
denom[,age:=gsub('5','5-',age)]
denom[age=='04',age:='0-4'] #tidy rest
denom[age=='65-',age:='65plus']
denom[age=='45-5-4',age:='45-54']
denom[age=='5-5-64',age:='55-64']
denom[,unique(age)] #check
denom$age <- factor(denom$age,levels=agz3,ordered=TRUE)
denom$sex <- factor(denom$sex,ordered=FALSE)
denom[,variable:=NULL]
names(denom)[2] <- 'pop'
summary(denom)
fn <- here('disaggregation/output/denom.Rdata')
save(file=fn,denom)



## hi30
h30splt <- incsplit[iso3 %in% hbc,
                    .(iso3,inc,newrel,age_group,sex)]
h30splt$age <- gsub('_','-',h30splt$age_group)
h30splt$age <- factor(h30splt$age,levels=agz3,ordered=TRUE)
h30splt$sex <- factor(h30splt$sex,ordered=FALSE)
h30splt[,age_group:=NULL]
h30splt <- merge(h30splt,hbcsh,by='iso3')
h30splt[iso3=='MOZ' & !age %in% c('0-4','5-14'),
        newrel:=NA] #NOTE this year - as has been split evenly in dx plt

h30splt <- merge(h30splt,denom,by=c('iso3','sex','age'),
             all.x=TRUE,all.y = FALSE)

fn <- here('disaggregation/reportoutput/h30splt.rda')
save(file=fn,h30splt)


## regional
whosh <- data.table(g.whoregion=whoz,name=whozt)
regsplt <- merge(regsplt,whosh,by='g.whoregion')
regdenom <- merge(denom,unique(est[,.(iso3,g.whoregion)]),by='iso3')
regdenom <- regdenom[,.(pop=sum(pop)),by=.(g.whoregion,sex,age)]
regsplt <- merge(regsplt,regdenom,by=c('g.whoregion','sex','age'))

fn <- here('disaggregation/reportoutput/regsplt.rda')
save(file=fn,regsplt)

## global
globsplt <- copy(glob)
globsplt$age <- factor(globsplt$age,levels=agz3,ordered=TRUE)
globsplt$sex <- factor(globsplt$sex,ordered=FALSE)
globdenom <- regsplt[,.(pop=sum(pop)),by=c('sex','age')]
globsplt <- merge(globsplt,globdenom,by=c('sex','age'))
fn <- here('disaggregation/reportoutput/globsplt.rda')
save(file=fn,globsplt)



## data preparation for analysis.
rm(list=ls())
## ------------libraries--------
library(data.table)
library(here)

## -------functions--------
oddit <- function(x) x/(1-x)
ioddit <- function(x) x/(1+x)
logit <- function(x) log(oddit(x))
ilogit <- function(x) ioddit(exp(x))

estyr <- 2021                           #NOTE

## needs hi/lo from est
## --- data ----
load(here('disaggregation/output/incsplits/data/tbimp.Rdata')) #notifications
load(here('disaggregation/R/02children/data/K.Rdata'))
load(here('data/est.rda')) #
load(here('data/pop.rda'))
load(here('data/unaids.rda'))
load(here('disaggregation/R/02children/data/BB.Rdata'))
BBD <- copy(BB)


## ART estimates
A <- fread(here('disaggregation/R/02children/mortinput/Treatment cascade_People living with HIV receiving ART (%)_Population_ All ages.csv'),
           header=TRUE)

## --------work----------
K <- merge(K,est[year==estyr][,.(iso3,inc.num=inc*pop/1e5)]
          ,by='iso3',all.x=TRUE,all.y = FALSE)
B <- K[,.(iso3,mid=inc.num*propu15b,var=(inc.num*propu15b.sd)^2)] #bug here (nb inc var included in prop.sd)
B[,inc.num.014:=mid]
B[,inc.num.014.lo:=pmax(0,mid - 1.96*sqrt(var))]
B[,inc.num.014.hi:=mid + 1.96*sqrt(var)]

## sample from the estimates assuming gamma (mean=k scale, var = k scale^2)
B$S <- B$var/B$mid
B$k <- B$mid/B$S

## ## inspect if reasonable representation of burden
## who <- 200+15                                #choose test country
## png(paste0('./test',B$country[who],'.png'))
## hist(rgamma(1e4,shape=B$k[who],scale=B$S[who]),main=B$country[who])
## abline(v=B$mid[who],col=2);
## abline(v=B$inc.num.014.hi[who],lty=2,col=2);
## abline(v=B$inc.num.014.lo[who],lty=2,col=2);
## dev.off()


## merge with notification data
N <- tb
nmz <- grep('14',names(N),value=TRUE)
nmz <- c(nmz,grep('04',names(N),value=TRUE))
N <- N[,c('country','iso3','year','g.whoregion',nmz),with=FALSE] #restsrict to relevant
dim(N)
N <- N[!is.na(newrel.m04)]
dim(N)
N <- N[,.SD[year==max(year)],by=iso3]   #most recent year with data
dim(N)

## all kids
N[,notifs:=newrel.m014 + newrel.f014 ]   #0-14 notifications, new system
N[,sum(is.na(notifs))]                   #7 use most recent
addon <- N$newrel.sexunk014
addon[is.na(addon)] <- 0                #most are NA
N$notifs <- N$notifs + addon            #add on those that aren't

## young
N[,notifsY:=newrel.m04 + newrel.f04]
addon <- N$newrel.sexunk04
addon[is.na(addon)] <- 0                #most are NA
N$notifsY <- N$notifsY + addon            #add on those that aren't

## older
N[,notifsO:=newrel.m514 + newrel.f514]
addon <- N$newrel.sexunk514
addon[is.na(addon)] <- 0                #most are NA
N$notifsO <- N$notifsO + addon            #add on those that aren't

## inspect
N[,.(notifs,notifsY,notifsO)]
N[,sum(is.na(notifs))]
N[,sum(is.na(notifsY))]
N[,sum(is.na(notifsO))]

N[is.na(notifs),notifs:=notifsY+notifsO]
N[,sum(is.na(notifs))]

## combind paediatric estimates and notification data
BB <- merge(B,
            N[,.(notifs,notifsY,notifsO,g.whoregion,iso3)],
            by='iso3') #merge with notifications
BB$id <- 1:nrow(BB)                     #for later

## merge in overall population data
BB <- merge(BB,
            pop[year==estyr][,.(iso3,e.pop.num,e.pop.014 = e.pop.m014 + e.pop.f014)],
            by='iso3')

## merge in overall incidence estimates
BB <- merge(BB,
            est[year==estyr][,.(iso3,year,inc,inc.lo,inc.hi)],
            by='iso3')
BB[,inc.num:=e.pop.num*inc/1e5]
BB[,inc.num.lo:=e.pop.num*inc.lo/1e5]
BB[,inc.num.hi:=e.pop.num*inc.hi/1e5]
BB[is.na(inc.num.hi),inc.num.hi:=inc.num*1.1] #
BB[is.na(inc.num.lo),inc.num.lo:=inc.num*0.9]


## ---- HIV data
H <- unaids[year==estyr,.(iso3,HIVmid=h014,
                         HIVlo=h014.lo, HIVhi=h014.hi)]
BBA <- merge(BB,H,by='iso3',all.x=TRUE)
BB <- BBA


## ART from AIDSinfo in horrible format (although slightly improved)
A <- A[,.(country=Country,ART=`2020`)]
A[ART=='...',ART:=NA]
A[grep('>',ART),ART:='97.5']
A[,ART:=as.numeric(ART)]
A <- A[country!='Global']
AA <- merge(A,unique(pop[,.(iso3,country)]),by='country')
missed <- setdiff(A$country,AA$country) #some missed by name merge
setkey(A,country)
(tmp <- A[missed])
## tmp[country=="C\303\264te d'Ivoire",iso3:='CIV']
tmp[country=="United Kingdom",iso3:='GBR']
tmp[country=="United States",iso3:='USA']
AA <- rbind(AA,tmp)
AA[,country:=NULL]
BB <- merge(BB,AA,by='iso3',all.x=TRUE)

## convert ART/HIV NAs to regional averages (countries with no HIV data ~ no HIV)
tmp <- BB[,.(tmp=mean(ART,na.rm=TRUE)),by=g.whoregion]
BB <- merge(BB,tmp,by='g.whoregion',all.x=TRUE,all.y=FALSE)
## BB[is.na(ART),.(iso3,HIVmid,ART)]
BB[is.na(ART),ART:=round(tmp)]
BB[,tmp:=NULL]

## same for HIV
tmp <- BB[,.(tmp=mean(HIVmid/e.pop.014,na.rm=TRUE),
             tmp.se=mean((HIVhi-HIVlo)/(2*e.pop.014),na.rm=TRUE)),
          by=g.whoregion]
BB <- merge(BB,tmp,by='g.whoregion',all.x=TRUE,all.y=FALSE)
## BB[is.na(HIVmid),.(iso3,HIVmid,HIVlo,HIVhi)]
BB[is.na(HIVmid),HIVmid:=as.integer(tmp*e.pop.014)]
BB[is.na(HIVlo),HIVlo:=as.integer(pmax(tmp-2*tmp.se,0)*e.pop.014)]
BB[is.na(HIVhi),HIVhi:=as.integer((tmp+2*tmp.se)*e.pop.014)]
BB[,c('tmp','tmp.se'):=NULL]
## BB[,.(iso3,HIVmid,HIVlo,HIVhi,e.pop.num,tmp,tmp.se)]

## proportion ART/HIV
## make ART uncertainty derive from HIV denominator uncertainty!! NB
BB$artmid <- BB$ART/100
ses <- (BB$HIVhi-BB$HIVlo) / (BB$HIVmid+1e-6) 
BB$artlo <- ilogit( logit(BB$artmid) * (1+ses/2))
BB$arthi <- ilogit( logit(BB$artmid) * (1-ses/2))


## HIV
BB$hivmid <- BB$HIVmid/BB$e.pop.014
BB$hivlo <- BB$HIVlo/BB$e.pop.014
BB$hivhi <- BB$HIVhi/BB$e.pop.014


## distributions
BB$hS <- ((BB$hivhi-BB$hivlo)/3.92)^2/BB$hivmid
BB$hk <- BB$hivmid/BB$hS
bad <- !is.finite(BB$hk) | !is.finite(BB$hS)
BB$hk[bad] <- BB$hS[bad] <- NA
BB$aS <- ((BB$arthi-BB$artlo)/3.92)^2/BB$artmid
BB$ak <- BB$artmid/BB$aS
bad <- !is.finite(BB$ak) | !is.finite(BB$aS) 
BB$ak[bad] <- BB$aS[bad] <- NA



## inspect if reasonable representation of burden
who <- which(BB$iso3=='BWA')                                #choose test country
png(here('disaggregation/R/02children/plots',
         paste0('test_HIV_',BB$iso3[who],'.png')))
hist(rgamma(1e4,shape=BB$hk[who],scale=BB$hS[who]),main=BB$country[who])
abline(v=BB$hivmid[who],col=2);abline(v=(BB$hivlo)[who],lty=2);
abline(v=(BB$hivhi)[who],lty=2);
dev.off()

png(here('disaggregation/R/02children/plots',
         paste0('test_ART_',BB$iso3[who],'.png')))
hist(rgamma(1e4,shape=BB$ak[who],scale=BB$aS[who]),main=BB$country[who])
abline(v=BB$artmid[who],col=2);abline(v=(BB$artlo)[who],lty=2);
abline(v=(BB$arthi)[who],lty=2);
dev.off()

## age splits
AF <- BBD[,.(iso3,aa,ab)]
BA <- merge(BB,AF,by='iso3',all.x = TRUE)
BA$aa[is.na(BA$aa)] <- mean(BA$aa,na.rm=TRUE)
BA$ab[is.na(BA$ab)] <- mean(BA$ab,na.rm=TRUE)
BB <- BA


## Avoid issues with remaining few NAs
## cf previous version:
## NAs in notifs x3
## hS hk, aS, ak
BB[is.na(notifs),notifs:=0]
BB[is.na(notifsY),notifsY:=0]
BB[is.na(notifsO),notifsO:=0]
BB[!is.finite(hS),hS:=0]
BB[!is.finite(hk),hk:=0]
BB[!is.finite(aS),aS:=0]
BB[!is.finite(ak),ak:=0]

## BB for mortality
MBB <- BB

## save/load here!
save(MBB,file=here('disaggregation/R/02children/mortinput/MBB.Rdata'))

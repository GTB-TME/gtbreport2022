## experiments to help with CV work
rm(list=ls())
## packages
library(here)
library(glue)
library(data.table)
## utilities
gh <- function(x) glue(here(x))
## load data
load(gh('drtb/data/D4S.Rdata'))
load(gh('drtb/data/isoidx5.Rdata'))
load(here('drtb/data/regkey.Rdata'))
load(here('drtb/data/RPD.Rdata'))
## data prep
ckey <- data.table(iso3=isoidx,cid=1:length(isoidx))
RPD <- merge(RPD[year>=2000],regkey[,.(iso3,g_whoregion)])
RPD[,max(year),by=type]
maxy <- max(RPD$year)


## country data counts
RPD[,cdc:=.N,by=iso3] #country data count
CC <- unique(RPD[,.(iso3,cdc)])
CC[,table(cdc)]
plot(CC[,table(cdc)])
CC <- CC[order(cdc)]
CC[,ccdc:=cumsum(cdc)]
CC[,pcdc:=cumsum(cdc)/sum(cdc)]
CC[,ncns:=nrow(CC):1]
CC <- merge(CC,regkey[,.(iso3,g.whoregion)],by='iso3')

## potential cut points
## OLD versions:
## CC[pcdc>1-0.8] #89 countries, 8 dp (C)
## CC[cdc>=10] #77 countries, 10 dp ~75% data points
## CC[pcdc>1-0.75] #77 countries, 10 dp ~75% data points (B)
CC[pcdc>1-0.5] #97 countries, 18 data points (A)
CC[pcdc>1-0.5,table(g.whoregion)]

## AFR AMR EMR EUR SEA WPR 
## 23  13  13  22   7  19 


pdf(gh('drtb/plots/cCV.pdf'))
plot(CC$ncns,1-CC$pcdc,
     ylab='Prop data points included',
     xlab = 'No. countries included in CV study' )
grid()
## abline(h=0.75,col=4,lty=2);abline(v=77,col=4,lty=2)
abline(h=0.5,col=6,lty=2);abline(v=97,col=6,lty=2)
dev.off()

## create & save worklists

c.omit.vA <- CC[pcdc>1-0.5,iso3]
c.omit.vB <- CC[pcdc>0,iso3] #NOTE all
length(c.omit.vB)            #198

cat(c.omit.vA,file=gh('drtb/data/c.omit.vA.txt'))
cat(c.omit.vB,file=gh('drtb/data/c.omit.vB.txt'))


## NOTE seems this year many more issues with reporting
## drop the MR drop analysis
## interested in countries with data in most recent year:
RPD[,dy:=maxy-year]
mrc <- RPD[dy==0,unique(iso3)]
(norecent <- setdiff(c.omit.vA,mrc)) #these countries don't have dy=0
setdiff(c.omit.vA,norecent)
c.omit.vAt <- intersect(c.omit.vA,mrc)
cat(c.omit.vAt,file=gh('data/c.omit.vAt.txt'))

## simple version of INH incidnence based on Philippe's previous approach to RR/MDR
rm(list = ls())

## libraries
library(here)
library(glue)
library(data.table)
library(metafor)
gh <- function(x)glue(here(x))
ssum <- function(x,...) sqrt(sum(x^2,...)) #sqrt sum sqrs

## utilities &  helper functions for binomial CIs, incidence calculator
source(gh('R/fun.R'))
source(gh('drtb/R/utils/CIhelpers.R'))
source(here('drtb/R/utils/incdr.R'))

## globals
yr <- 2021

## main data
load(here('data/est.rda')) #estimates
load(here('data/tb.rda'))  #notifications
load(file=here('drtb/data/rhofnr.Rdata')) #f & r data - see 5update_fr.R

## data from 1prepdata.R
load(gh('drtb/data/LB.Rdata')) #surveiLLance
load(gh('drtb/data/YB.Rdata')) #surveY

## restrict to that which is rated for use in H-computations
LB <- LB[Hrated=='yes']
YB <- YB[Hrated=='yes']

## for each country, restrict to the most recent data within last 5 years
LB[,mry:=max(year),by=iso3]
YB[,mry:=max(year),by=iso3]
LB <- LB[year==mry & yr-year<=5]
YB <- YB[year==mry & yr-year<=5]
nrow(LB) #165
nrow(YB) #35
intersect(LB$iso3,YB$iso3) #"LBN" "SWZ" "TGO" "TJK" "IND" "LKA"

## use H regardless of R TODO do we want MDR if missing?
YBR <- YB[!is.na(dst_rlt_hr_pcnt),
          .(id,iso3,year,patients,
            prop.H=dst_rlt_hr_pcnt/1e2,
            prop.H.lo=dst_rlt_hr_pcnt_lo/1e2,
            prop.H.hi=dst_rlt_hr_pcnt_hi/1e2)]
nrow(YBR) #27


## "dr_h_nr"    "dr_r_nh"    "mdr"        "dst_rlt"    "dst_rlt_hr"
## dst_rlt_hr NOTE other exists too
LB[!is.na(dr_h_nr)] #5 records with H & R separately, all lack dst_rlt_hr
LB[!is.na(dr_h_nr),dst_rlt_hr:=dr_h_nr + mdr] #aggregate H-resistance to create dst_rlt_hr

## create binomial confidence intervals
LB[,c('prop.H','prop.H.lo','prop.H.hi'):=MLH(dst_rlt_hr,dst_rlt)]
LB[,c('prop.H','prop.H.lo','prop.H.hi'):=.(prop.H/1e2,prop.H.lo/1e2,prop.H.hi/1e2)]
LBR <- LB[,.(id,iso3,year,patients,prop.H,prop.H.lo,prop.H.hi)]


## any in both now?
(both <- intersect(LBR$iso3,YBR$iso3)) #"SWZ" "TJK" "IND" "LKA"
YBR[iso3 %in% both]
LBR[iso3 %in% both]

## real overlaps are:
## SWZ/new 2018 vs 2017 - quite different
## TJK/new 2017 vs 2021
## IND/ret 2016 vs 2019 - quite different
## LKA/ret 2018 vs 2016
## SWZ/ret 2018 vs 2017
## TJK/ret 2017 vs 2021 - a bit different

## NOTE decision - use survey data for things in both
## TODO check
LBR <- LBR[!iso3 %in% both]

## combined dataset
HP <- rbind(YBR,LBR) #186 rows, 128 countries


## regional meta-analyses - following PGs example
HP[,prop.H.sd:=(prop.H.hi-prop.H.lo)/3.92]
tmp <- est[year==yr,.(iso3,g.mdr)]
nr <- nrow(tmp)
tmp <- rbind(tmp,tmp)
tmp[,patients:=c(rep('new',nr),rep('ret',nr))]
HP <- merge(HP,tmp,by=c('iso3','patients'),all.y=TRUE)

## new
fit.new <- rma(
  yi = prop.H,
  sei = prop.H.sd,
  data = HP[patients=='new' & !is.na(prop.H)],
  mods =  ~ g.mdr - 1
)

## retreatment
fit.ret <- rma(
  yi = prop.H,
  sei = prop.H.sd,
  data = HP[patients=='ret' & !is.na(prop.H)],
  mods =  ~ g.mdr - 1
)

reg <- names(table(HP$g.mdr))

## check the order of regions in reg is the same as in imp*
all.equal(reg, gsub('g.mdr', '', dimnames(fit.new$beta)[[1]]))
all.equal(reg, gsub('g.mdr', '', dimnames(fit.ret$beta)[[1]]))


imp.new <- fit.new$b[, 1]
imp.ret <- fit.ret$b[, 1]
imp.new.sd <- fit.new$se
imp.ret.sd <- fit.ret$se
imp.new.ui <- vlohi(unlist(imp.new), imp.new.sd)
imp.ret.ui <- vlohi(unlist(imp.ret), imp.ret.sd)


tmp <- copy(HP)
## not elegant, should set up a proper merge
for (i in 1:length(reg)) {
  sel <- HP$g.mdr == reg[i] & HP$patients=='new' & is.na(HP$prop.H)
  HP[sel, prop.H := imp.new[i]]
  HP[sel, prop.H.sd := imp.new.sd[i]]
  HP[sel, prop.H.lo := imp.new.ui[1, i]]
  HP[sel, prop.H.hi := imp.new.ui[2, i]]
  sel <- HP$g.mdr == reg[i]  & HP$patients=='ret' & is.na(HP$prop.H)
  HP[sel, prop.H := imp.ret[i]]
  HP[sel, prop.H.sd := imp.ret.sd[i]]
  HP[sel, prop.H.lo := imp.ret.ui[1, i]]
  HP[sel, prop.H.hi := imp.ret.ui[2, i]]
}


## TODO checking
HP

save(HP,file=gh('drtb/data/HP.Rdata'))


## incidence calculations
HPI <- dcast(HP[,.(iso3,patients,prop.H,prop.H.sd)],
             iso3 ~ patients,
             value.var = c('prop.H','prop.H.sd')) #reshap
HPI <- merge(HPI,rhofnr,by=c('iso3'))             #merge in f / r
HPI <- merge(HPI,est[year==yr,.(iso3,inc.num,inc.num.sd=(inc.hi.num-inc.lo.num)/3.92)],
             by=c('iso3'),all.x = TRUE,all.y=FALSE) #merge in incidence estimates

## all H incidence using (inc * ((1 - f) * pn * ((1 - r) + r * rho) + f * pr))
HPI[,c('inc.H','inc.H.sd'):=incdr(inc.num,
                                  inc.num.sd,
                                  f,
                                  f.sd,
                                  pn = prop.H_new,
                                  pn.sd = prop.H.sd_new,
                                  r,
                                  r.sd,
                                  rho = rho,
                                  rho.sd = rho.sd,
                                  pr = prop.H_ret,
                                  pr.sd = prop.H.sd_ret,
                                  TRUE),
    by=.(iso3)]



## ----------- aggregation etc -------
KOH <- HPI[,.(iso3,
               e_inc_hr_num=inc.H,
               e_inc_hr_num_lo=pmax(0,inc.H-1.96*inc.H.sd),
               e_inc_hr_num_hi=inc.H+1.96*inc.H.sd
               )]
attr(KOH, "timestamp") <- Sys.Date() #set date

## save all time version
fn <- gh('drtb/outdata/KOH.rda')
save(KOH,file=fn)

## save out as CSV and Rdata
dt <- gsub('-','_',Sys.Date())
fn <- gh('drtb/dboutput/HR_country_{dt}.csv')
fwrite(KOH,file=fn)

## === global summaries
## -- totals
glob <- HPI[,.(inc=sum(inc.num),
               inc.sd=ssum(inc.num.sd),
               inc.H=sum(inc.H),
               inc.H.sd=ssum(inc.H.sd)
               )]

glob[,c('inc.H.lo','inc.H.hi'):=
        .(pmax(0,inc.H-1.96*inc.H.sd),inc.H+1.96*inc.H.sd)
     ]

## --- proportions
glob[,c('prop.H'):=
        .(inc.H/inc)]
glob[,prop.H.sd:=prop.H * sqrt((inc.H.sd/inc.H)^2+(inc.sd/inc)^2)]
glob[,c('prop.H.lo','prop.H.hi'):=
        .(pmax(0,prop.H-1.96*prop.H.sd),
          pmin(1,prop.H+1.96*prop.H.sd)
          )]

glob.H <- glob

save(glob.H,file=gh('drtb/data/glob.H.Rdata'))


fn <- gh('drtb/dboutput/HR_glob_{dt}.csv')
fwrite(glob.H,file=fn)

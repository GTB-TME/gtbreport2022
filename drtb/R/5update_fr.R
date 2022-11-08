## adapting PG's analysis to create the f & r parms for the DR incidence formula
rm(list = ls())

## libraries
library(here)
library(glue)
library(data.table)
library(metafor)
gh <- function(x)glue(here(x))

## utilities
source(here('R/fun.R'))
vlohi <- Vectorize(lohi, c('ev', 'sd'))

## data
load(here('data/est.rda')) #estimates
load(here('data/tb.rda'))  #notifications
load(here('data/rel.rda')) #dataset on relapses
load(here('data/dre.rda')) #most recent DRS data

## globals
yr <- 2021
newversion <- TRUE #NOTE change this flag to either run old analysis, or new

## generate/load rr object
if(!newversion){
  ## run PGs approach (modified) to create rr object
  source(here('drtb/R/utils/PGprops.R'))
} else {
  ## load chosen estimates
  mdl <- scan(here('drtb/R/utils/modelchoice.txt'),what='char') #NOTE depends choice
  load(file=gh('drtb/outdata/KO.{mdl}.Rdata'))
  ## reformat
  rr <- dcast(KO[year==yr],
              iso3+year+type+g_whoregion ~ patients,
              value.var = c('RR.mid','RR.lo','RR.hi'))
  oldnmz <- c("g_whoregion",
              "RR.mid_new", "RR.mid_ret",
              "RR.lo_new", "RR.lo_ret",
              "RR.hi_new", "RR.hi_ret")
  newnmz <- c("g.whoregion",
              "prop.rr.new", "prop.rr.ret",
              "prop.rr.new.lo", "prop.rr.ret.lo",
              "prop.rr.new.hi", "prop.rr.ret.hi") 
  setnames(rr,oldnmz,newnmz)
  ## convert from percent to proportion
  xcols <- newnmz[-1]
  rr[,(xcols):=lapply(.SD,function(x)x/1e2),.SDcols=xcols]
  rr[,c('prop.rr.new.sd','prop.rr.ret.sd'):=
        .((prop.rr.new.hi-prop.rr.new.lo)/3.92,(prop.rr.ret.hi-prop.rr.ret.lo)/3.92)]

}

## rrtb[,length(unique(iso3))]
## est[,length(unique(iso3))]

## setdiff(est[,(unique(iso3))],rrtb[,unique((iso3))])
## setdiff(rrtb[,(unique(iso3))],rr[,unique((iso3))])


## relevant bits of notification data
rrtb <- tb[year > yr - 5, .(iso3, year, c.newinc, c.ret, ret.nrel)]
setkey(rrtb, iso3, year)
rrtb[, ret.rel := c.ret - ret.nrel] #relapse
rrtb[!is.na(ret.rel), test.ispos(ret.rel)] #check

## --- f ---
rrtb[, f := ret.nrel / c.newinc] #definition of f

## safeties
rrtb[c.newinc == 0, f := 1]
rrtb[c.newinc < c.ret, f := 1]
rrtb[!is.na(f), test.isbinom(f)]

## --- r ---
rrtb[f < 1,
     r :=
       pmin((ret.rel / c.newinc) / (1 - f), #defition
            1)                              #safety
     ]
rrtb[f == 1, r := 0] #safety

## check
rrtb[!is.na(r), test.isbinom(r)]

## collect f & r values by country
out <- rrtb[, .(
  r = mean(r, na.rm = TRUE),
  r.sd = sd(r, na.rm = TRUE),
  f = mean(f, na.rm = TRUE),
  f.sd = sd(f, na.rm = TRUE)
), by = iso3]
summary(out) #very few NAs

## checks
out[!is.na(r), test.ispos(r)]
out[!is.na(r.sd), test.ispos(r.sd)]
out[!is.na(f), test.ispos(f)]
out[!is.na(f.sd), test.ispos(f.sd)]

## where SD is missing, assume 20% of mean
out[is.na(r.sd) & !is.na(r), r.sd := r * .2]
out[is.na(f.sd) & !is.na(f), f.sd := f * .2]

## uniformize NAs
out[!is.finite(r), r := NA]
out[!is.finite(r.sd), r.sd := NA]
out[!is.finite(f), f := NA]
out[!is.finite(f.sd), f.sd := NA]
summary(out) #very few NAs

## imputation - pool r and f by mdr region
out <- merge(out, est[year == yr, .(iso3, g.mdr)], by = 'iso3')

fit.r <- rma(
  yi = r,
  sei = r.sd,
  data = out[r.sd > 0],
  mods =  ~ g.mdr - 1
)
fit.f <- rma(
  yi = f,
  sei = f.sd,
  data = out[f.sd > 0],
  mods =  ~ g.mdr - 1
)

imp.r <- fit.r$b[, 1]
imp.f <- fit.f$b[, 1]
imp.r.sd <- fit.r$se
imp.f.sd <- fit.f$se

## check the order of regions in reg is the same as in imp*
reg <- names(table(est$g.mdr))
all.equal(reg, gsub('g.mdr', '', dimnames(fit.r$beta)[[1]]))
all.equal(reg, gsub('g.mdr', '', dimnames(fit.f$beta)[[1]]))

## out[, g.mdr := NULL]
if('g.mdr' %in% names(rr)) rr[,g.mdr:=NULL]
rr <- merge(rr, out, by = 'iso3', all.x = TRUE)

## add in imputed f/r - regional MAs
table(!is.na(rr$r))
table(!is.na(rr$f))
for (i in 1:length(reg)) {
  sel <- rr$g.mdr == reg[i] & (is.na(rr$r) | rr$r == 0)
  rr[sel, r := imp.r[i]]
  rr[sel, r.sd := imp.r.sd[i]]
  sel <- rr$g.mdr == reg[i] & (is.na(rr$f) | rr$f == 0)
  rr[sel, f := imp.f[i]]
  rr[sel, f.sd := imp.f.sd[i]]
}

## check!
table(!is.na(rr$r))
table(!is.na(rr$f))

## --- rho ---
## $\rho$ and its SD
## (see rel.R)
rho <- 5.5 # reverting to last year's after USAID discussions
rho.sd <- (6.8 - 4.4) / 3.92
rho.cv <- rho.sd / rho

## use the relapse dataset
rel2 <- rel[level=='national',
            .(rho.rr=last(ratio.rel.new),
              rho.sd.rr=last(ratio.rel.new.sd)),
            by=iso3]

(rel3 <- rel2[rho.rr<=rho & rho.rr>1]) #15 countries with their own sensible value
rr <- merge(rr, rel3, by='iso3', all.x=TRUE) #use these for preference

## but check that these results are sensible & modify otherwise
rr[!is.na(rho.rr),
   rho.rr := ifelse(prop.rr.ret / prop.rr.new > rho.rr,
                    rho.rr,
                    prop.rr.ret / prop.rr.new)]

## for the countries not in rel use rho, but with same sense-corrections
rr[is.na(rho.rr), rho.rr := ifelse(prop.rr.ret / prop.rr.new > rho, rho, prop.rr.ret / prop.rr.new)]
rr[is.na(rho.sd.rr), rho.sd.rr := rho.rr * rho.cv]

## keep what we need:
rhofnr <- rr[,.(iso3,rho=rho.rr,rho.sd=rho.sd.rr,f,f.sd,r,r.sd)] #NOTE rename
summary(rhofnr)

if(!newversion){
  rhofnr.old <- copy(rhofnr)
  save(rhofnr.old,file=here('drtb/data/rhofnr.old.Rdata'))
} else {
  library(ggplot2)
  library(ggrepel)
  save(rhofnr,file=here('drtb/data/rhofnr.Rdata'))
  load(file=here('drtb/data/rhofnr.old.Rdata'))

  ## comparison test with 'old' version
  names(rhofnr.old)[2:ncol(rhofnr.old)] <- paste0(names(rhofnr.old)[2:ncol(rhofnr.old)],".old")
  CF <- merge(rhofnr,rhofnr.old,by='iso3')
  CFM <- melt(CF,id='iso3')
  CFM[,vn:=ifelse(grepl('old',variable),'old','new')]
  CFM[,variable:=gsub("\\.old","",variable)]
  CFM <- dcast(CFM,iso3+variable~vn)
  CFM[abs(new/old-1)>0.1]
  CFM[,summary(abs(new/old-1))]

  ## plot
  ggplot(CFM[!grepl('rho',variable)],aes(old,new,label=iso3))+
    geom_point()+
    geom_text_repel()+
    facet_wrap(~variable)+
    geom_abline(intercept = 0,slope=1,col=2)

  CFM[variable=='rho' & abs(new/old-1)>0.1] #TODO check

}


## TODO need to check whether use of specific year in new version causes trouble & fix if so

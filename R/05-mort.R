#' ---
#' title: Mortality
#' author: Philippe Glaziou
#' date: 18/07/2022
#' output:
#'    html_document:
#'      mode: selfcontained
#'      toc: true
#'      toc_depth: 3
#'      toc_float: true
#'      number_sections: true
#'      theme: flatly
#'      highlight: zenburn
#'      df_print: paged
#'      code_folding: hide
#' ---

#' # Preamble
#' (Last updated: `r Sys.Date()`)
#'
# Mortality HIV-neg
#
# Load libraries and data
#
library(data.table)
library(imputeTS)
library(propagate)
library(here)


load(here('data/cty.rda'))
load(here('data/vr.rda'))
load(here('data/old.rda'))
load(here('data/model.rda'))
load(here('data/extra.rda'))
load(here('data/unaids.rda'))


# !!! update the following line as needed !!!
est <- fread('csv/est_04inc_2022-09-27.csv')
setkey(est, iso3, year)


m <- 1e5
yr <- 2021

source(here('R/fun.R'))


# vectorized lohi
#
vlohi <- Vectorize(lohi, c('ev', 'sd'))


# import vr
dim(est)
est <-
  merge(est, vr[, .(
    iso3,
    year,
    vr.keep = keep.vr,
    vr.garbage = garbage,
    vr.coverage,
    vr.quality = codqual,
    vr.env = env,
    # ghe.env, ghe.env.lo, ghe.env.hi,
    vr.mort.nh = tb.adj * m / pop,
    vr.raw = tb * m / pop,
    vr.mort.nh.sd = tb.adj.sd * m / pop
  )],
  by = c('iso3', 'year'), all.x = TRUE)
dim(est)


# incorporate old values of mortality (2000 - 2019)
est <-
  merge(est, old[, .(iso3,
                     year,
                     mort.nh,
                     mort.nh.sd,
                     mort.h,
                     mort.h.sd,
                     mort,
                     mort.sd,
                     old.source.mort = source.mort)], by = c('iso3', 'year'), all.x = TRUE)
dim(est)



# RUS 2021
sel <- est$iso3=='RUS' & est$year==2021
est[sel, mort.nh := vr.mort.nh]
est[sel, mort.nh.sd := vr.mort.nh.sd]
est[sel, source.mort := 'VR']


# VR updates since the June 2021 snapshot of the mortality database
sel <- !is.na(est$vr.mort.nh) & est$vr.keep==TRUE & est$old.source.mort!='IHME'
table(sel)
est[sel, mort.nh := vr.mort.nh]
est[sel, mort.nh.sd := vr.mort.nh.sd]
est[sel, source.mort := 'VR']
est[is.na(source.mort), source.mort := old.source.mort]




# check missing values: should only be for 2021
sum(is.na(est$mort.nh) &
      est$year < yr) == 0   # TRUE: only year==yr inc values are missing

# check VR updates where usable
est[year < yr &
      vr.keep == TRUE, sum(vr.mort.nh / mort.nh > 1.1, na.rm = T)]

est[year < yr &
      vr.keep == TRUE & vr.mort.nh / mort.nh > 1.1, .(iso3, year, mort.nh, vr.mort.nh)]
# all in ZAF, we ignore due to major miscoding issues btw HIV and TB causes


est[year < yr &
      vr.keep == TRUE, sum(vr.mort.nh / mort.nh < .9, na.rm = T)]





# update indirect estimates, stored in e.mort.*
out1 <-
  est[, {
    tmp = inc2mort(inc, inc.sd, imp.newinc, tbhiv, tbhiv.sd, noHIV =
                     T)$prop

    list(mort.nh = tmp[2],
         mort.nh.sd = tmp[4])
  },
  by = .(iso3, year)]

out2 <-
  est[, {
    tmp = inc2mort(inc, inc.sd, imp.newinc, tbhiv, tbhiv.sd, noHIV =
                     F)$prop

    list(mort.h = tmp[2],
         mort.h.sd = tmp[4])
  },
  by = .(iso3, year)]

est[, e.mort.nh := out1$mort.nh]
est[, e.mort.nh.sd := out1$mort.nh.sd]
est[, e.mort.h := out2$mort.h]
est[, e.mort.h.sd := out2$mort.h.sd]
est[, e.mort := e.mort.h + e.mort.nh]
est[, e.mort.sd := sqrt(e.mort.h.sd ^ 2 + e.mort.nh.sd ^ 2)]



# complete RUS 2021 (add HIV+ and totals)
sel <- est$iso3=='RUS' & est$year %in% 2020:2021
est[sel, mort.h := e.mort.h]
est[sel, mort.h.sd := e.mort.h.sd]
est[sel, mort := mort.nh + mort.h]
est[sel, mort.sd := sqrt(mort.h.sd ^ 2 + mort.nh.sd ^ 2)]
est[sel, source.mort := 'VR']



# incorporate "focus" countries (Nim)
#

# HIV+ and HIV-
B2 <- copy(est) # backup point

(md.lst <- unique(model$iso3))

# exclude some countries from Nim's models:
# - Georgia (comms 2022)
# - Ethiopia, DPRK, South Africa because shortfalls in line with pre-2020 general decline
#   not necessarily due to covid-related drops in detection & treatment
# - CHN and RUS (model will not used for incidence)
md.exclude <- c('GEO', 'ETH', 'PRK', 'ZAF', 'CHN', 'RUS')


dim(est)
est <-
  merge(est, model[scenario == 'COVID' &
                     hiv == 'a' &
                     measure == 'mort' &
                     iso3 %ni% md.exclude &
                     year %in% 2020:yr, .(
                       iso3,
                       year,
                       mort.md = best,
                       mort.md.lo = lo,
                       mort.md.hi = hi
                     )],
        by = c('iso3', 'year'), all.x = T)
dim(est)


sel <- est$year %in% 2020:yr & est$iso3 %in% md.lst & est$iso3 %ni% md.exclude
table(sel)
est[sel, mort := mort.md]
est[sel, mort.sd := (mort.md.hi - mort.md.lo ) / 3.92]
est[sel, mort.lo := mort.md.lo]
est[sel, mort.hi := mort.md.hi]
est[sel, sum(is.na(mort))]
est[sel, sum(is.na(mort.sd))]

est[sel, source.mort := "Country model"]
est[year==yr, table(source.mort)]
est[year==yr-1, table(source.mort)]



# HIV+ only
(md.h.lst <- unique(model[hiv=='pos', iso3]))

dim(est)
est <-
  merge(est, model[scenario == 'COVID' &
                     hiv == 'pos' &
                     measure == 'mort' &
                     iso3 %ni% md.exclude &
                     year %in% 2020:yr, .(
                       iso3,
                       year,
                       mort.h.md = best,
                       mort.h.md.lo = lo,
                       mort.h.md.hi = hi
                     )],
        by = c('iso3', 'year'), all.x = T)
dim(est)

sel <- est$year %in% 2020:yr & est$iso3 %in% md.h.lst & est$iso3 %ni% md.exclude
table(sel)
est[sel, mort.h := mort.h.md]
est[sel, mort.h.sd := (mort.h.md.hi - mort.h.md.lo ) / 3.92]
est[sel, mort.h.lo := mort.h.md.lo]
est[sel, mort.h.hi := mort.h.md.hi]
est[sel, sum(is.na(mort.h))]
est[sel, sum(is.na(mort.h.sd))]

est[sel, mort.nh := mort - mort.h]
est[sel & mort.sd<mort.h.sd,.(iso3,year,mort,mort.h,mort.nh,mort.sd,mort.h.sd)]
est[sel & mort.sd>mort.h.sd, mort.nh.sd := sqrt(mort.sd^2 - mort.h.sd^2)]
est[sel & mort.sd<mort.h.sd, mort.nh.sd := mort.sd * mort.nh/mort]
est[sel & mort.sd<mort.h.sd, mort.h.sd := mort.sd * mort.h/mort]



# HIV+ missing
(mh.lst <- setdiff(md.lst, md.h.lst))
sel <- est$iso3 %in% mh.lst  & est$iso3 %ni% md.exclude
table(sel)

est[sel, h.ratio := mort.h / mort]
est[sel & is.na(h.ratio), table(year)]
est[sel, h.ratio := imputeTS::na_locf(h.ratio), by=iso3]
est[sel & is.na(mort.h), mort.h := mort * h.ratio]
est[sel & is.na(mort.h.sd), mort.h.sd := mort.sd * h.ratio]
est[sel & is.na(mort.nh), mort.nh := mort * (1 - h.ratio)]
est[sel & is.na(mort.nh.sd), mort.nh.sd := mort.sd * (1 - h.ratio)]
est['VNM',.(iso3,year,mort,mort.nh,mort.h,h.ratio)]
est[, h.ratio := NULL]



# incorporate "extrapolated" countries (Nim's regional models)
# exclude high-income and GEO
#
B3 <- copy(est)

(extra.lst <- unique(extra$iso3))

dim(est)
est <-
  merge(est, extra[scenario == 'COVID' &
                     hiv == 'a' &
                     measure == 'mort' &
                     iso3 %ni% md.exclude &
                     year %in% 2020:yr, .(
                       iso3,
                       year,
                       mort.rmd = best,
                       mort.rmd.lo = lo,
                       mort.rmd.hi = hi
                     )],
        by = c('iso3', 'year'), all.x = T)
dim(est)

(exclude <- c('GEO', unique(est[g.income=='HIC', iso3])))
sel <-
  est$year %in% 2020:yr &
  est$iso3 %in% extra.lst &
  est$iso3 %ni% unique(c(exclude,  md.exclude))
table(sel)

est[sel, unique(iso3)] # list of extrapolated countries


est[sel, mort := mort.rmd]
est[sel, mort.sd := (mort.rmd.hi - mort.rmd.lo ) / 3.92]
est[sel, mort.lo := mort.rmd.lo]
est[sel, mort.hi := mort.rmd.hi]
est[sel, source.mort := "Regional model"]

est[sel, sum(is.na(mort))]
est[sel, sum(is.na(mort.sd))]

est[year==yr, table(source.mort)]
est[year==yr-1, table(source.mort)]



(rmd.h.lst <- unique(extra[hiv=='pos', iso3]))

dim(est)
est <-
  merge(est, extra[scenario == 'COVID' &
                     hiv == 'pos' &
                     measure == 'mort' &
                     iso3 %ni% md.exclude &
                     year %in% 2020:yr, .(
                       iso3,
                       year,
                       mort.h.rmd = best,
                       mort.h.rmd.lo = lo,
                       mort.h.rmd.hi = hi
                     )],
        by = c('iso3', 'year'), all.x = T)
dim(est)

sel <- est$year %in% 2020:yr & est$iso3 %in% rmd.h.lst & est$iso3 %ni% unique(c(exclude, md.exclude))
table(sel)
est[sel, mort.h := mort.h.rmd]
est[sel, mort.h.sd := (mort.h.rmd.hi - mort.h.rmd.lo ) / 3.92]
est[sel, mort.h.lo := mort.h.rmd.lo]
est[sel, mort.h.hi := mort.h.rmd.hi]

est[sel, sum(is.na(mort.h))]
est[sel, sum(is.na(mort.h.sd))]

est[, sum(is.na(mort)), by=year]
est[, sum(is.na(mort.sd)), by=year]

est[sel, mort.nh := mort - mort.h]
est[sel & mort.sd<mort.h.sd,.(iso3,year,mort,mort.h,mort.nh,mort.sd,mort.h.sd)]
est[sel & mort.sd>mort.h.sd, mort.nh.sd := sqrt(mort.sd^2 - mort.h.sd^2)]
# est[sel & mort.sd<mort.h.sd, mort.nh.sd := mort.sd * mort.nh/mort]



# HIV+ missing
(rmh.lst <- setdiff(extra.lst, rmd.h.lst))
sel <- est$iso3 %in% rmh.lst  & est$iso3 %ni% unique(c(exclude, md.exclude))
table(sel)
(est[sel, unique(iso3)])

est[sel, h.ratio := mort.h / mort]
est[sel & is.na(h.ratio), table(year)]
est[sel, h.ratio := imputeTS::na_locf(h.ratio), by=iso3]
est[sel & is.na(mort.h), mort.h := mort * h.ratio]
est[sel & is.na(mort.h.sd), mort.h.sd := mort.sd * h.ratio]
est[sel & is.na(mort.nh), mort.nh := mort * (1 - h.ratio)]
est[sel & is.na(mort.nh.sd), mort.nh.sd := mort.sd * (1 - h.ratio)]
est[, h.ratio := NULL]




# VR/IHME countries not in previous lists:
# use LOCF
B3b <- copy(est)

(est['FRA',.(iso3,year,mort,mort.nh,e.mort.nh,mort.nh.hat,source.mort)])
(est['CHN',.(iso3,year,mort,mort.nh,e.mort.nh,mort.nh.hat,source.mort)])
(est['ZAF',.(iso3,year,mort,mort.nh,e.mort.nh,mort.nh.hat,source.mort)])
(est['RUS',.(iso3,year,mort,mort.nh,e.mort.nh,mort.nh.hat,source.mort)])
est[, sum(is.na(mort.nh)), by=year]

(nh.miss <- est[year==yr & is.na(mort.nh), iso3])
(vr.lst <- intersect(nh.miss, est[source.mort %in% c('VR','IHME'), unique(iso3)]))

est[iso3 %in% setdiff(vr.lst, 'CHN') & year==2020, mort := NA]
est[iso3 %in% setdiff(vr.lst, 'CHN') & year==2020, mort.sd := NA]
est[iso3 %in% setdiff(vr.lst, 'CHN') & year==2020, mort.h := NA]
est[iso3 %in% setdiff(vr.lst, 'CHN') & year==2020, mort.h.sd := NA]
est[iso3 %in% setdiff(vr.lst, 'CHN') & year==2020, mort.nh := NA]
est[iso3 %in% setdiff(vr.lst, 'CHN') & year==2020, mort.nh.sd := NA]

est[iso3 %in% vr.lst, mort := imputeTS::na_locf(mort), by=iso3]
est[iso3 %in% vr.lst, mort.sd := imputeTS::na_locf(mort.sd), by=iso3]
est[iso3 %in% vr.lst, mort.h := imputeTS::na_locf(mort.h), by=iso3]
est[iso3 %in% vr.lst, mort.h.sd := imputeTS::na_locf(mort.h.sd), by=iso3]
est[iso3 %in% vr.lst, mort.nh := imputeTS::na_locf(mort.nh), by=iso3]
est[iso3 %in% vr.lst, mort.nh.sd := imputeTS::na_locf(mort.nh.sd), by=iso3]

est[iso3 %in% vr.lst & year == 2021, source.mort := "VR/IHME, extrapolated"]
est[iso3 %in% setdiff(vr.lst, 'CHN') & year==2020, source.mort := "VR/IHME, extrapolated"]

# TEMP: Bug fix for ZAF 2020 source because inherited "model" from last year so change
est[iso3 == 'ZAF' & year == 2020, source.mort := "VR/IHME, extrapolated", by=iso3]




# replace LOCF with mort.nh.hat (predicted using logistic regression) 
# in subset of countries with VR/IHME located in AFR & Latin AMR + GEO
# to capture declines
# (keep SDs from LOCF above)
vr.sub <- intersect(vr.lst, est[g.whoregion %in% c('AFR','AMR'), unique(iso3)])
(vr.sub <- setdiff(vr.sub, c('USA','CAN')))
(vr.sub <- c(vr.sub, 'GEO'))
sel <- est$iso3 %in% vr.sub & est$year > 2019
table(sel)

# check series before update
qplot(year, mort.nh, data=est[iso3 %in% vr.sub & year>2009], geom='line') +
  facet_wrap(~iso3, scales='free_y')


est[sel, mort.nh := mort.nh.hat]
est[sel, mort.h := mort.h.hat]
est[sel, mort := mort.h + mort.nh]
est[sel, mort.sd := sqrt(mort.nh.sd^2 + mort.h.sd^2)]
est[sel, source.mort := "VR/IHME, extrapolated"]


# check updated series
qplot(year, mort.nh, data=est[iso3 %in% vr.sub & year>2009], geom='line') +
  facet_wrap(~iso3, scales='free_y')
qplot(year, mort.h, data=est[iso3 %in% vr.sub & year>2009], geom='line') +
  facet_wrap(~iso3, scales='free_y')

est[, sum(is.na(mort.nh)), by=year]



# others = use indirect estimates
#
B4 <- copy(est)

(lst <- est[is.na(mort.nh), unique(iso3)])

sel <- est$iso3 %in% lst & est$year %in% 2020:yr
table(sel)
est[sel, table(year)]

est[sel, mort := e.mort]
est[sel, mort.nh := e.mort.nh]
est[sel, mort.h := e.mort.h]

est[sel, mort.sd := e.mort.sd]
est[sel, mort.nh.sd := e.mort.nh.sd]
est[sel, mort.h.sd := e.mort.h.sd]

est[sel, sum(is.na(mort))]
est[sel, sum(is.na(mort.nh))]
est[sel, sum(is.na(mort.h))]

est[sel, table(source.mort)]
est[sel, source.mort := 'Indirect']

est[, sum(is.na(mort)), by=year]
est[, sum(is.na(mort.h)), by=year]
est[, sum(is.na(mort.nh)), by=year]

est[, sum(is.na(mort.sd)), by=year]
est[, sum(is.na(mort.h.sd)), by=year]
est[, sum(is.na(mort.nh.sd)), by=year]

est[, sum(is.na(source.mort)), by=year]

est[year==2021, table(source.mort)]
est[year==2020, table(source.mort)]
est[year==2019, table(source.mort)]


# # GEO inelegant fix
# # HT: also did this for ETH, PRK and ZAF who were also removed from Nim's model
# (est['GEO',.(iso3,year,mort,mort.nh,e.mort.nh,mort.nh.hat,source.mort)])
# (est['ETH',.(iso3,year,mort,mort.nh,e.mort.nh,mort.nh.hat,source.mort)])
# (est['PRK',.(iso3,year,mort,mort.nh,e.mort.nh,mort.nh.hat,source.mort)])
# (est['ZAF',.(iso3,year,mort,mort.nh,e.mort.nh,mort.nh.hat,source.mort)])

# lst <- c('ETH','PRK','ZAF')
# est[iso3 %in% lst & year>2019, mort.nh := mort.nh.hat]
# est[iso3 %in% lst & year>2019, mort.h := mort.h.hat]
# est[iso3 %in% lst & year>2019, mort := mort.h + mort.nh]
# 
# # START OF HACK - - - - - - - - - - - - - - - - -
# # Use same mort.nh.sd/mort.nh ratio from 2019 for subsequent years for eacg country
# # Hazim's very crude night hack -- donethis way to avoid any inadvertant errors, but should
# # be done in a better way in future
# 
# # 0.05 GEO nh from    est[iso3=='GEO' & year == 2019, .(mort.nh.sd/mort.nh)]
# # 0.21 GEO h from     est[iso3=='GEO' & year == 2019, .(mort.h.sd/mort.h)]
# 
# est[iso3 == 'GEO' & year>2019, mort.nh.sd := mort.nh * 0.05]
# est[iso3 == 'GEO' & year>2019, mort.h.sd := mort.h * 0.21]
# 
# # 0.22 ETH nh from    est[iso3=='ETH' & year == 2019, .(mort.nh.sd/mort.nh)]
# # 0.18 ETH h  from    est[iso3=='ETH' & year == 2019, .(mort.h.sd/mort.h)]
# 
# est[iso3 == 'ETH' & year>2019, mort.nh.sd := mort.nh * 0.22]
# est[iso3 == 'ETH' & year>2019, mort.h.sd := mort.h * 0.18]
# 
# # 0.17 PRK nh from    est[iso3=='PRK' & year == 2019, .(mort.nh.sd/mort.nh)]
# # 0.25 PRK h  from    est[iso3=='PRK' & year == 2019, .(mort.h.sd/mort.h)]
# 
# est[iso3 == 'PRK' & year>2019, mort.nh.sd := mort.nh * 0.17]
# est[iso3 == 'PRK' & year>2019, mort.h.sd := mort.h * 0.25]
# 
# # 0.03 ZAF nh from    est[iso3=='ZAF' & year == 2019, .(mort.nh.sd/mort.nh)]
# # 0.38 ZAF h  from    est[iso3=='ZAF' & year == 2019, .(mort.h.sd/mort.h)]
# 
# est[iso3 == 'ZAF' & year>2019, mort.nh.sd := mort.nh * 0.03]
# est[iso3 == 'ZAF' & year>2019, mort.h.sd := mort.h * 0.38]
# 
# # END OF HACK - - - - - - - - - - - - - - - - -


# est[iso3 %in% md.exclude & year>2019, mort.sd := sqrt(mort.nh.sd^2 + mort.h.sd^2)]



# # replace CHN with VR when communicated
# # locf for 2021 in the meantime ;
# # reuse reported 2020 values
(est['CHN',.(iso3,year,mort,mort.nh,e.mort.nh,mort.nh.hat,source.mort)])

# sel <- est$iso3 == 'CHN' & est$year == 2020
# sel2 <- old$iso3 == 'CHN' & old$year == 2020
# est$mort[sel] <- old$mort[sel2]
# est$mort.sd[sel] <- old$mort.sd[sel2]
# est$mort.h[sel] <- old$mort.h[sel2]
# est$mort.h.sd[sel] <- old$mort.h.sd[sel2]
# est$mort.nh[sel] <- old$mort.nh[sel2]
# est$mort.nh.sd[sel] <- old$mort.nh.sd[sel2]
# 
# sel <- est$iso3 == 'CHN'
# est[sel & year == 2021, mort.nh := NA]
# est[sel & year == 2021, mort.nh.sd := NA]
# est[sel, mort.h := e.mort.h]
# est[sel, mort.h.sd := e.mort.h.sd]
# 
# est[sel, mort.nh := imputeTS::na_locf(mort.nh)]
# est[sel, mort.nh.sd := imputeTS::na_locf(mort.nh.sd)]
# est[sel, mort.h := imputeTS::na_locf(mort.h)]
# est[sel, mort.h.sd := imputeTS::na_locf(mort.h.sd)]
# est[sel, mort := mort.nh + mort.h]
# est[sel, mort.sd := sqrt(mort.nh.sd^2 + mort.h.sd^2)]
# 
# est[sel & year==2021, source.mort := 'Current trends']


# mort.h greater than mort.hiv?
sel <- est$mort.h >= est$mort.hiv & est$mort.h>0 & est$mort.hiv > 0
table(sel)
est[sel, unique(iso3)]
est[sel, table(year)]
est[sel, summary(mort.h/mort.hiv)]
est[!sel & mort.hiv>0, summary(mort.h/mort.hiv)]
# set cap at 60%
est[sel, mort.h := mort.hiv * 0.6]
est[sel, mort.h.sd := mort.h * 0.25]
est[sel, mort := mort.nh + mort.h]
est[sel, mort.sd := sqrt(mort.nh.sd^2 + mort.h.sd^2)]

est[is.na(mort.h.sd), ]
est[is.na(mort.h.sd), mort.h.sd := 0]


# mort.nh greater than mort?
est[mort.nh>mort, .(iso3,year,mort,mort.nh,mort.h)]
est[mort.nh>mort, mort.nh := mort - mort.h]

# checks
(sum(is.na(est$mort)) == 0)
(sum(is.na(est$mort.sd)) == 0)
(sum(is.na(est$mort.nh)) == 0)
(sum(is.na(est$mort.nh.sd)) == 0)
(sum(is.na(est$mort.h)) == 0)
(sum(is.na(est$mort.h.sd)) == 0)

est[, test.AgeB(mort, mort.nh)]
est[mort.hiv>0 & mort.h>0, test.AgeB(mort.hiv, mort.h)]






# add bounds and counts
#
B5 <- copy(est)

sel <- est$mort > 0 & est$mort.sd > 0
table(sel)

est[sel & (mort.sd/1e5)^2 > mort/1e5 * (1-mort/1e5), mort.sd := mort * .2]

out <- vlohi(est$mort[sel] / m, est$mort.sd[sel] / m)

est$mort.lo[sel] <- out[1, ] * m
est$mort.hi[sel] <- out[2, ] * m

sel <- est$mort.sd == 0 & !is.na(est$mort.sd)
table(sel)
est$mort.lo[sel] <- est$mort[sel]
est$mort.hi[sel] <- est$mort[sel]

sel <-
  (est$mort.lo > est$mort) |
  (est$mort.hi < est$mort) &
  (!is.na(est$mort) &
     !is.na(est$mort.lo) & !is.na(est$mort.hi))
table(sel)
# est[sel, .(iso3, year, inc, mort, mort.sd, mort.lo, mort.hi)]
# est[sel, mort.hi := mort + 1.96 * mort.sd]
est[mort.lo == 0, .(iso3, year, mort, mort.sd, mort.lo, mort.hi)]
est[mort > 0, test.bounds(mort, mort.lo, mort.hi)]

est[is.na(mort.lo), mort.lo := 0]
est[is.na(mort.hi), mort.hi := mort + 1.96 * mort.sd]



# HIV-neg
sel <- est$mort.nh > 0 & est$mort.nh.sd > 0
table(sel)
est[sel & (mort.nh.sd/1e5)^2 > mort.nh/1e5 * (1-mort.nh/1e5)]
est[sel & (mort.nh.sd/1e5)^2 > mort.nh/1e5 * (1-mort.nh/1e5), mort.nh.sd := mort.nh * .2]

out <- vlohi(est$mort.nh[sel] / m, est$mort.nh.sd[sel] / m)

est$mort.nh.lo[sel] <- out[1, ] * m
est$mort.nh.hi[sel] <- out[2, ] * m

sel <- est$mort.nh.sd == 0 & !is.na(est$mort.nh.sd)
table(sel)
est$mort.nh.lo[sel] <- est$mort.nh[sel]
est$mort.nh.hi[sel] <- est$mort.nh[sel]

sel <-
  (est$mort.nh.lo > est$mort.nh) |
  (est$mort.nh.hi < est$mort.nh) &
  (!is.na(est$mort.nh) &
     !is.na(est$mort.nh.lo) & !is.na(est$mort.nh.hi))
table(sel)
# est[sel, .(iso3, year, inc, mort.nh, mort.nh.sd, mort.nh.lo, mort.nh.hi)]
# est[sel, mort.nh.hi := mort.nh + 1.96 * mort.nh.sd]
est[mort.nh.lo == 0, .(iso3, year, mort.nh, mort.nh.sd, mort.nh.lo, mort.nh.hi)]
est[mort.nh > 0, test.bounds(mort.nh, mort.nh.lo, mort.nh.hi)]

est[is.na(mort.nh.lo), mort.nh.lo := 0]
est[is.na(mort.nh.hi), mort.nh.hi := mort.nh + 1.96 * mort.nh.sd]



# HIV-pos
sel <- est$mort.h > 0 & est$mort.h.sd > 0 & (est$mort.h.sd/m)^2 >= est$mort.h/m * (1 - est$mort.h/m)
table(sel)
est[sel, .(iso3,mort.h, mort.h.sd)]
est[sel, mort.h.sd := mort.h * .2]

sel <- est$mort.h > 0 & est$mort.h.sd > 0
table(sel)


out <- vlohi(est$mort.h[sel] / m, est$mort.h.sd[sel] / m)

est$mort.h.lo[sel] <- out[1, ] * m
est$mort.h.hi[sel] <- out[2, ] * m

sel <- est$mort.h.sd == 0 & !is.na(est$mort.h.sd)
table(sel)
est$mort.h.lo[sel] <- est$mort.h[sel]
est$mort.h.hi[sel] <- est$mort.h[sel]

sel <-
  (est$mort.h.lo > est$mort.h) |
  (est$mort.h.hi < est$mort.h) &
  (!is.na(est$mort.h) &
     !is.na(est$mort.h.lo) & !is.na(est$mort.h.hi))
table(sel)
# est[sel, .(iso3, year, inc, mort.h, mort.h.sd, mort.h.lo, mort.h.hi)]
# est[sel, mort.h.hi := mort.h + 1.96 * mort.h.sd]
est[mort.h.lo == 0, .(iso3, year, mort.h, mort.h.sd, mort.h.lo, mort.h.hi)]
est[mort.h > 0, test.bounds(mort.h, mort.h.lo, mort.h.hi)]

est[is.na(mort.h.lo), mort.h.lo := 0]
est[is.na(mort.h.hi), mort.h.hi := mort.h + 1.96 * mort.h.sd]

sel <- est$iso3=='MSR' & est$year==yr
est[sel, .(iso3,year,mort,mort.sd,mort.nh,mort.nh.sd,mort.h,mort.h.sd,mort.h.lo,mort.h.hi)]
est[sel, mort.sd := 0]
est[sel, mort.hi := 0]
est[sel, mort.nh.sd := 0]
est[sel, mort.nh.hi := 0]
est[sel, mort.h.sd := 0]
est[sel, mort.h.hi := 0]

est <- within(est, {
  mort.num <- mort * pop / m
  mort.lo.num <- mort.lo * pop / m
  mort.hi.num <- mort.hi * pop / m

  mort.nh.num <- mort.nh * pop / m
  mort.nh.lo.num <- mort.nh.lo * pop / m
  mort.nh.hi.num <- mort.nh.hi * pop / m

  mort.h.num <- mort.h * pop / m
  mort.h.lo.num <- mort.h.lo * pop / m
  mort.h.hi.num <- mort.h.hi * pop / m

  inc.num <- inc * pop / m
  inc.lo.num <- inc.lo * pop / m
  inc.hi.num <- inc.hi * pop / m
  inc.nh.num <- inc.nh * pop / m
  inc.nh.lo.num <- inc.nh.lo * pop / m
  inc.nh.hi.num <- inc.nh.hi * pop / m
  inc.h.num <- inc.h * pop / m
  inc.h.lo.num <- inc.h.lo * pop / m
  inc.h.hi.num <- inc.h.hi * pop / m

})


# checks
#
est[, test.bounds(mort, mort.lo, mort.hi)]
est[, test.bounds(mort.nh, mort.nh.lo, mort.nh.hi)]
est[, test.bounds(mort.h, mort.h.lo, mort.h.hi)]

est[, .(
  sums(mort.nh.num),
  sums(inc.num),
  sums(inc.nh.num),
  sums(mort.nh.num) / sums(inc.nh.num)
), by = year]
old[, .(
  sums(mort.nh.num),
  sums(inc.num),
  sums(inc.nh.num),
  sums(mort.nh.num) / sums(inc.nh.num)
), by = year]

wr <- unique(as.character(est$g.whoregion))




for (i in wr) {
  p <-
    qplot(
      year,
      mort,
      data = subset(est, g.whoregion == i),
      geom = 'line',
      colour = I('blue')
    ) +
    geom_ribbon(
      aes(year, ymin = mort.lo, ymax = mort.hi),
      fill = I('blue'),
      alpha = I(0.4)
    ) +
    geom_line(
      aes(year, mort),
      data = subset(old, g.whoregion == i),
      colour = I('red'),
      linetype = I(2)
    ) +
    facet_wrap(~ iso3, scales = 'free_y') + xlab('') + ylab('Rate per 100,000/year')
  suppressWarnings(print(p))

  suppressWarnings(ggsave(here(
    paste('output/checks/mort_', i, '_compare.pdf', sep = '')
  ),
  width = 14,
  height = 8))
}




for (i in wr) {
  p <-
    qplot(
      year,
      mort.h,
      data = subset(est, g.whoregion == i),
      geom = 'line',
      colour = I('blue')
    ) +
    geom_ribbon(
      aes(year, ymin = mort.h.lo, ymax = mort.h.hi),
      fill = I('blue'),
      alpha = I(0.4)
    ) +
    geom_line(
      aes(year, mort.h),
      data = subset(old, g.whoregion == i),
      colour = I('red'),
      linetype = I(2)
    ) +
    facet_wrap(~ iso3, scales = 'free_y') + xlab('') + ylab('Rate per 100,000/year')
  suppressWarnings(print(p))

  suppressWarnings(ggsave(here(
    paste('output/checks/mort.h_', i, '_compare.pdf', sep = '')
  ),
  width = 14,
  height = 8))
}



for (i in wr) {
  p <-
    qplot(
      year,
      mort.nh,
      data = subset(est, g.whoregion == i),
      geom = 'line',
      colour = I('blue')
    ) +
    geom_ribbon(
      aes(year, ymin = mort.nh.lo, ymax = mort.nh.hi),
      fill = I('blue'),
      alpha = I(0.4)
    ) +
    geom_line(
      aes(year, mort.nh),
      data = subset(old, g.whoregion == i),
      colour = I('red'),
      linetype = I(2)
    ) +
    facet_wrap(~ iso3, scales = 'free_y') + xlab('') + ylab('Rate per 100,000/year')
  suppressWarnings(print(p))

  suppressWarnings(ggsave(here(
    paste('output/checks/mort.nh_', i, '_compare.pdf', sep = '')
  ),
  width = 14,
  height = 8))
}


for (i in wr) {
  p <-
    qplot(
      year,
      mort.nh,
      data = subset(est, g.whoregion == i & year>=yr-5),
      geom = 'line',
      colour = I('blue')
    ) +
    geom_ribbon(
      aes(year, ymin = mort.nh.lo, ymax = mort.nh.hi),
      fill = I('blue'),
      alpha = I(0.4)
    ) +
    geom_line(
      aes(year, mort.nh),
      data = subset(old, g.whoregion == i & year>=yr-5),
      colour = I('red'),
      linetype = I(2)
    ) +
    facet_wrap(~ iso3, scales = 'free_y') + xlab('') + ylab('Rate per 100,000/year')
  suppressWarnings(print(p))

  suppressWarnings(ggsave(here(
    paste('output/checks/mort.nh_', i, '_5y.pdf', sep = '')
  ),
  width = 14,
  height = 8))
}

est[, table(source.mort)]
est[, sum(is.na(source.mort)), by=year]
est[is.na(source.mort) & year==2021, source.mort := "Current trends"]

# Bug fix suggested by Philippe
# (when mort.nh is the same as mort.nh.hat then the method is indirect using CFR, notification and incidence)
est[mort.nh == e.mort.nh & source.mort != "Indirect", source.mort := "Indirect"]

est[, sum(is.na(source.mort))==0]
est[, sum(is.na(source.inc))==0]






# check global aggregates
#
(est[, sum(mort.num), by=year])
(est[iso3 %in% md.lst & iso3 %ni% md.exclude, sum(mort.num), by=year])
(est[, sum(mort.h.num), by=year])
(est[, sum(mort.nh.num), by=year])



# late update from CHN and RUS
# est[iso3=='CHN' & year==yr, source.mort := 'VR']
# est[iso3=='RUS' & year==yr, source.mort := 'VR']


# save
#
save(est, file = here('data/est.rda'))
fwrite(est, file = here(paste0('csv/est_', Sys.Date(), '.csv')))


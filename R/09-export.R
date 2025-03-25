#' ---
#' title: Export files to GTB database
#' author: Philippe Glaziou
#' date: 19/07/2022
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---

#' # Preamble
#' (Last updated: `r Sys.Date()`)
#'
#' Export files to GTB database
#'

library(here)

load(here('data/est.rda'))
load(here('data/global.rda'))
load(here('data/regional.rda'))
load(here('data/hbc.rda'))
# load(here('data/rr.rda'))
# load(here('data/rrg.rda'))
load(here('data/vr.rda'))
load(here('data/cty.rda'))
load(here('data/att.rda'))
load(here('data/bpos.rda'))
load(here('data/model.rda'))


source(here('R/fun.R'))

yr <- 2021
vlohi <- Vectorize(lohi, c('ev', 'sd'))



# save est
#
if (!exists('country', est)) {
  est <-
    merge(est, cty[, .(iso3, country)], by = 'iso3', all.x = TRUE)
  save(est, file = here('data/est.rda'))
  fwrite(est, file = here(paste0('csv/est_', Sys.Date(), '.csv')))
}



# # Main estimates
#
export <- est[, list(
  iso3,
  g.whoregion,
  year,
  pop,
  inc,
  inc.lo,
  inc.hi,
  inc.sd,
  inc.h,
  inc.h.lo,
  inc.h.hi,
  inc.sd,
  tbhiv,
  tbhiv.lo,
  tbhiv.hi,
  tbhiv.sd,
  mort.nh,
  mort.nh.lo,
  mort.nh.hi,
  mort.nh.sd,
  mort.h,
  mort.h.lo,
  mort.h.hi,
  mort.h.sd,
  mort,
  mort.lo,
  mort.hi,
  mort.sd,
  inc.num,
  inc.lo.num,
  inc.hi.num,
  inc.h.num,
  inc.h.lo.num,
  inc.h.hi.num,
  mort.nh.num,
  mort.nh.lo.num,
  mort.nh.hi.num,
  mort.h.num,
  mort.h.lo.num,
  mort.h.hi.num,
  mort.num,
  mort.lo.num,
  mort.hi.num,
  source.inc,
  source.mort
)]

# blank out TBHIV and mortality estimates in PRK
# requested by NTP/PRK in Aug 2020
#
sel <- export$iso3 == 'PRK'
export[sel, tbhiv := NA]
export[sel, tbhiv.lo := NA]
export[sel, tbhiv.hi := NA]

export[sel, inc.h := NA]
export[sel, inc.h.lo := NA]
export[sel, inc.h.hi := NA]
export[sel, inc.h.num := NA]
export[sel, inc.h.lo.num := NA]
export[sel, inc.h.hi.num := NA]

export[sel, mort.h := NA]
export[sel, mort.h.lo := NA]
export[sel, mort.h.hi := NA]
export[sel, mort.h.num := NA]
export[sel, mort.h.lo.num := NA]
export[sel, mort.h.hi.num := NA]

export[sel, mort.nh := NA]
export[sel, mort.nh.lo := NA]
export[sel, mort.nh.hi := NA]
export[sel, mort.nh.num := NA]
export[sel, mort.nh.lo.num := NA]
export[sel, mort.nh.hi.num := NA]

export[sel, mort := NA]
export[sel, mort.lo := NA]
export[sel, mort.hi := NA]
export[sel, mort.num := NA]
export[sel, mort.lo.num := NA]
export[sel, mort.hi.num := NA]


fwrite(export, file = here(paste0('csv/db/db_est_country_', Sys.Date(), '.csv')))

# B+ Pulmonary TB incidence (2019)
#
fwrite(bpos, file = here(paste0('csv/db/db_bactpos_country_', Sys.Date(), '.csv')))


# # CFR
#
cfr <- est[year == yr,
           {
             tmp = divXY(mort, inc, mort.sd, inc.sd)

             list(cfr = tmp[[1]],
                  cfr.sd = tmp[[2]])
           }, by = iso3]

cfr$cfr[cfr$cfr > 1 & !is.na(cfr$cfr)] <- 1

sel <- cfr$cfr > 0 & cfr$cfr < 1 & !is.na(cfr$cfr)
out <- with(cfr[sel], vlohi(cfr, cfr.sd))
cfr$cfr.lo[sel] <- out[1,]
cfr$cfr.hi[sel] <- out[2,]

cfr[!sel]
cfr$cfr[!sel] <- NA
cfr$cfr.sd[!sel] <- NA
cfr$cfr.lo[!sel] <- cfr$cfr.hi[!sel] <- NA

fwrite(cfr, file = here(paste0('csv/db/db_cfr_', Sys.Date(), '.csv')))


# # Attributable cases
#
att2 <- att[, list(
  iso3,
  sex,
  inc.at.hiv.num,
  inc.at.hiv.lo.num,
  inc.at.hiv.hi.num,
  inc.at.dia.num,
  inc.at.dia.lo.num,
  inc.at.dia.hi.num,
  inc.at.alc.num,
  inc.at.alc.lo.num,
  inc.at.alc.hi.num,
  inc.at.smk.num,
  inc.at.smk.lo.num,
  inc.at.smk.hi.num,
  inc.at.und.num,
  inc.at.und.lo.num,
  inc.at.und.hi.num,
  sex
)]


# reshape to long
#
att3 <- melt(att2, id.vars = c(1,2, 4:5), measure.vars = 3)
att4 <- melt(att2, id.vars = c(1,2, 7:8), measure.vars = 6)
att5 <- melt(att2, id.vars = c(1,2, 10:11), measure.vars = 9)
att6 <- melt(att2, id.vars = c(1,2, 13:14), measure.vars = 12)
att7 <- melt(att2, id.vars = c(1,2, 16:17), measure.vars = 15)
for (i in list(att3, att4, att5, att6, att7)) setnames(i, c('iso3', 'sex', 'lo', 'hi', 'risk.factor', 'best'))
att3[, risk.factor := 'hiv']
att4[, risk.factor := 'diabetes']
att5[, risk.factor := 'alcohol']
att6[, risk.factor := 'smoking']
att7[, risk.factor := 'undernutrition']
rf <-
  rbind(att3, att4, att5, att6, att7)
rf[, measure := 'inc']
rf[, unit := 'num']
setkey(rf, iso3)

rf[, age := 'a']
rf[risk.factor %in% c('smoking', 'alcohol'), age := '15+']
rf[risk.factor %in% c('diabetes'), age := '18+']

rf <- merge(rf, cty[, .(iso3, g.whoregion)], by = 'iso3')

sel <- rf$risk.factor %in% c('hiv', 'undernutrition') & rf$sex %in% c('m','f')
table(sel)
rf <- rf[!sel]


fwrite(rf[, .(iso3, year = yr, risk.factor, measure, unit, sex, age, best, lo, hi)],
       file = here(paste0('csv/db/db_inc_risk_factor_country_', Sys.Date(), '.csv')))

save(rf, file = here('data/rf.rda'))


# aggregate rf
#
vglohi <- Vectorize(glohi, c('ev', 'sd'))

rf.regional <- rf[, .(best = as.integer(sums(best)),
                      sd = as.integer(sqrt(sums( ((hi - lo) / 3.92) ^ 2) ))),
                  by = .(g.whoregion, risk.factor, measure, unit, sex, age)]

out <- with(rf.regional, vglohi(best, sd))
rf.regional[, lo := as.integer(out[1,])]
rf.regional[, hi := as.integer(out[2,])]
rf.regional[, sd := NULL]

rf.global <- rf[, .(best = as.integer(sums(best)),
                    sd = as.integer(sqrt(sums( ((hi - lo) / 3.92) ^ 2) ))),
                by = .(risk.factor, measure, unit, sex, age)]

out <- with(rf.global, vglohi(best, sd))
rf.global[, lo := as.integer(out[1,])]
rf.global[, hi := as.integer(out[2,])]
rf.global[, sd := NULL]

rf.regional[, group_type := 'g_whoregion']
rf.regional[, group_name := g.whoregion]
rf.regional[, g.whoregion := NULL]
rf.regional[, year := yr]

rf.global[, group_type := 'global']
rf.global[, group_name := 'global']
rf.global[, year := yr]

cols <- c(
  'group_type',
  'group_name',
  'year',
  'risk.factor',
  'measure',
  'unit',
  'sex',
  'age',
  'best',
  'lo',
  'hi'
)

setcolorder(rf.global, cols)
setcolorder(rf.regional, cols)

save(rf.global, file = here('data/rf.global.rda'))
save(rf.regional, file = here('data/rf.regional.rda'))




# # Aggregates
#
fwrite(global, file = here(paste0('csv/db/db_est_global_est_', Sys.Date(), '.csv')))
fwrite(regional,
       file = here(paste0('csv/db/db_est_regional_est_', Sys.Date(), '.csv')))
fwrite(hbc, file = here(paste0('csv/db/db_est_hbc_est_', Sys.Date(), '.csv')))

fwrite(rf.global,
       file = here(paste0('csv/db/db_inc_risk_factor_global', Sys.Date(), '.csv')))
fwrite(rf.regional,
       file = here(paste0('csv/db/db_inc_risk_factor_regional', Sys.Date(), '.csv')))



# # RR
#
# exportrr <- rr[, .(
#   iso3,
#   year = yr,
#   source_new = source.new,
#   source_drs_year_new = year.new,
#   source_drs_all_areas_new = all.areas.covered.new,
#   e_rr_prop_new = prop.rr.new,
#   e_rr_prop_new_lo = prop.rr.new.lo,
#   e_rr_prop_new_hi = prop.rr.new.hi,
#   e_rr_prop_new_se = prop.rr.new.sd,
#   e_mdr_prop_rr_new = mdr.rr.new,
#   source_ret = source.ret,
#   source_drs_year_ret = year.ret,
#   source_drs_all_areas_ret = all.areas.covered.ret,
#   e_rr_prop_ret = prop.rr.ret,
#   e_rr_prop_ret_lo = prop.rr.ret.lo,
#   e_rr_prop_ret_hi = prop.rr.ret.hi,
#   e_rr_prop_ret_se = prop.rr.ret.sd,
#   e_mdr_prop_rr_ret = mdr.rr.ret,
#   e_inc_rr_num = inc.rr.num,
#   e_inc_rr_num_lo = inc.rr.lo.num,
#   e_inc_rr_num_hi = inc.rr.hi.num,
#   e_mdr_prop_rr = inc.mdr / inc.rr
# )]
#
# fwrite(exportrr, file = here(paste0('csv/db/db_dr_country_', Sys.Date(), '.csv')))
#
#
# # DR aggregates
# #
# xp1 <- rrg[, .(
#   region,
#   year = yr,
#   e_rr_prop_new = prop.rr.new,
#   e_rr_prop_new_lo = prop.rr.new.lo,
#   e_rr_prop_new_hi = prop.rr.new.hi,
#   e_rr_prop_new_se = (prop.rr.new.hi - prop.rr.new.lo) / 3.92,
#   e_rr_prop_ret = prop.rr.ret,
#   e_rr_prop_ret_lo = prop.rr.ret.lo,
#   e_rr_prop_ret_hi = prop.rr.ret.hi,
#   e_rr_prop_ret_se = (prop.rr.ret.hi - prop.rr.ret.lo) / 3.92,
#   e_inc_rr_num = inc.rr.num,
#   e_inc_rr_num_lo = inc.rr.num.lo,
#   e_inc_rr_num_hi = inc.rr.num.hi,
#   e_mdr_prop_rr = pctMDR / 100
# )]
#
# xp2 <- xp1[2:7]
# xp3 <- xp1[8]
# xp4 <- xp1[1]
#
# xp2[, group_type := 'g_whoregion']
# xp2[, group_name := region]
# xp2[, region := NULL]
#
# xp3[, group_type := 'global']
# xp3[, group_name := 'global']
# xp3[, region := NULL]
#
# xp <- rbind(xp2, xp3)
# nm <- names(xp)
# setcolorder(xp, c(nm[14:15], nm[1:13]))
#
# fwrite(xp, file = here(paste0('csv/db/db_dr_group_', Sys.Date(), '.csv')))


# dynamic model
#
fwrite(model, file = here(paste0('csv/db/db_model_', Sys.Date(), '.csv')))




# Export to spectrum
#
y0 <- 2000   # first year of estimates' series
ys <- yr - 5 # period to determine forecasts
m <-
  1e5     # constant to convert rates per 100000 in rates per individual


tospectrum <- est[, list(
  iso3,
  year,
  g.whoregion,
  e.pop.num = pop,
  newrel.num = newinc * pop / m,
  inc.num = inc * pop / m,
  inc.lo.num = inc.lo * pop / m,
  inc.hi.num = inc.hi * pop / m,
  inc.sd.num = inc.sd * pop / m,
  newrel = newinc,
  inc,
  inc.lo,
  inc.hi,
  inc.sd,
  source.inc,
  c.newinc,
  c.notified,
  tbhiv.routine,
  hivtest.coverage,
  tbhiv.surv,
  tbhiv.surv.lo,
  tbhiv.surv.hi,
  tbhiv.sentin,
  tbhiv.sentin.lo,
  tbhiv.sentin.hi,
  hivtest.p,
  hivtest.f,
  hivtest.pos.p,
  hivtest.pos.f,
  hiv.art.p,
  hiv.art.f,
  newrel.art,
  newrel.hivtest,
  newrel.hivpos,
  tbhiv.routine,
  tbhiv.routine.sd
)]

fwrite(tospectrum,
       file = here(paste0('csv/spectrum/tospectrum_', Sys.Date(), '.csv')))


# # VR data -> Carel
# load(here('Rdata/vr.Rdata'))
# tospectrum.vr <- vr[source.mort %in% c('Sample VR', 'VR', 'VA') & keep.vr & !is.na(tb.adj),
#                     .(iso3, year, tb.adj, tb.adj.sd)]
#
# write.csv(tospectrum.vr, file='output/tospectrumVR.csv', row.names=FALSE)

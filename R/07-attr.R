#' ---
#' title: Attributable cases by country (SDG)
#' author: Philippe Glaziou
#' date: 19/07/2022
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
# # Attributable cases by country, 2019 (SDG)
#
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(propagate))
library(here)

source(here('R/fun.R'))
m <- 1e5
yr <- 2021
yrsplit <- 2020

load(here('data/est.rda'))
load(here('data/sdg.rda'))
load(here('data/pop.rda'))



# check timestamps and use most recent
#

incs <- fread(here('disaggregation/output/incsplits/data/A_country_incidence_disaggregated_age_sex_rate_2021-07-15.csv'))

incs <- incs[, c(16, 24:32) := lapply(.SD, function(x)as.integer(gsub(' ','', x))), .SDcols = c(16, 24:32)]
incs[, V1 := NULL]
setkey(incs, iso3)
# save(incs, file = 'Rdata/incs.Rdata')

age <- fread(here('disaggregation/output/incsplits/data/db_estimates_country_2021-07-15.csv'))
setkey(age, iso3)
age[, analysed := NULL]
# save(age, file = 'Rdata/age.rdata')




# returns the last non missing value of a vector or NA
#
lastv <-
  function(x)
    return(ifelse(any(!is.na(x)), last(x[!is.na(x)]), last(x)))

# attributable fraction, propagating errors
#
paf2 <- function(pe, rr) {
  # @param pe proportion exposed: c(mean, sd)
  # @param rr relative risk: c(mean, sd)
  # @export
  require(propagate)

  DT <- cbind(pe, rr)
  EXPR <- expression(pe * (rr - 1) / (pe * (rr - 1) + 1))
  out <-
    propagate(
      expr = EXPR,
      data = DT,
      type = 'stat',
      do.sim = F,
      second.order = T
    )
  return(out)
}

# vectorized glohi
#
vglohi <- Vectorize(glohi, c('ev', 'sd'))


# # Risk Ratios (mean, sd)
#
rr.alc <-
  c(3.33, (5.19 - 2.14) / 3.92) # alcohol, Eur Resp Dis 2017
rr.dia <-
  c(1.5, (1.76 - 1.28) / 3.92) # diabetes, Trop Med Int Health 2018
rr.und <- c(3.2, 0.2 / 3.92) # under-nourishment, GTB 2018
rr.smk <- c(1.57, (2.1 - 1.18) / 3.92) # smoking, GTB 2019



# # Merge data sources
#
dta <-
  est[year == yr, .(iso3, inc, inc.lo, inc.hi, inc.sd, hiv, hiv.sd, irr, irr.sd)]

dta2 <-
  sdg[ind %in% c('alcohol',
                 'diabetes',
                 'undernutrition',
                 'smoking'), .(iso3, year, sex, ind, value, sd = (hi - lo)/3.92)]


dta3 <-
  reshape(
    dta2,
    idvar = c("iso3", "year", "sex"),
    timevar = "ind",
    direction = "wide"
  ) %>%
  filter(year <= yr)
names(dta3) <- gsub('value.', '', names(dta3))
setkey(dta3, iso3)




dta4 <- dta3[, .(
  diabetes = lastv(diabetes) / 100,
  alcohol = lastv(alcohol) / 100,
  undernutrition = lastv(undernutrition) / 100,
  smoking = lastv(smoking) / 100,
  diabetes.sd = lastv(sd.diabetes) / 100,
  alcohol.sd = lastv(sd.alcohol) / 100,
  undernutrition.sd = lastv(sd.undernutrition) / 100,
  smoking.sd = lastv(sd.smoking) / 100
),
by = .(iso3, sex)]


dta5 <-
  merge(dta4, incs[, .(iso3,
                       inc.m,
                       inc.m.sd,
                       inc.f,
                       inc.f.sd,
                       inc.15plus,
                       inc.15plus.sd,
                       inc.15plus.m,
                       inc.15plus.m.sd,
                       inc.15plus.f,
                       inc.15plus.f.sd)], by = 'iso3')
dim(dta4)
dim(dta5)

dta5 <-
  merge(dta5, pop[year==yrsplit, .(iso3,
                       pop = e.pop.num,
                       pop.m = e.pop.m04 + e.pop.m15plus,
                       pop.f = e.pop.f04 + e.pop.f15plus,
                       adults.m = e.pop.m15plus,
                       adults.f = e.pop.f15plus,
                       adults = e.pop.15plus)], by = 'iso3')
(dim(dta5))

# pop18 <- pop[year==yrsplit, .(iso3,
#   adults18.m = e.pop.m15plus - .3 * e.pop.m1524,
#   adults18.f = e.pop.f15plus - .3 * e.pop.f1524,
#   adults18 = e.pop.15plus - .3 * (e.pop.m1524 + e.pop.f1524)
# )]
#
# dta5 <- merge(dta5, pop18, by = 'iso3')


# inelegant reshaping
dta5[sex=='m', pop := pop.m]
dta5[sex=='f', pop := pop.f]
dta5[, pop.m := NULL]
dta5[, pop.f := NULL]

dta5[sex=='m', adults := adults.m]
dta5[sex=='f', adults := adults.f]
dta5[, adults.m := NULL]
dta5[, adults.f := NULL]

# dta5[sex=='m', adults18 := adults18.m]
# dta5[sex=='f', adults18 := adults18.f]
# dta5[, adults18.m := NULL]
# dta5[, adults18.f := NULL]

dta5[sex=='m', inc.15plus := inc.15plus.m]
dta5[sex=='f', inc.15plus := inc.15plus.f]
dta5[, inc.15plus.m := NULL]
dta5[, inc.15plus.f := NULL]

dta5[sex=='m', inc.15plus.sd := inc.15plus.m.sd]
dta5[sex=='f', inc.15plus.sd := inc.15plus.f.sd]
dta5[, inc.15plus.m.sd := NULL]
dta5[, inc.15plus.f.sd := NULL]


ac <- merge(dta5, dta, by = 'iso3', all.x = TRUE)

ac[sex=='m', inc := inc.m]
ac[sex=='f', inc := inc.f]
ac[, inc.m := NULL]
ac[, inc.f := NULL]

ac[sex=='m', inc.sd := inc.m.sd]
ac[sex=='f', inc.sd := inc.f.sd]
ac[, inc.m.sd := NULL]
ac[, inc.f.sd := NULL]

ac[sex %ni% 'a', hiv := NA]
ac[sex %ni% 'a', hiv.sd := NA]
ac[sex %ni% 'a', irr := NA]
ac[sex %ni% 'a', irr.sd := NA]



# # Attributable fractions
#

# HIV
#
out.hiv <- ac[, {
  tmp = paf2(c(hiv, hiv.sd), c(irr, irr.sd))$prop
  list(paf.hiv = tmp[2],
       paf.hiv.sd = tmp[4])
}, by = .(iso3, sex)]

out.hiv[paf.hiv<0, paf.hiv := 0]


# diabetes
#
out.dia <- ac[, {
  tmp = paf2(c(diabetes, diabetes.sd), rr.dia)$prop
  list(paf.dia = tmp[2],
       paf.dia.sd = tmp[4])
}, by = .(iso3, sex)]


# alcohol
#
out.alc <- ac[, {
  tmp = paf2(c(alcohol, alcohol.sd), rr.alc)$prop
  list(paf.alc = tmp[2],
       paf.alc.sd = tmp[4])
}, by = .(iso3, sex)]

out.alc[paf.alc<0, paf.alc := 0]


# smoking
#
out.smk <- ac[, {
  tmp = paf2(c(smoking, smoking.sd), rr.smk)$prop
  list(paf.smk = tmp[2],
       paf.smk.sd = tmp[4])
}, by = .(iso3, sex)]

out.smk[paf.smk<0, paf.smk := 0]


# undernutrition
#
out.und <- ac[, {
  tmp = paf2(c(undernutrition, undernutrition.sd), rr.und)$prop
  list(paf.und = tmp[2],
       paf.und.sd = tmp[4])
}, by = .(iso3, sex)]


ac <-
  cbind(
    ac,
    out.hiv[, -1],
    out.dia[, -1],
    out.alc[, -1],
    out.smk[, -1],
    out.und[, -1]
  )




# # Attributable incidence rate
#
inc.hiv <-
  ac[, prodXY(inc, paf.hiv, inc.sd, paf.hiv.sd), by = iso3]
setnames(inc.hiv, c('iso3', 'inc.at.hiv', 'inc.at.hiv.sd'))
# inc.hiv[, inc.at.hiv.sd := inc.at.hiv.sd]

inc.dia <-
  ac[, prodXY(inc.15plus, paf.dia, inc.15plus.sd, paf.dia.sd), by = iso3]
setnames(inc.dia, c('iso3', 'inc.at.dia', 'inc.at.dia.sd'))
# inc.dia[, inc.at.dia.sd := inc.at.dia.sd]

inc.alc <-
  ac[, prodXY(inc.15plus, paf.alc, inc.15plus.sd, paf.alc.sd), by = iso3]
setnames(inc.alc, c('iso3', 'inc.at.alc', 'inc.at.alc.sd'))
# inc.alc[, inc.at.alc.sd := inc.at.alc.sd]

inc.smk <-
  ac[, prodXY(inc.15plus, paf.smk, inc.15plus.sd, paf.smk.sd), by = iso3]
setnames(inc.smk, c('iso3', 'inc.at.smk', 'inc.at.smk.sd'))
# inc.smk[, inc.at.smk.sd := inc.at.smk.sd]

inc.und <-
  ac[, prodXY(inc, paf.und, inc.sd, paf.und.sd), by = iso3]
setnames(inc.und, c('iso3', 'inc.at.und', 'inc.at.und.sd'))
# inc.und[, inc.at.und.sd := inc.at.und.sd]

ac <-
  cbind(
    ac,
    inc.hiv[, -1],
    inc.dia[, -1],
    inc.alc[, -1],
    inc.smk[, -1],
    inc.und[, -1]
  )




# # Attributable cases
#
m <- 1e5

ac[, inc.at.hiv.num := inc.at.hiv * pop / m]
ac[, inc.at.dia.num := inc.at.dia * adults / m]
ac[, inc.at.alc.num := inc.at.alc * adults / m]
ac[, inc.at.smk.num := inc.at.smk * adults / m]
ac[, inc.at.und.num := inc.at.und * pop / m]




# bounds
#
sel <- !is.na(ac$inc.at.hiv)
out <- ac[sel, vglohi(inc.at.hiv.num, inc.at.hiv.sd * pop / m)]
ac[sel, inc.at.hiv.lo.num := out[1,]]
ac[sel, inc.at.hiv.hi.num := out[2,]]


sel <- !is.na(ac$inc.at.dia) & ac$inc.at.dia.sd > 0
out <- ac[sel, vglohi(inc.at.dia.num, inc.at.dia.sd * pop / m)]
ac[sel, inc.at.dia.lo.num := out[1,]]
ac[sel, inc.at.dia.hi.num := out[2,]]


sel <- !is.na(ac$inc.at.alc) & ac$inc.at.alc.sd > 0
out <- ac[sel, vglohi(inc.at.alc.num, inc.at.alc.sd * pop / m)]
ac[sel, inc.at.alc.lo.num := out[1,]]
ac[sel, inc.at.alc.hi.num := out[2,]]


sel <- !is.na(ac$inc.at.smk) & ac$inc.at.smk.sd > 0
out <- ac[sel, vglohi(inc.at.smk.num, inc.at.smk.sd * pop / m)]
ac[sel, inc.at.smk.lo.num := out[1,]]
ac[sel, inc.at.smk.hi.num := out[2,]]


sel <- !is.na(ac$inc.at.und) & ac$inc.at.und.sd > 0
out <- ac[sel, vglohi(inc.at.und.num, inc.at.und.sd * pop / m)]
ac[sel, inc.at.und.lo.num := out[1,]]
ac[sel, inc.at.und.hi.num := out[2,]]

ac[, 48:62 := lapply(.SD, function(x)
  as.integer(x)), .SDcols = 48:62]

ac[inc.at.smk.hi.num < inc.at.smk.num, inc.at.smk.hi.num := as.integer(3 * inc.at.smk.num)]



# checks
#
ac[inc.at.dia.num>0, test.bounds(inc.at.dia.num, inc.at.dia.lo.num, inc.at.dia.hi.num)]
ac[inc.at.alc.num>0, test.bounds(inc.at.alc.num, inc.at.alc.lo.num, inc.at.alc.hi.num)]
ac[inc.at.smk.num>0, test.bounds(inc.at.smk.num, inc.at.smk.lo.num, inc.at.smk.hi.num)]
ac[inc.at.hiv.num>0, test.bounds(inc.at.hiv.num, inc.at.hiv.lo.num, inc.at.hiv.hi.num)]
ac[inc.at.und.num>0, test.bounds(inc.at.und.num, inc.at.und.lo.num, inc.at.und.hi.num)]

att <- copy(ac)



# check aggregates
#
att[sex=='a', sums(inc.at.und.num)]
att[sex=='a', sums(inc.at.dia.num)]
att[sex=='a', sums(inc.at.smk.num)]
att[sex=='a', sums(inc.at.hiv.num)]
att[sex=='a', sums(inc.at.alc.num)]


rm(dta, dta2, dta3, dta4, dta5)



# save
#
save(att, file = here('data/att.rda'))
fwrite(att, file = here(paste0('csv/attributable_', Sys.Date(), '.csv')))





# compare with last year's aggregates
#
# load(here('../gtb2019/Rdata/ac.Rdata'))
# ac[, sums(inc.at.und.num)]
# ac[, sums(inc.at.dia.num)]
# ac[, sums(inc.at.smk.num)]
# ac[, sums(inc.at.hiv.num)]
# ac[, sums(inc.at.alc.num)]
# rm(ac)




#' ---
#' title: Aggregates GTB 2020
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
# # Aggregates of TB estimates
#
library(data.table)
library(here)


rm(list=ls())

load(here('data/est.rda'))
load(here('data/old.rda'))
load(here('data/tb.rda'))
load(here('data/cty.rda'))

source(here('R/fun.R'))

# vectorized lohi
#
vlohi <- Vectorize(lohi, c('ev', 'sd'))

m <- 1e5
yr <- 2020


# global incidence
#
global.inc <-
  est[, addXY(inc / m, r.sd = inc.sd / m, weights = pop), by = year]
setnames(
  global.inc,
  c(
    "r",
    "r.lo",
    "r.hi",
    "r.sd",
    "r.num",
    "r.lo.num",
    "r.hi.num",
    "pop"
  ),
  c(
    'inc',
    'inc.lo',
    'inc.hi',
    'inc.sd',
    'inc.num',
    'inc.lo.num',
    'inc.hi.num',
    'pop'
  )
)
global.inc <-
  cbind(global.inc[, -(2:5), with = FALSE] , global.inc[, 2:5, with = FALSE] * m)

# global incidence HIV-
#
global.inc.nh <-
  est[, addXY(inc.nh / m, r.sd = inc.nh.sd / m, weights = pop), by = year]
setnames(
  global.inc.nh,
  c(
    "r",
    "r.lo",
    "r.hi",
    "r.sd",
    "r.num",
    "r.lo.num",
    "r.hi.num",
    "pop"
  ),
  c(
    'inc.nh',
    'inc.nh.lo',
    'inc.nh.hi',
    'inc.nh.sd',
    'inc.nh.num',
    'inc.nh.lo.num',
    'inc.nh.hi.num',
    'pop'
  )
)
global.inc.nh <-
  cbind(global.inc.nh[, -(2:5), with = FALSE] , global.inc.nh[, 2:5, with =
                                                                FALSE] * m)

# global incidence HIV+
#
global.inc.h <-
  est[, addXY(inc.h / m, r.sd = inc.h.sd / m, weights = pop), by = year]
setnames(
  global.inc.h,
  c(
    "r",
    "r.lo",
    "r.hi",
    "r.sd",
    "r.num",
    "r.lo.num",
    "r.hi.num",
    "pop"
  ),
  c(
    'inc.h',
    'inc.h.lo',
    'inc.h.hi',
    'inc.h.sd',
    'inc.h.num',
    'inc.h.lo.num',
    'inc.h.hi.num',
    'pop'
  )
)
global.inc.h <-
  cbind(global.inc.h[, -(2:5), with = FALSE] , global.inc.h[, 2:5, with =
                                                              FALSE] * m)


# global mortality HIV-neg
#
global.mort.nh <-
  est[, addXY(mort.nh / m, r.sd = mort.nh.sd / m, weights = pop), by = year]
setnames(
  global.mort.nh,
  c(
    "r",
    "r.lo",
    "r.hi",
    "r.sd",
    "r.num",
    "r.lo.num",
    "r.hi.num",
    "pop"
  ),
  c(
    'mort.nh',
    'mort.nh.lo',
    'mort.nh.hi',
    'mort.nh.sd',
    'mort.nh.num',
    'mort.nh.lo.num',
    'mort.nh.hi.num',
    'pop'
  )
)
global.mort.nh <-
  cbind(global.mort.nh[, -(2:5), with = FALSE] , global.mort.nh[, 2:5, with =
                                                                  FALSE] * m)


# global mortality HIV-pos
#
global.mort.h <-
  est[, addXY(mort.h / m, r.sd = mort.h.sd / m, weights = pop), by = year]
setnames(
  global.mort.h,
  c(
    "r",
    "r.lo",
    "r.hi",
    "r.sd",
    "r.num",
    "r.lo.num",
    "r.hi.num",
    "pop"
  ),
  c(
    'mort.h',
    'mort.h.lo',
    'mort.h.hi',
    'mort.h.sd',
    'mort.h.num',
    'mort.h.lo.num',
    'mort.h.hi.num',
    'pop'
  )
)
global.mort.h <-
  cbind(global.mort.h[, -(2:5), with = FALSE] , global.mort.h[, 2:5, with =
                                                                FALSE] * m)


# global total mortality
#
global.mort <-
  est[, addXY(mort / m, r.sd = mort.sd / m, weights = pop), by = year]
setnames(
  global.mort,
  c(
    "r",
    "r.lo",
    "r.hi",
    "r.sd",
    "r.num",
    "r.lo.num",
    "r.hi.num",
    "pop"
  ),
  c(
    'mort',
    'mort.lo',
    'mort.hi',
    'mort.sd',
    'mort.num',
    'mort.lo.num',
    'mort.hi.num',
    'pop'
  )
)
global.mort <-
  cbind(global.mort[, -(2:5), with = FALSE] , global.mort[, 2:5, with = FALSE] * m)


# put the whole global thing together
#
rg <- c(2:4, 6:9)
global <-
  cbind(
    global.inc,
    global.inc.nh[, rg, with = FALSE],
    global.inc.h[, rg, with = FALSE],
    global.mort.nh[, rg, with = FALSE],
    global.mort.h[, rg, with = FALSE],
    global.mort[, rg, with = FALSE]
  )

# add tbhiv
#
out <- global[, divXY(inc.h, inc, inc.h.sd, inc.sd)]
out2 <- vlohi(out[[1]], out[[2]])
global$tbhiv <- out[[1]]
global$tbhiv.lo <- out2[1,]
global$tbhiv.hi <- out2[2,]
global[, tbhiv.sd := out[[2]]]


# add c.newinc
#
out <- tb[year > 1999, sum(c.newinc, na.rm = TRUE), by = year]

global$c.newinc <- out$V1
global$newinc <- global$c.newinc * m / global$pop


# add CFR
#
out <- global[, divXY(mort, inc, mort.sd, inc.sd)]
out2 <- vlohi(out[[1]], out[[2]])
global$cfr <- out[[1]]
global$cfr.lo <- out2[1,]
global$cfr.hi <- out2[2,]
global[, cfr.sd := out[[2]]]


save(global, file = here('data/global.rda'))
fwrite(global, file = here(paste0('csv/global_', Sys.Date(), '.csv')))




# Aggregates by WHO region
#
# regional incidence
#
regional.inc <-
  est[, addXY(inc / m, r.sd = inc.sd / m, weights = pop), by = c("g.whoregion", "year")]
setnames(
  regional.inc,
  c(
    "r",
    "r.lo",
    "r.hi",
    "r.sd",
    "r.num",
    "r.lo.num",
    "r.hi.num",
    "pop"
  ),
  c(
    'inc',
    'inc.lo',
    'inc.hi',
    'inc.sd',
    'inc.num',
    'inc.lo.num',
    'inc.hi.num',
    'pop'
  )
)
regional.inc <-
  cbind(regional.inc[, -(3:6), with = FALSE] , regional.inc[, 3:6, with =
                                                              FALSE] * m)

regional.inc.nh <-
  est[, addXY(inc.nh / m, r.sd = inc.nh.sd / m, weights = pop), by = c("g.whoregion", "year")]
setnames(
  regional.inc.nh,
  c(
    "r",
    "r.lo",
    "r.hi",
    "r.sd",
    "r.num",
    "r.lo.num",
    "r.hi.num",
    "pop"
  ),
  c(
    'inc.nh',
    'inc.nh.lo',
    'inc.nh.hi',
    'inc.nh.sd',
    'inc.nh.num',
    'inc.nh.lo.num',
    'inc.nh.hi.num',
    'pop'
  )
)
regional.inc.nh <-
  cbind(regional.inc.nh[, -(3:6), with = FALSE] , regional.inc.nh[, 3:6, with =
                                                                    FALSE] * m)

regional.inc.h <-
  est[, addXY(inc.h / m, r.sd = inc.h.sd / m, weights = pop), by = c("g.whoregion", "year")]
setnames(
  regional.inc.h,
  c(
    "r",
    "r.lo",
    "r.hi",
    "r.sd",
    "r.num",
    "r.lo.num",
    "r.hi.num",
    "pop"
  ),
  c(
    'inc.h',
    'inc.h.lo',
    'inc.h.hi',
    'inc.h.sd',
    'inc.h.num',
    'inc.h.lo.num',
    'inc.h.hi.num',
    'pop'
  )
)
regional.inc.h <-
  cbind(regional.inc.h[, -(3:6), with = FALSE] , regional.inc.h[, 3:6, with =
                                                                  FALSE] * m)


# regional mortality HIV-neg
#
regional.mort.nh <-
  est[, addXY(mort.nh / m, r.sd = mort.nh.sd / m, weights = pop), by = c("g.whoregion", "year")]
setnames(
  regional.mort.nh,
  c(
    "r",
    "r.lo",
    "r.hi",
    "r.sd",
    "r.num",
    "r.lo.num",
    "r.hi.num",
    "pop"
  ),
  c(
    'mort.nh',
    'mort.nh.lo',
    'mort.nh.hi',
    'mort.nh.sd',
    'mort.nh.num',
    'mort.nh.lo.num',
    'mort.nh.hi.num',
    'pop'
  )
)
regional.mort.nh <-
  cbind(regional.mort.nh[, -(3:6), with = FALSE] , regional.mort.nh[, 3:6, with =
                                                                      FALSE] * m)


# regional mortality HIV-pos
#
regional.mort.h <-
  est[, addXY(mort.h / m, r.sd = mort.h.sd / m, weights = pop), by = c("g.whoregion", "year")]
setnames(
  regional.mort.h,
  c(
    "r",
    "r.lo",
    "r.hi",
    "r.sd",
    "r.num",
    "r.lo.num",
    "r.hi.num",
    "pop"
  ),
  c(
    'mort.h',
    'mort.h.lo',
    'mort.h.hi',
    'mort.h.sd',
    'mort.h.num',
    'mort.h.lo.num',
    'mort.h.hi.num',
    'pop'
  )
)
regional.mort.h <-
  cbind(regional.mort.h[, -(3:6), with = FALSE] , regional.mort.h[, 3:6, with =
                                                                    FALSE] * m)


# regional total mortality
#
regional.mort <-
  est[, addXY(mort / m, r.sd = mort.sd / m, weights = pop), by = c("g.whoregion", "year")]
setnames(
  regional.mort,
  c(
    "r",
    "r.lo",
    "r.hi",
    "r.sd",
    "r.num",
    "r.lo.num",
    "r.hi.num",
    "pop"
  ),
  c(
    'mort',
    'mort.lo',
    'mort.hi',
    'mort.sd',
    'mort.num',
    'mort.lo.num',
    'mort.hi.num',
    'pop'
  )
)
regional.mort <-
  cbind(regional.mort[, -(3:6), with = FALSE] , regional.mort[, 3:6, with =
                                                                FALSE] * m)


# put the whole regional thing together
#
rg <- c(3:5, 7:10)
regional <-
  cbind(
    regional.inc,
    regional.inc.nh[, rg, with = FALSE],
    regional.inc.h[, rg, with = FALSE],
    regional.mort.nh[, rg, with = FALSE],
    regional.mort.h[, rg, with = FALSE],
    regional.mort[, rg, with = FALSE]
  )

# add tbhiv
#
out <- regional[, divXY(inc.h, inc, inc.h.sd, inc.sd)]
out2 <- vlohi(out$mean, out$sd)
regional$tbhiv <- out[[1]]
regional$tbhiv.lo <- out2[1,]
regional$tbhiv.hi <- out2[2,]
regional$tbhiv.sd <- out[[2]]




# add c.newinc
#
out <-
  tb[year > 1999, sum(c.newinc, na.rm = TRUE), by = c("g.whoregion", "year")]
regional$c.newinc <- out$V1
regional$newinc <- regional$c.newinc * m / regional$pop


# add CFR
#
out <- regional[, divXY(mort, inc, mort.sd, inc.sd)]
out2 <- vlohi(out[[1]], out[[2]])
regional$cfr <- out[[1]]
regional$cfr.lo <- out2[1,]
regional$cfr.hi <- out2[2,]
regional[, cfr.sd := out[[2]]]

save(regional, file = here('data/regional.rda'))
fwrite(regional, file = here(paste0('csv/regional_', Sys.Date(), '.csv')))





# # Aggregates in HBCs
#
select <- est$g.hbc == TRUE

# hbc incidence
#
hbc.inc <-
  est[select][, addXY(inc / m, r.sd = inc.sd / m, weights = pop), by = year]
setnames(
  hbc.inc,
  c(
    "r",
    "r.lo",
    "r.hi",
    "r.sd",
    "r.num",
    "r.lo.num",
    "r.hi.num",
    "pop"
  ),
  c(
    'inc',
    'inc.lo',
    'inc.hi',
    'inc.sd',
    'inc.num',
    'inc.lo.num',
    'inc.hi.num',
    'pop'
  )
)
hbc.inc <-
  cbind(hbc.inc[, -(2:5), with = FALSE] , hbc.inc[, 2:5, with = FALSE] * m)

# hbc incidence HIV+
#
hbc.inc.h <-
  est[select][, addXY(inc.h / m, r.sd = inc.h.sd / m, weights = pop), by =
                year]
setnames(
  hbc.inc.h,
  c(
    "r",
    "r.lo",
    "r.hi",
    "r.sd",
    "r.num",
    "r.lo.num",
    "r.hi.num",
    "pop"
  ),
  c(
    'inc.h',
    'inc.h.lo',
    'inc.h.hi',
    'inc.h.sd',
    'inc.h.num',
    'inc.h.lo.num',
    'inc.h.hi.num',
    'pop'
  )
)
hbc.inc.h <-
  cbind(hbc.inc.h[, -(2:5), with = FALSE] , hbc.inc.h[, 2:5, with = FALSE] * m)


# hbc mortality HIV-neg
#
hbc.mort.nh <-
  est[select][, addXY(mort.nh / m, r.sd = mort.nh.sd / m, weights = pop), by =
                year]
setnames(
  hbc.mort.nh,
  c(
    "r",
    "r.lo",
    "r.hi",
    "r.sd",
    "r.num",
    "r.lo.num",
    "r.hi.num",
    "pop"
  ),
  c(
    'mort.nh',
    'mort.nh.lo',
    'mort.nh.hi',
    'mort.nh.sd',
    'mort.nh.num',
    'mort.nh.lo.num',
    'mort.nh.hi.num',
    'pop'
  )
)
hbc.mort.nh <-
  cbind(hbc.mort.nh[, -(2:5), with = FALSE] , hbc.mort.nh[, 2:5, with = FALSE] * m)


# hbc mortality HIV-pos
#
hbc.mort.h <-
  est[select][, addXY(mort.h / m, r.sd = mort.h.sd / m, weights = pop), by =
                year]
setnames(
  hbc.mort.h,
  c(
    "r",
    "r.lo",
    "r.hi",
    "r.sd",
    "r.num",
    "r.lo.num",
    "r.hi.num",
    "pop"
  ),
  c(
    'mort.h',
    'mort.h.lo',
    'mort.h.hi',
    'mort.h.sd',
    'mort.h.num',
    'mort.h.lo.num',
    'mort.h.hi.num',
    'pop'
  )
)
hbc.mort.h <-
  cbind(hbc.mort.h[, -(2:5), with = FALSE] , hbc.mort.h[, 2:5, with = FALSE] * m)


# hbc total mortality
#
hbc.mort <-
  est[select][, addXY(mort / m, r.sd = mort.sd / m, weights = pop), by =
                year]
setnames(
  hbc.mort,
  c(
    "r",
    "r.lo",
    "r.hi",
    "r.sd",
    "r.num",
    "r.lo.num",
    "r.hi.num",
    "pop"
  ),
  c(
    'mort',
    'mort.lo',
    'mort.hi',
    'mort.sd',
    'mort.num',
    'mort.lo.num',
    'mort.hi.num',
    'pop'
  )
)
hbc.mort <-
  cbind(hbc.mort[, -(2:5), with = FALSE] , hbc.mort[, 2:5, with = FALSE] * m)


# put the whole hbc thing together
#
rg <- c(2:4, 6:9)
hbc <- cbind(hbc.inc, hbc.inc.h[, rg, with = FALSE],
             hbc.mort.nh[, rg, with = FALSE], hbc.mort.h[, rg, with = FALSE],
             hbc.mort[, rg, with = FALSE])

# add tbhiv
#
out <- hbc[, divXY(inc.h, inc, inc.h.sd, inc.sd)]
out2 <- vlohi(out[[1]], out[[2]])
hbc$tbhiv <- out[[1]]
hbc$tbhiv.lo <- out2[1,]
hbc$tbhiv.hi <- out2[2,]
hbc$tbhiv.sd <- out[[2]]


# add c.newinc
#
out <- est[g.hbc == TRUE, sum(c.newinc, na.rm = TRUE), by = c("year")]
hbc$c.newinc <- out$V1
hbc$newinc <- hbc$c.newinc * m / hbc$pop

save(hbc, file = here('data/hbc.rda'))
fwrite(hbc, file = here(paste0('csv/hbc_', Sys.Date(), '.csv')))


#' ---
#' title: incidence B+ pulmonary
#' author: Philippe Glaziou
#' date: 19-07-2022
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

library(data.table)
library(here)

yr <- 2021
m <- 1e5

load(here('data/est.rda'))
load(here('data/tb.rda'))

source(here('R/fun.R'))

# vectorized lohi
#
vlohi <- Vectorize(lohi, c('ev', 'sd'))


# weighted average b+ by income group and year
#
tb <- merge(tb, est[year == yr, .(iso3, g.income)], by = 'iso3')
tb[year > 2013, conf.denom := rowSums(cbind(new.ep, ret.rel.ep), na.rm = TRUE)]
out <-
  tb[year > 2017 &
       !is.na(g.income) &
       c.newinc > 1000, .(conf = weighted.mean(conf, w = conf.denom, na.rm = T)), by =
       .(g.income, year)]
out.sd <-
  tb[year > 2017 &
       !is.na(g.income) &
       c.newinc > 1000, .(conf.sd = weighted.sd(conf, w = conf.denom, na.rm = T)), by =
       .(g.income, year)]

tb[year > 2017 &
     !is.na(g.income) &
     c.newinc > 1000, .(min(conf)), by = .(g.income, year)]
tb[year > 2017 &
     !is.na(g.income) &
     c.newinc > 1000, .(max(conf)), by = .(g.income, year)]


# use HIC values of conf to estimate incidence B+
#
hic <- last(out$conf)
hic.sd <- last(out.sd$conf.sd)

est2 <- merge(est[year==yr, .(iso3,pop,g.whoregion,g.income,inc,inc.sd)], tb[year==yr, .(iso3,ep)], by='iso3')

# impute missing ep
#
out <- est2[, .(ep.imp = mean(ep, na.rm=T)), by=g.whoregion]
est3 <- merge(est2, out, by='g.whoregion')
est3[is.na(ep), ep := ep.imp]

bpos <- est3[, {
  tmp = prodXY(inc * (1 - ep), hic, inc.sd * (1 - ep), hic.sd)

  list (inc.bc = tmp[[1]],
        inc.bc.sd = tmp[[2]])
},
by = .(iso3)]

sel <- bpos$inc.bc > 0
out <- vlohi(bpos$inc.bc[sel] / m, bpos$inc.bc.sd[sel] / m)

bpos$inc.bc.lo[sel] <- out[1, ] * m
bpos$inc.bc.hi[sel] <- out[2, ] * m

bpos[!sel]
bpos[sel, test.bounds(inc.bc, inc.bc.lo, inc.bc.hi)]
bpos <- merge(bpos, est[year == yr, .(iso3, inc)], by = 'iso3')
bpos[, test.AgeB(inc, inc.bc)]
bpos[, inc := NULL]



save(bpos, file = here('data/bpos.rda'))
fwrite(bpos, file = here(paste0('csv/bpos', Sys.Date(), '.csv')))

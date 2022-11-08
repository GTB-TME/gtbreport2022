#' ---
#' title: Estimated lives saved
#' author: Philippe Glaziou
#' date: 10/08/2022
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---

#' (Last updated: `r Sys.Date()`)
#'

load(here('data/global.rda'))
load(here('data/regional.rda'))
load(here('data/est.rda'))

m <- 1e5
M <- 1e6


# TODO: fix inefficient code



# vectorized betaop
betaop <-
  function(ev1,
           ev2,
           se1,
           se2,
           op = "*",
           distr = F,
           nsim = 500000) {
    # defaults to product for a ratio, op='/'
    # ev1, ev2 = expected values
    # se1, se2 = standard errors
    # dist = returns distribution in a vector of size nsim

    if (is.na(ev1) ||
        is.na(ev2) ||
        is.na(se1) || is.na(se2))
      return (c(
        ev = NA,
        se = NA,
        lo = NA,
        hi = NA
      ))

    stopifnot(ev1 <= 1 & ev2 <= 1)
    stopifnot(se1 <= 0.5 & se2 <= 0.5)
    stopifnot(ev1 >= 0 & ev2 >= 0)
    stopifnot(se1 >= 0 & se2 >= 0)
    stopifnot(op %in% c('*', '/', '+', '-'))

    get.beta <- function(ev, sd){
      S = (ev * (1 - ev) / sd ^ 2) - 1
      a = S * ev
      b = S * (1 - ev)
      return(c(a = a, b = b))
    }

    par1 <- get.beta(ev1, se1)
    par2 <- get.beta(ev2, se2)

    out <-
      get(op)(rbeta(nsim, par1[1], par1[2]), rbeta(nsim, par2[1], par2[2]))
    if (distr)
      return(out)
    else
      return(c(
        ev = mean(out),
        se = sd(out),
        lo = quantile(
          out,
          prob = 0.025,
          names = FALSE,
          na.rm = TRUE
        ),
        hi = quantile(
          out,
          prob = 0.975,
          names = FALSE,
          na.rm = TRUE
        )
      ))
  }
vbetaop <- Vectorize(betaop, c('ev1', 'ev2', 'se1', 'se2'))



lsaved <- function(dta) {
  # CFRs untreated
  # HIV negative not on TB treatment  0.43 (0.28 - 0.53), assume ~beta
  cfrn <- 0.43
  cfrn.sd <- (0.53 - 0.28) / 4
  # HIV positive not on ART, not on TB treatment  0.78 (0.65 - 0.94), assume ~beta
  cfrp <- 0.78
  cfrp.sd <- (0.94 - 0.65) / 4

  get.beta <- function(ev, sd){
    S = (ev * (1 - ev) / sd ^ 2) - 1
    a = S * ev
    b = S * (1 - ev)
    return(c(a = a, b = b))
  }
  # lives saved HIV-neg, HIV-pos and total by row
  for (i in 1:dim(dta)[1]) {
    # TODO: optimize this
    #
    # counterfactual (rates)
    cfn <-
      betaop(
        dta$inc.nh[i] / m,
        cfrn,
        dta$inc.nh.sd[i] / m,
        cfrn.sd,
        op = '*',
        dist = T,
        nsim = M
      )
    cfp <-
      betaop(
        dta$inc.h[i] / m,
        cfrp,
        dta$inc.h.sd[i] / m,
        cfrp.sd,
        op = '*',
        dist = T,
        nsim = M
      )

    # factual (rates)
    parn <- get.beta(dta$mort.nh[i] / m, dta$mort.nh.sd[i] / m)
    fn <- rbeta(M, parn[1], parn[2])
    if (dta$mort.h[i] > 0) {
      parp <- get.beta(dta$mort.h[i] / m, dta$mort.h.sd[i] / m)
      fp <- rbeta(M, parp[1], parp[2])
    } else
      fp <- 0

    # difference counterfactual minus factual (absolute numbers)
    savedn <- pmax((cfn - fn), 0) * dta$pop[i]
    savedp <- pmax((cfp - fp), 0) * dta$pop[i]
    saved <- savedn + savedp

    dta$savedn[i] <- mean(savedn)
    dta$savedn.sd[i] <- sd(savedn)
    dta$savedn.lo[i] <- quantile(savedn, probs = 0.025, na.rm=TRUE)
    dta$savedn.hi[i] <- quantile(savedn, probs = 0.975, na.rm=TRUE)

    dta$savedp[i] <- mean(savedp)
    dta$savedp.sd[i] <- sd(savedp)
    dta$savedp.lo[i] <- quantile(savedp, probs = 0.025, na.rm=TRUE)
    dta$savedp.hi[i] <- quantile(savedp, probs = 0.975, na.rm=TRUE)

    dta$saved[i] <- mean(saved)
    dta$saved.sd[i] <- sd(saved)
    dta$saved.lo[i] <- quantile(saved, probs = 0.025, na.rm=TRUE)
    dta$saved.hi[i] <- quantile(saved, probs = 0.975, na.rm=TRUE)
  }

  return(dta)
}



clsaved <- function(dta) {
  # cumulative lives saved
  savedn <- sum(dta$savedn)
  savedn.sd <- sqrt(sum(dta$savedn.sd ^ 2))
  savedn.lo <- savedn + qnorm(0.025) * savedn.sd
  savedn.hi <- savedn + qnorm(0.975) * savedn.sd

  savedp <- sum(dta$savedp)
  savedp.sd <- sqrt(sum(dta$savedp.sd ^ 2))
  savedp.lo <- savedp + qnorm(0.025) * savedp.sd
  savedp.hi <- savedp + qnorm(0.975) * savedp.sd

  saved <- sum(dta$saved)
  saved.sd <- sqrt(sum(dta$saved.sd ^ 2))
  saved.lo <- saved + qnorm(0.025) * saved.sd
  saved.hi <- saved + qnorm(0.975) * saved.sd
  dta <- data.table(
    savedn,
    savedn.lo,
    savedn.hi,
    savedn.sd,
    savedp,
    savedp.lo,
    savedp.hi,
    savedp.sd,
    saved,
    saved.lo,
    saved.hi,
    saved.sd
  )
  return(dta)
}


# vectorized lohi
vlohi <- Vectorize(lohi, c('ev', 'sd'))




#' # Globally (2000-2018)
#'
set.seed(123)
global.lsaved <- lsaved(global)
global.lsaved <-
  global.lsaved[, list(
    year = as.character(year),
    savedn,
    savedn.lo,
    savedn.hi,
    savedn.sd,
    savedp,
    savedp.lo,
    savedp.hi,
    savedp.sd,
    saved,
    saved.lo,
    saved.hi,
    saved.sd
  )]
global.clsaved <- clsaved(global.lsaved)
global.clsaved$year <- 'Cumulative'

saved.global <- rbind(global.lsaved, global.clsaved, use.names = TRUE)

saved.global.print <- saved.global[, list(
  year,
  saved.hivneg = signif(savedn /
                          M, 3),
  saved.hivneg.ui = paste(
    "(",
    as.character(signif(savedn.lo / M, 3)),
    "-",
    as.character(signif(savedn.hi /
                          M, 3)),
    ")",
    sep = ''
  ),
  saved.hivpos = signif(savedp /
                          M, 3),
  saved.hivpos.ui = paste(
    "(",
    as.character(signif(savedp.lo / M, 3)),
    "-",
    as.character(signif(savedp.hi /
                          M, 3)),
    ")",
    sep = ''
  ),
  saved = signif(saved / M, 3),
  saved.ui = paste(
    "(",
    as.character(signif(saved.lo / M, 3)),
    "-",
    as.character(signif(saved.hi /
                          M, 3)),
    ")",
    sep = ''
  )
)]
(saved.global.print)
fwrite(saved.global.print,
       file = paste0(here('output/globalSaved'), Sys.Date(), '.csv'))





#' HIV+
global.lsaved <- lsaved(global[year > 2004])
global.lsaved <-
  global.lsaved[, list(
    year = as.character(year),
    savedn,
    savedn.lo,
    savedn.hi,
    savedn.sd,
    savedp,
    savedp.lo,
    savedp.hi,
    savedp.sd,
    saved,
    saved.lo,
    saved.hi,
    saved.sd
  )]
global.clsaved <- clsaved(global.lsaved)
global.clsaved$year <- 'Cumulative'

hsaved.global <-
  rbind(global.lsaved, global.clsaved, use.names = TRUE)

hsaved.global.print <- hsaved.global[, list(
  year,
  saved.hivpos = signif(savedp /
                          M, 2),
  saved.hivpos.ui = paste(
    "(",
    as.character(signif(savedp.lo / M, 3)),
    "-",
    as.character(signif(savedp.hi /
                          M, 2)),
    ")",
    sep = ''
  )
)]
(hsaved.global.print)
fwrite(hsaved.global.print,
       file = paste0(here('output/HIVposGlobalSaved'), Sys.Date(), '.csv'))





#' by WHO region (2000-2019)
#'
set.seed(101)
saved.global <- lsaved(global)
csaved.global <- clsaved(saved.global)
regional.lsaved <- lsaved(regional)
regional.clsaved <-
  regional.lsaved[, clsaved(.SD), keyby = g.whoregion]
regional.clsaved$year <- 'Cumulative'
saved.regional <-
  rbind(regional.clsaved,
        csaved.global,
        use.names = TRUE,
        fill = T)
names(saved.regional)[1] <- 'Region'
saved.regional.print <-
  saved.regional[, list(
    Region = as.character(Region),
    saved.hivneg = signif(savedn /
                            M, 3),
    saved.hivneg.ui = paste(
      "(",
      as.character(signif(savedn.lo / M, 3)),
      "-",
      as.character(signif(savedn.hi /
                            M, 3)),
      ")",
      sep = ''
    ),
    saved.hivpos = signif(savedp /
                            M, 3),
    saved.hivpos.ui = paste(
      "(",
      as.character(signif(savedp.lo / M, 3)),
      "-",
      as.character(signif(savedp.hi /
                            M, 3)),
      ")",
      sep = ''
    ),
    saved = signif(saved / M, 3),
    saved.ui = paste(
      "(",
      as.character(signif(saved.lo / M, 3)),
      "-",
      as.character(signif(saved.hi /
                            M, 3)),
      ")",
      sep = ''
    )
  )]
saved.regional.print$Region[7] <- 'Global'
(saved.regional.print)

fwrite(
  saved.regional.print,
  file = paste0(here('output/RegionalLivesSaved_'), Sys.Date(), '.csv')
)





#' by country - request from GF
#'
#' BGD since 2000
cty.lsaved <- function(iso = 'BGD',
                       start = 2000,
                       csv = FALSE) {
  th <- 1000
  bgd <- lsaved(est[iso3 == iso & year >= start])
  bgd.cum <- clsaved(bgd)
  bgd.saved <-
    bgd[, list(
      year = as.character(year),
      savedn,
      savedn.lo,
      savedn.hi,
      savedn.sd,
      savedp,
      savedp.lo,
      savedp.hi,
      savedp.sd,
      saved,
      saved.lo,
      saved.hi,
      saved.sd
    )]
  bgd.cum$year <- 'Cumulative'
  bgd.saved <- rbind(bgd.saved, bgd.cum, use.names = T)
  bgd.saved <- bgd.saved[, list(
    year,
    saved.hivneg = signif(savedn / th, 3),
    saved.hivneg.ui = paste(
      "(",
      as.character(signif(savedn.lo / th, 3)),
      "-",
      as.character(signif(savedn.hi /
                            th, 3)),
      ")",
      sep = ''
    ),
    saved.hivpos = signif(savedp / th, 3),
    saved.hivpos.ui = paste(
      "(",
      as.character(signif(savedp.lo / th, 3)),
      "-",
      as.character(signif(savedp.hi /
                            th, 3)),
      ")",
      sep = ''
    ),
    saved = signif(saved / th, 3),
    saved.ui = paste(
      "(",
      as.character(signif(saved.lo / th, 3)),
      "-",
      as.character(signif(saved.hi /
                            th, 3)),
      ")",
      sep = ''
    )
  )]

  if (csv)
    fwrite(bgd.saved, file = paste(here('output/'), iso, '_livesSaved.csv', sep = ''))
  return(bgd.saved)
}


#---- temp code
sel <- is.na(est$inc.h) | est$inc.h == 0
est$inc.h[sel] <- 1e-5
sel <- is.na(est$inc.h.sd) | est$inc.h.sd == 0
est$inc.h.sd[sel] <- 1e-5

sel <- is.na(est$inc.nh) | est$inc.nh == 0
est$inc.nh[sel] <- 1e-5
sel <- is.na(est$inc.nh.sd) | est$inc.nh.sd == 0
est$inc.nh.sd[sel] <- 1e-5

sel <- is.na(est$mort.h) | est$mort.h == 0
est$mort.h[sel] <- 1e-6
sel <- is.na(est$mort.h.sd) | est$mort.h.sd == 0
est$mort.h.sd[sel] <- 1e-6

sel <- is.na(est$mort.nh) | est$mort.nh == 0
est$mort.nh[sel] <- 1e-6
sel <- is.na(est$mort.nh.sd) | est$mort.nh.sd == 0
est$mort.nh.sd[sel] <- 1e-6
#---- end temp code


cty.lsaved2 <- function(dta, start = 2000) {
  th <- 1000
  bgd <- lsaved(dta[year >= start])
  bgd.cum <- clsaved(bgd)
  bgd.saved <-
    bgd[, list(
      year = as.character(year),
      savedn,
      savedn.lo,
      savedn.hi,
      savedn.sd,
      savedp,
      savedp.lo,
      savedp.hi,
      savedp.sd,
      saved,
      saved.lo,
      saved.hi,
      saved.sd
    )]
  bgd.cum$year <- 'Cumulative'

  bgd.saved <- rbind(bgd.saved, bgd.cum, use.names = T)
  bgd.saved <- bgd.saved[, list(
    year,
    saved.hivneg = signif(savedn / th, 3),
    saved.hivneg.lo = signif(savedn.lo / th, 3),
    saved.hivneg.hi = signif(savedn.hi / th, 3),
    saved.hivpos = signif(savedp / th, 3),
    saved.hivpos.lo = signif(savedp.lo / th, 3),
    saved.hivpos.hi = signif(savedp.hi / th, 3),
    saved = signif(saved / th, 3),
    saved.lo = signif(saved.lo / th, 3),
    saved.hi = signif(saved.hi / th, 3)
  )]

  return(bgd.saved)
}


#' this takes time...
#'
sel <- est$mort.nh.sd/m > est$mort.nh/m * (1 - est$mort.nh/m)
table(sel)  # SD too large for a beta to fit

knitr::kable(est[sel, .(iso3, year, inc, mort.nh, mort.nh.sd)])
est[sel, mort.nh.sd := mort.nh * 0.2]

sel <- est$inc.nh.sd/m > est$inc.nh/m * (1 - est$inc.nh/m)
table(sel)
est[sel, inc.nh.sd := inc.nh * 0.2]

sel <- est$inc.h.sd/m > est$inc.h/m * (1 - est$inc.h/m)
table(sel)
est[sel, inc.h.sd := inc.h * 0.2]

knitr::kable(est[sel, .(iso3, year, inc, inc.h, inc.h.sd)])
est[sel, inc.h.sd := inc.h * 0.2]

sel <- est$inc.sd/m > est$inc/m * (1 - est$inc/m)
table(sel)

out <- est[, cty.lsaved2(.SD, start = 2000), by = iso3]

fwrite(out, file = paste0(here('csv/lsaved_'), Sys.Date(), '.csv'))


#' reload clean est
#'
load(here('data/est.rda'))


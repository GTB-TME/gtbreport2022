#' ---
#' title: functions
#' author: Philippe Glaziou
#' date: 2019-05-22
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---

#' # Preamble
#'
#' (Last updated: `r Sys.Date()`)
#'
#' Utility functions
#'
#' # load deps and data
library('data.table')
library('ggplot2')

theme_set(theme_bw())
devAskNewPage(ask = FALSE)



isocty <- function(iso){
  #' get country name from iso3 code
  #'
  #' @param iso ISO3 country code.
  #'
  #' @export
  cty[iso, country]
}


iso <- function (ct = "China") {
  #' get iso code from country name
  #'
  #' @param ct country name.
  #'
  #' @export
  suppressWarnings({
    country <-
      as.character(cty$country[grep(ct, cty$country, ignore.case = TRUE)])
    iso3 <- as.character(cty$iso3[cty$country %in% country])
    out <- cbind(iso3, country)
    return(out)
  })
}



#' prime?
#'
is.prime <- function(n){
  #' is prime?
  #'
  #' @param n integer.
  #'
  #' @export
  stopifnot(is.integer(n))

  n == 2L || all(n %% 2L:ceiling(sqrt(n)) != 0)
}


sums <- function(x, ..., na.rm = T) {
  #' sum() defaulting to na.rm=T except when all elements are NA
  #'
  #' @param x vector to be summed.
  #'
  #' @export
  stopifnot(is.numeric(x))

  if (all(is.na(x)))
    as.numeric(NA)
  else
    base::sum(x, ..., na.rm = na.rm)
}


#' certification
#

test.bounds <- function(best, lo, hi) {
  stopifnot(sum(is.na(best)) == 0)
  stopifnot(sum(is.na(lo)) == 0)
  stopifnot(sum(is.na(hi)) == 0)
  stopifnot(best >= lo)
  stopifnot(best <= hi)
  invisible(T)
}

test.ispos <- function(var) {
  stopifnot(sum(is.na(var)) == 0)
  stopifnot(var >= 0)
  invisible(T)
}

test.AgB <- function(A, B) {
  stopifnot(sum(is.na(A)) == 0)
  stopifnot(sum(is.na(B)) == 0)
  stopifnot(A > B)
  invisible(T)
}

test.AgeB <- function(A, B) {
  stopifnot(sum(is.na(A)) == 0)
  stopifnot(sum(is.na(B)) == 0)
  stopifnot(A >= B)
  invisible(T)
}

test.isbinom <- function(var) {
  test.AgeB(1, var) && test.AgeB(var, 0)
}


#' weighted sd (adapted from Hmisc)
#'
weighted.sd <-
  function(x,
           w = NULL,
           normwt = FALSE,
           na.rm = TRUE,
           method = c("unbiased",
                      "ML"))
  {
    method <- match.arg(method)
    if (!length(w)) {
      if (na.rm)
        x <- x[!is.na(x)]
      return(var(x))
    }
    if (na.rm) {
      s <- !is.na(x + w)
      x <- x[s]
      w <- w[s]
    }
    if (normwt)
      w <- w * length(x) / sum(w)
    if (normwt || method == "ML")
      return(as.numeric(stats::cov.wt(cbind(x), w, method = method)$cov))
    sw <- sum(w)
    if (sw <= 1)
      warning("only one effective observation; variance estimate undefined")
    xbar <- sum(w * x) / sw
    sqrt(sum(w * ((x - xbar) ^ 2)) / (sw - 1))
  }



#' product of two random variables X and Y using Taylor expansion
#'
prodXY <- function(X, Y, sdX, sdY, covXY = 0) {
  eXY <- X * Y + covXY
  sdXY <-
    sqrt(X^2 * sdY^2 + Y^2 * sdX^2 + 2 * X * Y * covXY + sdX^2 * sdY^2 + covXY ^
    2)
  return(list("mean" = eXY, "sd" = sdXY))
}


#' ratio of two random variables using Taylor expansion
#'
divXY <- function(X, Y, sdX, sdY, covXY = 0) {
  eXY <- X / Y - covXY / Y ^ 2 + X * sdY^2 / Y ^ 3
  sdXY <-
    sqrt(sdX^2 / Y ^ 2 - 2 * X * covXY / Y ^ 3 + X ^ 2 * sdY^2 / Y ^ 4)
  return(list("mean" = eXY, "sd" = sdXY))
}


#' EV and variance of the log of a RV
#' based on Taylor expansion
#'
logX <- function(X, sdX) {
  #' @param X = EV.
  #' @param varX = variance.
  #' @export
  E.logX <- log(X) - sdX^2 / (2 * X ^ 2)
  sd.logX <- sqrt(sdX^2 / X ^ 2)
  return(list("mean" = E.logX, "sd" = sd.logX))
}

#' logit and invlogit
#'
logit <- function (x)
  log (x / (1 - x))
invlogit <- function (x)
  1 / (1 + exp(-x))



#' returns Beta shape and scale params using the method of moments
#'
get.beta <- function(ev, sd) {
  #' @param ev expected value.
  #' @param sd standard deviation.
  #' @export
  stopifnot(ev > 0 & ev < 1)
  stopifnot(sd > 0)

  S = (ev * (1 - ev) / sd ^ 2) - 1
  if (S < 0)
    stop('Not distributed Beta: sd^2 >= ev*(1-ev)')

  a = S * ev
  b = S * (1 - ev)
  return(c(a = a, b = b))
}


#' returns gamma shape and scale params using the method of moments
#'
get.gamma <- function(ev, sd) {
  #' @param ev expected value.
  #' @param sd standard deviation.
  #' @export
  stopifnot(sd > 0)

  v = ev / sd ^ 2
  r = ev ^ 2 / sd ^ 2
  return(c(r = r, v = v))
}


#' generate low and high bounds assuming Beta distribution
#'
lohi <- function(ev, sd) {
  #' @param ev expected value.
  #' @param sd standard deviation.
  #' @export
  stopifnot(ev > 0 & ev < 1)
  stopifnot(sd > 0)

  par <- get.beta(ev, sd)
  lo <- qbeta(0.025, par[1], par[2])
  hi <- qbeta(0.975, par[1], par[2])
  return(c(lo = lo, hi = hi))
}

#' generate low and high bounds assuming Gamma distribution
#'
glohi <- function(ev, sd) {
  #' @param ev expected value.
  #' @param sd standard deviation.
  #' @export
  stopifnot(sd > 0)

  par <- get.gamma(ev, sd)
  lo <- qgamma(0.025, par[1], par[2])
  hi <- qgamma(0.975, par[1], par[2])
  return(c(lo = lo, hi = hi))
}


#' aggregated rates
#'
addXY <-
  function (r,
            r.lo,
            r.hi,
            r.sd,
            weights = 1,
            method = "beta")
  {
    #' @param r rates
    #' @param r.lo low bounds
    #' @param r.hi high bounds
    #' @param r.sd SDs
    #' @param weights weights
    #' @param method distribution of rates, defaults to beta else normal
    #' @export
    if (is.null(r) || length(r) == 0)
      stop("Error: r must contain at least one value")
    if (sum(r < 0 & !is.na(r) & method == "beta"))
      stop("Error: r must be positive with method 'beta'")
    if (sum(r > 1 & !is.na(r) & method == "beta"))
      stop("Error: r must be between 0 and 1 with method 'beta'")
    if (missing(r.sd))
      r.sd <- (r.hi - r.lo) / 4
    if (missing(r.lo) & !missing(r.sd))
      r.lo <- numeric()
    if (missing(r.hi) & !missing(r.sd))
      r.hi <- numeric()
    if (sum(r.lo < 0 & !is.na(r.lo) & method == "beta"))
      stop("Error: r.lo must be positive with method 'beta'")
    if (sum(r.lo > 1 & !is.na(r.lo) & method == "beta"))
      stop("Error: r.lo must be between 0 and 1 with method 'beta'")
    if (sum(r.hi < 0 & !is.na(r.hi) & method == "beta"))
      stop("Error: r.hi must be positive with method 'beta'")
    if (sum(r.hi > 1 & !is.na(r.hi) & method == "beta"))
      stop("Error: r.hi must be between 0 and 1 with method 'beta'")
    if (sum(r.sd > 1 & !is.na(r.sd) & method == "beta"))
      stop("Error: sd must be between 0 and 1 with method 'beta'")
    if (sum(r[!is.na(r) & is.na(r.sd)]))
      stop("Error: some values for r are supplied without uncertainty")
    if (sum(r.sd < 0 & !is.null(r.sd) & !is.na(r.sd)))
      stop("Error: sd must be positive")
    if (!is.null(r.sd))
      v <- r.sd ^ 2
    else
      v <- ((r.hi - r.lo) / 4) ^ 2
    sw <- ifelse(length(weights) > 1, sum(weights[!is.na(r)],
                                          na.rm = TRUE), 1)
    out.m <- sum(r * weights, na.rm = TRUE) / sw
    out.v <-
      ifelse(length(weights) > 1,
             sum(v[!is.na(r)] * weights[!is.na(r)] ^ 2,
                 na.rm = TRUE) / sw ^ 2,
             sum(v))
    if (method == "beta") {
      S <- (out.m * (1 - out.m) / out.v) - 1
      a <- S * out.m
      b <- S * (1 - out.m)
      lo <- qbeta(0.025, a, b)
      hi <- qbeta(0.975, a, b)
    }
    else {
      lo <- qnorm(0.025, out.m, sqrt(out.v))
      hi <- qnorm(0.975, out.m, sqrt(out.v))
    }
    if (all(weights == 1))
      return(data.frame(
        r = out.m,
        r.lo = lo,
        r.hi = hi,
        r.sd = sqrt(out.v)
      ))
    else
      return(
        data.frame(
          r = out.m,
          r.lo = lo,
          r.hi = hi,
          r.sd = sqrt(out.v),
          r.num = out.m * sw,
          r.lo.num = lo * sw,
          r.hi.num = hi *
            sw,
          pop = sw
        )
      )
  }


#' **Statistical ensemble**. The rate $R$ obtained using method i is assumed
#' distributed Beta with shape and scale parameters $\alpha_i+1$ and
#' $\beta_i+1$, respectively, and determined using the method of moments:
#' $R_i \sim B(\alpha_i+1,\beta_i+1)$ so that
#'
#'
#' $\textrm{Prob}(x = \textrm{TB}) = \int_{0}^{1} x B(\alpha_i, \beta_i)\, dx = \frac{\alpha_i + 1}{\alpha_i + \beta_i + 2}$
#'
#' The combined probability is then expressed as
#'
#' $\textrm{Prob}(x = \textrm{TB}) = \frac{\sum{\alpha_i} + 1}{\sum{\alpha_i} + \sum{\beta_i} + 2}$
#'
#' $\textrm{Var} = \frac{(\sum{\alpha_i+1})(\sum{\beta_i+1})}{(\sum{\alpha}+\sum{\beta+2})^2(\sum{\alpha}+\sum{\beta+3})}$
#'
#'
ensbeta <- function(xi, xi.sd) {
  stopifnot(xi < 1 & xi.sd < 1)
  stopifnot(xi > 0 & xi.sd > 0)

  get.beta <- function(ev, sd) {
    S = (ev * (1 - ev) / sd ^ 2) - 1
    a = S * ev
    b = S * (1 - ev)
    return(c(a = a, b = b))
  }

  vget.beta <- Vectorize(get.beta, c('ev', 'sd'))

  w <- vget.beta(xi, xi.sd) - 1
  a <- sum(w[1, ])
  b <- sum(w[2, ])
  pw <- list(c = a + 1, d = b + 1)
  k <- pw$c + pw$d
  post.ev <- pw$c / k
  post.lo <- qbeta(0.025, pw$c, pw$d)
  post.hi <- qbeta(0.975, pw$c, pw$d)
  post.sd <- sqrt(pw$c * pw$d / (k ^ 2 * (k + 1)))
  return(
    list(
      post.param = c(shape = pw$c, scale = pw$d),
      post.ev = post.ev,
      post.lo = post.lo,
      post.hi = post.hi,
      post.sd = post.sd
    )
  )
}


#' Weighted quantile
#'
#' Function copied from **spatstat** package.
#'
#' @param x Vector of values.
#' @param w Vector of weights.
#' @param probs Vector of probabilities.
#' @param na.rm Ignore missing data?
#' @export
weighted.quantile <- function(x, w, probs=seq(0,1,0.25), na.rm=TRUE) {
  x <- as.numeric(as.vector(x))
  w <- as.numeric(as.vector(w))
  if(anyNA(x) || anyNA(w)) {
    ok <- !(is.na(x) | is.na(w))
    x <- x[ok]
    w <- w[ok]
  }
  stopifnot(all(w >= 0))
  if(all(w == 0)) stop("All weights are zero", call.=FALSE)
  #'
  oo <- order(x)
  x <- x[oo]
  w <- w[oo]
  Fx <- cumsum(w)/sum(w)
  #'
  result <- numeric(length(probs))
  for(i in seq_along(result)) {
    p <- probs[i]
    lefties <- which(Fx <= p)
    if(length(lefties) == 0) {
      result[i] <- x[1]
    } else {
      left <- max(lefties)
      result[i] <- x[left]
      if(Fx[left] < p && left < length(x)) {
        right <- left+1
        y <- x[left] + (x[right]-x[left]) * (p-Fx[left])/(Fx[right]-Fx[left])
        if(is.finite(y)) result[i] <- y
      }
    }
  }
  names(result) <- paste0(format(100 * probs, trim = TRUE), "%")
  return(result)
}


#' Weighted median
#'
#' Function copied from **spatstat** package.
#'
#' @param x Vector of values.
#' @param w Vector of weights.
#' @param na.rm Ignore missing data?
#' @export
weighted.median <- function(x, w, na.rm=TRUE) {
  unname(weighted.quantile(x, probs=0.5, w=w, na.rm=na.rm))
}


#' combine plots - source: R cookbook
#'
multiplot <- function(..., plotlist = NULL, cols) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # Make the panel
  plotCols = cols                       # Number of columns of plots
  plotRows = ceiling(numPlots / plotCols) # Number of rows needed, calculated from # of cols

  # Set up the page
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
  vplayout <- function(x, y)
    viewport(layout.pos.row = x, layout.pos.col = y)

  # Make each plot, in the correct location
  for (i in 1:numPlots) {
    curRow = ceiling(i / plotCols)
    curCol = (i - 1) %% plotCols + 1
    print(plots[[i]], vp = vplayout(curRow, curCol))
  }
}




#' function to plot mortality trends
#'
mplot <- function(iso,
                  dta = est,
                  hiv = FALSE,
                  ylog = TRUE) {
  #' @param iso iso3 country code.
  #' @param dta dataset.
  #' @param hiv hiv lines?
  #' @param ylog log scale on y-axis?
  #' @export
  if (missing(iso) || !is.character(iso) || nchar(iso) != 3)
    stop("Error: iso must be a valid ISO3 country code")
  if (is.na(match('cty', objects(
    envir = parent.frame(), all.names = TRUE
  ))))
    stop("Error: cty must be loaded")
  iso <- toupper(iso)
  isocty <- function(iso)
    cty[iso, country]

  p <-
    qplot(
      year,
      mort.nh,
      data = subset(dta, iso3 == iso),
      geom = 'line',
      main = paste('TB mortality in', as.character(isocty(iso)))
    ) +
#    geom_point(aes(year, vr.raw), shape = I(4), size = I(3)) +
    geom_ribbon(
      aes(year, ymin = mort.nh.lo, ymax = mort.nh.hi),
      fill = I('blue'),
      outline.type = 'full',
      alpha = I(0.3)
    ) +
    xlab('') +
    ylab('Rate per 100,000/year') +
    theme_bw(base_size = 16)

  q <- p + geom_line(aes(year, mort.h), colour = I('red')) +
    geom_ribbon(
      aes(year, ymin = mort.h.lo, ymax = mort.h.hi),
      fill = I('red'),
      outline.type = 'full',
      alpha = I(0.3)
    )
  if (ylog == TRUE) {
    p <-
      p + coord_trans(y = "log10") + ylab('Rate per 100,000/year (log scale)')
    q <-
      q + coord_trans(y = "log10") + ylab('Rate per 100,000/year (log scale)')
  }

  if (hiv)
    return(q)
  else
    return(p)
}



#' function to plot incidence trends
#'
iplot <- function(iso,
                  dta = est,
                  hiv = FALSE,
                  ylog = TRUE) {
  #' @param iso iso3 country code.
  #' @param dta dataset.
  #' @param hiv hiv lines?
  #' @param ylog log scale on y-axis?
  #' @export
  if (missing(iso) || !is.character(iso) || nchar(iso) != 3)
    stop("Error: iso must be a valid ISO3 country code")
  if (is.na(match('cty', objects(
    envir = parent.frame(), all.names = TRUE
  ))))
    stop("Error: cty must be loaded")
  iso <- toupper(iso)
  isocty <- function(iso)
    cty[iso, country]

  p <- qplot(
    year,
    inc,
    data = subset(dta, iso3 == iso),
    geom = 'line',
    main = paste('TB incidence in', as.character(isocty(iso)))
  ) +
    geom_ribbon(
      aes(year, ymin = inc.lo, ymax = inc.hi),
      fill = I('blue'),
      outline.type = 'full',
      alpha = I(0.3)
    ) +
    geom_line(aes(year, newinc)) +
    xlab('') + ylab('Rate per 100,000/year') +
    theme_bw(base_size = 16)

  q <- p + geom_line(aes(year, inc.h), colour = I('red')) +
    geom_ribbon(
      aes(year, ymin = inc.h.lo, ymax = inc.h.hi),
      fill = I('red'),
      outline.type = 'full',
      alpha = I(0.3)
    )

  if (ylog == TRUE) {
    p <-
      p + coord_trans(y = "log10") + ylab('Rate per 100,000/year (log scale)')
    q <-
      q + coord_trans(y = "log10") + ylab('Rate per 100,000/year (log scale)')
  }

  if (hiv)
    return(q)
  else
    return(p)
}


#' incidence comparison plot
#'
cplot <- function(iso, new, old, ylog = F) {
  #' @param iso iso3 country code.
  #' @param new dataset.
  #' @param old comparison dataset.
  #' @param ylog log scale on y-axis?
  #' @export
  if (missing(iso) || !is.character(iso) || nchar(iso) != 3)
    stop("Error: iso must be a valid ISO3 country code")
  if (is.na(match('cty', objects(
    envir = parent.frame(), all.names = TRUE
  ))))
    stop("Error: cty must be loaded")
  iso <- toupper(iso)
  isocty <- function(iso)
    cty[iso, country]

  p <-
    qplot(
      year,
      inc,
      data = subset(new, iso3 == iso),
      geom = 'line',
      colour = I('grey90'),
      main = paste('TB incidence in', as.character(isocty(iso)))
    ) +
    geom_ribbon(
      aes(year, ymin = inc.lo, ymax = inc.hi),
      fill = I('blue'),
      outline.type = 'full',
      alpha = I(.4)
    ) +
    geom_line(
      aes(year, inc),
      data = subset(old, iso3 == iso),
      colour = I('yellow'),
      linetype = I(2)
    ) +
    geom_ribbon(
      aes(year, ymin = inc.lo, ymax = inc.hi),
      data = subset(old, iso3 == iso),
      outline.type = 'full',
      fill = I('yellow'),
      alpha = I(.4)
    ) +
    geom_line(aes(year, newinc)) + xlab('') + ylab('Incidence rate per 100k/yr')
  if (ylog)
    return(p + coord_trans(y = 'log10'))
  else
    (return(p))
}





#' function to plot tbhiv trends
#'
hplot <- function(iso, toplot = T) {
  #' @param iso iso3 country code.
  #' @param toplot returns plot or data.
  #' @export
  if (missing(iso) || !is.character(iso) || nchar(iso) != 3)
    stop("Error: iso must be a valid ISO3 country code")
  if (is.na(match('est', objects(
    envir = parent.frame(), all.names = TRUE
  ))) ||
  is.na(match('tb', objects(
    envir = parent.frame(), all.names = TRUE
  ))) ||
  is.na(match('sty', objects(
    envir = parent.frame(), all.names = TRUE
  ))) ||
  is.na(match('cty', objects(
    envir = parent.frame(), all.names = TRUE
  ))))
  stop("Error: est, tb, sty and cty must be loaded")
  iso <- toupper(iso)
  isocty <- function(iso)
    cty[iso, country]

  sources <- length(table(est$source.tbhiv[est$iso3 == iso]))

  dta <-
    merge(est[iso3 == iso, list(iso3, year, tbhiv, tbhiv.lo, tbhiv.hi)],
          tb[iso3 == iso, list(
            iso3,
            year,
            tot.newrel,
            hivtest.f,
            hivtest.p,
            hivtest.pos.f,
            hivtest.pos.p,
            hiv.art.p,
            hiv.art.f
          )],
          by = c('iso3', 'year'), all.x = TRUE)
  dta <-
    merge(dta, sty[, list(
      iso3,
      year,
      tbhiv.surv.prev,
      tbhiv.surv.cil,
      tbhiv.surv.ciu,
      tbhiv.sentin.prev,
      tbhiv.sentin.cil,
      tbhiv.sentin.ciu
    )],
    by = c('iso3', 'year'), all.x = TRUE)

  dta <- within(dta, {
    test.coverage <- 100 * hivtest.f / tot.newrel
    test.pos <- 100 * hivtest.pos.f / hivtest.f
    art <- 100 * hiv.art.f / hivtest.pos.f
  })

  lastyr <- max(dta$year[dta$iso3 == iso])
  dta$test.coverage[dta$year == lastyr] <-
    dta$hivtest.p[dta$year == lastyr] * 100 / dta$tot.newrel[dta$year == lastyr]
  dta$test.pos[dta$year == lastyr] <-
    dta$hivtest.pos.p[dta$year == lastyr] * 100 / dta$tot.newrel[dta$year ==
                                                                   lastyr]
  dta$art[dta$year == lastyr] <-
    dta$hiv.art.p[dta$year == lastyr] * 100 / dta$hivtest.pos.p[dta$year ==
                                                                  lastyr]


  p <- qplot(
    year,
    tbhiv * 100,
    data = dta,
    geom = 'line',
    main = paste('HIV prevalence in TB in', as.character(isocty(iso)))
  ) +
    scale_linetype_discrete(name = "Data Source") +
    geom_ribbon(
      aes(year, ymin = tbhiv.lo * 100, ymax = tbhiv.hi * 100),
      fill = I('red'),
      outline.type = 'full',
      alpha = I(0.3)
    ) +
    expand_limits(y = 0) +
    xlab('') + ylab('Percent') +
    geom_point(
      aes(year, tbhiv.surv.prev),
      colour = I('blue'),
      shape = I(16),
      size = I(3)
    ) +
    geom_linerange(aes(year, ymin = tbhiv.surv.cil, ymax = tbhiv.surv.ciu),
                   colour = I('blue')) +
    geom_point(
      aes(year, tbhiv.sentin.prev),
      colour = I('darkgreen'),
      shape = I(2),
      size = I(4)
    ) +
    geom_linerange(
      aes(year, ymin = tbhiv.sentin.cil, ymax = tbhiv.sentin.ciu),
      colour = I('darkgreen'),
      size = I(3)
    ) +
    geom_point(aes(year, test.pos),
               shape = I(4),
               size = I(3)) +
    theme_bw(base_size = 16)
  if (toplot)
    return(p)
  else
    (return(dta))
}


#' the above xplot functions can be combined with multiplot, e.g.:
#' multiplot(iplot('CHN', hiv=T), mplot('CHN'), pplot('CHN'), hplot('CHN'), cols=2)
#'
hivplot <- function(iso, start = 1990, toplot = T) {
  if (missing(iso) || !is.character(iso) || nchar(iso) != 3)
    stop("Error: iso must be a valid ISO3 country code")
  if (is.na(match('unaids', objects(
    envir = parent.frame(), all.names = TRUE
  ))) ||
  is.na(match('cty', objects(
    envir = parent.frame(), all.names = TRUE
  ))))
    stop("Error: unaids and cty must be loaded")
  iso <- toupper(iso)
  isocty <- function(iso)
    cty[iso, country]


  p <-
    qplot(
      year,
      hiv * 100,
      data = unaids[iso3 == iso & year >= start],
      geom = 'line',
      main = paste('HIV prevalence (all ages) in', as.character(isocty(iso)))
    ) +
    geom_ribbon(
      aes(year, ymin = hiv.lo * 100, ymax = hiv.hi * 100),
      fill = I('red'),
      outline.type = 'full',
      alpha = I(0.3)
    ) +
    expand_limits(y = 0) +
    xlab('') + ylab('Percent') +
    theme_bw(base_size = 16)
  if (toplot)
    return(p)
  else
    (return(dta))
}







#' population pyramids
#'
pyramid <- function(iso, yr = 2010) {
  #' @param iso iso3 country code.
  #' @param yr ref year for demographics.
  #' @export
  if (missing(iso) || !is.character(iso) || nchar(iso) != 3)
    stop("Error: iso must be a valid ISO3 country code")
  if (is.na(match('pop', objects(
    envir = parent.frame(), all.names = TRUE
  ))) ||
  is.na(match('cty', objects(
    envir = parent.frame(), all.names = TRUE
  ))))
    stop("Error: pop and cty must be loaded")

  library(plotrix)

  iso <- toupper(iso)

  dta <- as.data.frame(pop[iso3 == iso & year == yr])
  dta <- within(dta, {
    m04 <- e.pop.m04 / e.pop.num
    m514 <- e.pop.m514 / e.pop.num
    m1524 <- e.pop.m1524 / e.pop.num
    m2534 <- e.pop.m2534 / e.pop.num
    m3544 <- e.pop.m3544 / e.pop.num
    m4554 <- e.pop.m4554 / e.pop.num
    m5564 <- e.pop.m5564 / e.pop.num
    m65 <- e.pop.m65 / e.pop.num
    f04 <- e.pop.f04 / e.pop.num
    f514 <- e.pop.f514 / e.pop.num
    f1524 <- e.pop.f1524 / e.pop.num
    f2534 <- e.pop.f2534 / e.pop.num
    f3544 <- e.pop.f3544 / e.pop.num
    f4554 <- e.pop.f4554 / e.pop.num
    f5564 <- e.pop.f5564 / e.pop.num
    f65 <- e.pop.f65 / e.pop.num
  })
  country <- as.character(dta$country[1])
  mpop <-
    with(dta, c(m04, m514, m1524, m2534, m3544, m4554, m5564, m65) * 100)
  fpop <-
    with(dta, c(f04, f514, f1524, f2534, f3544, f4554, f5564, f65) * 100)
  agelab <-
    c("0-4",
      "5-14",
      "15-24",
      "25-34",
      "35-44",
      "45-54",
      "55-64",
      "65+")
  mcol = 'blue'
  fcol = 'red'

  pyramid.plot(
    mpop,
    fpop,
    labels = agelab,
    main = paste("Population pyramid in ", country, " (", yr, ")", sep =
                   ""),
    lxcol = mcol,
    rxcol = fcol,
    gap = 1.2,
    show.values = FALSE
  )
}





#' incidence, prevalence, mortality
#'
prev2inc <- function(prev,
                     prev.sd,
                     prevk = NA,
                     newinc,
                     tbhiv,
                     tbhiv.sd,
                     rtbhiv = NA,
                     rtbhiv.sd = NA)
{
  #' function to derive incidence from prevalence
  #'
  #' @param prev prevalence per capita.
  #' @param prev.sd standard deviation of prevalence.
  #' @param prevk prevalence on tx (known cases).
  #' @param newinc detection rate (new+relapse) per capita.
  #' @param tbhiv proportion HIV+ among incident cases.
  #' @param tbhiv.sd SD of proportion HIV+.
  #' @param rtbhiv HIV+ rate ratio (prevalent/incident).
  #' @param rtbhiv.sd SD of HIV+ rate ratio.
  #' @export
  require(propagate) # use second-order Taylor expansion about moments

  stopifnot(prev > 0 & prev.sd > 0)
  stopifnot((prevk > 0 &
               prevk <= prev) | is.na(prevk))
  stopifnot(newinc >= 0)
  stopifnot(tbhiv > 0 & tbhiv < 1)
  stopifnot(tbhiv.sd > 0 & tbhiv.sd < 1)
  stopifnot((rtbhiv > 0 &
               rtbhiv < 1 &
               rtbhiv.sd > 0 &
               rtbhiv.sd < 1) | is.na(rtbhiv)) & !is.na(rtbhiv.sd)

  # durations
  Du.nh <-
    c(2.5, sqrt(3 / 4))           # not detected, HIV-neg ~U(1,4)
  Dn.nh <-
    c(1.1, sqrt(1.8 ^ 2 / 12))    # detected, HIV-neg ~U(0.2,2)
  Du.h <-
    c(0.105, sqrt(0.19 ^ 2 / 12)) # not detected, HIV-pos ~U(0.01,0.2)
  Dn.h <-
    c(0.505, sqrt(0.99 ^ 2 / 12)) # detected, HIV-pos ~U(0.01,1)

  Pr <- c(prev, prev.sd)
  Prk <- c(prevk, prev.sd * prevk / prev)
  Ni <- c(newinc, 0)
  H <- c(tbhiv, tbhiv.sd)
  r <- c(rtbhiv, rtbhiv.sd)
  if (is.na(rtbhiv) | is.na(rtbhiv.sd))
    r <- c(1, 0)
  known <- !is.na(prevk)

  if (known) {
    # $I = Sum_{i=u,n} \frac{P_i}{D_i}$
    # P denotes prevalence; D duration; n detected; u undetected
    p <- c(0.025, sqrt(0.05^2/12)) # self cures or dies before tx
    DT <- cbind(Pr, Prk, r, H, p, Dn.nh, Dn.h)
    EXPR <-
      expression((Pr - Prk) /((1 - p) * (r * H * Dn.h + (1 - r * H) * Dn.nh)))
  } else if (Pr[1] > Ni[1] * ((1 - H[1]) * Dn.nh[1] + H[1] * Dn.h[1])) {
    # I = P_u/D_u + newinc, P_u > 0
    DT <-
      cbind(Pr, r, H, Ni, Du.nh, Dn.nh, Du.h, Dn.h)
    EXPR <-
      expression((Pr - Ni * ((1 - H) * Dn.nh + H * Dn.h)) / (r * H * Du.h + (1 -
                                                                               r * H) * Du.nh) + Ni)

  } else
    stop('prev too low compared with newinc')

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


inc2prev <- function(inc,
                     inc.sd,
                     newinc,
                     tbhiv,
                     tbhiv.sd)
{
  #' function to derive incidence from prevalence
  #'
  #' @param inc incidence per capita.
  #' @param inc.sd standard deviation of incidence.
  #' @param newinc detection rate (new+relapse) per capita.
  #' @param tbhiv proportion HIV+ among incident cases.
  #' @param tbhiv.sd SD of proportion HIV+.
  #' @export
  require(propagate) # use second-order Taylor expansion about moments

  stopifnot(inc >= 0 & inc.sd >= 0)
  stopifnot(newinc >= 0)
  stopifnot((tbhiv >= 0 &
               tbhiv <= 1) | is.na(tbhiv))
  stopifnot((tbhiv.sd >= 0 &
               tbhiv.sd < 1) | is.na(tbhiv.sd))

  if (inc <= newinc) {
    warning('inc <= newinc', call. = TRUE, immediate. = TRUE)
    inc <- newinc
  }

  if (is.na(tbhiv))
    tbhiv <- 0
  if (is.na(tbhiv.sd))
    tbhiv.sd <- 0

  # durations
  Du.nh <-
    c(2.5, sqrt(3 / 4))        # not detected, HIV-neg ~U(1,4)
  Dn.nh <-
    c(1.1, sqrt(1.8 ^ 2 / 12))   # detected, HIV-neg ~U(0.2,2)
  Du.h <-
    c(0.105, sqrt(0.19 ^ 2 / 12)) # not detected, HIV-pos ~U(0.01,0.2)
  Dn.h <-
    c(0.505, sqrt(0.99 ^ 2 / 12)) # detected, HIV-pos ~U(0.01,1)

  I <- c(max(inc, newinc), inc.sd)
  Ni <- c(newinc, 0)
  H <- c(tbhiv, tbhiv.sd)

  DT <- cbind(I, H, Ni, Du.nh, Dn.nh, Du.h, Dn.h)
  EXPR <-
    expression(Ni * ((1 - H) * Dn.nh + H * Dn.h) + (I - Ni) * (Du.nh * (1 -
                                                                          H) + Du.h * H))

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

i2p <- function(iso, yr) {
  #' inc2prev wrapper, returns prevalence with uncertainty bounds
  #'
  #' @param iso ISO3 code
  #' @param yr year
  #' @export
  res <- unname(with(
    est[iso3==iso & year==yr],
    inc2prev(
      inc = inc,
      inc.sd = inc.sd,
      newinc = newinc,
      tbhiv = tbhiv,
      tbhiv.sd = tbhiv.sd
    )$prop
  ))
  out <- lohi(res[2] / 1e5, res[4] / 1e5)
  return(c(prev = res[2],
           out[1] * 1e5,
           out[2] * 1e5))
}

inc2mort <- function(inc,
                     inc.sd,
                     newinc,
                     tbhiv,
                     tbhiv.sd,
                     art = NA,
                     art0 = 0,
                     noHIV = TRUE)
{
  #' function to derive mortality from incidence
  #'
  #' @param inc incidence rate.
  #' @param inc.sd SD of incidence rate.
  #' @param newinc detection rate.
  #' @param tbhiv proportion HIV+.
  #' @param tbhiv.sd SD of proportion HIV+.
  #' @param art ART coverage in detected TB.
  #' @param art0 ART coverage in undetected TB.
  #' @param HIV-neg (noHIV=TRUE) or HIV-pos (noHIV=FALSE).
  #' @export
  require(propagate) # use second-order Taylor expansion about moments

  stopifnot(inc >= 0 & inc.sd >= 0)
  stopifnot(newinc >= 0)
  stopifnot((tbhiv >= 0 & tbhiv <= 1) | is.na(tbhiv))
  stopifnot((tbhiv.sd >= 0 & tbhiv.sd < 1) | is.na(tbhiv))
  stopifnot((art >= 0 & art <= 1) | is.na(art))
  stopifnot((art0 >= 0 & art0 <= 1) | is.na(art0))
  stopifnot(noHIV == TRUE | noHIV == FALSE)

  if (is.na(tbhiv))
    tbhiv <- 0
  if (is.na(tbhiv.sd))
    tbhiv.sd <- 0
  if (!noHIV)
    if (is.na(art))
      art <- 0
  if (!noHIV)
    if (is.na(art0))
      art0 <- 0
  if (inc < newinc)
    inc <- newinc

  # CFRs HIV-negative
  CFRu <- c(0.43, (0.53 - 0.28) / 3.92) # undetected
  CFRd <- c(0.03, 0.07 / 3.92)          # detected

  # CFRs HIV-positive
  # noART, untx 0.78 (0.65-0.94)
  CFR.noARTu <- c(0.78, (0.94 - 0.65) / 3.92)

  # noART, tx 0.09 (0.03-0.15)
  CFR.noARTd <- c(0.09, (0.15 - 0.03) / 3.92)

  # ART1, untx 0.62 (0.39-0.86) less than one year
  ART1.u.cfr <- 0.62
  ART1.u.cfr.sd <- (0.86 - 0.39) / 3.92

  # ART1, tx 0.06 (0.01-0.13)
  ART1.n.cfr <- 0.06
  ART1.n.cfr.sd <- (0.13 - 0.01) / 3.92

  # ART2, untx 0.49 (0.31-0.70) more than one year
  ART2.u.cfr <- 0.49
  ART2.u.cfr.sd <- (0.7 - 0.31) / 3.92

  # ART2, tx 0.04 (0.00-0.10)
  ART2.n.cfr <- 0.04
  ART2.n.cfr.sd <- 0.1 / 3.92

  # ART unknown duration, take unweighted average
  CFR.ARTu <-
    c(mean(ART1.u.cfr, ART2.u.cfr),
      mean(ART1.u.cfr.sd, ART2.u.cfr.sd))
  CFR.ARTd <-
    c(mean(ART1.n.cfr, ART2.n.cfr),
      mean(ART1.n.cfr.sd, ART2.n.cfr.sd))

  I <- c(inc, inc.sd)

  # split error about I between treated N and untreated U
  U <-
    c(max(inc - newinc, 0), inc.sd * max(inc - newinc, 0) / inc) # untreated
  N <- c(newinc, inc.sd * newinc / inc) # treated

  H <- c(tbhiv, tbhiv.sd)
  ARTd <- c(art, 0)
  ARTu <- c(art0, 0)

  if (noHIV) {
    DT <- cbind(CFRu, CFRd, U, H, N)
    EXPR <- expression((U * CFRu + N * CFRd) * (1 - H))
  } else {
    DT <-
      cbind(CFR.noARTu,
            CFR.noARTd,
            CFR.ARTu,
            CFR.ARTd,
            U,
            H,
            N,
            ARTd,
            ARTu)
    EXPR <- expression((
      U * (ARTu * CFR.ARTu + (1 - ARTu) * CFR.noARTu) +
        N * (ARTd * CFR.ARTd + (1 - ARTd) * CFR.noARTd)
    ) * H)
  }

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


h2t <- function(hiv = 0,
                hiv.sd = 0,
                irr,
                irr.sd,
                sim = FALSE)
{
  #' function to derive the prevalence of HIV in TB from the IRR
  #'
  #' @param hiv proportion HIV+.
  #' @param hiv.sd SD of proportion HIV+.
  #' @param irr incidence rate ratio (HIV+/HIV-).
  #' @param irr.sd SD of IRR.
  #' @param sim do sim?
  #' @export
  require(propagate) # use second-order Taylor expansion about moments

  stopifnot((hiv >= 0 &
               hiv <= 1 &
               hiv.sd >= 0 &
               hiv.sd <= 1) | is.na(hiv) & is.na(hiv.sd))
  if (is.na(irr))
    irr <- 6
  if (is.na(irr.sd))
    irr.sd <- 1.2

  h <- c(hiv, hiv.sd)
  r <- c(irr, irr.sd)
  DT <- cbind(h, r)
  EXPR <- expression(h * r / (1 + h * (r - 1)))

  out <-
    propagate(
      expr = EXPR,
      data = DT,
      type = 'stat',
      do.sim = sim,
      second.order = T
    )

  return(out)
}



#' Number formatter according to GTB rounding rules
#'
#' Formats vectors of numbers (<2 billion)
#'
#' GTB rounding convention:
#'
#' - 0 is written as "0" (output "0")
#'
#' - values under 0.1 are written "<0.1" ("<0.1")
#'
#' - from 0.1 to under 1 are rounded to 2
#'
#'   significant figures (0.NN)
#'
#' - from 1 to under 10 are rounded to 2 significant figures ("N.N")
#'
#' - 10 to under 100 are rounded to 2 significant figures ("NN")
#'
#' - 100 to under 1000 are rounded to 3 significant figures ("NNN")
#'
#' - 1000 upwards are rounded to 3 significant figures ("N NN0 000")
#'
#' - data that are not reported, but could be are represented
#'   as empty cells and should be accompanied by a footnote.
#'
#' - data that cannot be calculated, either because of missing
#'   data, data was not requested, or any other reason are represented
#'   with a dash.
#'
#' When the number represents thousands, show numbers between 0.01 and 0.1:
#'
#' - values under 0.01 are written "<0.01"
#'
#' - values between 0.01 and under 0.1 are rounded to
#'   2 significant figures ("0.0NN")
#'
#' @param x Vector of numbers
#' @examples
#' ftb(348838)
#' ftb(c(0.0359, 0.00036))
#'
#' @export
#'
ftb <- Vectorize(function(x) {
  # formatter according to GTB rounding rules
  # https://docs.google.com/document/d/1cu_syknBiF3scAX9d7hEfN0LZwdG40i8ttN6yua2xTQ/edit
  #' @param x vector of values
  #' @export
  stopifnot(!is.character(x))
  stopifnot(x < 2e9)

  if (!is.na(x)) {
    smallpos <- x > 0 & x < 0.01
    one2ten <- x >= 1 & x < 10
    zero2one <- x >= 0.1 & x < 1

    dg <- ifelse(abs(x) > 0.01 & abs(x) < 100, 2, 3)
    x2 <- signif(x, dg)

    trailing.0 <- x2 == round2(x) & one2ten == TRUE
    trailing0 <- x2 * 10 == round2(x * 10) & zero2one == TRUE & x2 < 1

    x2 <-
      format(
        x2,
        digits = dg,
        nsmall = 0L,
        big.mark = " ",
        justify = 'right',
        drop0trailing = T,
        scientific = F
      )
    if (smallpos)
      x2 <- '<0.01'
    if (trailing.0)
      x2 <- paste0(x2, '.0')
    if (trailing0)
      x2 <- paste0(x2, '0')
  } else
    x2 <- '-'
  return(x2)
}, 'x')



#' always round 0.5 up
#'
round2 <- function(x, digits = 0) sign(x) * trunc(abs(x) * 10^digits + 0.5) / 10^digits



#' Converts integers into English words,
#'
int2word <- function(x) {
  #' Converts integers into English words,
  #' Based on a function by John Fox:
  #' http://tolstoy.newcastle.edu.au/R/help/05/04/2715.html
  #' @param x an integer
  #' @examples
  #' int2word(324513)
  #' int2word(-3)
  #'
  #' @export
  if (x < 0){
    x <- abs(x)
    neg <- TRUE
  } else {
    neg <- FALSE
  }
  if (x == 0) {
    print("zero")
  } else{
    helper <- function(x) {
      digits <- rev(strsplit(as.character(x), "")[[1]])
      nDigits <- length(digits)
      if (nDigits == 1)
        as.vector(ones[digits])
      else if (nDigits == 2)
        if (x <= 19)
          as.vector(teens[digits[1]])
      else
        trim(paste(tens[digits[2]],
                   Recall(as.numeric(digits[1]))))
      else if (nDigits == 3)
        trim(paste(ones[digits[3]], "hundred and",
                   Recall(makeNumber(digits[2:1]))))
      else {
        nSuffix <- ((nDigits + 2) %/% 3) - 1
        if (nSuffix > length(suffixes))
          stop(paste(x, "is too large!"))
        trim(paste(
          Recall(makeNumber(digits[nDigits:(3 * nSuffix + 1)])),
          suffixes[nSuffix],
          "," ,
          Recall(makeNumber(digits[(3 * nSuffix):1]))
        ))
      }
    }
    trim <- function(text) {
      # Tidy leading/trailing whitespace, space before comma
      text = gsub("^\ ", "", gsub("\ *$", "", gsub("\ ,", ",", text)))
      # Clear any trailing " and"
      text = gsub(" and$", "", text)
      # Clear any trailing comma
      gsub("\ *,$", "", text)
    }
    makeNumber <-
      function(...)
        as.numeric(paste(..., collapse = ""))
    #Disable scientific notation
    opts <- options(scipen = 100)
    on.exit(options(opts))
    ones <-
      c("",
        "one",
        "two",
        "three",
        "four",
        "five",
        "six",
        "seven",
        "eight",
        "nine")
    names(ones) <- 0:9
    teens <-
      c(
        "ten",
        "eleven",
        "twelve",
        "thirteen",
        "fourteen",
        "fifteen",
        "sixteen",
        "seventeen",
        "eighteen",
        "nineteen"
      )
    names(teens) <- 0:9
    tens <-
      c("twenty",
        "thirty",
        "forty",
        "fifty",
        "sixty",
        "seventy",
        "eighty",
        "ninety")
    names(tens) <- 2:9
    x <- round(x)
    suffixes <- c("thousand", "million", "billion", "trillion")
    if (length(x) > 1)
      return(trim(sapply(x, helper)))
    if (neg == TRUE)
      return(paste('minus', helper(x)))
    helper(x)
  }

}




#' prevalence adjusted for se and sp
#'
aprev <- function(tpos, n, se, se.sd, sp, sp.sd) {
  #' bayesian post-hoc prevalence accounting for sensitivity and specifity
  #'
  #' @param tpos number testing positive.
  #' @param n number tested.
  #' @param se sensitivity.
  #' @param se.sd SD of se.
  #' @param sp specificity.
  #' @param sp.sd SD of sp.
  #' @export
  require(rjags)

  stopifnot(tpos >= 0 | n > 0)
  stopifnot(!is.na(se) |
              !is.na(se.sd) | !is.na(sp) | !is.na(sp.sd))
  stopifnot(se > 0 & se < 1)
  stopifnot(se.sd > 0 & se.sd < 1)
  stopifnot(sp > 0 & sp < 1)
  stopifnot(sp.sd > 0 & sp.sd < 1)


  se.par <- get.beta(se, se.sd)
  sp.par <- get.beta(sp, sp.sd)

  # model specs
  m <- "model {
  tpos ~ dbin(theta, n)
  theta <- se*phi + (1-sp)*(1-phi)
  se ~ dbeta(se1, se2)
  sp ~ dbeta(sp1, sp2)
  phi ~ dbeta(1, 1)
}"

  out <- jags.model(
    textConnection(m),
    data = list(
      tpos = tpos,
      n = n,
      se1 = se.par[1],
      se2 = se.par[2],
      sp1 = sp.par[1],
      sp2 = sp.par[2]
    ),
    n.chains = 3
  )

  update(out, 10000)

  mcmc <-
    coda.samples(out, variable.names = c("phi"), n.iter = 20000)
  omcmc <- summary(mcmc)
  post <- omcmc$statistics[1]
  post.lo <- omcmc$quantiles[1]
  post.hi <- omcmc$quantiles[5]
  post.sd <- omcmc$statistics[2]

  return(list(
    mcmc = mcmc,
    res = c(
      post = post,
      post.lo = post.lo,
      post.hi = post.hi,
      post.sd = post.sd
    )
  ))
}


#' shortcut to aprev
#'
trueprev <- function(testpos, se, sp)(testpos + sp - 1) / (se + sp -1)



#' Robust SE of glm objects similar to Stata
#'
robustse <- function(x, coef = "logit") {
  #' Robust standard error for glm objects
  #' returns coefficients as logit (default), OR, prob
  #' adapted from http://stackoverflow.com/questions/27367974/
  #'
  #' @paramx glm object.
  #'

  suppressMessages(require(lmtest))
  suppressMessages(require(sandwich))

  # calculate SE's
  sandwich1 <- function(object, ...)
    sandwich(object) *
    nobs(object) / (nobs(object) - 1)

  # apply to the variance-covariance matrix
  mod1 <- coeftest(x, vcov = sandwich1)


  if (coef == "logit") {
    return(mod1) # return logit with robust SE's
  } else if (coef == "or") {
    mod1[, 1] <- exp(mod1[, 1]) # return odd ratios with robust SE's
    mod1[, 2] <- mod1[, 1] * mod1[, 2]
    return(mod1)
  } else {
    mod1[, 1] <- (mod1[, 1] / 4) # return probabilites with robust SE's
    mod1[, 2] <- mod1[, 2] / 4
    return(mod1)
  }
}




#' Miscellenia utilities
#'
#' recode -- from car
#' example:
#' x<-rep(1:3,3)
#' recode(x, "c(1,2)='A'; else='B'")
#' [1] "A" "A" "B" "A" "A" "B" "A" "A" "B"
#' recode(x, "1:2='A'; 3='B'")
#' [1] "A" "A" "B" "A" "A" "B" "A" "A" "B"

recode <- function (var, recodes, as.factor.result, levels)
{
  recode.list <- rev(strsplit(recodes, ";")[[1]])
  is.fac <- is.factor(var)
  if (missing(as.factor.result))
    as.factor.result <- is.fac
  if (is.fac)
    var <- as.character(var)
  result <- var
  if (is.numeric(var)) {
    lo <- min(var, na.rm = TRUE)
    hi <- max(var, na.rm = TRUE)
  }
  for (term in recode.list) {
    if (0 < length(grep(":", term))) {
      range <- strsplit(strsplit(term, "=")[[1]][1], ":")
      low <- eval(parse(text = range[[1]][1]))
      high <- eval(parse(text = range[[1]][2]))
      target <- eval(parse(text = strsplit(term, "=")[[1]][2]))
      result[(var >= low) & (var <= high)] <- target
    }
    else if (0 < length(grep("else", term))) {
      target <- eval(parse(text = strsplit(term, "=")[[1]][2]))
      result[1:length(var)] <- target
    }
    else {
      set <- eval(parse(text = strsplit(term, "=")[[1]][1]))
      target <- eval(parse(text = strsplit(term, "=")[[1]][2]))
      for (val in set) {
        if (is.na(val))
          result[is.na(var)] <- target
        else
          result[var == val] <- target
      }
    }
  }
  if (as.factor.result) {
    result <- if (!missing(levels))
      factor(result, levels = levels)
    else
      as.factor(result)
  }
  else if (!is.numeric(result)) {
    result.valid <- na.omit(result)
    opt <- options(warn = -1)
    result.valid <- as.numeric(result.valid)
    options(opt)
    if (!any(is.na(result.valid)))
      result <- as.numeric(result)
  }
  result
}



#' Compute sample size to observe at least N events
#' translated from C code (sampsize utility)
#' http://sampsize.sourceforge.net
#'
small.size <- function(obs = 10,
                       pr = 0.5,
                       level = 0.9) {
  mid <- obs
  cubinom <- bot <- 0

  while (cubinom < level) {
    cubinom <- pbeta(pr, obs, mid - obs + 1)
    mid <- mid * 2
  }


  repeat {
    if (cubinom >= level)
      top <- mid
    else
      bot <- mid

    mid <- (bot + top) / 2
    cubinom <- pbeta(pr, obs, mid - obs + 1)

    if (abs(bot - mid) <= 0.5)
      break
  }
  return (ceiling(top))
}




#' Compute sample size based on hypergeometric distribution
#'
n.size <- function(p0 = 0.05,
                   N = 1000,
                   d = 1,
                   alpha = 0.05) {
  s <- N
  for (n in N:1) {
    m <- N - n
    k <- trunc(p0 * N)
    if (dhyper(d, n, m, k) > alpha)
      break
    s <- n
  }
  return(s)
}




#' Survey size -- adapted from epicalc
#'
ssize <-
  function (p,
            delta = "auto",
            popsize = NULL,
            deff = 1,
            alpha = 0.05)
  {
    q <- 1 - p
    pq <- cbind(p, q)
    minpq <- apply(pq, 1, min)
    if (any(delta == "auto")) {
      delta <- ifelse(minpq >= 0.3, 0.1, ifelse(minpq >= 0.1,
                                                0.05, minpq / 2))
    }
    if (any(p >= 1) | any(delta >= 1) | any(popsize < 2))
      stop("p and delta both must be < 1, popsize must be >=2")
    else {
      n1 <- qnorm(1 - alpha / 2) ^ 2 * p * (1 - p) / delta ^ 2
      if (!is.null(popsize)) {
        n1 = n1 / (1 + n1 / popsize)
      }
      if (deff != 1) {
        n1 = n1 * deff
      }
    }
    deff1 <- deff
    if (deff == 1)
      deff1 <- NULL
    tab <- cbind(p, popsize, deff1, delta, alpha, round(n1))
    colnames(tab)[colnames(tab) == "deff1"] <- "deff"
    colnames(tab)[ncol(tab)] <- "n"
    return(as.data.frame(tab))
  }




#' Bayesian sample size, binomial
#' translated from Fulvia Mecatti's Mathematica code
#' in the Lime book
#'
bsize <- function(M = 2000,
                  m = 2000,
                  lo = 10000,
                  hi = 100000,
                  a = 16,
                  b = 7504,
                  d = 0.00053) {
  # Vectorize rbinom() to avoid creating an inner loop
  vrbinom <- Vectorize(rbinom, "prob")

  # generate M random integers between lo and hi
  N <- as.integer(runif(M, min = lo, max = hi))

  # outer loop
  fnab <- numeric(M)
  for (j in 1:M) {
    # generate m random values, distributed Beta(a, b)
    p <- rbeta(m, a, b)

    # inner loop vectorized to improve efficiency
    # generate m random values, distributed Bin(N_j, p_i)
    # stored in an M X m matrix
    x <- vrbinom(m, size = N[j], prob = p)
    k <- 1 / (x + a) + 1 / (N[j] + b - x)
    h <- (N[j] + b - x) / (x + a) + (x + a) / (N[j] + b - x)
    nab <- N[j] + a + b
    dj <- 2 / (nab * sqrt(k)) *
      (
        1.96 - (1.96 ^ 3 + 3 * 1.96) * (h - 1) / (4 * nab) +
          1.96 * h / (2 * nab) +
          5 * (1.96 ^ 3 + 3 * 1.96) * (h - 2) / (18 * nab) -
          1.96 * (h - 2) / nab
      )
    fnab[j] <- mean(dj)
  }

  # OLS: 1 / f(N, a, b) = alpha1 + alpha2 * N
  alpha <- coef(lm(I(1 / fnab ^ 2) ~ N))

  size <- (1 / (4 * d ^ 2) - alpha[1]) / alpha[2]
  names(size) <- 'N'

  return(ceiling(size))
}


#' Frequentist sample size, lime book eq 9.1 and 9.2
#'
ssize1 <-
  function(p1 = 0.002,
           p2 = 0.0014,
           m = 600,
           k1 = 0.3,
           k2 = 0.1,
           beta = 0.8) {
    zb <- qnorm(beta)
    N <-
      (zb * sqrt(p2 * (1 - p2) + (m - 1) * (p2 * k2) ^ 2) + 1.65 * sqrt(p1 *
                                                                          (1 - p1) + (m - 1) * (p1 * k1) ^ 2)) / (p1 - p2)
    N <- ceiling(N ^ 2)
    return(N)
  }



ssize2 <-
  function(pii = 0.0014,
           p1 = 0.002,
           m1 = 600,
           m2 = 600,
           k1 = 0.3,
           k2 = 0.1,
           N1 = 50000,
           beta = 0.8) {
    zb <- qnorm(beta)
    num <-
      (1.65 + zb) ^ 2 * (pii * (1 - pii) + (m2 - 1) * pii ^ 2 * k2 ^
                           2) * N1
    den <-
      (pii - p1) ^ 2 * N1 - (1.65 + zb) ^ 2 * (p1 * (1 - p1) + (m1 - 1) * p1 ^
                                                 2 * k1 ^ 2)
    N <- ceiling(num / den)
    if (N < 0)
      N <- 'undefined'
    return(N)
  }


#' sampling design effects utils
#' k from DEFF
#'
deff2k <- function(deff = 2,
                   pii = 0.002,
                   m = 600) {
  k <- sqrt((deff - 1) * (1 - pii) / ((m - 1) * pii))
  return(k)
}


#' DEFF from k
#'
k2deff <- function(k = 1,
                   pii = 0.002,
                   m = 600) {
  deff <- 1 + (m - 1) * k ^ 2 * pii / (1 - pii)
  return(deff)
}






#' Sample size McNemar test
#' adapted from: Lehr RG. Drug Information Journal 2001;35:1227-1233.
#'
mcn <- function(p1,
                p2,
                alpha = 0.05,
                power = 0.8,
                two.tailed = TRUE) {
  stopifnot(p2 < p1)
  stopifnot(p1 <= 1 & p1 >= 0)
  stopifnot(p2 <= 1 & p2 >= 0)
  stopifnot(alpha > 0 & alpha < 1)
  stopifnot(power > 0 & power < 1)

  p <- (p1 + p2) / 2
  q <- 1 - p
  d <- abs(p2 - p1)
  r <- sqrt(p2 * (1 - p1) / (p1 * (1 - p2)))
  k <- ifelse(two.tailed, 2, 1)
  m = 2 * (qnorm(1 - alpha / k) + (qnorm(power))) ^ 2
  psi <- p1 - p2
  n.lo <- m * p * q * (1 - r) / d ^ 2
  n.hi <- m * p * q / d ^ 2
  n.connor <-
    (qnorm(1 - alpha / k) * sqrt(psi) + qnorm(power) * sqrt(psi - d ^ 2)) ^
    2 / d ^ 2
  return(list(
    p1 = p1,
    p2 = p2,
    psi = psi,
    r = r,
    n.lo = ceiling(n.lo),
    n.connor = ceiling(n.connor),
    n.hi = ceiling(n.hi)
  ))
}








#' binomial CI
#' inspired by stata cii command
#'
cii <- function (size, x, precision, alpha = 0.05)
{
  success <- x
  if (missing(size)) {
    success1 <- success
    if (min(success, na.rm = TRUE) != 0 |
        max(success, na.rm = TRUE) !=
        1) {
      stop("This is not a binary vector.")
    }
    success <- length(na.omit(success1)[na.omit(success1) >
                                          0])
    size <- length(na.omit(success1))
  }
  reverse <- rep(FALSE, length(success))
  reverse[success / size > 0.5] <- TRUE
  success[reverse] <- size[reverse] - success[reverse]
  if (missing(precision)) {
    precision <- success / size / 10000
  }
  precision[success == 0 |
              success == size] <- 0.01 / size[success ==
                                                0 |
                                                success == size]
  probab <- success / size
  success1 <- success
  success1[success > 0] <- success[success > 0] - 1
  for (i in 1:length(success)) {
    while (pbinom(success1[i], size[i], probab[i], lower.tail = FALSE) >
           alpha / 2) {
      probab[i] <- probab[i] - precision[i]
    }
  }
  estimate <- success / size
  se <- sqrt(estimate * (1 - estimate) / size)
  ll <- probab
  probab <- success / size
  for (i in 1:length(success)) {
    while (pbinom(success[i], size[i], probab[i], lower.tail = TRUE) >
           alpha / 2) {
      probab[i] <- probab[i] + precision[i]
    }
  }
  ul <- probab
  data.frame.a <- data.frame(
    events = success,
    total = size,
    prob = estimate,
    se = se,
    ll = ll,
    ul = ul
  )
  data.frame.a[reverse, ] <- data.frame(
    events = size[reverse] -
      success[reverse],
    total = size[reverse],
    prob = 1 -
      estimate[reverse],
    se = se[reverse],
    ll = 1 - ul[reverse],
    ul = 1 - ll[reverse]
  )
  names(data.frame.a)[5] <- paste("lower", 100 * (1 -
                                                    alpha), "ci", sep = "")
  names(data.frame.a)[6] <- paste("upper", 100 * (1 -
                                                    alpha), "ci", sep = "")
  if (nrow(data.frame.a) == 1) {
    rownames(data.frame.a) <- ""
  }
  data.frame.a
}






#' PPV, PPN, PAF
#'
ppv <-
  function(prev = 0.15,
           se = 0.95,
           sp = 0.98)
    se * prev / (se * prev + (1 - sp) * (1 - prev))
ppn <-
  function(prev = 0.15,
           se = 0.95,
           sp = 0.98)
    sp * (1 - prev) / (sp * (1 - prev) + (1 - se) * prev)
paf <- function(pe, rr)
  pe * (rr - 1) / (pe * (rr - 1) + 1)

tprev <- function(x, n, se, sp) (x/n + sp - 1) / (sp + se - 1)



#' capture recapture 3 lists, loglinear
#'
capture <-
  function (A,
            B,
            C,
            AB,
            AC,
            BC,
            ABC,
            deps = 'saturated',
            level = 0.05) {
    #' capture recapture 3 lists, loglinear
    #'
    #' @param A count in A only (not in B, not in C)
    #' @param B count in B only
    #' @param C count in C only
    #' @param AB count in A and B not in C
    #' @param AC count in A and C not in B
    #' @param BC count in B and C not in A
    #' @param ABC count in A and B and C
    #' @param deps model selection, one of: 'AIC','independent','AB','AC','BC','ABAC','ABBC','ACBC','saturated'
    #' @param level confidence level
    #' @export
    stopifnot(A >= 0)
    stopifnot(B >= 0)
    stopifnot(C >= 0)
    stopifnot(AB >= 0)
    stopifnot(AC >= 0)
    stopifnot(BC >= 0)
    stopifnot(ABC >= 0)
    stopifnot(level > 0 & level < 1)

    library(MASS)

    lista <- c(0, 1, 1, 0, 0, 1, 0)
    listb <- c(1, 0, 1, 0, 1, 0, 0)
    listc <- c(1, 1, 0, 1, 0, 0, 0)

    frq <- c(A, B, C, AB, AC, BC, ABC)

    dta <- data.frame(frq, lista, listb, listc)
    if (deps == 'AIC')
      fit <- stepAIC(
        glm(frq ~ lista * listb * listc, family = poisson, data = dta),
        direction = ('backward'),
        trace = 0
      )
    else if (deps == 'independent')
      fit <-
      glm(frq ~ lista + listb + listc, family = poisson, data = dta)
    else if (deps == 'AB')
      fit <-
      glm(frq ~ lista + listb + listc + lista:listb,
          family = poisson,
          data = dta)
    else if (deps == 'AC')
      fit <-
      glm(frq ~ lista + listb + listc + lista:listc,
          family = poisson,
          data = dta)
    else if (deps == 'BC')
      fit <-
      glm(frq ~ lista + listb + listc + lista:listb,
          family = poisson,
          data = dta)
    else if (deps == 'ABAC')
      fit <-
      glm(
        frq ~ lista + listb + listc + lista:listb + lista:listc,
        family = poisson,
        data = dta
      )
    else if (deps == 'ABBC')
      fit <-
      glm(
        frq ~ lista + listb + listc + lista:listb + listb:listc,
        family = poisson,
        data = dta
      )
    else if (deps == 'ACBC')
      fit <-
      glm(
        frq ~ lista + listb + listc + lista:listc + listb:listc,
        family = poisson,
        data = dta
      )
    else if (deps == 'saturated')
      fit <-
      glm(
        frq ~  lista + listb + listc + lista:listb + lista:listc + listb:listc,
        family = poisson,
        data = dta
      )

    n <- sum(frq)
    x <- coef (fit)[1]
    ci.fit <- confint(fit)
    low <- ci.fit[1, 1]
    high <- ci.fit[1, 2]

    n.point <- round(n + exp(x))
    n.low <- floor(n + exp(low))
    n.high <- ceiling(n + exp(high))

    return(list(
      fit = summary(fit),
      estimated = c(n = n.point, lo = n.low, hi = n.high),
      listed = c(
        A = sum(A, AB, AC, ABC),
        B = sum(B, AB, BC, ABC),
        C = sum(C, AC, BC, ABC)
      ),
      inventory = n,
      added = c(n = n.point - n, lo = n.low - n, hi = n.high - n)
    ))
  }





#' misc utils
#'
lp <- function(x)
  #' long print of data.table x
  #'
  #' @param x data.table
  #' @export
  print(x, nrow = Inf)


#' not in operator
#'
`%ni%` <- Negate(`%in%`)




#' lives saved (GF method, counterfactual of no treatment)
#'
cty.lsaved <- function(iso = 'CHN',
                       start = 1995,
                       csv = FALSE) {
  m <- 1e5
  M <- 1e6

  lsaved <- function(dta) {
    # CFRs untreated
    # HIV negative not on TB treatment  0.43 (0.28 - 0.53), assume ~beta
    cfrn <- 0.43
    cfrn.se <- (0.53 - 0.28) / 4
    # HIV positive not on ART, not on TB treatment  0.78 (0.65 - 0.94), assume ~beta
    cfrp <- 0.78
    cfrp.se <- (0.94 - 0.65) / 4

    # lives saved HIV-neg, HIV-pos and total by row
    for (i in 1:dim(dta)[1]) {
      # ugly loop, will optimize later
      # counterfactual (rates)
      cfn <-
        betaop(
          dta$inc.nh[i] / m,
          cfrn,
          dta$inc.nh.se[i] / m,
          cfrn.se,
          op = '*',
          dist = T,
          nsim = M
        )
      cfp <-
        betaop(
          dta$inc.h[i] / m,
          cfrp,
          dta$inc.h.se[i] / m,
          cfrp.se,
          op = '*',
          dist = T,
          nsim = M
        )

      # factual (rates)
      parn <- get.beta(dta$mort.nh[i] / m, dta$mort.nh.se[i] / m)
      fn <- rbeta(M, parn[1], parn[2])
      parp <- get.beta(dta$mort.h[i] / m, dta$mort.h.se[i] / m)
      fp <- rbeta(M, parp[1], parp[2])

      # difference counterfactual minus factual (absolute numbers)
      savedn <- pmax((cfn - fn), 0) * dta$e.pop.num[i]
      savedp <- pmax((cfp - fp), 0) * dta$e.pop.num[i]
      saved <- savedn + savedp

      dta$savedn[i] <- mean(savedn)
      dta$savedn.se[i] <- sd(savedn)
      dta$savedn.lo[i] <- quantile(savedn, probs = 0.025)
      dta$savedn.hi[i] <- quantile(savedn, probs = 0.975)

      dta$savedp[i] <- mean(savedp)
      dta$savedp.se[i] <- sd(savedp)
      dta$savedp.lo[i] <- quantile(savedp, probs = 0.025)
      dta$savedp.hi[i] <- quantile(savedp, probs = 0.975)

      dta$saved[i] <- mean(saved)
      dta$saved.se[i] <- sd(saved)
      dta$saved.lo[i] <- quantile(saved, probs = 0.025)
      dta$saved.hi[i] <- quantile(saved, probs = 0.975)
    }

    return(dta)
  }

  clsaved <- function(dta) {
    # cumulative lives saved
    savedn <- sum(dta$savedn)
    savedn.se <- sqrt(sum(dta$savedn.se ^ 2))
    savedn.lo <- savedn - 1.96 * savedn.se
    savedn.hi <- savedn + 1.96 * savedn.se

    savedp <- sum(dta$savedp)
    savedp.se <- sqrt(sum(dta$savedp.se ^ 2))
    savedp.lo <- savedp - 1.96 * savedp.se
    savedp.hi <- savedp + 1.96 * savedp.se

    saved <- sum(dta$saved)
    saved.se <- sqrt(sum(dta$saved.se ^ 2))
    saved.lo <- saved - 1.96 * saved.se
    saved.hi <- saved + 1.96 * saved.se
    dta <- data.table(
      savedn,
      savedn.lo,
      savedn.hi,
      savedn.se,
      savedp,
      savedp.lo,
      savedp.hi,
      savedp.se,
      saved,
      saved.lo,
      saved.hi,
      saved.se
    )
    return(dta)
  }

  th <- 1000
  bgd <- lsaved(est[iso3 == iso & year >= start])
  bgd.cum <- clsaved(bgd)
  bgd.saved <-
    bgd[, list(
      year = as.character(year),
      savedn,
      savedn.lo,
      savedn.hi,
      savedn.se,
      savedp,
      savedp.lo,
      savedp.hi,
      savedp.se,
      saved,
      saved.lo,
      saved.hi,
      saved.se
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
    write.csv(
      bgd.saved,
      file = paste('output/', iso, '_livesSaved.csv', sep = ''),
      row.names = FALSE
    )
  return(bgd.saved)
}




#' install from github behind WHO firewall
#'
install_zip <- function(path,
                        pkg = sub("(-[^-]+)?\\.[^.]+$", "", basename(path))) {
  dir1 <- tempfile()
  dir2 <- file.path(dir1, pkg)
  dir.create(dir2, recursive = TRUE)
  on.exit(unlink(dir1, recursive = TRUE, force = TRUE))
  suppressMessages(unzip(path, exdir = dir2,
                         unzip = getOption("unzip")))
  temp_contents <- dir(dir2)
  if (length(temp_contents) == 1L) {
    dir3 <- file.path(dir2, temp_contents)
    if (file.info(dir3, extra_cols = FALSE)[["isdir"]]) {
      dir2 <- file.path(dir2, pkg)
      file.rename(dir3, dir2)
    }
  }
  install.packages(dir2, repos = NULL, type = "source")
}


#' map colour names
#'
col <- function() {
  d = data.frame(
    c = colors(),
    y = seq(0, length(colors()) - 1) %% 66,
    x = seq(0, length(colors()) - 1) %/% 66
  )
  p <- ggplot() +
    scale_x_continuous(name = "",
                       breaks = NULL,
                       expand = c(0, 0)) +
    scale_y_continuous(name = "",
                       breaks = NULL,
                       expand = c(0, 0)) +
    scale_fill_identity() +
    geom_rect(
      data = d,
      mapping = aes(
        xmin = x,
        xmax = x + 1,
        ymin = y,
        ymax = y + 1
      ),
      fill = "white"
    ) +
    geom_rect(
      data = d,
      mapping = aes(
        xmin = x + 0.05,
        xmax = x + 0.95,
        ymin = y + 0.5,
        ymax = y + 1,
        fill = c
      )
    ) +
    geom_text(
      data = d,
      mapping = aes(
        x = x + 0.5,
        y = y + 0.5,
        label = c
      ),
      colour = "black",
      hjust = 0.5,
      vjust = 1,
      size = 3
    )
  print(p)
}



#' ---
#' title: Reformat UNAIDS estimates
#' author: Philippe Glaziou
#' date: 2022-06-16
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

#' Last updated: `r Sys.Date()`
#'
#' # Reformat UNAIDS estimates
#'
#' This script converts UNAIDS HIV incidence and mortality estimates (Excel file) to wide format rda
#' and dated csv files which can then be used by later scripts.
#'
#' Output:
#'
#' - unaids:  (as dated csv file and as rda file)
#'


# Load libraries and data

library(data.table)
library(bit64)
library(readxl)
library(here)

load('data/cty.rda')

#' vectorize lohi
#'
source(here('R/fun.R'))
vlohi <- Vectorize(lohi, c('ev', 'sd'))
m <- 1e5

dic <- fread('input/UNAIDS/dic.csv')

typz <- c(rep('numeric', 1),
          rep('text', 2),
          rep('numeric', 60))

suppressWarnings(unaids <-
                   as.data.table(
                     read_excel(
                       'input/UNAIDS/HIV2022Estimates_AllCountries_7July22.xlsx',
                       sheet = 1,
                       skip = 2,
                       col_types = typz
                     )
                   ))

setnames(unaids, dic$varnames)

sel <- unaids$country %ni% c('Global') & nchar(unaids$iso3, keepNA=F) == 3
dim(unaids); dim(unaids[sel])
unaids <- unaids[sel, ]

setkey(unaids, iso3, year)
unaids[, country := NULL]

unaids[, pop := pop014 + pop15]
unaids[, hiv.num := h014 + h15]

(unaids[is.na(iso3)])

unaids <- within(unaids, {
  hiv15 <- h15 / pop15
  hiv15.lo <- h15.lo / pop15
  hiv15.hi <- h15.hi / pop15
  hiv014 <- h014 / pop
  hiv014.lo <- h014.lo / pop
  hiv014.hi <- h014.hi / pop
  inchiv <- (inch014 + inch15) / pop
  inchiv.num <- inch014 + inch15
  d <- rowSums(cbind(d15, d014), na.rm = T)
})

unaids <- within(unaids, {
  hiv <- hiv.num / pop
  hiv.sd <-
    sqrt(((h014.hi - h014.lo) / 3.92) ^ 2 + ((h15.hi - h15.lo) /
                                               3.92) ^ 2) / pop
  inchiv.sd <-
    sqrt(((inch014.hi - inch014.lo) / 3.92) ^ 2 + ((inch15.hi - inch15.lo) /
                                                     3.92) ^ 2) / pop
  mort.hiv <-
    d / pop
  mort.hiv.sd <-
    sqrt(((d15.hi - d15.lo) / 3.92) ^ 2 + ((d014.hi - d014.lo) /
                                             3.92) ^ 2) / pop
})

sel <-
  unaids$hiv > 0 &
  !is.na(unaids$hiv) &
  unaids$hiv.sd > 0 &
  unaids$hiv.sd^2 < unaids$hiv * (1 - unaids$hiv)
(table(sel))
out <- vlohi(unaids$hiv[sel], unaids$hiv.sd[sel])
unaids[sel, hiv.lo := out[1,]]
unaids[sel, hiv.hi := out[2,]]
unaids[hiv>hiv.hi, hiv.hi := hiv + 1.96 * hiv.sd]
unaids[sel, test.bounds(hiv, hiv.lo, hiv.hi)]


sel <-
  unaids$mort.hiv > 0 &
  unaids$mort.hiv.sd > 0 & !is.na(unaids$mort.hiv) &
  unaids$mort.hiv.sd^2 < unaids$mort.hiv * (1 - unaids$mort.hiv)
(table(sel))
out <- vlohi(unaids$mort.hiv[sel], unaids$mort.hiv.sd[sel])
unaids$mort.hiv.lo[sel] <- out[1,]
unaids$mort.hiv.hi[sel] <- out[2,]
unaids[mort.hiv > mort.hiv.hi, mort.hiv.hi := mort.hiv + 1.96 * mort.hiv.sd]
unaids[sel, test.bounds(mort.hiv, mort.hiv.lo, mort.hiv.hi)]

(unaids$mort.hiv[!sel])
(unaids$mort.hiv.sd[!sel])
unaids[!sel, mort.hiv.lo := 0]
unaids[!sel, mort.hiv.hi := mort.hiv + mort.hiv.sd * 1.96]

unaids <- within(unaids, {
  mort.hiv <- mort.hiv * m
  mort.hiv.lo <- mort.hiv.lo * m
  mort.hiv.hi <- mort.hiv.hi * m
  mort.hiv.sd <- mort.hiv.sd * m
})

unaids <- within(unaids, {
  hiv.lo.num <- hiv.lo * pop
  hiv.hi.num <- hiv.hi * pop
  mort.hiv.num <- d
  mort.hiv.lo.num <- mort.hiv.lo * pop / m
  mort.hiv.hi.num <- mort.hiv.hi * pop / m
})


#' check aggregates
#'
(unaids[, sums(hiv.num), by=year])


#' save
#'
save(unaids, file = 'data/unaids.rda')
fwrite(unaids, file = paste0('csv/unaids_', Sys.Date(), '.csv'))

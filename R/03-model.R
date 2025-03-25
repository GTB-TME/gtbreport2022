#' ---
#' title: Reformat model output
#' author: Philippe Glaziou
#' date: 2022/07/18
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
#' # Reformat output from Nim's models
#'
#' This script converts the output from Nim’s model (Excel files) to long format rda
#' and dated csv files which can then be used by later scripts.
#'
#' Output:
#'
#' - model:   Annual estimates for countries where dynamic modelling was used
#' - model2:  Monthly estimates for countries where dynamic modelling was used
#' - extra:   Annual estimates for countries where a regional model was used ('extrapolated')
#'



# Load libraries and data
library(data.table)
library(readxl)
library(here)

load('data/old.rda')
load('data/pop.rda')

source('R/fun.R')


# check that the latest model outputs are used
# xlfn <- 'Nim/Focal_countries.xlsx'
# xlfn <- 'Nim/Focal_countries_19Jul22.xlsx'
# xlfn <- 'Nim/Focal_countries_28Jul22.xlsx'
# xlfn <- 'Nim/Focal_countries_29Jul22.xlsx'
xlfn <- 'Nim/Focal_countries_05Sept22.xlsx'
xlfn2 <- 'Nim/Extrapolated_countries_07Aug22.xlsx'





# Import "focus" countries, annual
whiv <- rep(c('BRA','COL','KEN','LSO','PNG','THA','ZAF','ZWE'), each = 6)

# format xls data
fmt <- function(M, snames) {
  setnames(M, new = snames)
  M[, index := rep(1:(nrow(M) / 6), each = 6)]
  M[, index2 := rep(1:(nrow(M) / 3), each = 3)]

  M[, iso3 := iso3[1], by=index]
  M[is.na(iso3), iso3 := whiv]

  M[iso3 %in% unique(whiv), scenario := scenario[1], by = index]
  M[iso3 %ni% unique(whiv), scenario := scenario[1], by = index2]

  M[, hiv := hiv[1], by = index2]
  M[is.na(hiv), hiv := 'All']
  M[, hiv := gsub('All','a',hiv)]
  M[hiv != 'a', hiv := 'pos']

  M[, index := NULL]
  M[, index2 := NULL]

  setkey(M, iso3)
  return(M)
}


# incidence, annual
M1 <-
  as.data.table(read_excel(xlfn,
                           sheet = 2,
                           skip = 0))
snames <- c(
  'iso3',
  'scenario',
  'hiv',
  'stat',
  '2019','2020','2021','2022','2023','2024'
)

M1 <- fmt(M1, snames)


# mortality, annual
M2 <-
  as.data.table(read_excel(xlfn,
                           sheet = 4,
                           skip = 0))

M2 <- fmt(M2, snames)






# reshape to long
fmt2 <- function(m) {
  .m <- melt(m, id.vars = c(1:4))
  .m2 <- dcast(.m, ... ~ stat)
  return(.m2)
}


m1 <- fmt2(M1)
m2 <- fmt2(M2)
m2[, hi := as.numeric(hi)]
m2[, lo := as.numeric(lo)]
m2[, md := as.numeric(md)]



# combine
setnames(m1, old = 'variable', new = 'year')
m1[, measure := 'inc']
setnames(m2, old = 'variable', new = 'year')
m2[, measure := 'mort']
model <-
  rbind(m1[, .(
    iso3,
    scenario,
    hiv,
    measure,
    year,
    best = md,
    lo,
    hi
  )],
  m2[, .(
    iso3,
    scenario,
    hiv,
    measure,
    year,
    best = md,
    lo,
    hi
  )])



model[, year := as.integer(as.character(year))]

setkey(model, iso3, year)

fwrite(model, file = here(paste0("csv/model_", Sys.Date(), '.csv')))
save(model, file = here('data/model.rda'))




# import "extrapolated" countries
whiv2 <- rep(c('BWA','NAM','SWZ','SYC'), each = 6)

# format xl data
fmt3 <- function(M, snames) {
  setnames(M, new = snames)
  M[, index := rep(1:(nrow(M) / 6), each = 6)]
  M[, index2 := rep(1:(nrow(M) / 3), each = 3)]

  M[, iso3 := iso3[1], by=index]
  M[is.na(iso3), iso3 := whiv2]

  M[iso3 %in% unique(whiv2), scenario := scenario[1], by = index]
  M[iso3 %ni% unique(whiv2), scenario := scenario[1], by = index2]

  M[, hiv := hiv[1], by = index2]
  M[is.na(hiv), hiv := 'All']
  M[, hiv := gsub('All','a',hiv)]
  M[hiv != 'a', hiv := 'pos']

  M[, index := NULL]
  M[, index2 := NULL]
  M[, region := NULL]

  setkey(M, iso3)
  return(M)
}


# extrapolated incidence, annual
E1 <-
  as.data.table(read_excel(xlfn2,
                           sheet = 2,
                           skip = 0))
snames <- c(
  'region',
  'iso3',
  'scenario',
  'hiv',
  'stat',
  '2019','2020','2021','2022','2023','2024'
)

E1 <- fmt3(E1, snames)


# extrapolated mortality, annual
E2 <-
  as.data.table(read_excel(xlfn2,
                           sheet = 4,
                           skip = 0))

E2 <- fmt3(E2, snames)


# reshape to long
e1 <- fmt2(E1)
e2 <- fmt2(E2)


# combine
setnames(e1, old = 'variable', new = 'year')
e1[, measure := 'inc']
setnames(e2, old = 'variable', new = 'year')
e2[, measure := 'mort']
extra <-
  rbind(e1[, .(
    iso3,
    scenario,
    hiv,
    measure,
    year,
    best = md,
    lo,
    hi
  )],
  e2[, .(
    iso3,
    scenario,
    hiv,
    measure,
    year,
    best = md,
    lo,
    hi
  )])

extra[scenario=='With COVID', scenario := 'COVID']
extra[, year := as.integer(as.character(year))]

setkey(extra, iso3, year)

fwrite(extra, file = here(paste0("csv/extra_", Sys.Date(), '.csv')))
save(extra, file = here('data/extra.rda'))





# Import monthly estimates in "focus" countries
#
# Import "focus" countries, annual
whiv <- rep(c('BRA','COL','KEN','LSO','PNG','THA','ZAF','ZWE'), each = 6)

# format xl data
mfmt <- function(M, snames) {
  setnames(M, new = snames)
  M[, index := rep(1:(nrow(M) / 6), each = 6)]
  M[, index2 := rep(1:(nrow(M) / 3), each = 3)]

  M[, iso3 := iso3[1], by=index]
  M[is.na(iso3), iso3 := whiv]

  M[iso3 %in% unique(whiv), scenario := scenario[1], by = index]
  M[iso3 %ni% unique(whiv), scenario := scenario[1], by = index2]

  M[, hiv := hiv[1], by = index2]
  M[is.na(hiv), hiv := 'All']
  M[, hiv := gsub('All','a',hiv)]
  M[hiv != 'a', hiv := 'pos']

  M[, index := NULL]
  M[, index2 := NULL]

  setkey(M, iso3)
  return(M)
}


# incidence, monthly
mM1 <-
  as.data.table(read_excel(xlfn,
                           sheet = 1,
                           skip = 0))
snames <- c(
  'iso3',
  'scenario',
  'hiv',
  'stat',
  as.character(as.Date(as.integer(names(
    mM1
  )[5:76]), origin = "1899-12-30"))
)

mM1 <- mfmt(mM1, snames)


# mortality, annual
mM2 <-
  as.data.table(read_excel(xlfn,
                           sheet = 3,
                           skip = 0))

mM2 <- mfmt(mM2, snames)






# reshape to long
fmt2 <- function(m) {
  .m <- melt(m, id.vars = c(1:4))
  .m2 <- dcast(.m, ... ~ stat)
  return(.m2)
}


mm1 <- fmt2(mM1)
mm2 <- fmt2(mM2)
mm2[, hi := as.numeric(hi)]
mm2[, lo := as.numeric(lo)]
mm2[, md := as.numeric(md)]



# combine
setnames(mm1, old = 'variable', new = 'date')
mm1[, measure := 'inc']
setnames(mm2, old = 'variable', new = 'date')
mm2[, measure := 'mort']
model2 <-
  rbind(mm1[, .(
    iso3,
    scenario,
    hiv,
    measure,
    date,
    best = md,
    lo,
    hi
  )],
  mm2[, .(
    iso3,
    scenario,
    hiv,
    measure,
    date,
    best = md,
    lo,
    hi
  )])

model2[, year := year(as.Date(date))]
setkey(model2, iso3, year)

fwrite(model2, file = here(paste0("csv/model2_", Sys.Date(), '.csv')))
save(model2, file = here('data/model2.rda'))



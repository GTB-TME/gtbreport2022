#' ---
#' title: init GTB2022
#' author: Philippe Glaziou
#' date: 2022-05-17
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---

#' # Description
#'
#' This script is run once only at the beginning while setting up the data for the repository.
#' It uses the final set of estimates from the 2021 global TB report as a starting point.
#'
#' Output:
#' - old:    A copy incidence and mortality estimates from last year's repository
#' - ovr:    A copy of VR data from last year's repository
#'
#' Also outputs these files, though not sure why because first two
#' are downloaded again from the database (script ~/R/01-odbc.R) and the final one
#' is a subset of ovr:
#'
#' - cty:    Country names and iso codes
#' - pop:    Population estimates
#' - vrcov:  A subset of last year's VR data (coverage and quality)
#'

library(data.table)

load('../gtb2021/data/est.rda')
load('../gtb2021/data/cty.rda')
load('../gtb2021/data/pop.rda')
load('../gtb2021/data/vr.rda')
old <- copy(est)
ovr <- copy(vr)
vrcov <- vr[,.(iso3, year, vr.coverage, codqual)]


save(old, file = 'data/old.rda')
save(cty, file = 'data/cty.rda')
save(pop, file = 'data/pop.rda')
save(vrcov, file = 'data/vrcov.rda')
save(ovr, file = 'data/ovr.rda')

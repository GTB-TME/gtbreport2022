

# Set the report year ----
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
report_year <- 2022

# Set whether or not to include objects with estimates in the output ----
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
show_estimates <- TRUE

# Kill any attempt at using factors, unless we explicitly want them!
options(stringsAsFactors=FALSE)


# Load packages ----
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

library(tidyverse)
library(dplyr)
library(tidyr)
library(readxl)
library(magrittr)
library(scales)
library(stringr)
#devtools::install_github('glaziou/gtbreport')
library(gtbreport)
library(metafor)
library(here)


# Load functions ----
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


# Load data if not already loaded----
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Function to read in timestamped CSV file and to undo Philippe's conversion of underscores to dots
get_timestamped_csv <- function(csv_filename_stub, timestamp) {
  df <- read.csv(here(paste0('./csv/', csv_filename_stub, '_', timestamp,'.csv')))
  names(df) <- gsub('[.]', '_', names(df))
  return(df)
}

data_loc <- "./ch6_data/"
csv_datestamp <- '2022-08-29'
csv_datestamp2 <- '2022-10-04' # for sdg data updated on 4 Oct
csv_datestamp3 <- '2022-11-08' # for risk factor updated on 8 Nov
csv_estimate_datestamp <- '2022-09-13'

sdg    <- get_timestamped_csv("sdg",csv_datestamp2)
sdgdef    <- get_timestamped_csv("sdgdef",csv_datestamp)
grpmbr    <- get_timestamped_csv("grpmbr",csv_datestamp)
tb    <- get_timestamped_csv("tb",csv_datestamp)
est    <- get_timestamped_csv("./db/db_est_country",csv_estimate_datestamp)
dic    <- get_timestamped_csv("dic",csv_datestamp)
pop    <- get_timestamped_csv("pop",csv_datestamp)
rf_global    <- get_timestamped_csv("./db/db_inc_risk_factor_global",csv_datestamp3)
rf_country    <- get_timestamped_csv("./db/db_inc_risk_factor_country",csv_datestamp3)

# load(here::here("./report/ch6_data/att.Rdata"))
# load(here::here("./data/att.rda"))
# load(here::here("./report/ch6_data/att_agg.Rdata"))

hgf <- read.csv(here::here("./csv/TBreport_SDG381382.csv"))

# Get Pete's aggregate incidence estimates by age group and sex
load(here('disaggregation/dboutput/db_estimates_country.Rdata'))
load(here('disaggregation/dboutput/db_estimates_group.Rdata'))

#The 30 high TB burden countries – 3 new entries and 3 exits
#• Gabon, Mongolia and Uganda added.
#• Cambodia, the Russian Federation and Zimbabwe removed.

# list iso3 - country
list_iso3_country <-
  grpmbr %>% filter(group_type == 'g_whoregion') %>%
  select(iso3,country)

# HBCs exit 2021 = watchlist countries
list_watch_list <-
  list_iso3_country %>%
  filter(iso3 %in% c('KHM','RUS','ZWE'))

list_hbcs <-
  grpmbr %>%
  filter(group_type == "g_hb_tb" & group_name == 1) %>%
  select(iso3,country) %>%
  arrange(country)

list_hbcs_plus_wl <- list_hbcs %>% add_row(list_watch_list) %>%   arrange(country)

iso3_hbc <- list_hbcs$iso3
iso3_hbc_plus_wl <- list_hbcs_plus_wl$iso3

iso3_income <- grpmbr %>%
  filter(group_type == "g_income") %>% select(iso3,group_name) %>% rename(g_income=2)

list_hbcs_income <-
  list_hbcs %>%
  left_join(iso3_income)

list_hbcs_plus_wl_income <-
  list_hbcs_plus_wl %>%
  left_join(iso3_income)

dummydf <-
  list_hbcs_income %>%
  mutate(year=0,value=0) %>%
  select(iso3,year,value,country,g_income)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Common data preparation script for chapters (sections) 3 and 4
# Hazim Timimi, July 2022
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Load data packages ----
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

library(stringr)
library(dplyr)
library(tidyr)
library(here)
library(readr) # to save csv
library(magrittr) # to use tee pipe
library(data.table)
library(jsonlite) # for provisinal monthly notification in India

# Set the report year ----
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
report_year <- 2022

# Set whether or not to include objects with estimates in the output ----
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

show_estimates <- TRUE

# Set datestamps of CSV files to use ----
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csv_datestamp <- '2022-08-29'
csv_datestamp2 <- '2022-09-19' # tentative for live saved csv file as of 10 Aug
csv_estimate_datestamp <- '2022-09-27'

csv_tpt_fix_datestamp <- '2022-09-30' # Contains late fixes for TPT in PLHIV from GAM, nothing else


# Load TB data (CSV,not rda because the latter are datatables ...) ----
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Function to read in timestamped CSV file and to undo Philippe's conversion of underscores to dots
get_timestamped_csv <- function(csv_filename_stub, timestamp) {
  df <- read.csv(here(paste0('csv/', csv_filename_stub, '_', timestamp,'.csv')))
  names(df) <- gsub('[.]', '_', names(df))
  return(df)
}

notification <- get_timestamped_csv('tb', csv_tpt_fix_datestamp)
data_collection <- get_timestamped_csv('datacoll', csv_datestamp)
strategy <- get_timestamped_csv('sty', csv_datestamp)
country_group_membership <- get_timestamped_csv('grpmbr', csv_datestamp)
TBHIV_for_aggregates <- get_timestamped_csv('agg', csv_datestamp)
dr_surveillance <- get_timestamped_csv('dr', csv_datestamp)
outcomes  <- get_timestamped_csv('tx', csv_datestamp)

estimates_population <- get_timestamped_csv('pop', csv_datestamp)


# Tweak dr_surveillance: remove records that aren't explicitly covering a whole country
dr_surveillance <- filter(dr_surveillance, all_areas_covered==1 | is.na(all_areas_covered))

# set e_pop_num to numeric to avoid integer overflow error
estimates_population$e_pop_num <- as.numeric(estimates_population$e_pop_num)


# Fix lists of the three sets of 30 high burden countries (used to filter records for some figures)
hbc30 <- country_group_membership %>%
  filter(group_type == "g_hb_tb" & group_name == 1) %>%
  select(iso3,group_type)

hbmdr30 <- country_group_membership %>%
  filter(group_type == "g_hb_mdr" & group_name == 1) %>%
  select(iso3,group_type)

hbtbhiv30 <- country_group_membership %>%
  filter(group_type == "g_hb_tbhiv" & group_name == 1) %>%
  select(iso3,group_type)

# list iso3 - country
list_iso3_country <-
  country_group_membership %>% filter(group_type == 'g_whoregion') %>%
  select(iso3,country)

# HBCs exit 2021 = watchlist countries
list_watch_list <-
  list_iso3_country %>%
  filter(iso3 %in% c('KHM','RUS','ZWE'))

list_hbcs <-
  country_group_membership %>%
  filter(group_type == "g_hb_tb" & group_name == 1) %>%
  select(iso3,country) %>%
  arrange(country)

list_hbcs_plus_wl <- list_hbcs %>% add_row(list_watch_list) %>%   arrange(country)


# Get global and regional estimates directly from the files that get imported into the database
if(show_estimates) {

  est_country  <- read.csv(here(paste0('csv/db/db_est_country_', csv_estimate_datestamp, '.csv')))
  est_regional <- read.csv(here(paste0('csv/db/db_est_regional_est_', csv_estimate_datestamp, '.csv')))
  est_global <- read.csv(here(paste0('csv/db/db_est_global_est_', csv_estimate_datestamp, '.csv')))

  est_dr_country <- get_timestamped_csv('dr.est', csv_datestamp)
  est_dr_group <- get_timestamped_csv('aggdr', csv_datestamp)

  lives_saved <- read.csv(here(paste0('output/RegionalLivesSaved_', csv_datestamp2, '.csv')))

  # Get Pete's aggregate incidence estimates by age group and sex
  load(here('disaggregation/dboutput/db_estimates_country.Rdata'))
  load(here('disaggregation/dboutput/db_estimates_group.Rdata'))

}

# Create a set of WHO region short names to use in figures and tables
who_region_shortnames <- get_timestamped_csv('grp', csv_datestamp) %>%
  filter(group_type=="g_whoregion") %>%
  select(g_whoregion = group_name, entity= group_description) %>%
  # Remove the WHO bit at the beginning of regional names
  mutate(entity = gsub("(WHO )|(WHO/PAHO )", "", entity)) %>%
  # Change the order based on the bizarre method chosen by WHO Press ...
  mutate(entity = factor(entity,
                         levels = c("African Region", "Region of the Americas", "South-East Asia Region",
                                    "European Region", "Eastern Mediterranean Region", "Western Pacific Region")))


# Add simple function to convert NA to zero
NZ <- function(x) {
  ifelse(is.na(x),
         0,
         x)
}


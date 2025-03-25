# - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - -
# Data preparation script for ch2-3.rmd on drug resistance surveillance
# Anna Dean, adapted from Olga Tosas Auguet,
# with some extra tweaks by Hazim Timimi to switch to using the CSV files
# within the repository and to cut use of out intermediary files. Made
# it use dataframes numbered according to the figure numbers, as already
# done for other chapters
#
# July 2021 - August 2022
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


# Load data packages ----
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

library(stringr)
library(dplyr)
library(tidyr)
library(here)


# Set the report year ----
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
report_year <- 2022

# Set whether or not to include objects with estimates in the output ----
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

show_estimates <- TRUE

# Set datestamps of CSV files to use ----
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csv_datestamp <- '2022-08-29'
csv_estimate_datestamp <- '2022-09-14'

# Need the next one to calculation proportions in the text
csv_inc_estimate_datestamp <- '2022-09-27'


# Load TB data (CSV,not rda because the latter are datatables ...) ----
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Function to read in timestamped CSV file and to undo Philippe's conversion of underscores to dots
get_timestamped_csv <- function(csv_filename_stub, timestamp) {
  df <- read.csv(here(paste0('csv/', csv_filename_stub, '_', timestamp,'.csv')))
  names(df) <- gsub('[.]', '_', names(df))
  return(df)
}

# Philippe saved view_DRS_most_recent_for_estimation as dre.csv
most_recent_for_estimation  <- get_timestamped_csv('dre', csv_datestamp)
dr_surveillance <- get_timestamped_csv('dr', csv_datestamp)
drs <- get_timestamped_csv('drs', csv_datestamp)
drs_exclusions <- get_timestamped_csv('drsexcl', csv_datestamp)

pop <- get_timestamped_csv('pop', csv_datestamp)

# Tweak dr_surveillance: remove records that aren't explicitly covering a whole country
dr_surveillance <- filter(dr_surveillance, all_areas_covered==1 | is.na(all_areas_covered))

# Get global and regional DR-TB estimates directly from the files that get imported into the database
if(show_estimates) {

  est_dr_country <- read.csv(here(paste0('drtb/dboutput/db_dr_country_', csv_estimate_datestamp, '.csv')))
  est_dr_group <- read.csv(here(paste0('drtb/dboutput/db_dr_group_', csv_estimate_datestamp, '.csv')))

  # Get incidence estimates to calculate proportion in the text
  est_country  <- read.csv(here(paste0('csv/db/db_est_country_', csv_inc_estimate_datestamp, '.csv')))

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


# Get list of the 30 high MDR-TB burden countries (used to filter records for some figures)
country_group_membership <- get_timestamped_csv('grpmbr', csv_datestamp)
hbmdr30 <- country_group_membership %>%
  filter(group_type == "g_hb_mdr" & group_name == 1) %>%
  select(country, iso3)

# Get list of 43 countries in both high-TB and high-MDR lists for calculating stats in the text
hbtb_hbmdr <- country_group_membership %>%
  filter(group_type %in% c("g_hb_mdr", "g_hb_tb") & group_name == 1) %>%
  select(iso3) %>%
  unique()

# Load the country names
country_names <- get_timestamped_csv('cty', csv_datestamp) %>%
  select(iso3, country)




# pre-XDR estimates ----
# These were produced by Anna but not included in the database, so will hard-code them here for now :-(
#
est_rr_fdr_global <- data.frame(best = 20, lo = 15.5, hi = 25.7)





#---------------- Fig 2.3.1 and 2.3.3 Panel plot of RR-TB incidence estimates by WHO region and globally since 2015  -------------#
# Late renumbering from 2.3.4 to 2.3.1 (global) and 2.3.3 (regional)

f2_3_01_data <- est_dr_group %>%
  filter(group_type == 'global' & year>=2015) %>%
  mutate(entity = 'Global') %>%
  select(year,
         entity,
         e_inc_rr_num,
         e_inc_rr_num_lo,
         e_inc_rr_num_hi)

# Summary dataset for simple quoting of numbers in the text
f2_3_01_txt <- f2_3_01_data %>%
  arrange(year) %>%
  filter(year >= report_year - 2) %>%
  # Calculate % change between the last two years
  mutate(previous = lag(e_inc_rr_num)) %>%
  mutate(pct_diff = abs(e_inc_rr_num - previous)*100/previous)



f2_3_03_data <- est_dr_group %>%
  filter(group_type == 'g_whoregion' & year>=2015) %>%
  select(year,
         g_whoregion = group_name,
         e_inc_rr_num,
         e_inc_rr_num_lo,
         e_inc_rr_num_hi) %>%

  # merge with regional names
  inner_join(who_region_shortnames, by = "g_whoregion") %>%
  select(-g_whoregion) %>%

  # Set the entity order for plotting
  mutate(entity = factor(entity,
                         levels = c("African Region", "Region of the Americas", "South-East Asia Region",
                                    "European Region", "Eastern Mediterranean Region", "Western Pacific Region")))





# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Create dataframe for figure 2.3.2 ----
# (Panel plot of global proportion of TB cases with MDR/RR-TB)
# Late renumbering from 2.3.5
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

f2_3_02_data <- est_dr_group %>%
  filter(group_name=="global") %>%
  mutate(pct_new_best    = e_rr_prop_new * 100,
         pct_new_lo = e_rr_prop_new_lo * 100,
         pct_new_hi = e_rr_prop_new_hi * 100,
         pct_ret_best    = e_rr_prop_ret * 100,
         pct_ret_lo = e_rr_prop_ret_lo * 100,
         pct_ret_hi = e_rr_prop_ret_hi * 100) %>%
  select(year, starts_with("pct_")) %>%

  # Switch to a long format
  pivot_longer(cols = starts_with("pct_"),
               names_to = c("pct", "case_type", "val"),
               names_sep = "_") %>%
  select(-pct) %>%

  # Switch back to wide but keep case type as identifier
  pivot_wider(names_from = val,
              values_from = value) %>%

  # Change case type to a factor with descriptive names
  mutate(case_type = factor(case_type,
                            levels = c("new", "ret"),
                            labels = c("New TB cases", "Previously treated TB cases")))


# Summary dataset for simple quoting of numbers in the text
f2_3_02_txt <- f2_3_02_data %>%
  filter(year %in% c(2015, report_year - 1)) %>%
  arrange(year, case_type)


#---------------- Fig 2.3.4 Panel plot of RR-TB incidence estimates in the 30 MDR high burden countries since 2015  -------------#
# Late renumbering from 2.3.6

if(show_estimates){

  f2_3_04_data <- est_dr_country %>%
    inner_join(hbmdr30, by = "iso3")  %>%
    select(country,
           year,
           e_inc_rr_num,
           e_inc_rr_num_lo,
           e_inc_rr_num_hi) %>%
    arrange(country)

}




#---------------- Fig 2.3.5 Bubble map of estimated incidence of MDR/RR-TB for countries with at least 1000 incident cases  -------------#
# Late renumbering from 2.3.3 to 2.3.5


f2_3_05_data <- est_dr_country %>%
  filter(year == report_year - 1 & e_inc_rr_num >= 1000) %>%
  select(iso3,
         size = e_inc_rr_num
  )

# Summary dataset for simple quoting of numbers in the text
f2_3_05_txt <- f2_3_05_data %>%
  arrange(desc(size)) %>%
  inner_join(country_names, by="iso3") %>%
  # pick the top 3
  head(3) %>%
  # Calculate proportion of 2021 global incidence
  mutate(pct = size * 100 / f2_3_01_data[f2_3_01_data$year==2021, "e_inc_rr_num"]) %>%
  mutate(country = ifelse(country == "Russian Federation", "the Russian Federation", country)) %>%
  select(-size)



#-------------- Fig 2.3.6 - Percentage of new TB cases with MDR/RR-TB ----------------#
# Percentages were based on the most recent data point for countries with representative data from 2007 to 2022.
# Changed in 2022 to instead use the estimates produced by Pete

# Late renumbering from 2.3.1 to 2.3.6

# Create dataframe for figure 2.3.6 ----

f2_3_06_data <- est_dr_country %>%
  filter(year == report_year - 1) %>%
  # COnvert proportion to %
  mutate(var = e_rr_prop_new * 100) %>%
  select(iso3, var) %>%

  # Assign the categories for the map
  mutate(var = cut(
    var,
    c(0, 3.0, 6.0, 12.0, 20, Inf),
    c('0-2.9', '3-5.9', '6-11.9', '12-19.9','\u226520'),
    right = FALSE
  ))



# Summary dataset for simple quoting of numbers in the text
# Text before 2.3.6 refers to regional averages

f2_3_06_txt <- est_dr_group %>%
  filter(year==2021) %>%
  arrange(desc(e_rr_prop_new)) %>%
  inner_join(who_region_shortnames, by = c("group_name" = "g_whoregion")) %>%
  select(g_whoregion = group_name, entity,  year, e_rr_prop_new)



#---------------- Fig 2.3.7 - Percentage of previously treated TB cases with MDR/RR-TB  -------------#
# Percentages are based on the most recent data point for countries with representative data from 2007 to 2022.
# Changed in 2022 to instead use the estimates produced by Pete
#
# Late renumbering from 2.3.2 to 2.3.7

# Create dataframe for figure 2.3.6 ----

f2_3_07_data <- est_dr_country %>%
  filter(year == report_year - 1) %>%
  # Convert proportion to %
  mutate(var = e_rr_prop_ret * 100) %>%
  select(iso3, var) %>%

  # Assign the categories for the map
  mutate(var = cut(
    var,
    c(0, 6.0, 12.0, 30.0, 50, Inf),
    c('0-5.9', '6-11.9', '12-29.9', '30-49.9','\u226550'),
    right = FALSE
  ))


# Summary dataset for simple quoting of numbers in the text
# Text before 2.3.7 refers to regional averages

f2_3_07_txt <- est_dr_group %>%
  filter(year==2021) %>%
  arrange(desc(e_rr_prop_ret)) %>%
  inner_join(who_region_shortnames, by = c("group_name" = "g_whoregion")) %>%
  select(g_whoregion = group_name, entity,  year, e_rr_prop_ret)




#------------------ Data for historical maps - Figs 2.3.8, 2.3.9, 2.3.10 ----------------#
# HT: Have not re-written this yet in a way that removes the need for the many intermediate
# dataframes (n2, n3, ...etc)
# Late renumbering from 2.3.7 to 2.3.8  etc.

Country_Template <- most_recent_for_estimation %>%
  select(country,iso3, year_new, source_new, all_areas_covered_new, surv_quality_new)



##we will report data points excluding sub-national data, and removing unreliable records as informed by the "view_dr_surveillance_exclusions" file
##we only use A, N data for new cases.
## the only sub-national data we want to keep is that from CAF, PNG, BRZ and PRK because these data are considered by WHO when producing global estimates.
##to count available data points for RR-TB, we will check that data is indeed available to calculate prevalence of RR-TB among new cases, by using the same code that is used for the
#prevalence maps. Note that if the RR-TB prevalence can only be calculated for new and previously treated cases combined (e.g. new plus prev treated [prevalence only provided for "unk" category]),
#then, this also DOES NOT COUNT as a data point of RR-TB among new cases.


n2 <- drs %>% select(country,iso3, g_whoregion, year, drs_record_type, surv_quality, all_areas_covered, dst_rlt_new,
                  mdr_new, mdr_new_pcnt, dr_r_nh_new, dr_r_nh_new_pcnt,
                  xpert_new, xpert_dr_r_new,xpert_dr_r_new_pcnt,
                  r_rlt_new, rr_new, rr_new_pcnt)

#HT: next line gives an error and I couldn't figure out from dplyr documentation how to make it work
# so using ifelse instead to do the same recoding
#n2$drs_record_type <- recode(n2$drs_record_type, "pre-2010 TME data collection" = "Surveillance")

n2 <- n2 %>%
  mutate(drs_record_type = ifelse(drs_record_type == "pre-2010 TME data collection", "Surveillance", drs_record_type))

#select data for the past 25 years. For the 2022 Global TB report, this would correspond to data from 1996-2022 (which is the same as >1995)

n3<- n2 %>%
  filter(year>=(report_year-26))%>%
  droplevels()


#check which records with surv_quality "A" or "N" actually had no data to estimate the proportion of RR-TB, and should consequently be excluded


n4<-n3 %>%
  mutate(RR_NEW_GR_2022 = case_when((surv_quality =="A" | surv_quality =="N") & year >=2017 & drs_record_type == "Surveillance" ~ round((rr_new / r_rlt_new)*100,1),
                                    (surv_quality =="A" | surv_quality =="N") & between(year, 2012, 2016) & drs_record_type == "Surveillance" ~ round(((coalesce(dr_r_nh_new,0) + coalesce(mdr_new,0) + coalesce(xpert_dr_r_new,0)) / (coalesce(dst_rlt_new,0) + coalesce(xpert_new,0)))*100,1),
                                    (surv_quality =="A" | surv_quality =="N") & year <=2011 & drs_record_type == "Surveillance" ~ round(((coalesce(dr_r_nh_new,0) + coalesce(mdr_new,0)) / dst_rlt_new)*100,1),
                                    (surv_quality =="A" | surv_quality =="N") & year >=2018 & drs_record_type == "Survey" ~ round(rr_new_pcnt,1),
                                    (surv_quality =="A" | surv_quality =="N") & year <=2012 & drs_record_type == "Survey" ~ round(coalesce(dr_r_nh_new_pcnt,0) + coalesce(mdr_new_pcnt,0),1),
                                    (surv_quality =="A" | surv_quality =="N") & between(year, 2013, 2017) & drs_record_type == "Survey" &  !is.na(xpert_dr_r_new_pcnt) & is.na(rr_new_pcnt) ~ round(xpert_dr_r_new_pcnt,1),
                                    (surv_quality =="A" | surv_quality =="N") & between(year, 2013, 2017) & drs_record_type == "Survey" &  is.na(xpert_dr_r_new_pcnt) & is.na(rr_new_pcnt)~ round(coalesce(dr_r_nh_new_pcnt,0) + coalesce(mdr_new_pcnt,0),1),
                                    (surv_quality =="A" | surv_quality =="N") & between(year, 2013, 2017) & drs_record_type == "Survey" &  !is.na(rr_new_pcnt) ~ round(rr_new_pcnt,1)))

# replace prevalence with NA where we have summed values that were coerced to 0 but were all NA in truth

n5<-n4 %>%
  mutate(RR_NEW_GR_2022 = replace(RR_NEW_GR_2022, !is.na(surv_quality) & year <=2012 & drs_record_type == "Survey" & is.na(dr_r_nh_new_pcnt) & is.na(mdr_new_pcnt), NA))

n6<-n5 %>%
  mutate(RR_NEW_GR_2022 = replace(RR_NEW_GR_2022, !is.na(surv_quality) & between(year, 2013, 2017) & drs_record_type == "Survey" &  is.na(xpert_dr_r_new_pcnt) & is.na(dr_r_nh_new_pcnt) & is.na(mdr_new_pcnt) & is.na(rr_new_pcnt), NA))



#exclude records where the percentage of RR-TB among new cases cannot be calculated. This is a combination of data where surv_quality is not "A" and is not "N", and cases where survey quality is "A" or "N" BUT
#there was actually no data or denominator data was zero:

n9<- n6 %>%
  filter(!is.na(RR_NEW_GR_2022))%>%
  droplevels()


#now we must exclude sub-national data, except for sub-national data from CAF, BRA, PRK, PNG

n10<- n9 %>%
  filter(all_areas_covered == 1 | (all_areas_covered == 0 & (iso3 == "CAF" | iso3 == "BRA" | iso3 == "PRK" | iso3 == "PNG")) )%>%
  droplevels()


#we must also exclude surveillance data from new cases, if this is in "view_dr_surveillance_exclusions"

ex2<- drs_exclusions %>%
  filter(rm_new == 1)%>%
  droplevels()%>%
  select(iso3,year)%>%
  mutate(drs_record_type = "Surveillance")%>%
  mutate(exclude = "exclude")

n11<-merge(n10,ex2,by=c("iso3","year","drs_record_type"),all.x = TRUE, all.y = FALSE)


n12<- n11 %>%
  filter(is.na(exclude))%>%
  droplevels()%>%
  select(-exclude)

#COUNT NUMBER OF DATAPOINTS PER COUNTRY
#A difference this year is that if for a particular year both survey and surveillance data are available, then this is counted as two data points.

n_count<-n12 %>% group_by(iso3) %>% tally() %>%
  rename(data_points=n)

#GET LATEST YEAR OF DATA AND LATEST SOURCE OF DATA
#these data include more years than that considered in view_drs_most_recent_for_estimation (2007-2022), and hence is taken from view_TME_master_drs (but is restricted to the past 25 years [1996-2022])
#now we select the most recent year of data along with the drs_record_type; if both survey and surveillance data are available for the last year, there will be 2 data points for the latest year of data. To
#avoid duplicate records, we spread drs_record_type data into two columns. In these cases, because surveillance data of quality A or N is prioritized over survey data (unless listed in the exclusions list), this will be
#reflected in the final selection of source in cases where data is too old to be in view_drs_most_recent_for_estimation (data taken from the latter is already de-duplicated).

n13<-n12 %>% select(country,iso3,year,drs_record_type, all_areas_covered) %>%
  spread(drs_record_type, drs_record_type) %>%
  group_by(country) %>% top_n(1, year)

#Now we merge n_count and n13 with our Country_Template data, which is based on view_drs_most_recent for estimation. This also serves to validate the data from n13

Data_Figs_7_8_9a<-merge(Country_Template,n_count,by=c("iso3"),all.x = TRUE, all.y = FALSE)

Data_Figs_7_8_9_b<-merge(Data_Figs_7_8_9a,n13,by=c("iso3"),all.x = TRUE, all.y = FALSE)


Data_Figs_7_8_9_c<-Data_Figs_7_8_9_b %>%
  mutate(most_recent_year = case_when(!is.na(year_new)  ~ year_new,
                                      TRUE ~ year)) %>%
  mutate(most_recent_source = case_when(!is.na(year_new)  ~ source_new,
                                        is.na(year_new)& !is.na(Surveillance) ~ Surveillance,
                                        TRUE ~ Survey))%>%
  rename(country=country.x) %>%
  select(country, iso3, all_areas_covered,data_points, most_recent_year, most_recent_source)

#for countries where there is no previous data, flag those that are planning or ongoing surveys this year. This part of the code
#must be reviewed every year

f2_3_7_and_8_and_9_data <- Data_Figs_7_8_9_c %>%
  mutate(most_recent_year = replace(most_recent_year, iso3 == "BDI", "Ongoing in 2022")) %>%
  mutate(most_recent_year = replace(most_recent_year,is.na(most_recent_year), "No data")) %>%
  mutate(most_recent_source = replace(most_recent_source,is.na(most_recent_source), "No data")) %>%
  mutate(all_areas_covered = case_when(all_areas_covered ==1  ~ "Nationwide",
                                       all_areas_covered ==0  ~ "Subnational",
                                       TRUE ~ "No data")) %>%
  mutate(most_recent_year_ranked = case_when(most_recent_year == "No data" | most_recent_year == "Ongoing in 2022"~ most_recent_year,
                                             between(most_recent_year, 1996, 2000) ~ "1996 - 2000",
                                             between(most_recent_year, 2001, 2005) ~ "2001 - 2005",
                                             between(most_recent_year, 2006, 2010) ~ "2006 - 2010",
                                             between(most_recent_year, 2011, 2015) ~ "2011 - 2015",
                                             between(most_recent_year, 2016, 2022) ~ "2016 - 2022"))

# Tidy up
rm(n2, n3, n4, n5, n6, n9, n10, n11, n12, n13, n_count)



# Summary dataset for simple quoting of numbers in the text
# Identify the countries with survey or surveillance measurements in the past 15 years
f2_3_8910_txt <- f2_3_7_and_8_and_9_data %>%
  filter(!is.na(data_points) & most_recent_year >= report_year - 15) %>%
  select(iso3, most_recent_source)

# Proportion of population accounted for by these countries
pop_drs <- pop %>%
  filter(year == 2021) %>%
  inner_join(f2_3_8910_txt, by = "iso3") %>%
  summarise_at("e_pop_num", sum, na.rm=TRUE)

pop_tot <-pop %>%
  filter(year == 2021) %>%
  summarise_at("e_pop_num", sum, na.rm=TRUE)

pop_drs_pct = pop_drs * 100 /pop_tot
pop_drs_pct <- pop_drs_pct %>% select(pct = e_pop_num)
rm(pop_drs, pop_tot, pop)

# Proportion of TB incidence for by these countries
inc_drs <- est_country %>%
  filter(year == 2021) %>%
  inner_join(f2_3_8910_txt, by = "iso3") %>%
  summarise_at("inc.num", sum, na.rm=TRUE)

inc_tot <- est_country %>%
  filter(year == 2021) %>%
  summarise_at("inc.num", sum, na.rm=TRUE)

inc_drs_pct = inc_drs * 100 / inc_tot
inc_drs_pct <- inc_drs_pct %>% select(pct = inc.num)
rm(inc_drs, inc_tot, est_country)

# Create the individual dataframes for each of the figures

f2_3_08_data <- f2_3_7_and_8_and_9_data %>%
  select(iso3, most_recent_source) %>%
  rename(var=most_recent_source)

f2_3_08_data$var <- factor(f2_3_08_data$var, levels = c("Surveillance", "Survey"))


f2_3_09_data <- f2_3_7_and_8_and_9_data %>%
  select(iso3, most_recent_year_ranked) %>%
  rename(var=most_recent_year_ranked)

f2_3_09_data$var <- factor(f2_3_09_data$var,
                          levels = c('1996 - 2000',
                                     '2001 - 2005',
                                     '2006 - 2010',
                                     '2011 - 2015',
                                     '2016 - 2022',
                                     'Ongoing in 2022'))


f2_3_10_data <- f2_3_7_and_8_and_9_data %>%
  select(iso3, data_points) %>%
  rename(var=data_points)


f2_3_10_data$var <- case_when(f2_3_10_data$var<= 1 ~ '1',
                             f2_3_10_data$var== 2 ~ '2',
                             between(f2_3_10_data$var, 3, 5) ~ '3-5',
                             between(f2_3_10_data$var, 6, 10) ~ '6-10',
                             between(f2_3_10_data$var, 11, 15) ~ '11-15',
                             f2_3_10_data$var >= 16 ~ "\u226516",
                             is.na(f2_3_10_data$var)~ 'No data')

f2_3_10_data$var <- factor(f2_3_10_data$var, levels = c('1',
                                                      '2',
                                                      '3-5',
                                                      '6-10',
                                                      '11-15',
                                                      '\u226516'))


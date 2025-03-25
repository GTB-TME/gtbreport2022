# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data preparation script for ch3-4.rmd
# Hazim Timimi, July 2022
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


# Load chapter 3 and 4 packages, settings and data
source(here::here('report/ch3_ch4_load_data.r'))


# Function to calculate outcome percentages for plotting as stacked bars ----
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

calculate_outcomes_pct <- function(df, prefix_string) {

  # Simplify outcome variable names so can calculate %
  # Remove the prefix) so data frame has columns called coh, succ, fail, died, lost and c_neval

  names(df) <- str_replace(names(df), prefix_string, "")

  df <- mutate(df,
               `Treatment success` = ifelse(NZ(coh) > 0,
                                            succ * 100 / coh,
                                            NA),
               Failure = ifelse(NZ(coh) > 0,
                                fail * 100 / coh,
                                NA),
               Died = ifelse(NZ(coh) > 0,
                             died * 100 / coh,
                             NA),
               `Lost to follow-up` = ifelse(NZ(coh) > 0,
                                            lost * 100 / coh,
                                            NA),
               `Not evaluated` = ifelse(NZ(coh) > 0,
                                        c_neval * 100 / coh,
                                        NA))

  return(df)

}




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.4.1 (previously 3.4.3) ----
# (Bar chart of numbers treated for MDR-TB each year since 2015)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

f3.4.1_data <- filter(notification, year >=2015) %>%

  select(year,
         unconf_rrmdr_tx,
         conf_rrmdr_tx,
         unconf_rr_nfqr_tx,
         conf_rr_nfqr_tx,
         conf_rr_fqr_tx,
         rrmdr_014_tx) %>%

  group_by(year) %>%
  summarise(across(unconf_rrmdr_tx:rrmdr_014_tx, sum, na.rm=TRUE)) %>%
  ungroup() %>%

  # Calculate total treated (all ages and over-15s). Set all ages to NA when we have numbers for children
  # and set adults to NA when we don;t have children. That way we can have a stacked bar chart that
  # handles the early years without children and the latter with children
  mutate(rrmdr_all_tx = ifelse(rrmdr_014_tx > 0,
                               NA,
                               unconf_rrmdr_tx + conf_rrmdr_tx + unconf_rr_nfqr_tx + conf_rr_nfqr_tx + conf_rr_fqr_tx)) %>%
  mutate(rrmdr_15plus_tx = ifelse(rrmdr_014_tx > 0,
                                  unconf_rrmdr_tx + conf_rrmdr_tx + unconf_rr_nfqr_tx + conf_rr_nfqr_tx + conf_rr_fqr_tx - rrmdr_014_tx,
                                  NA)) %>%
  mutate(rrmdr_014_tx = ifelse(rrmdr_014_tx == 0,
                               NA,
                               rrmdr_014_tx)) %>%

  # restrict to essential variables
  select(year,
         rrmdr_all_tx,
         rrmdr_014_tx,
         rrmdr_15plus_tx) %>%

  # pivot to long format
  pivot_longer(cols = starts_with("rrmdr_"),
               names_to = "age_group",
               names_pattern = "rrmdr_(.*)_tx",
               values_to = "how_many")

f3.4.1_txt <- f3.4.1_data %>%
  group_by(year) %>%
  summarise(across(how_many, sum, na.rm=TRUE)) %>%
  filter(year >=2019) %>%
  pivot_wider(names_from = year, values_from = how_many) %>%
  mutate(pct_increase_2021 = (`2021`/`2020`-1)*100) %>%
  rename(rr_2019=`2019`,
         rr_2020=`2020`,
         rr_2021=`2021`)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.4.2 ----
# (Panel of line charts of MDR/RR-TB cases detected, MDR/RR-TB put on treatment
# since 2009 for 30 countries)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

f3.4.2_data <- filter(notification, year >= 2010 & iso3 %in% hbmdr30$iso3) %>%
  select(country,
         iso3,
         year,
         rapid_dx_dr_r,
         conf_mdr,
         conf_rrmdr,
         conf_rr_nfqr,
         conf_rr_fqr,
         conf_mdr_tx,
         unconf_mdr_tx,
         conf_rrmdr_tx,
         unconf_rrmdr_tx,
         conf_rr_nfqr_tx,
         unconf_rr_nfqr_tx,
         conf_rr_fqr_tx) %>%
  mutate(conf_rr_nfqr=ifelse(conf_rr_nfqr==130124,2938,conf_rr_nfqr)) %>% # tentative correction for Nigeria
  rowwise() %>%
  # Derive total number detected and total enrolled on treatment
  mutate(rr_detected = ifelse(year < 2014,
                              sum(across(rapid_dx_dr_r:conf_mdr), na.rm = TRUE),
                              # the next three are mutually exclusive so can be added
                              sum(across(conf_rrmdr:conf_rr_fqr), na.rm = TRUE)),
         # treatment variables are in mutually exclusive sets so again can be added
         rr_treated = sum(across(conf_mdr_tx:conf_rr_fqr_tx), na.rm = TRUE)) %>%
  ungroup() %>%

  # drop unneeded variables
  select(country,
         iso3,
         year,
         rr_detected,
         rr_treated)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.4.3 ----
# (Irwin's doughnuts on completion of UNHLM targets)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Summary dataset for simple quoting of numbers
f3.4.3_txt <- filter(notification, year >= 2018) %>%
  summarise(across(c(conf_rrmdr_tx,
                     unconf_rrmdr_tx,
                     conf_rr_nfqr_tx,
                     unconf_rr_nfqr_tx,
                     conf_rr_fqr_tx,
                     rrmdr_014_tx), sum, na.rm=TRUE)) %>%

  # Derive total enrolled on MDR treatment
  rowwise() %>%
  mutate(rr_treated = sum(across(conf_rrmdr_tx:conf_rr_fqr_tx), na.rm = TRUE)) %>%
  ungroup() %>%
  select(-conf_rrmdr_tx,
         -unconf_rrmdr_tx,
         -conf_rr_nfqr_tx,
         -unconf_rr_nfqr_tx,
         -conf_rr_fqr_tx) %>%

  # Calculate percentage complete for each UNHLM 2018-2022 target
  mutate(rr_treated_pct  = rr_treated * 100 / 1.5e6,  # target: 1.5 million treated for drug-resistant TB
         rrmdr_014_tx_pct = rrmdr_014_tx * 100 / 115e3 # target: 115 thousand children treated for drug-resistant TB
  )

f3.4.3_data <- f3.4.3_txt %>%

  select(contains("_pct"))  %>%

  pivot_longer(cols = contains("_pct"),
               names_to = "target_completion")

# S
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.4.4 ----
# (Forest plot of TB treatment coverage in 30 countries, regionally and globally)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# notification
f3.4.4_data <- filter(notification, year >= 2015 ) %>%
  select(country,
         g_whoregion,
         iso3,
         year,
         rapid_dx_dr_r,
         conf_mdr,
         conf_rrmdr,
         conf_rr_nfqr,
         conf_rr_fqr,
         conf_mdr_tx,
         unconf_mdr_tx,
         conf_rrmdr_tx,
         unconf_rrmdr_tx,
         conf_rr_nfqr_tx,
         unconf_rr_nfqr_tx,
         conf_rr_fqr_tx) %>%
  mutate(conf_rr_nfqr=ifelse(conf_rr_nfqr==130124,2938,conf_rr_nfqr)) %>% # tentative correction for Nigeria

  rowwise() %>%
  # Derive total number detected and total enrolled on treatment
  mutate(rr_detected = ifelse(year < 2014,
                              sum(across(rapid_dx_dr_r:conf_mdr), na.rm = TRUE),
                              # the next three are mutually exclusive so can be added
                              sum(across(conf_rrmdr:conf_rr_fqr), na.rm = TRUE)),
         # treatment variables are in mutually exclusive sets so again can be added
         rr_treated = sum(across(conf_mdr_tx:conf_rr_fqr_tx), na.rm = TRUE)) %>%

  # calculate regional aggregates
  group_by(g_whoregion,year) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%

  # merge with regional names
  inner_join(who_region_shortnames, by = "g_whoregion") %>%
  arrange(entity) %>%
  select(-g_whoregion) %>%

  ungroup() %>%

  # drop unneeded variables
  select(entity,
         year,
         rr_detected,
         rr_treated)

# Add global summary to the regional summary
f3.4.4_global <- f3.4.4_data %>%
  group_by(year) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  mutate(entity="Global") %>%
  ungroup()

f3.4.4_data <- rbind(f3.4.4_data, f3.4.4_global)

f3.4.4_drest <- est_dr_group %>%
  mutate(g_whoregion=group_name) %>%
  # merge with regional names
  inner_join(who_region_shortnames, by = "g_whoregion") %>%
  arrange(entity) %>%
  select(-g_whoregion) %>%
  select(entity,
         year,
         e_inc_rr_num,
         e_inc_rr_num_lo,
         e_inc_rr_num_hi)

# Add global summary to the regional summary
f3.4.4_drest_global <- est_dr_group %>%
  mutate(entity=group_name) %>%
  subset(entity=="global") %>%
  mutate(entity="Global") %>%
  select(entity,
         year,
         e_inc_rr_num,
         e_inc_rr_num_lo,
         e_inc_rr_num_hi)

f3.4.4_drest <- rbind(f3.4.4_drest, f3.4.4_drest_global)

# merge notification and estimates
f3.4.4_data <- left_join(f3.4.4_data,f3.4.4_drest, by=c("entity","year"))
f3.4.4_data <- f3.4.4_data %>%
  mutate(entity = factor(entity,
                         levels = c("Global","African Region", "Region of the Americas", "South-East Asia Region",
                                    "European Region", "Eastern Mediterranean Region", "Western Pacific Region")))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.4.5 ----
# (Forest plot of TB treatment coverage for MDR/RR-TB in 30 countries, regionally and globally)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if(show_estimates){

  # A. Countries
  # - - - - - - - -
  coverage_inc_country <- filter(est_dr_country, year==report_year-1) %>%
    select(year,
           iso3,
           e_inc_rr_num,
           e_inc_rr_num_lo,
           e_inc_rr_num_hi) %>%

    # restrict to high burden countries
    inner_join(hbmdr30, by = "iso3")

  coverage_country <- filter(notification, year==report_year-1) %>%
    select(entity = country,
           iso3,
           unconf_rr_nfqr_tx,
           conf_rr_nfqr_tx,
           conf_rr_fqr_tx)  %>%
    mutate(rr_tx = ifelse(is.na(unconf_rr_nfqr_tx) & is.na(conf_rr_nfqr_tx) & is.na(conf_rr_fqr_tx),
                          NA,
                          NZ(unconf_rr_nfqr_tx) + NZ(conf_rr_nfqr_tx) + NZ(conf_rr_fqr_tx))) %>%
    inner_join(coverage_inc_country, by = "iso3") %>%
    select(-iso3) %>%
    mutate(c_rr_coverage = rr_tx * 100 / e_inc_rr_num,
           c_rr_coverage_lo = rr_tx * 100  / e_inc_rr_num_hi,
           c_rr_coverage_hi = rr_tx * 100  / e_inc_rr_num_lo,
           # highlight countries with no data
           entity = ifelse(is.na(rr_tx), paste0(entity, "*"), entity )) %>%
    select(entity,
           c_rr_coverage,
           c_rr_coverage_lo,
           c_rr_coverage_hi) %>%
    arrange(desc(c_rr_coverage))

  # B. Regions
  # - - - - - - - -
  coverage_inc_region <- filter(est_dr_group, year==report_year-1 & group_type=="g_whoregion") %>%
    select(g_whoregion = group_name,
           e_inc_rr_num,
           e_inc_rr_num_lo,
           e_inc_rr_num_hi)

  coverage_region <- filter(notification, year==report_year-1) %>%
    select(g_whoregion,
           unconf_rr_nfqr_tx,
           conf_rr_nfqr_tx,
           conf_rr_fqr_tx)  %>%
    # calculate regional aggregates
    group_by(g_whoregion) %>%
    summarise(across(unconf_rr_nfqr_tx:conf_rr_fqr_tx, sum, na.rm=TRUE)) %>%
    ungroup() %>%
    mutate(rr_tx = unconf_rr_nfqr_tx + conf_rr_nfqr_tx + conf_rr_fqr_tx) %>%

    # merge with incidence estimates
    inner_join(coverage_inc_region, by = "g_whoregion") %>%

    # Calculate coverage
    mutate(c_rr_coverage = rr_tx * 100 / e_inc_rr_num,
           c_rr_coverage_lo = rr_tx * 100  / e_inc_rr_num_hi,
           c_rr_coverage_hi = rr_tx * 100  / e_inc_rr_num_lo) %>%

    # merge with regional names and simplify to match structure of country table
    inner_join(who_region_shortnames, by = "g_whoregion") %>%
    select(entity,
           c_rr_coverage,
           c_rr_coverage_lo,
           c_rr_coverage_hi) %>%
    arrange(desc(c_rr_coverage))

  # C. Global
  # - - - - - - - -
  coverage_inc_global <- filter(est_dr_group, year==report_year-1 & group_type=="global") %>%
    select(e_inc_rr_num,
           e_inc_rr_num_lo,
           e_inc_rr_num_hi) %>%
    mutate(entity="Global")

  coverage_global <- filter(notification, year==report_year-1) %>%
    select(unconf_rr_nfqr_tx,
           conf_rr_nfqr_tx,
           conf_rr_fqr_tx)  %>%
    # calculate global aggregate
    summarise(across(unconf_rr_nfqr_tx:conf_rr_fqr_tx, sum, na.rm=TRUE)) %>%
    mutate(rr_tx = unconf_rr_nfqr_tx + conf_rr_nfqr_tx + conf_rr_fqr_tx) %>%
    mutate(entity="Global") %>%

    inner_join(coverage_inc_global, by="entity") %>%

    # Calculate coverage
    mutate(c_rr_coverage = rr_tx * 100 / e_inc_rr_num,
           c_rr_coverage_lo = rr_tx * 100  / e_inc_rr_num_hi,
           c_rr_coverage_hi = rr_tx * 100  / e_inc_rr_num_lo) %>%
    select(entity,
           c_rr_coverage,
           c_rr_coverage_lo,
           c_rr_coverage_hi)

  # D. Bring them all together
  # - - - - - - - - - - - - -

  # Create dummy records so can see a horizontal line in the output to separate countries, regions and global parts
  coverage_dummy1 <- data.frame(entity = " ", c_rr_coverage = NA, c_rr_coverage_lo = 0, c_rr_coverage_hi = 100)
  coverage_dummy2 <- data.frame(entity = "  ", c_rr_coverage = NA, c_rr_coverage_lo = 0, c_rr_coverage_hi = 100)


  # Create combined dataframe in order of countries then regional and global estimates
  f3.4.5_data <- rbind(coverage_country,
                       coverage_dummy1,
                       coverage_region,
                       coverage_dummy2,
                       coverage_global) %>%

    # The dataframe is in the order I want, so make entity an ordered factor based on
    # what I already have. That way ggplot will not reorder by entity name
    # But I need to reverse order for plotting

    mutate(entity = factor(entity,
                           levels = rev(entity)))

  # Simple dataframe of numbers for section text
  f3.4.5_txt <- coverage_global %>%
    select(c_rr_coverage)

  f3.4.5_txt_list_region <- coverage_region %>%
    slice(1:2)

  f3.4.5_txt_list_hi <- f3.4.5_data %>%
    slice(1:30) %>%
    filter(c_rr_coverage > 40) %>%
    select(entity)

  f3.4.5_txt_list_lo <- f3.4.5_data %>%
    slice(1:30) %>%
    filter(c_rr_coverage < 12) %>%
    select(entity)

  # remove the temporary dataframes
  rm(list=ls(pattern = "^coverage"))

}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.4.06  ----
# (Bubble map of difference between MDR/RR-TB patients starting treatment and
#  estimated MDR/RR-TB incidence for 10 countries)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if(show_estimates){

  f3.4.6_data  <- filter(est_dr_country, year == report_year - 1) %>%
    select(iso3,
           year,
           e_inc_rr_num) %>%

    # Link to notifications
    inner_join(notification, by = c("iso3", "year")) %>%
    select(iso3,
           country,
           e_inc_rr_num,
           conf_rr_nfqr_tx,
           unconf_rr_nfqr_tx,
           conf_rr_fqr_tx) %>%

    # Filter out countries that have not reported anything yet for the latest year
    filter(!is.na(conf_rr_nfqr_tx) | !is.na(unconf_rr_nfqr_tx) | !is.na(conf_rr_fqr_tx)) %>%

    # Calculate the gap and use that for the bubble sizes
    mutate(size = e_inc_rr_num - (NZ(conf_rr_nfqr_tx) + NZ(unconf_rr_nfqr_tx) + NZ(conf_rr_fqr_tx))) %>%

    # limit to the top 10 by size of gap
    top_n(10, size) %>%

    # sort in descending order so can list the country names in the figure footnote
    arrange(desc(size))

  # Summary number of gaps for the section text
  # Get global incidence
  f3.4.6_txt <- filter(est_dr_group, year == report_year-1 & group_type == "global") %>%
    select(e_inc_rr_num)

  # Add global enrolments
  f3.4.6_txt <- filter(notification, year == report_year-1) %>%
    select(year,
           conf_rr_nfqr_tx,
           unconf_rr_nfqr_tx,
           conf_rr_fqr_tx) %>%
    # calculate global aggregate
    group_by(year) %>%
    summarise(across(conf_rr_nfqr_tx:conf_rr_fqr_tx, sum, na.rm=TRUE)) %>%
    ungroup() %>%
    mutate(enrolled = conf_rr_nfqr_tx + unconf_rr_nfqr_tx + conf_rr_fqr_tx) %>%
    cbind(f3.4.6_txt) %>%

    # Calculate the global gap and drop the other variables
    mutate(gap = e_inc_rr_num - enrolled) %>%
    select(gap)

  # Calculate % of global gap contributed by the top 10 countries
  f3.4.6_txt <- f3.4.6_data %>%
    summarise(gap_top_ten = sum(size)) %>%
    cbind(f3.4.6_txt) %>%
    mutate(pct_gap_top_ten = gap_top_ten * 100 / gap)

}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.4.7  ----
# (horizontal bar chart showing TB treatment outcomes globally by year since 2012 for MDR/RR-TB)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

f3.4.7_data <- outcomes %>%

  filter(between(year, 2012, report_year - 3)) %>%
  select(iso2,
         year,
         contains("mdr_")) %>%

  # calculate global aggregates
  group_by(year) %>%
  summarise(across(mdr_coh:c_mdr_neval, sum, na.rm=TRUE)) %>%
  ungroup() %>%

  # Calculate outcome proportions for plotting as stacked bars
  calculate_outcomes_pct("mdr_") %>%

  # Drop the actual numbers and keep percentages, plu sother unwanted variables
  select(-coh, -succ, -fail, -died, -lost, -c_neval, -cur, -cmplt) %>%

  # Add tx group type
  mutate(subgroup = "MDR/RR-TB cases") %>%

  # flip into long format
  pivot_longer(cols = `Treatment success`:`Not evaluated`,
               names_to = "outcome")

# Summary for section text
f3.4.7_txt <- filter(f3.4.7_data, year %in% c(2012, report_year-4,  report_year-3) & outcome=="Treatment success") %>%
  select(-outcome) %>%
  pivot_wider(names_from = year,
              names_prefix = "c_tsr_",
              values_from = value)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.4.8 ----
# (Horizontal bar chart showing TB treatment outcomes in MDR/RR-TB cases for WHO regions and globally")
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Code similar to others so maybe make a reusable function?

# A. Regional aggregates
# - - - - - - - - - - - -
txout_region  <- filter(outcomes, year==report_year - 3) %>%

  select(iso2,
         g_whoregion,
         contains("mdr_")) %>%
  select(-mdr_cur, -mdr_cmplt) %>%

  group_by(g_whoregion) %>%
  summarise(across(mdr_coh:c_mdr_neval, sum, na.rm=TRUE)) %>%
  ungroup() %>%

  # merge with regional names
  inner_join(who_region_shortnames, by = "g_whoregion") %>%
  select(-g_whoregion) %>%

  # Calculate outcome proportions for plotting as stacked bars
  calculate_outcomes_pct("mdr_") %>%

  # Sort regions in descending order of success rate
  arrange(desc(`Treatment success`))


# B. Global aggregates
# - - - - - - - - - - - -
txout_global  <- filter(outcomes, year==report_year - 3) %>%

  select(iso2,
         contains("mdr_")) %>%
  select(-mdr_cur, -mdr_cmplt) %>%

  summarise(across(mdr_coh:c_mdr_neval, sum, na.rm=TRUE)) %>%
  ungroup()  %>%
  mutate(entity="Global")  %>%

  # Calculate outcome proportions for plotting as stacked bars
  calculate_outcomes_pct("mdr_")


# Create a dummy record a gap in the output to separate countries, regions and global parts
txout_dummy <- data.frame(entity = " ", coh = NA, succ = NA, fail = NA,
                          died = NA, lost = NA, c_neval = NA,
                          Failure = NA, Died = NA)

# Had to use mutate to create the next 3 fields because data.frame converted spaces to dots. Grrr
txout_dummy <- txout_dummy %>%
  mutate(`Treatment success` = NA,
         `Lost to follow-up` = NA,
         `Not evaluated` = NA)


# Create combined table in order of countries then regional and global estimates
f3.4.8_data <- rbind(txout_region, txout_dummy, txout_global) %>%

  # Keep record of current order (in reverse) so plot comes out as we want it
  mutate(entity = factor(entity, levels=rev(entity))) %>%

  # Drop the actual numbers and keep percentages
  select(-coh, -succ, -fail, -died, -lost, -c_neval) %>%

  # Flip into long mode for stacked bar plotting
  pivot_longer(cols = `Treatment success`:`Not evaluated`,
               names_to = "outcome")

# Summary for section text
f3.4.8_txt <- filter(txout_region, entity %in% c("African Region", "European Region")) %>%
  select(entity,
         c_tsr = `Treatment success`) %>%
  pivot_wider(names_from = entity,
              names_prefix = "c_tsr_",
              values_from = c_tsr) %>%
  # deal with spaces in variable names
  select(c_tsr_AFR = `c_tsr_African Region`,
         c_tsr_EUR = `c_tsr_European Region`)

# tidy up
rm(list = ls(pattern = "^txout_"))



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.4.9 ----
# (World map showing which countries used bedaquiline in drug-resistant TB treatment regimens)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

f3.4.9_data <- filter(notification, year == report_year-1) %>%
  select(country,
         iso3,
         mdrxdr_bdq_used) %>%

  # Change the "No data" option 3 to NA to avoid weird effects in the map legend
  mutate(mdrxdr_bdq_used = ifelse(mdrxdr_bdq_used==3, NA, mdrxdr_bdq_used)) %>%

  # Assign the categories for the map
  mutate(var = factor(mdrxdr_bdq_used,
                      levels = c(1, 0),
                      labels = c("Used", "Not used")))

f3.4.9_txt <- filter(notification, year %in% c(2013, report_year-2, report_year-1)) %>%

    select(country,
         iso3,
         year,
         mdrxdr_bdq_used) %>%

  # Change the "No data" option 3 to NA to avoid weird effects in the map legend
  mutate(mdrxdr_bdq_used = ifelse(mdrxdr_bdq_used==3, NA, mdrxdr_bdq_used)) %>%

  group_by(year) %>%
  summarise(across(mdrxdr_bdq_used, sum, na.rm=TRUE)) %>%
  ungroup() %>%

  pivot_wider(names_from = year,
              names_prefix = "bdq_used_",
              values_from = mdrxdr_bdq_used)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.4.10  ----
# (World map showing which countries used all-oral longer MDR-TB treatment regimens)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

f3.4.10_data <- filter(notification, year == report_year-1) %>%
  select(country,
         iso3,
         mdrxdr_alloral_used) %>%

  # Change the "No data" option 3 to NA to avoid weird effects in the map legend
  mutate(mdrxdr_alloral_used = ifelse(mdrxdr_alloral_used==3, NA, mdrxdr_alloral_used)) %>%

  # Assign the categories for the map
  mutate(var = factor(mdrxdr_alloral_used,
                      levels = c(1, 0),
                      labels = c("Used", "Not used")))



f3.4.10_txt <- filter(notification, year %in% c(2019, report_year-2, report_year-1)) %>%

  select(country,
         iso3,
         year,
         mdrxdr_alloral_used) %>%

  # Change the "No data" option 3 to NA to avoid weird effects in the map legend
  mutate(mdrxdr_alloral_used = ifelse(mdrxdr_alloral_used==3, NA, mdrxdr_alloral_used)) %>%

  group_by(year) %>%
  summarise(across(mdrxdr_alloral_used, sum, na.rm=TRUE)) %>%
  ungroup() %>%

  pivot_wider(names_from = year,
              names_prefix = "oral_used_",
              values_from = mdrxdr_alloral_used)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.4.11  ----
# (World map showing which countries used shorter MDR-TB treatment regimens)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

f3.4.11_data <- filter(notification, year == report_year-1) %>%
  select(country,
         iso3,
         mdr_alloral_short_used) %>%

  # Change the "No data" option 3 to NA to avoid weird effects in the map legend
  mutate(mdr_alloral_short_used = ifelse(mdr_alloral_short_used==3, NA, mdr_alloral_short_used)) %>%


  # Assign the categories for the map
  mutate(var = factor(mdr_alloral_short_used,
                      levels = c(1, 0),
                      labels = c("Used", "Not used")))


f3.4.11_txt <-  filter(notification, year %in% c(2020, report_year-2, report_year-1)) %>%

  select(country,
         iso3,
         year,
         mdr_alloral_short_used) %>%

  # Change the "No data" option 3 to NA to avoid weird effects in the map legend
  mutate(mdr_alloral_short_used = ifelse(mdr_alloral_short_used==3, NA, mdr_alloral_short_used)) %>%

  group_by(year) %>%
  summarise(across(mdr_alloral_short_used, sum, na.rm=TRUE)) %>%
  ungroup() %>%

  pivot_wider(names_from = year,
              names_prefix = "oral_short_used_",
              values_from = mdr_alloral_short_used)


# Summary for section text from the last 3 figures' data
f3.4.9_11_txt <- filter(f3.4.9_data, mdrxdr_bdq_used == 1) %>%
  summarise(bdq = n())

f3.4.9_11_txt <- filter(f3.4.10_data, mdrxdr_alloral_used == 1) %>%
  summarise(alloral_long = n()) %>%
  cbind(f3.4.9_11_txt)

f3.4.9_11_txt <- filter(f3.4.11_data, mdr_alloral_short_used == 1) %>%
  summarise(alloral_short = n()) %>%
  cbind(f3.4.9_11_txt)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.4.12  ----
# (World map showing proportion of drug-resistant TB patients followed up for adverse events)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

f3.4.12_data <- filter(notification, year==report_year-1) %>%

  select(country,
         iso3,
         unconf_rr_nfqr_tx,
         conf_rr_nfqr_tx,
         conf_rr_fqr_tx,
         mdr_tx_adsm) %>%

  # Calculate percent with active follow up
  mutate(adsm_pct = ifelse(NZ(unconf_rr_nfqr_tx) + NZ(conf_rr_nfqr_tx) + NZ(conf_rr_fqr_tx) > 0,
                           mdr_tx_adsm * 100 / (NZ(conf_rr_nfqr_tx) + NZ(conf_rr_nfqr_tx) + NZ(conf_rr_fqr_tx)),
                           NA)) %>%

  # Assign the categories for the map
  mutate(var = cut(adsm_pct,
                   c(0, 25, 50, 75, Inf),
                   c('0-24', '25-49', '50-74', '\u226575'),
                   right=FALSE)) %>%

  # get rid of extra variables
  select(country,
         iso3,
         adsm_pct,
         var)



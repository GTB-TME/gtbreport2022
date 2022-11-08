# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data preparation script for ch3-2.rmd
# Hazim Timimi, July 2022
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


# Load chapter 3 and 4 packages, settings and data
source(here::here('report/ch3_ch4_load_data.r'))

# WB income group list
wb_incomelist <- country_group_membership %>%
  filter(group_type == "g_income") %>% select(iso3,group_name) %>% rename(income=2)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: bac confirmation ----
# Determine numerators and denominators for bacteriological confirmation and
# then use this dataframe for figures 3.2.1, 3.2.2 and 3.3.3
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


bacconf_data <- notification %>%
  filter(year >= 2000) %>%
  select(iso3,
         country,
         year,
         g_whoregion,
         # old variables pre-2013
         new_sp,
         new_sn,
         new_su,
         # new variables
         new_labconf, new_clindx,
         ret_rel_labconf, ret_rel_clindx) %>%

  #calculate % of pulmonary cases with bac confirmation
  rowwise() %>%
  # a bit tricky for years before 2013, so do for new only by smear only
  mutate(bacconf_pct_numerator = ifelse(year < 2013 & g_whoregion != 'EUR',
                                        # old variables, do for new only outside EUR
                                        new_sp,
                                        # new variables
                                        sum(c_across(contains("labconf")), na.rm = TRUE)),
         bacconf_pct_denominator = ifelse(year < 2013 & g_whoregion != 'EUR',
                                          # old variables, do for new only outside EUR
                                          sum(c_across(new_sp:new_su), na.rm = TRUE),
                                          # new variables
                                          sum(c_across(new_labconf:ret_rel_clindx), na.rm = TRUE))) %>%

  # Adjust calculation for EUR pre-2013 (applies to years 2002 - 2012)
  mutate(bacconf_pct_numerator = ifelse(between(year, 2002, 2012) & g_whoregion == 'EUR',
                                        # old variables, but using new_labconf
                                        new_labconf,
                                        # otherwise keep calculation from previous step
                                        bacconf_pct_numerator),
         bacconf_pct_denominator = ifelse(between(year, 2002, 2012) & g_whoregion == 'EUR',
                                          # old variables
                                          sum(c_across(new_sp:new_su), na.rm = TRUE),
                                          # otherwise keep calculation from previous step
                                          bacconf_pct_denominator)) %>%

  # Finally deal with EUR 2000 and 2001 numerator
  mutate(bacconf_pct_numerator = ifelse(between(year, 2000, 2001) & g_whoregion == 'EUR',
                                        # old variables
                                        new_sp,
                                        # otherwise keep calculation from previous step
                                        bacconf_pct_numerator),
         bacconf_pct_denominator = ifelse(between(year, 2000, 2001) & g_whoregion == 'EUR',
                                          # old variables
                                          sum(c_across(new_sp:new_su), na.rm = TRUE),
                                          # otherwise keep calculation from previous step
                                          bacconf_pct_denominator)) %>%

  ungroup() %>%

  # reduce to needed variables
  select(country,
         iso3,
         year,
         g_whoregion,
         bacconf_pct_numerator,
         bacconf_pct_denominator)




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.2.1 ----
# (Panel plot of TB cases with bacteriological confirmation by WHO region and globally since 2000)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# by income group for main texts
bacconf_data_income <- bacconf_data %>%
  mutate(bacconf_pct = bacconf_pct_numerator * 100 / bacconf_pct_denominator) %>%
  inner_join(wb_incomelist, by = "iso3") %>%
  filter(year==2021) %>%
  group_by(income) %>%
  summarise(across(bacconf_pct, median, na.rm = TRUE)) 
  

# Calculate aggregates
bacconf_data_regional <- bacconf_data %>%
  group_by(year, g_whoregion) %>%
  summarise(across(bacconf_pct_numerator:bacconf_pct_denominator, sum, na.rm = TRUE)) %>%

  # merge with regional names
  inner_join(who_region_shortnames, by = "g_whoregion") %>%
  select(-g_whoregion) %>%
  ungroup()

bacconf_data_global <- bacconf_data %>%
  group_by(year) %>%
  summarise(across(bacconf_pct_numerator:bacconf_pct_denominator, sum, na.rm = TRUE)) %>%
  mutate(entity = 'Global')

# Add global to the regional aggregates
f3.2.1_data <- rbind(bacconf_data_regional, bacconf_data_global) %>%

  # Calculate the percentages
  mutate(bacconf_pct = bacconf_pct_numerator * 100 / bacconf_pct_denominator) %>%

  # Change the order of the entities
  mutate(entity = factor(entity,
                         levels = c("Global", "African Region", "Region of the Americas", "South-East Asia Region",
                                    "European Region", "Eastern Mediterranean Region", "Western Pacific Region")))



# summary dataset for quoting numbers in the text on bac confirmation in pulmonary TB
f3.2.1_txt <- filter(notification, year == (report_year - 1)) %>%
  select(c_newinc,
         new_labconf, new_clindx,
         ret_rel_labconf, ret_rel_clindx) %>%
  summarise(across(c_newinc:ret_rel_clindx, sum, na.rm=TRUE)) %>%
  mutate(pulm = new_labconf + new_clindx + ret_rel_labconf + ret_rel_clindx) %>%
  mutate(pulm_pct = pulm * 100/ c_newinc) %>%
  select(c_newinc, pulm, pulm_pct)

# Calculate global % bac conf for the last two years and percent change and add to the summary
f3.2.1_txt <- filter(f3.2.1_data, entity == "Global" & year >= report_year-2) %>%
  select(year, bacconf_pct) %>%
  pivot_wider(names_from = year,
              names_prefix = "bc_pct_",
              values_from = bacconf_pct) %>%
  cbind(f3.2.1_txt)

# Add regional max and min values for 2020
f3.2.1_txt <- filter(f3.2.1_data, year==report_year-1 & entity %in% c("Region of the Americas", "Western Pacific Region")) %>%
  select(entity, bacconf_pct) %>%
  pivot_wider(names_from = entity,
              names_prefix = "bc_pct_",
              values_from = bacconf_pct) %>%
  # handle spaces in column names
  select(bc_pct_AMR = `bc_pct_Region of the Americas`,
         bc_pct_WPR = `bc_pct_Western Pacific Region`) %>%
  cbind(f3.2.1_txt)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.2.2 ----
# (World map showing percent of TB cases with bacteriological confirmation)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

f3.2.2_data <- filter(bacconf_data, year>=report_year-2) %>%
  
  # Calculate the percentages
  mutate(bacconf_pct = ifelse(bacconf_pct_denominator > 0,
                              bacconf_pct_numerator * 100 / bacconf_pct_denominator,
                              NA)) %>%
  
  # Assign the categories for the map
  mutate(var = cut(bacconf_pct,
                   c(0, 50, 65, 80, Inf),
                   c('0-49', '50-64', '65-79', '\u226580'),
                   right=FALSE))

# Find the countries with empty data for latest year and see if there are data for the previous year
bacconf_prev_year_data <- f3.2.2_data %>%
  filter(year == report_year - 1 & is.na(bacconf_pct)) %>%
  select(iso3) %>%
  inner_join(filter(f3.2.2_data, year == report_year - 2), by = "iso3") %>%
  filter(!is.na(bacconf_pct))

# Now combine into one data frame, with previous data used if latest year's data are not available
f3.2.2_data <- f3.2.2_data %>%
  filter(year == report_year - 1) %>%
  anti_join(bacconf_prev_year_data, by= "iso3") %>%
  rbind(bacconf_prev_year_data)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.2.3 ----
# (Box plot showing percent of TB cases tested with rapid diagnostics)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Calculate aggregates
f3.2.3_data <- bacconf_data %>%
  mutate(bacconf_pct=bacconf_pct_numerator/bacconf_pct_denominator) %>%
  left_join(wb_incomelist, by = "iso3") %>%  #NA: (no TB) Anguilla, Cook Islands, Montserrat, Niue, Tokelau, Wallis and Futuna, (with TB cases) Venezuela 
  mutate(income=factor(income,labels=c("High-income","Low-income","Lower-middle-income","Upper-middle-income"))) %>%
  mutate(income=factor(income,levels=c("Low-income","Lower-middle-income","Upper-middle-income","High-income")))


# Just uses the conf dataframe
# create summary for the text
# f3.2.3_txt <- filter(conf, year==report_year-1 & g.income=="HIC") %>%
#   mutate(median_HIC = conf.med * 100) %>%
#   select(median_HIC)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.2.4 ----
# (Panel plot of TB cases with bacteriological confirmation for 30 countries since 2000)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

f3.2.4_data <- bacconf_data %>%
  inner_join(hbc30, by = "iso3") %>%

  # Calculate the percentages
  mutate(bacconf_pct = ifelse(bacconf_pct_denominator > 0,
                              bacconf_pct_numerator * 100 / bacconf_pct_denominator,
                              NA)) %>%

  # get rid of extra variables
  select(country,
         year,
         bacconf_pct)

# summary datasets for quoting numbers in the text
f3.2.4_txt_MOZ <- filter(f3.2.4_data, year==report_year-1 & country=="Mozambique")
f3.2.4_txt_list_hi <- filter(f3.2.4_data, year==report_year-1 & bacconf_pct > 73) %>%
  select(country)
f3.2.4_txt_list_lo <- filter(f3.2.4_data, year==report_year-1 & bacconf_pct < 45) %>%
  select(country)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.2.5 ----
# (World map showing percent of TB cases tested with rapid diagnostics)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

f3.2.5_data <- notification %>%
  filter(year  >= report_year - 2) %>%
  select(iso3,
         country,
         year,
         c_newinc,
         rdx_data_available,
         newinc_rdx,
         rdxsurvey_newinc,
         rdxsurvey_newinc_rdx) %>%

  # calculate the percentage for each country depending on data availability
  mutate(wrd_pcnt_nu = ifelse(rdx_data_available == 60 & NZ(c_newinc) > 0,
                              newinc_rdx,
                              ifelse(rdx_data_available == 61 & NZ(rdxsurvey_newinc) > 0,
                                     rdxsurvey_newinc_rdx,
                                     NA)),
         wrd_pcnt_de = ifelse(rdx_data_available == 60 & NZ(c_newinc) > 0,
                              c_newinc,
                              ifelse(rdx_data_available == 61 & NZ(rdxsurvey_newinc) > 0,
                                     rdxsurvey_newinc,
                                     NA)),
         wrd_pcnt =ifelse(rdx_data_available == 60 & NZ(c_newinc) > 0,
                          newinc_rdx * 100 / c_newinc,
                          ifelse(rdx_data_available == 61 & NZ(rdxsurvey_newinc) > 0,
                                 rdxsurvey_newinc_rdx * 100 / rdxsurvey_newinc,
                                 NA))) %>%

  # Assign the categories for the map
  mutate(var = cut(wrd_pcnt,
                   c(0, 25, 50, 75, 90, Inf),
                   c('<25', '25-49', '50-75', '76-90','\u226590'),
                   right=FALSE)) %>%

  # get rid of extra variables
  select(country,
         iso3,
         year,
         wrd_pcnt_nu,
         wrd_pcnt_de,
         wrd_pcnt,
         var)


# Find the countries with empty data for latest year and see if there are data for the previous year
wrd_prev_year_data <- f3.2.5_data %>%
  filter(year == report_year - 1 & is.na(wrd_pcnt)) %>%
  select(iso3) %>%
  inner_join(filter(f3.2.5_data, year == report_year - 2), by = "iso3") %>%
  filter(!is.na(wrd_pcnt))

# Now combine into one dataframe, with previous data used if latest year's data are not available
  
f3.2.5_data <- f3.2.5_data %>%
  filter(year == report_year - 1) %>%
  anti_join(wrd_prev_year_data, by= "iso3") %>%
  rbind(wrd_prev_year_data) %>%
  left_join(hbc30, by= "iso3") %>%
  left_join(hbtbhiv30, by= "iso3") %>%
  left_join(hbmdr30, by= "iso3") %>%
  rename(hbc30=group_type.x,hbtbhiv30=group_type.y,hbmdr30=group_type) %T>%
  write_csv(here("./report/ch3_data/f3.2.5_data.csv"))


# summary dataset for quoting numbers in the text
f3.2.5_txt <- filter(notification, year  >= report_year - 3) %>%
  select(year,
         c_newinc,
         newinc_rdx) %>%
  group_by(year) %>%
  summarise(across(c_newinc:newinc_rdx, sum, na.rm=TRUE)) %>%
  ungroup() %>%
  mutate(wrd_pct = newinc_rdx * 100/c_newinc) %>%
  select(-c_newinc) %>%
  pivot_wider(names_from = year,
              values_from = c(newinc_rdx, wrd_pct))

# Add info about testing at least half of new cases in the high burden countries
f3.2.5_txt <- filter(notification,
                     year  >= report_year - 3 &
                       newinc_rdx >= 0.5 * c_newinc &
                       (iso3 %in% hbc30$iso3 | iso3 %in% hbtbhiv30$iso3 | iso3 %in% hbmdr30$iso3)) %>%
  select(year, country) %>%
  group_by(year) %>%
  summarise(hbcs = n()) %>%
  ungroup() %>%
  pivot_wider(names_from = year,
              names_prefix = "hbcs_",
              values_from = hbcs) %>%
  cbind(f3.2.5_txt) 

f3.2.5_txt_list <- f3.2.5_data %>%
  filter(!is.na(hbc30), wrd_pcnt>90) 

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.2.6 ----
# (World map showing proportion of diagnostic sites with WRDs)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

f3.2.6_data <- strategy %>%
  filter(year == report_year - 1) %>%
  select(iso3,
         country,
         year,
         m_wrd,
         dx_test_sites) %>%
  # Calculate the proportion
  mutate(wrd_pct = ifelse(dx_test_sites > 0,
                          m_wrd * 100 / dx_test_sites,
                              NA)) %>%
  # Assign the categories for the map
  mutate(var = cut(wrd_pct,
                   c(0, 25, 50, 75, 90, Inf),
                   c('<25', '25-49', '50-75', '76-90','\u226590'),
                   right=FALSE)) %>%
  # get rid of extra variables
  select(country,
         iso3,
         year,
         m_wrd,
         dx_test_sites,
         wrd_pct,
         var) %>%
  left_join(hbc30, by= "iso3") %>%
  left_join(hbtbhiv30, by= "iso3") %>%
  left_join(hbmdr30, by= "iso3") %>%
  rename(hbc30=group_type.x,hbtbhiv30=group_type.y,hbmdr30=group_type) %T>%
  write_csv(here("./report/ch3_data/f3.2.6_data.csv"))

  filter(f3.2.6_data, !is.na(hbc30), wrd_pct>50) %>%
  nrow()

f3.2.6_txt <- f3.2.6_data %>%
  summarise(median = median(wrd_pct, na.rm=TRUE),
            quantile = list(quantile(wrd_pct, probs = seq(.25, 0.75, by = .5), na.rm = TRUE))) %>%
  unnest_wider(quantile) %>%
  ungroup() %>%
  mutate(entity = "Global") %>%
  select(entity,
         median,
         q1=`25%`,
         q3=`75%`) %>%
  mutate(hbc_wrd50over = filter(f3.2.6_data, !is.na(hbc30), wrd_pct>=50) %>% nrow())

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.2.7 ----
# (World map showing number of WRD tests per capita as an indicator of the case finding effort)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# wrd_percapita is being used also for 3.2.9
wrd_percapita <- strategy %>%
  filter(year == report_year - 1) %>%
  select(country,iso3,year,g_whoregion,m_wrd_tests_performed) 

pop <- estimates_population %>%
  filter(year == report_year - 1) %>%
  select(iso3, year, pop)

inc <- filter(est_country, year == report_year - 1)  %>%
  select(iso3,
         year,
         inc,
         e_inc_num = inc.num) 

not <- filter(notification, year == report_year - 1)  %>%
  select(iso3,
         year,
         c_newinc) 


wrd_percapita <- wrd_percapita %>%
  inner_join(.,pop,by=c("iso3","year")) %>%
  inner_join(.,inc,by=c("iso3","year")) %>%
  inner_join(.,not,by=c("iso3","year")) 
  
wrd_percapita <- wrd_percapita %>%
  # Calculate the proportion
  mutate(wrd_ppp           = ifelse(m_wrd_tests_performed > 0,
                          m_wrd_tests_performed *100000 / pop,
                          NA)) %>%

  mutate(wrd_per_notif     = ifelse(m_wrd_tests_performed > 0,
                                m_wrd_tests_performed / c_newinc,
                                NA)) %>%
  
  # Assign the categories for the map
  mutate(var_ppp = cut(wrd_ppp,
                   c(0, 200, 400, 800, 1600, Inf),
                   c('<200', '200-399', '400-799', '800-1599','\u22651600'),
                   right=FALSE)) %>%
  mutate(var_notif = cut(wrd_per_notif,
                         c(0, 2, 4, 6, 8, 10, Inf),
                         c('<2', '2-3.9', '4-5.9', '6-7.9', '8-9.9', '\u226510'),
                         right=FALSE)) 

f3.2.7_data <- wrd_percapita %>%
# get rid of extra variables
  select(country,
         iso3,
         year,
         g_whoregion, 
         wrd_ppp,
         var_ppp
         # wrd_per_notif,
         # var_notif
         ) %>%
  left_join(hbc30, by= "iso3") %>%
  left_join(hbtbhiv30, by= "iso3") %>%
  left_join(hbmdr30, by= "iso3") %>%
  rename(hbc30=group_type.x,hbtbhiv30=group_type.y,hbmdr30=group_type) %T>%
  write_csv(here("./report/ch3_data/f3.2.7_data.csv"))


f3.2.7_txt <- f3.2.7_data %>%
  summarise(median = median(wrd_ppp, na.rm=TRUE),
            quantile = list(quantile(wrd_ppp, probs = seq(.25, 0.75, by = .5), na.rm = TRUE))) %>%
  unnest_wider(quantile) %>%
  ungroup() %>%
  mutate(entity = "Global") %>%
  select(entity,
         median,
         q1=`25%`,
         q3=`75%`) %>%
  mutate(hbc_wrd1000over = filter(f3.2.7_data, !is.na(hbc30), wrd_ppp>=1000) %>% nrow())

f3.2.7_txt_list <- f3.2.7_data %>%
  filter(wrd_ppp>=1000, !is.na(hbc30)) %>%
  select(country) %>%
  arrange()


# sample ploting
f3.2.7b_country <- wrd_percapita %>%
  filter(iso3 %in% hbc30$iso3) %>%
  select(entity = country,
         median = wrd_ppp)%>%
  mutate(q1=NA,
         q3=NA) %>%
  arrange(desc(median))

# f3.2.7b_region <- wrd_percapita %>%
#   group_by(g_whoregion) %>%
#   summarise(median = median(wrd_ppp, na.rm=TRUE),
#             quantile = list(quantile(wrd_ppp, probs = seq(.25, 0.75, by = .5), na.rm = TRUE))) %>%
#   unnest_wider(quantile) %>%
#   ungroup() %>%
#   
#   # merge with regional names and simplify to match structure of country table
#   inner_join(who_region_shortnames, by = "g_whoregion") %>%
#   select(entity,
#          median,
#          q1=`25%`,
#          q3=`75%`)
# 
# f3.2.7b_global <- wrd_percapita %>%
#   summarise(median = median(wrd_ppp, na.rm=TRUE),
#             quantile = list(quantile(wrd_ppp, probs = seq(.25, 0.75, by = .5), na.rm = TRUE))) %>%
#   unnest_wider(quantile) %>%
#   ungroup() %>%
#   mutate(entity = "Global") %>%
#   # merge with regional names and simplify to match structure of country table
#   select(entity,
#          median,
#          q1=`25%`,
#          q3=`75%`)

f3.2.7b_txt <- f3.2.7b_country %>% filter(median<2) %>% nrow()

f3.2.7b_txt_list <- f3.2.7b_country %>% slice(1)


# Bring them all together
# - - - - - - - - - - - - -

# Create dummy records so can see a horizontal line in the output to separate countries, regions and global parts
# wrd_dummy1 <- data.frame(entity = " ", median = NA, q1 = 0, q3 = 100)
# wrd_dummy2 <- data.frame(entity = "  ", median = NA, q1 = 0, q3 = 100)


# Create combined dataframe in order of countries then regional and global estimates
f3.2.7b_data <- rbind(f3.2.7b_country#,
                      # wrd_dummy1,
                      # f3.2.7b_region,
                      # wrd_dummy2,
                      # f3.2.7b_global
                      ) %>%
  
  # The dataframe is in the order I want, so make entity an ordered factor based on
  # what I already have. That way ggplot will not reorder by entity name
  # But I need to reverse order for plotting
  
  mutate(entity = factor(entity,
                         levels = rev(entity))) 


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.2.8 ----
# (World map showing proportion of positive WHO-recommended rapid diagnostic tests)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

wrd_country <- strategy %>%
  filter(year == report_year - 1) %>%
  select(iso3,
         country,
         year,
         g_whoregion,
         m_wrd_tests_positive,
         m_wrd_tests_performed) %>%
  
  # Calculate the proportion
  mutate(pos_wrd_pct = ifelse(m_wrd_tests_performed > 0,
                          m_wrd_tests_positive * 100 / m_wrd_tests_performed,
                          NA))

# for mapping
# wrd_country <- wrd_country %>%
#   # Assign the categories for the map
#   mutate(var = cut(pos_wrd_pct,
#                    c(0, 8, 14, 20, 30, Inf),
#                    c('<8', '8-14', '15-19', '20-29','\u226530'),
#                    right=FALSE)) %>%
#   # get rid of extra variables
#   select(country,
#          iso3,
#          year,
#          pos_wrd_pct,
#          var) 

# for HBC plotting
f3.2.8_country <- wrd_country %>%
  filter(iso3 %in% hbc30$iso3) %>%
  select(entity = country,
         median = pos_wrd_pct) %>%
  mutate(q1=NA,
         q3=NA) %>%
  arrange(desc(median))

f3.2.8_region <- wrd_country %>%
  group_by(g_whoregion) %>%
  summarise(median = median(pos_wrd_pct, na.rm=TRUE),
            quantile = list(quantile(pos_wrd_pct, probs = seq(.25, 0.75, by = .5), na.rm = TRUE))) %>%
  unnest_wider(quantile) %>%
  ungroup() %>%

  # merge with regional names and simplify to match structure of country table
  inner_join(who_region_shortnames, by = "g_whoregion") %>%
  select(entity,
         median,
         q1=`25%`,
         q3=`75%`) %>%
  arrange(desc(median))

f3.2.8_global <- wrd_country %>%
  summarise(median = median(pos_wrd_pct, na.rm=TRUE),
            quantile = list(quantile(pos_wrd_pct, probs = seq(.25, 0.75, by = .5), na.rm = TRUE))) %>%
  unnest_wider(quantile) %>%
  ungroup() %>%
  mutate(entity = "Global") %>%
  # merge with regional names and simplify to match structure of country table
  select(entity,
         median,
         q1=`25%`,
         q3=`75%`) 
  
# Bring them all together
# - - - - - - - - - - - - -

# Create dummy records so can see a horizontal line in the output to separate countries, regions and global parts
wrd_dummy1 <- data.frame(entity = " ", median = NA, q1 = 0, q3 = 100)
wrd_dummy2 <- data.frame(entity = "  ", median = NA, q1 = 0, q3 = 100)


# Create combined dataframe in order of countries then regional and global estimates
f3.2.8_data <- rbind(f3.2.8_country,
                      wrd_dummy1,
                      f3.2.8_region,
                      wrd_dummy2,
                      f3.2.8_global) %>%
  
  # The dataframe is in the order I want, so make entity an ordered factor based on
  # what I already have. That way ggplot will not reorder by entity name
  # But I need to reverse order for plotting
  
  mutate(entity = factor(entity,
                         levels = rev(entity)))  %T>%
  write_csv(here("./report/ch3_data/f3.2.8_data.csv"))

# textS
f3.2.8_txt <- f3.2.8_data %>%
  filter(entity == "Global")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.2.9 ----
# (World map showing Number of WRD per notification of TB cases as an indicator of all TB cases are detected and notified)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

f3.2.9_data <- wrd_percapita %>%
  # get rid of extra variables
  select(country,
         iso3,
         year,
         g_whoregion,
         wrd_per_notif,
         var_notif
  ) %>%
  left_join(hbc30, by= "iso3") %>%
  left_join(hbtbhiv30, by= "iso3") %>%
  left_join(hbmdr30, by= "iso3") %>%
  rename(hbc30=group_type.x,hbtbhiv30=group_type.y,hbmdr30=group_type) %T>%
  write_csv(here("./report/ch3_data/f3.2.9_data.csv"))


f3.2.9_text <- f3.2.9_data %>%
  summarise(median = median(wrd_per_notif, na.rm=TRUE),
            quantile = list(quantile(wrd_per_notif, probs = seq(.25, 0.75, by = .5), na.rm = TRUE))) %>%
  unnest_wider(quantile) %>%
  ungroup() %>%
  mutate(entity = "Global") %>%
  select(entity,
         median,
         q1=`25%`,
         q3=`75%`) %>%
  mutate(hbc_wrd10over = filter(f3.2.9_data, !is.na(hbc30), wrd_per_notif>=7) %>% nrow())

# sample ploting
f3.2.9b_country <- wrd_percapita %>%
  filter(iso3 %in% hbc30$iso3) %>%
  select(entity = country,
         median = wrd_per_notif)%>%
  mutate(q1=NA,
         q3=NA) %>%
  arrange(desc(median))

# f3.2.9b_region <- wrd_percapita %>%
#   group_by(g_whoregion) %>%
#   summarise(median = median(wrd_per_notif, na.rm=TRUE),
#             quantile = list(quantile(wrd_per_notif, probs = seq(.25, 0.75, by = .5), na.rm = TRUE))) %>%
#   unnest_wider(quantile) %>%
#   ungroup() %>%
# 
#   # merge with regional names and simplify to match structure of country table
#   inner_join(who_region_shortnames, by = "g_whoregion") %>%
#   select(entity,
#          median,
#          q1=`25%`,
#          q3=`75%`)
# 
# f3.2.9b_global <- wrd_percapita %>%
#   summarise(median = median(wrd_per_notif, na.rm=TRUE),
#             quantile = list(quantile(wrd_per_notif, probs = seq(.25, 0.75, by = .5), na.rm = TRUE))) %>%
#   unnest_wider(quantile) %>%
#   ungroup() %>%
#   mutate(entity = "Global") %>%
#   # merge with regional names and simplify to match structure of country table
#   select(entity,
#          median,
#          q1=`25%`,
#          q3=`75%`)

f3.2.9b_txt <- f3.2.9b_country %>% filter(median<2) %>% nrow()

f3.2.9b_txt_list <- f3.2.9b_country %>% slice(1)


# Bring them all together
# - - - - - - - - - - - - -

# Create dummy records so can see a horizontal line in the output to separate countries, regions and global parts
# wrd_dummy1 <- data.frame(entity = " ", median = NA, q1 = 0, q3 = 100)
# wrd_dummy2 <- data.frame(entity = "  ", median = NA, q1 = 0, q3 = 100)


# Create combined dataframe in order of countries then regional and global estimates
f3.2.9b_data <- rbind(f3.2.9b_country#,
                      # wrd_dummy1,
                      # f3.2.9b_region,
                      # wrd_dummy2,
                      # f3.2.9b_global
                      ) %>%

  # The dataframe is in the order I want, so make entity an ordered factor based on
  # what I already have. That way ggplot will not reorder by entity name
  # But I need to reverse order for plotting

  mutate(entity = factor(entity,
                         levels = rev(entity)))  %T>%
  write_csv(here("./report/ch3_data/f3.2.9b_data.csv"))

# texts!
f3.2.9_txt <- f3.2.9_data %>%
  filter(wrd_per_notif<2, !is.na(hbc30)) 

f3.2.9_txt_list <- f3.2.9_data %>%
  filter(wrd_per_notif>7, !is.na(hbc30)) 


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.2.10 ----
# (Panel plot of TB cases with known HIV status by WHO region and globally since 2004)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Calculate regional aggregates
f3.2.10_data <- TBHIV_for_aggregates %>%
  filter(year >= 2004) %>%
  select(g_whoregion,
         year,
         hivtest_pct_numerator,
         hivtest_pct_denominator) %>%
  group_by(year, g_whoregion) %>%
  summarise(across(hivtest_pct_numerator:hivtest_pct_denominator, sum, na.rm = TRUE)) %>%
  # merge with regional names
  inner_join(who_region_shortnames, by = "g_whoregion") %>%
  select(-g_whoregion) %>%
  ungroup()

# Calculate Global aggregates
hivstatus_global <- TBHIV_for_aggregates %>%
  filter(year >= 2004) %>%
  select(year,
         hivtest_pct_numerator,
         hivtest_pct_denominator) %>%
  group_by(year) %>%
  summarise(across(hivtest_pct_numerator:hivtest_pct_denominator, sum, na.rm = TRUE)) %>%
  mutate(entity = "Global")

# COmbine regional with global
f3.2.10_data <- rbind(f3.2.10_data, hivstatus_global) %>%

  # Calculate % with known HIV status
  mutate(hivstatus_pct = hivtest_pct_numerator * 100 / hivtest_pct_denominator) %>%

  # Change the order of the entities
  mutate(entity = factor(entity,
                         levels = c("Global", "African Region", "Region of the Americas", "South-East Asia Region",
                                    "European Region", "Eastern Mediterranean Region", "Western Pacific Region")))

# summary dataset for quoting numbers in the text
f3.2.10_txt <- filter(f3.2.10_data, year>=report_year-2 & entity %in% c("African Region", "European Region", "Global")) %>%
  select(year, entity, hivstatus_pct) %>%
  pivot_wider(names_from = c(entity, year),
              values_from = hivstatus_pct) %>%
  # handle spaces in column names
  mutate(AFR_2021 = `African Region_2021`,
         EUR_2021 = `European Region_2021`)

# Add the numbers that are HIV-positive
f3.2.10_txt <- filter(TBHIV_for_aggregates, year == report_year-1) %>%
  summarise(across(c(hivtest_pos_pct_numerator, hivtest_pos_pct_denominator), sum, na.rm = TRUE)) %>%
  mutate(hivtest_pos_pct = hivtest_pos_pct_numerator * 100 / hivtest_pos_pct_denominator) %>%
  cbind(f3.2.10_txt)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.2.11 ----
# (World map showing percent of TB cases with known HIV status)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

f3.2.11_data <- TBHIV_for_aggregates %>%
  filter(year >= report_year - 2) %>%
  select(iso3,
         country,
         year,
         hivtest_pct_denominator,
         hivtest_pct_numerator) %>%

  # Calculate % with known HIV status
  mutate(hivstatus_pct = ifelse(hivtest_pct_denominator > 0,
                                hivtest_pct_numerator * 110 / hivtest_pct_denominator,
                                NA))  %>%

  # Assign the categories for the map
  mutate(var = cut(hivstatus_pct,
                   c(0, 50, 76, 90, Inf),
                   c('0-49', '50-75', '76-89', "\u226590"),
                   right=FALSE)) %>%

  # get rid of extra variables
  select(country,
         iso3,
         year,
         hivstatus_pct,
         var)

# Find the countries with empty data for latest year and see if there are data for the previous year
hivstatus_prev_year_data <- f3.2.11_data %>%
  filter(year == report_year - 1 & is.na(hivstatus_pct)) %>%
  select(iso3) %>%
  inner_join(filter(f3.2.11_data, year == report_year - 2)) %>%
  filter(!is.na(hivstatus_pct))

# Now combine into one dataframe, with previous data used if latest year's data are not available
f3.2.11_data <- f3.2.11_data %>%
  filter(year == report_year - 1) %>%
  anti_join(hivstatus_prev_year_data, by= "iso3") %>%
  rbind(hivstatus_prev_year_data)

# Summary data for text
f3.2.11_txt <- filter(f3.2.11_data, hivstatus_pct >= 90) %>%
  summarise(over_90 = n())


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.2.12 ----
# (Panel plot of TB cases tested for susceptibility to rifampicin by WHO region and globally since 2009)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# The calculations here are a bit messy for years prior to 2018. It is much simpler starting in 2018 when
# we switched to using routine DR surveillance records only. Prior to that we mixed and matched between
# that and the DR-TB detection section of the data collection form. Keeping this for consistency with
# previously published reports

# 1. Get DR-TB detection data
dst_notif_data <- notification %>%
  filter(year >= 2009) %>%
  select(year,
         iso3,
         g_whoregion,
         new_labconf,
         c_ret,
         rdst_new,
         rdst_ret) %>%

  rowwise() %>%
  mutate(

    # numerator
    dst_notif_num = ifelse(is.na(rdst_new) & is.na(rdst_ret),
                           NA,
                           sum(c_across(rdst_new:rdst_ret), na.rm = TRUE)),

    # denominator is a bit of a fudge: new_labconf + c_ret
    dst_notif_denom = ifelse(is.na(new_labconf) & is.na(c_ret),
                             NA,
                             sum(c_across(new_labconf:c_ret), na.rm = TRUE))) %>%
  ungroup()


# 2. Get routine DR surveillance data
# (numerator and denominator variables are different according to the year, but
#  don't use 2015, 2016 numerator data)
dst_drs_data <- dr_surveillance %>%
  filter(year >= 2009) %>%
  select(year,
         iso3,
         dst_rlt_new,
         dst_rlt_ret,
         pulm_labconf_new,
         pulm_labconf_ret,
         r_rlt_new,
         r_rlt_ret) %>%

  rowwise() %>%
  mutate(
    # numerator
    dst_drs_num = ifelse(year < 2015,
                         ifelse(is.na(dst_rlt_new) & is.na(dst_rlt_ret),
                                NA,
                                sum(c_across(dst_rlt_new:dst_rlt_ret), na.rm = TRUE)),
                         ifelse(year >= 2017,
                                sum(c_across(r_rlt_new:r_rlt_ret), na.rm = TRUE),
                                NA)),
    # denominator
    dst_drs_denom = ifelse(year >= 2017,
                           sum(c_across(pulm_labconf_new:pulm_labconf_ret), na.rm = TRUE),
                            NA)
  ) %>%
  ungroup()


# Link the two data sets

dst_data <- dst_notif_data %>%
  left_join(dst_drs_data, by = c("year", "iso3")) %>%

  # To calculate the percentage DST coverage we need to identify the greater of the two numerators
  # Note the exception made for South Africa in 2017


  mutate(
    dst_num = ifelse(year == 2017 & iso3 == "ZAF",
                     dst_notif_num,
                     ifelse(year >= 2017,
                            dst_drs_num,
                           ifelse(NZ(dst_drs_num) >= NZ(dst_notif_num),
                                  dst_drs_num,
                                  dst_notif_num))),

    dst_denom = ifelse(year == 2017 & iso3 == "ZAF",
                       dst_notif_denom,
                       ifelse(year >= 2017,
                              dst_drs_denom,
                              dst_notif_denom))) %>%

  # Set numerator to NA if the denominator is NA for a country-year
  mutate(dst_num = ifelse(is.na(dst_denom), NA, dst_num)) %>%

  # Drop unwanted variables
  select(iso3,
         year,
         g_whoregion,
         dst_num,
         dst_denom) %>%

  # Drop rows with empty numerators and denominator
  filter(!is.na(dst_num) & !is.na(dst_denom))


f3.2.12_data <- dst_data %>%
  group_by(g_whoregion, year) %>%
  summarise(across(dst_num:dst_denom, sum, na.rm = TRUE)) %>%

  # merge with regional names
  inner_join(who_region_shortnames, by = "g_whoregion") %>%
  ungroup() %>%
  select(-g_whoregion)

dst_global <- dst_data %>%
  group_by(year) %>%
  summarise(across(dst_num:dst_denom, sum, na.rm = TRUE)) %>%
  mutate(entity = "Global") %>%
  ungroup()

# Phew! Bring it all together now
# COmbine regional with global
f3.2.12_data <- rbind(f3.2.12_data, dst_global) %>%

  # Calculate % tested for rifampicin resistance
  mutate(dst_pcnt = dst_num * 100 / dst_denom) %>%

  # Change the order of the entities
  mutate(entity = factor(entity,
                         levels = c("Global", "African Region", "Region of the Americas", "South-East Asia Region",
                                    "European Region", "Eastern Mediterranean Region", "Western Pacific Region")))


# summary data for text
f3.2.12_txt <- filter(f3.2.12_data, year >= report_year-3 & entity %in% c("European Region", "Global")) %>%
  select(year, entity, dst_pcnt) %>%
  pivot_wider(names_from = c(entity, year),
              names_prefix = "dst_pct_",
              values_from = dst_pcnt) %>%
  # Handle spaces in variable name
  mutate(dst_pct_EUR_2021 = `dst_pct_European Region_2021`)

# Add numbers of drug-resistant cases detected
f3.2.12_txt <- filter(notification, year >= report_year - 2) %>%
  select(year,
         c_newinc,
         conf_rr_nfqr,
         conf_rr_fqr,
         conf_rrmdr) %>%
  group_by(year) %>%
  summarise(across(c_newinc:conf_rrmdr, sum, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(dr_tb = conf_rr_nfqr + conf_rr_fqr + conf_rrmdr) %>%
  select(-conf_rrmdr) %>%
  pivot_wider(names_from = year,
              values_from = c(c_newinc,
                              conf_rr_nfqr,
                              conf_rr_fqr,
                              dr_tb)) %>%
  select(-conf_rr_nfqr_2020, -conf_rr_fqr_2020) %>%
  mutate(dr_tb_change_pct = abs(dr_tb_2021 - dr_tb_2020) * 100 / dr_tb_2020,
         c_newinc_change_pct= abs(c_newinc_2021 - c_newinc_2020) * 100 / c_newinc_2020 ) %>%
  select(-c_newinc_2021, -c_newinc_2020) %>%
  cbind(f3.2.12_txt)






# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.2.13 ----
# (World map showing percent of TB cases tested for susceptibility to rifampicin)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

f3.2.13_data <- dr_surveillance %>%
  filter(year >= report_year - 2) %>%
  select(year,
         country,
         iso3,
         pulm_labconf_new,
         pulm_labconf_ret,
         r_rlt_new,
         r_rlt_ret) %>%

  # Calculate coverage of DST percentages
  mutate(
    dst_pct = ifelse((NZ(pulm_labconf_new) + NZ(pulm_labconf_ret)) == 0 |
                       is.na(r_rlt_new) & is.na(r_rlt_ret), NA,
                     (NZ(r_rlt_new) + NZ(r_rlt_ret)) * 100 /
                       (NZ(pulm_labconf_new) + NZ(pulm_labconf_ret)))
  ) %>%

  # Assign the categories for the map
  mutate(var = cut(dst_pct,
                   c(0, 20, 50, 80, Inf),
                   c('0-19.9', '20-49.9', '50-79.9', '\u226580'),
                   right=FALSE)) %>%

  # get rid of extra variables
  select(country,
         iso3,
         year,
         dst_pct,
         var)


# Find the countries with empty data for latest year and see if there are data for the previous year
dst_prev_year_data <- f3.2.13_data %>%
  filter(year == report_year - 1 & is.na(dst_pct)) %>%
  select(iso3) %>%
  inner_join(filter(f3.2.13_data, year == report_year - 2), by = "iso3") %>%
  filter(!is.na(dst_pct))

# Now combine into one dataframe, with previous data used if latest year's data are not available
f3.2.13_data <- f3.2.13_data %>%
  filter(year == report_year - 1) %>%
  anti_join(dst_prev_year_data, by= "iso3") %>%
  rbind(dst_prev_year_data)


# summary numbers for the text
f3.2.13_txt <- filter(f3.2.13_data, dst_pct >= 80 & iso3 %in% hbmdr30$iso3) %>%
  arrange(country) %>%
  select(country)




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.2.14 ----
# (Panel plot of RR-TB cases tested for susceptibility to fluoroquinolones by WHO region and globally since 2015)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

f3.2.14_data <- dr_surveillance %>%
  filter(year >= 2015) %>%
  select(iso3,
         g_whoregion,
         year,
         # denominator changed in 2017 from mdr to rr
         mdr_new,
         mdr_ret,
         xpert_dr_r_new,
         xpert_dr_r_ret,
         rr_new,
         rr_ret,
         # numerator changed in 2017 from mdr to rr and in 2019 from sld to fq
         mdr_dst_rlt,
         rr_dst_rlt,
         rr_dst_rlt_fq) %>%

  group_by(year, g_whoregion) %>%
  summarise(across(mdr_new:rr_dst_rlt_fq, sum, na.rm = TRUE)) %>%
  ungroup() %>%

  # Calculate the numerators and denominators depending on the year
  mutate(fqdst_pct_denominator = ifelse(year < 2017,
                                        mdr_new + mdr_ret + xpert_dr_r_new + xpert_dr_r_ret,
                                        NA),
         fqdst_pct_numerator = ifelse(year < 2017,
                                      mdr_dst_rlt,
                                      NA)) %>%

  mutate(fqdst_pct_denominator = ifelse(year >= 2017,
                                        rr_new + rr_ret,
                                        fqdst_pct_denominator),
         fqdst_pct_numerator = ifelse(year %in% c(2017, 2018),
                                      rr_dst_rlt,
                                      fqdst_pct_numerator)) %>%

  mutate(fqdst_pct_numerator = ifelse(year >= 2019,
                                      rr_dst_rlt_fq,
                                      fqdst_pct_numerator)) %>%

  # merge with regional names
  inner_join(who_region_shortnames, by = "g_whoregion") %>%

  # get rid of extra variables
  select(entity,
         year,
         fqdst_pct_numerator,
         fqdst_pct_denominator)

# Calculate global aggregaes
fqdst_global <- f3.2.14_data %>%
  group_by(year) %>%
  summarise(across(fqdst_pct_numerator:fqdst_pct_denominator, sum, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(entity = 'Global')

# Add global to the regional aggregates
f3.2.14_data <- rbind(f3.2.14_data, fqdst_global) %>%

  # Calculate the percentages
  mutate(fqdst_pct = fqdst_pct_numerator * 100 / fqdst_pct_denominator) %>%

  # Change the order of the entities
  mutate(entity = factor(entity,
                         levels = c("Global", "African Region", "Region of the Americas", "South-East Asia Region",
                                    "European Region", "Eastern Mediterranean Region", "Western Pacific Region")))




f3.2.14_txt <- filter(notification, year >= (report_year - 2)) %>%
  select(year,
         g_whoregion,
         c_notified,
         c_newinc,
         new_labconf, new_clindx, new_ep,
         ret_rel_labconf, ret_rel_clindx, ret_rel_ep,
         newrel_hivpos,
         conf_rr_nfqr,
         conf_rr_fqr) %>%
  
  # calculate regional aggregates
  group_by(year, g_whoregion) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  
  # merge with regional names
  inner_join(who_region_shortnames, by = "g_whoregion") %>%
  arrange(entity) %>%
  ungroup() %>%
  select(-g_whoregion)


# Add global summary to the regional summary
f3.2.14_txt <- f3.2.14_txt %>%
  group_by(year) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  mutate(entity="Global") %>%
  mutate(rr_nfqr_fqr = conf_rr_nfqr+conf_rr_fqr) %>%
  mutate(rr_nfqr_fqr_p = lag(rr_nfqr_fqr))%>%
  mutate(rr_nfqr_pct_dif = (rr_nfqr_fqr - rr_nfqr_fqr_p)*100/rr_nfqr_fqr_p) %>%
  mutate(c_newinc_p   = lag(c_newinc)) %>%
  mutate(c_newinc_pct_dif = (c_newinc - c_newinc_p)*100/c_newinc_p) %>%
  select(entity, year, 
         conf_rr_nfqr,
         conf_rr_fqr,
         rr_nfqr_fqr, rr_nfqr_pct_dif, c_newinc_pct_dif) %>%
  pivot_wider(names_from = year,
              values_from = conf_rr_nfqr:c_newinc_pct_dif) 



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.2.15 ----
# (World map showing percent of MDR/RR-TB cases tested for susceptibility to fluoroquinolones)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

f3.2.15_data <- dr_surveillance %>%
  filter(year >= report_year - 1) %>%
  select(iso3,
         country,
         year,
         rr_new,
         rr_ret,
         rr_dst_rlt_fq) %>%
  
  # Calculate percentage RR cases with 2nd line DST
  mutate(fqdst_pct = ifelse( (NZ(rr_new) + NZ(rr_ret)) > 0,
                             rr_dst_rlt_fq * 100 / (NZ(rr_new) + NZ(rr_ret)),
                             NA)) %>%
  
  # Assign the categories for the map
  mutate(var = cut(fqdst_pct,
                   c(0, 20, 50, 80, Inf),
                   c('0-19.9', '20-49.9', '50-79.9', '\u226580'),
                   right=FALSE)) %>%
  
  # get rid of extra variables
  select(country,
         iso3,
         year,
         fqdst_pct,
         var)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.2.16 ----
# (World map showing percent of pre-XDR-TB cases tested for susceptibility to bedaquiline)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

f3.2.16_data <- dr_surveillance %>%
  filter(year >= report_year - 1) %>%
  select(iso3,
         country,
         year,
         rr_fqr,
         rr_fqr_bdqr_lzdr,rr_fqr_bdqs_lzdr,rr_fqr_bdqu_lzdr,
         rr_fqr_bdqr_lzds,rr_fqr_bdqs_lzds,rr_fqr_bdqu_lzds,
         rr_fqr_bdqr_lzdu,rr_fqr_bdqs_lzdu,rr_fqr_bdqu_lzdu) %>%
  
  # Calculate percentage RR cases with 2nd line DST
  mutate(bddst_pct = ifelse( NZ(rr_fqr) > 0,
                             (NZ(rr_fqr_bdqr_lzdr)+NZ(rr_fqr_bdqs_lzdr)+NZ(rr_fqr_bdqr_lzds)+NZ(rr_fqr_bdqs_lzds)+NZ(rr_fqr_bdqr_lzdu)+NZ(rr_fqr_bdqs_lzdu)) * 100 / NZ(rr_fqr),
                             NA)) %>%
  
  # Assign the categories for the map
  mutate(var = cut(bddst_pct,
                   c(0, 20, 50, 80, Inf),
                   c('0-19.9', '20-49.9', '50-79.9', '\u226580'),
                   right=FALSE)) %>%
  
  # get rid of extra variables
  select(country,
         iso3,
         year,
         bddst_pct,
         var)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.2.17 ----
# (World map showing percent of pre-XDR-TB cases tested for susceptibility to linezolid)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

f3.2.17_data <- dr_surveillance %>%
  filter(year >= report_year - 1) %>%
  select(iso3,
         country,
         year,
         rr_fqr,
         rr_fqr_bdqr_lzdr,rr_fqr_bdqs_lzdr,rr_fqr_bdqu_lzdr,
         rr_fqr_bdqr_lzds,rr_fqr_bdqs_lzds,rr_fqr_bdqu_lzds,
         rr_fqr_bdqr_lzdu,rr_fqr_bdqs_lzdu,rr_fqr_bdqu_lzdu) %>%
  
  # Calculate percentage RR cases with 2nd line DST
  mutate(lzdst_pct = ifelse( NZ(rr_fqr) > 0,
                             (NZ(rr_fqr_bdqr_lzdr)+NZ(rr_fqr_bdqs_lzdr)+NZ(rr_fqr_bdqr_lzds)+NZ(rr_fqr_bdqs_lzds)+NZ(rr_fqr_bdqu_lzdr)+NZ(rr_fqr_bdqu_lzds)) * 100 / NZ(rr_fqr),
                             NA)) %>%
  
  # Assign the categories for the map
  mutate(var = cut(lzdst_pct,
                   c(0, 20, 50, 80, Inf),
                   c('0-19.9', '20-49.9', '50-79.9', '\u226580'),
                   right=FALSE)) %>%
  
  # get rid of extra variables
  select(country,
         iso3,
         year,
         lzdst_pct,
         var)




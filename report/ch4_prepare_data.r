# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data preparation script for ch4.rmd
# Hazim Timimi, July 2022
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


# Load chapter 3 and 4 packages, settings and data
source(here::here('report/ch3_ch4_load_data.r'))


# Also load the GHO data package
library(ghost)

# Also load the ltbi estimates
estimates_ltbi <- get_timestamped_csv('ltbi', csv_datestamp)


# Set whether or not to show rifapentine and BCG figures
show_rifapentine <- TRUE
show_bcg <- TRUE

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Create dataframe for figure 4.1 ----
# (Bar chart showing numbers provided with TB preventive treatment each year since 2015)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

f4.1_data <- filter(notification, year %in% seq(2015, report_year - 1)) %>%
  select(iso2,
         year,
         hiv_ipt_reg_all,
         hiv_ipt,
         # These next ones introduced dcyear 2021 by GAM
         hiv_all_tpt,
         hiv_new_tpt,
         hiv_elig_all_tpt,
         hiv_elig_new_tpt) %>%

  # Use the first non-empty tpt/ipt variable (not sure this is the best approach as it has become so
  # darn complicated in 2021 ...)
  mutate(hiv_tpt = coalesce(hiv_ipt_reg_all, hiv_ipt, hiv_all_tpt, hiv_elig_all_tpt, hiv_new_tpt, hiv_elig_new_tpt)) %>%

  select(iso2, year, hiv_tpt) %>%

  # Join to tpt data for household contacts in the strategy view
  inner_join( select(strategy,
                     iso2, year,
                     # next one added 2016 dcyear
                     newinc_con04_prevtx,
                     # next one used 2018 dcyear only
                     newinc_con5plus_prevtx,
                     # next one added 2019 dcyear
                     newinc_con_prevtx), by=c("iso2", "year")) %>%

  # Calculate the "5 and over" fraction of tpt for household contacts
  mutate(prevtx_5plus = ifelse(NZ(newinc_con_prevtx) > 0 & NZ(newinc_con04_prevtx) > 0,
                               newinc_con_prevtx - newinc_con04_prevtx,
                               newinc_con_prevtx)) %>%

  # Convert negative prevtx_5plus caused by weird combination of carry overs to zero
  mutate(prevtx_5plus = ifelse(NZ(prevtx_5plus) < 0 , 0, prevtx_5plus)) %>%

  # deal with 2017 variable
  mutate(prevtx_5plus = ifelse(year == 2017 ,
                               newinc_con5plus_prevtx,
                               prevtx_5plus)) %>%

  # Keep variables for HIV, contacts < 5 and contacts 5 plus
  select(iso2, year,
         hiv_tpt,
         house_con04_tpt = newinc_con04_prevtx,
         house_con5plus_tpt = prevtx_5plus) %>%

  # Calculate the global totals by year ready for the plot
  group_by(year) %>%
  summarise_at(vars(-iso2), sum, na.rm=TRUE) %>%
  ungroup() %>%

  # Finally, switch to a long format ready for plotting
  pivot_longer(cols = hiv_tpt:house_con5plus_tpt,
               names_to = "TPT_category",
               values_to = "how_many")

# Create summary stats for the text
f4.1_txt <- f4.1_data %>%
  group_by(year) %>%
  summarise(tot_tpt = sum(how_many)) %>%
  pivot_wider(names_from = year,
              names_prefix = "tot_tpt_",
              values_from = tot_tpt) %>%
  mutate(delta_19_20 = tot_tpt_2020 - tot_tpt_2019,
         pct_21_20 = abs(tot_tpt_2021 - tot_tpt_2020) * 100 / tot_tpt_2020,
         tot_tpt_18_21 = tot_tpt_2018 + tot_tpt_2019 + tot_tpt_2020 + tot_tpt_2021) %>%
  mutate(pct_tpt_target = tot_tpt_18_21 * 100 / 30e6)

# Add stats just for household contacts aged 5 and over
f4.1_txt <- filter(f4.1_data, TPT_category == "house_con5plus_tpt" & year > 2019) %>%
  group_by(year) %>%
  summarise(tot_con_tpt = sum(how_many)) %>%
  pivot_wider(names_from = year,
              names_prefix = "tot_con_tpt_",
              values_from = tot_con_tpt) %>%
  mutate(con_21_20_pct = abs(tot_con_tpt_2021 - tot_con_tpt_2020) * 100 / tot_con_tpt_2020) %>%
  cbind(f4.1_txt)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 4.2 ----
# (Irwin's doughnuts -- % completion of UNHLM targets)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

f4.2_data <- filter(f4.1_data, year>=2018) %>%
  group_by(TPT_category) %>%
  summarise(tot=sum(how_many)) %>%
  ungroup()

# Add total for all categories
f4.2_data <- f4.2_data %>%
  summarise(tot = sum(tot)) %>%
  mutate(TPT_category = "all_tpt") %>%
  rbind(f4.2_data)

# Calculate percent of targets reached
f4.2_targets <- data.frame(TPT_category = c('all_tpt', 'hiv_tpt', 'house_con04_tpt', 'house_con5plus_tpt'),
                           target = c(30e6, 6e6, 4e6, 20e6))

f4.2_data <- f4.2_data %>%
  inner_join(f4.2_targets, by = "TPT_category") %>%
  mutate(target_completion = tot * 100 / target)


f4.2_txt <- filter(f4.1_data, year>=2018) %>%
  group_by(TPT_category) %>%
  summarise(tot=sum(how_many)) %>%
  pivot_wider(names_from = TPT_category,
              values_from = tot) %>%
  mutate(all_con = house_con04_tpt + house_con5plus_tpt) %>%
  # Calculate percent of target
  mutate(con04_target = house_con04_tpt * 100 / 4e6,
         con5plus_target = house_con5plus_tpt * 100 / 20e6,
         all_con_target = all_con * 100 / 24e6)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Create dataframe for figure 4.3 ----
# (Map showing countries using rifapentine for TB preventive treatment)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Suppress this part in July 2022 for later resurrection

if(show_rifapentine){

  # Load external file sent by Sanofi to Dennis on use of rifapentime
  f4.3_data <- read.csv(here('csv/rifapentine_2022-07-25.csv')) %>%
    # simplify file
    select(iso3,
           description = var) %>%

    #re-create the var variable as used in the 2021 report
    mutate(var = case_when(description=="Current or past use" ~ 1,
                           description=="Used in trials only" ~ 2)) %>%

    # Assign the categories for the map
    mutate(var = factor(var,
                        levels = c(1, 2),
                        labels = c("Current or past use", "Used in trials only"))) %>%

    # remove the description field
    select(-description)


  # Summary for the text
  f4.3_txt <- filter(f4.3_data, var=="Current or past use") %>%
    summarise(countries_used = n())

  # Get number of countries and people treated from the GTB database
  rif_nums <- strategy %>%
    filter(year == report_year - 1 & tpt_short_regimens_used == 1) %>%
    select(iso3, tpt_short_regimens_used,
           tpt_1hp,
           tpt_3hp,
           tpt_3hr,
           tpt_4r,
           prevtx_short_rifamycin)

  # Stats for countries reporting numbers of people
  f4.3_txt$people_tx_2021 <- filter(rif_nums, prevtx_short_rifamycin > 0) %>%
    select(prevtx_short_rifamycin) %>%
    sum()

  f4.3_txt$country_tx_2021 <- filter(rif_nums, prevtx_short_rifamycin > 0) %>%
    nrow()


  # Stats for countries reporting regimen type (some didn;t report numbers of people)
  f4.3_txt$rifapentine <- filter(rif_nums, tpt_1hp == 1 | tpt_3hp == 1) %>%
    nrow()

  f4.3_txt$rifampicin <- filter(rif_nums, tpt_3hr == 1 | tpt_4r == 1) %>%
    nrow()

}






# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Create dataframe for figure 4.4  ----
# (Map showing evaluation for TB disease and TB infection among household contacts of confirmed pulmonary TB cases)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

f4.4_data <- filter(strategy, year == report_year - 1) %>%
  select(iso3,
         country,
         newinc_con,
         newinc_con_screen,
         rev_newinc_con,
         rev_newinc_con_screen) %>%

  # Calculate the proportion screened
  mutate(screened_pct = ifelse(!is.na(newinc_con) & newinc_con > 0,
                               newinc_con_screen * 100 / newinc_con,
                               NA)) %>%

  # Add any countries that reported from review of patient records
  mutate(screened_pct = ifelse(is.na(screened_pct) & !is.na(rev_newinc_con) & NZ(rev_newinc_con) > 0,
                               rev_newinc_con_screen * 100 / rev_newinc_con,
                               screened_pct))  %>%

  # Assign the categories for the map
  mutate(var = cut(screened_pct,
                   c(0, 25, 50, 75, Inf),
                   c('0-24', '25-49', '50-74', '\u226575'),
                   right = FALSE))


# Summary for the text
f4.4_txt <- filter(strategy, year >= report_year - 2 &
                     !is.na(newinc_con) & !is.na(newinc_con_screen)) %>%
  select(iso3,
         year,
         newinc_con,
         newinc_con_screen) %>%
  group_by(year) %>%
  summarise(across(newinc_con:newinc_con_screen, sum, na.rm=TRUE)) %>%
  ungroup() %>%
  mutate(screened_pct = newinc_con_screen * 100 / newinc_con) %>%
  pivot_wider(names_from = year,
              values_from = newinc_con:screened_pct) %>%
  mutate(change_con_21_20_pct = abs(newinc_con_2021 - newinc_con_2020) * 100 / newinc_con_2020,
         change_screen_21_20_pct = abs(newinc_con_screen_2021 - newinc_con_screen_2020) * 100 / newinc_con_screen_2020)





# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Create dataframe for figure 4.5  ----
# (Map showing percentage of household contacts aged under 5 years provided with TB preventive treatment)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

f4.5_data <- filter(estimates_ltbi, year == report_year - 1) %>%

  select(iso3,
         country,
         e_prevtx_kids_pct,
         e_prevtx_eligible,
         newinc_con04_prevtx)  %>%

  # Assign the categories for the map
  mutate(var = cut(e_prevtx_kids_pct,
                   c(0, 25, 50, 90, Inf),
                   c('0-24', '25-49', '50-89', '\u226590'),
                   right = FALSE))






# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Create dataframe for figure 4.6  ----
# (Panel plot showing percentage completion vs number contacts started TPT by WHO region)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

f4.6_data <- filter(strategy, year == report_year - 2) %>%

  select(country,
         g_whoregion,
         newinc_con_prevtx,
         newinc_con_prevtx_cmplt) %>%

  # Calculate completion rate
  mutate(pct_completed = ifelse(!is.na(newinc_con_prevtx_cmplt) & NZ(newinc_con_prevtx) > 0,
                                newinc_con_prevtx_cmplt * 100 /newinc_con_prevtx ,
                                NA)) %>%

  # merge with regional names
  inner_join(who_region_shortnames, by = "g_whoregion") %>%
  select(-g_whoregion) %>%

  # filter out empty lines
  filter(!is.na(pct_completed))


# Summary for the text
f4.6_txt <- f4.6_data %>%
  summarise(countries = n(),
            median = median(pct_completed),
            q1 = unname(quantile(pct_completed, 0.25)),
            q3 = unname(quantile(pct_completed, 0.75)))



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Create dataframe for figure 4.5 (OLD -- FIGURE DROPPED )----
# (Map showing provision of TB preventive treatment among PLHIV who started antiretroviral treatment)
# Keeping code as we use summary stats in the text
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

f4.5_old_data <- filter(notification, year == report_year - 1) %>%
  select(country,
         iso3,
         hiv_elig_new_tpt,
         hiv_elig_new,
         hiv_new_tpt,
         hiv_new) %>%

  # Calculate % coverage among newly enrolled on ART. Prioritise the variables that are not restricted
  # to "eligible" PLHIV

  mutate(coverage = ifelse(!is.na(hiv_new_tpt) & NZ(hiv_new) > 0,
                           hiv_new_tpt * 100 / hiv_new,
                           ifelse(!is.na(hiv_elig_new_tpt) & NZ(hiv_elig_new) > 0,
                                  hiv_elig_new_tpt * 100 / hiv_elig_new,
                                  NA))) %>%

  # Assign the categories for the map
  mutate(var = cut(coverage,
                   c(0, 25, 50, 75, Inf),
                   c('0-24', '25-49', '50-74', '\u226575'),
                   right = FALSE))

# Summary for the text
f4.5_old_txt <- filter(f4.5_old_data,
                       !is.na(hiv_elig_new_tpt) & !is.na(hiv_elig_new) & iso3 %in% hbtbhiv30$iso3) %>%
  summarise(countries = n(),
            median = median(coverage),
            q1 = unname(quantile(coverage, 0.25)),
            q3 = unname(quantile(coverage, 0.75)))

rm(f4.5_old_data)


# Add summary for the text on TPT completion in PLHIV (the figure we were going to use was dropped)
completion_txt <- filter(notification, year == report_year - 2 &
                           hiv_all_tpt_completed > 0 & hiv_all_tpt_started > 0) %>%
  select(iso3,
         hiv_all_tpt_completed,
         hiv_all_tpt_started)  %>%
  mutate(completion = hiv_all_tpt_completed * 100 / hiv_all_tpt_started) %>%
  summarise(countries = n(),
            median = median(completion),
            q1 = unname(quantile(completion, 0.25)),
            q3 = unname(quantile(completion, 0.75)))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Create dataframe for figure 4.7  ----
# (Panel plots showing numbers of people living with HIV provided with TB preventive treatment each year since 2005 by WHO region and globally)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

f4.7_data <- filter(notification, year %in% seq(2005, report_year - 1)) %>%
  select(iso2,
         g_whoregion,
         year,
         hiv_ipt_reg_all,
         hiv_ipt,
         # These next ones introduced dcyear 2021 by GAM
         hiv_all_tpt,
         hiv_new_tpt,
         hiv_elig_all_tpt,
         hiv_elig_new_tpt) %>%

  # Create calculated variables for TPT among all enrolled and TPT among newly enrolled
  # filling in gaps for missing data
  # The choice for 2020 GAM data is rather murky ...
  mutate(hiv_tpt_all = coalesce(hiv_ipt_reg_all, hiv_ipt, hiv_all_tpt, hiv_elig_all_tpt, hiv_new_tpt, hiv_elig_new_tpt),
         hiv_tpt_new = coalesce(hiv_ipt, hiv_ipt_reg_all, hiv_new_tpt, hiv_elig_new_tpt, hiv_all_tpt, hiv_elig_all_tpt)) %>%

  # Calculate regional aggregates
  group_by(year, g_whoregion) %>%
  summarise_at(vars(hiv_tpt_all:hiv_tpt_new),
               sum,
               na.rm = TRUE) %>%
  ungroup() %>%

  # merge with regional names
  inner_join(who_region_shortnames, by = "g_whoregion") %>%
  select(-g_whoregion)

# Create global aggregates
f4.7_data_global <- f4.7_data %>%
  group_by(year) %>%
  summarise_at(vars(hiv_tpt_all:hiv_tpt_new),
               sum,
               na.rm = TRUE) %>%
  mutate(entity = 'Global')

# Add global to the regional aggregates
f4.7_data <- rbind(f4.7_data, f4.7_data_global)

# Only want hiv_tpt_all for years after 2016
f4.7_data <- f4.7_data %>%
  mutate(hiv_tpt_all = ifelse(year < 2017,
                              NA,
                              hiv_tpt_all))

# Change the entity order
f4.7_data$entity <- factor(f4.7_data$entity,
                           levels = c("Global", "African Region", "Region of the Americas", "South-East Asia Region",
                                      "European Region", "Eastern Mediterranean Region", "Western Pacific Region"))

# Summary for quoting in the text
f4.7_txt <- filter(f4.7_data, entity=="Global" & year >= 2018) %>%
  select(-hiv_tpt_new, -entity) %>%
  pivot_wider(names_from = year,
              names_prefix = "hiv_tpt_",
              values_from = hiv_tpt_all) %>%

  mutate(pct_21_20 = abs(hiv_tpt_2021 - hiv_tpt_2020) * 100 / hiv_tpt_2020,
         hiv_tpt_18_21 = hiv_tpt_2018 + hiv_tpt_2019 + hiv_tpt_2020  + hiv_tpt_2021)

# Add cumulative total since 2005
f4.7_txt <- filter(f4.7_data, entity=="Global" & year %in% 2005:2017) %>%
  summarise(hiv_tpt_05_17 = sum(hiv_tpt_new)) %>%
  cbind(f4.7_txt) %>%
  mutate(hiv_tpt_05_21 = hiv_tpt_05_17 + hiv_tpt_18_21) %>%
  mutate(hiv_tpt_05_21_pct = hiv_tpt_05_21 * 100 / 37.7e6)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Create dataframe for figure 4.8  ----
# (Tree map showing the provision of TB preventive treatment highlighting countries reported >= 200 000 people)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

f4.8_data_all <- filter(notification, year == report_year - 1) %>%
  select(iso2,
         country,
         hiv_ipt_reg_all,
         hiv_ipt,
         # These next ones introduced dcyear 2021 by GAM
         hiv_all_tpt,
         hiv_new_tpt,
         hiv_elig_all_tpt,
         hiv_elig_new_tpt) %>%

  # Use the first non-empty tpt/ipt variable (not sure this is the best approach as it has become so
  # darn complicated in 2021 ...)
  mutate(hiv_tpt = coalesce(hiv_ipt_reg_all, hiv_ipt, hiv_all_tpt, hiv_elig_all_tpt, hiv_new_tpt, hiv_elig_new_tpt)) %>%
  select(iso2, country, hiv_tpt) %>%

  # Calculate the proportion of global total for each country
  mutate(proportion = hiv_tpt * 100 / sum(hiv_tpt, na.rm = TRUE)) %>%

  # Sort  in descending order of TPT provision
  arrange(desc(proportion))

# Now cut data to size for the treemap
f4.8_data <- filter(f4.8_data_all, hiv_tpt < 2e5) %>%
  summarise(hiv_tpt=sum(hiv_tpt)) %>%
  mutate(entity = "(Rest of the world)")

f4.8_data <- filter(f4.8_data_all, hiv_tpt >= 2e5) %>%
  select(hiv_tpt, entity = country) %>%
  rbind(f4.8_data)



# Summary for quoting in the text: Capture the global total for the section text before filtering for the countries reporting
# at least 200k
f4.8_txt <- f4.8_data_all %>%
  summarise(hiv_tpt_glob = sum(hiv_tpt, na.rm=TRUE))

# Add proportion accounted for by the countries reporting >= 200k  to the summary for quoting in the text
f4.8_txt <- f4.8_data_all %>%
  filter(hiv_tpt >= 2e5) %>%
  summarise(prop_top = sum(proportion)) %>%
  cbind(f4.8_txt)

# Make list of the top  countries sorted alphabetically
f4.8_txt_list <- f4.8_data_all %>%
  filter(hiv_tpt >= 2e5) %>%
  select(country) %>%
  arrange(country)

rm(f4.8_data_all)

# Add info about coverage of tpt among newly enrolled in 2021
f4.8_txt_coverage <- filter(notification, year == report_year - 1 & !is.na(hiv_new_tpt) & !is.na(hiv_new) ) %>%
  inner_join(hbtbhiv30, by = "iso3") %>%
  select(iso2,
         hiv_new_tpt,
         hiv_new) %>%
  mutate(coverage = hiv_new_tpt * 100 / hiv_new) %>%
  summarise(countries = n(),
            median = median(coverage),
            q1 = unname(quantile(coverage, 0.25)),
            q3 = unname(quantile(coverage, 0.75)))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Create dataframe for figure 4.9  ----
# (Panel plot showing percentage completion vs number PLHIV started TPT by WHO region)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

f4.9_data <- filter(notification, year == report_year - 2) %>%

  select(country,
         g_whoregion,
         hiv_all_tpt_started,
         hiv_all_tpt_completed) %>%

  # Calculate completion rate
  mutate(pct_completed = ifelse(!is.na(hiv_all_tpt_completed) & NZ(hiv_all_tpt_started) > 0,
                                hiv_all_tpt_completed * 100 /hiv_all_tpt_started ,
                                NA)) %>%

  # merge with regional names
  inner_join(who_region_shortnames, by = "g_whoregion") %>%
  select(-g_whoregion) %>%

  # filter out empty lines
  filter(!is.na(pct_completed))

f4.9_txt <- f4.9_data %>%
  summarise(countries = n(),
            median = median(pct_completed),
            q1 = unname(quantile(pct_completed, 0.25)),
            q3 = unname(quantile(pct_completed, 0.75)))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Create dataframe for figure 4.10 ----
# (Map showing use of diagnostic tests for TB infection)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


f4.10_data <- strategy %>%
  select(iso3, g_whoregion, year, infection_tests_used) %>%
  filter(year == report_year -1) %>%
  mutate(var = factor(infection_tests_used,
                      levels = c(177, 175, 176, 194, 3),
                      labels = c("IGRA and TST",
                                 "Interferon Gamma Release Assays (IGRA)",
                                 "Tuberculin Skin Test (TST)",
                                 "Not used",
                                 "Don’t know")))



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Create dataframe for figure 4.11 ----
# (Map showing ratio of TB notification rates among health care workers to those among the adult population)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

f4.11_data <- filter(strategy, year == report_year - 1) %>%

  select(iso3,
         country,
         hcw_tb_infected,
         hcw_tot)

# Get the total adult population aged between 15 and 64 (exclude those aged 65 and above)
pop_adults <- filter(estimates_population, year == report_year - 1) %>%

  select(iso3,
         e_pop_f15plus,
         e_pop_f65,
         e_pop_m15plus,
         e_pop_m65) %>%

  mutate(e_pop_adult = e_pop_f15plus + e_pop_m15plus - e_pop_f65 -  e_pop_m65 ) %>%

  select(iso3,
         e_pop_adult)

# Get the total notifications among adults aged between 15 and 64 (exclude those aged 65 and above)
notif_adults <- filter(notification, year == report_year - 1) %>%

  select(iso3,
         newrel_f15plus,
         newrel_f65,
         newrel_m15plus,
         newrel_m65) %>%

  mutate(newrel_adult = newrel_f15plus + newrel_m15plus - NZ(newrel_f65) -  NZ(newrel_m65) ) %>%

  select(iso3,
         newrel_adult) %>%

  # Join to the adult population
  inner_join(pop_adults, by = "iso3")

f4.11_data <- f4.11_data %>%

  inner_join(notif_adults, by = "iso3") %>%

  # Calculate notification rate ratio
  # Use as.numeric() to avoid integer overflow
  mutate(nrr = ifelse(NZ(hcw_tot) > 0 & NZ(newrel_adult) > 0,
                      (as.numeric(hcw_tb_infected) * as.numeric(e_pop_adult))
                      /
                        (as.numeric(hcw_tot) * as.numeric(newrel_adult)),
                      NA)) %>%

  # in previous years I had filtered out countries with fewer than 100 health care workers
  # as the rate ratios jumped around a lot but because these are very small countries they
  # don;t show up in the maps so won't bother anymore

  # Assign the categories for the map
  mutate(var = cut(nrr,
                   c(0, 1, 2, 3, Inf),
                   c('0-0.9', '1-1.9', '2-2.9', '\u22653'),
                   right=FALSE))

# Summary for the text
f4.11_txt <- filter(f4.11_data, hcw_tb_infected>0) %>%
  summarise(tot_hcw_tb = sum(hcw_tb_infected, na.rm=TRUE),
            countries_hcw_tb = n())

# Add number with nrr more than one when number of TB among hcw is 5 or more
f4.11_txt <- filter(f4.11_data, nrr > 1 & hcw_tb_infected >= 5) %>%
  summarise(countries_nrr = n()) %>%
  cbind(f4.11_txt)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Create dataframe for figure 4.12 ----
# (Map showing BCG immunisation policies)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Suppress this part in July 2022 for later resurrection

if(show_bcg){

  # # Downloaded the BCG indicator for 2019-2021 from the GHO on 2022-08-08 using this code:
  # # -------- start of BCG download code -----
  # # indicator code WHS4_543 is the BCG immunisation coverage among 1-year-olds
  # # https://www.who.int/data/gho/data/indicators/indicator-details/GHO/bcg-immunization-coverage-among-1-year-olds-(-)
  # bcg_gho <- ghost::gho_data("WHS4_543") %>%
  #   filter(TimeDim >= 2015 & SpatialDimType=="COUNTRY") %>%
  #   mutate(bcg_coverage = as.integer(Value)) %>%
  #   select(iso3 = SpatialDim,
  #          year = TimeDim,
  #          bcg_coverage)
  #
  # # Save to the CSV folder. Note that GHO web page says data were last updated on 2022-07-15
  # write.csv(bcg_gho, here(paste0("csv/bcg_gho_", Sys.Date(), '.csv')), row.names = FALSE)
  # #
  # # Also downloaded aggregates (global and by WHO region)
  #
  #   bcg_gho_agg <- ghost::gho_data("WHS4_543") %>%
  #     filter(TimeDim >= 2015 & SpatialDimType=="REGION") %>%
  #     mutate(bcg_coverage = as.integer(Value),
  #            group_type = ifelse(SpatialDim == "GLOBAL", "global", "g_whoregion"),
  #            group_name = ifelse(SpatialDim == "GLOBAL", "global", substr(SpatialDim, 1, 3)) ) %>%
  #     select(group_type,
  #            group_name,
  #            year = TimeDim,
  #            bcg_coverage)
  #
  #     # # Save to the CSV folder. Note that GHO web page says data were last updated on 2022-07-15
  #     write.csv(bcg_gho_agg, here(paste0("csv/bcg_gho_agg_", Sys.Date(), '.csv')), row.names = FALSE)
  #
  # # -------- end of BCG download code -----


  # Third option is to look at aggregate values
  f4.12_data <- read.csv(here('csv/bcg_gho_agg_2022-08-10.csv'))

  # merge with regional names
  f4.12_reg <- f4.12_data %>%
    inner_join(who_region_shortnames, by = c("group_name" = "g_whoregion")) %>%
    select(entity,
           year,
           bcg_coverage)

  # Or to avoid factor problems
  f4.12_glob <- f4.12_data %>%
    filter(group_name == "global") %>%
    mutate(entity = "Global") %>%
    select(entity,
           year,
           bcg_coverage)

  # Combine
  f4.12_data <- rbind(f4.12_reg, f4.12_glob)

  # Change the entity order
  f4.12_data$entity <- factor(f4.12_data$entity,
                              levels = c("Global", "African Region", "Region of the Americas", "South-East Asia Region",
                                         "European Region", "Eastern Mediterranean Region", "Western Pacific Region"))

  f4.12_txt <- f4.12_glob %>%
    filter(year %in% c(2019, 2021)) %>%
    pivot_wider(names_from = year,
                names_prefix = "bcg_cov_",
                values_from = bcg_coverage)


  # Dennis wants to show number of reporting countries

  # Load the country names and the BCG by country files
  member_states <- get_timestamped_csv('cty', csv_datestamp) %>%
    filter(g_whostatus == "M") %>%
    select(iso3, g_whostatus, g_whoregion)

  f4.12_txt_reps <- read.csv(here('csv/bcg_gho_2022-08-10.csv')) %>%
    filter(year==2021) %>%
    right_join(member_states, by="iso3") %>%
    mutate(reported = ifelse(is.na(bcg_coverage), 0, 1)) %>%
    select(g_whoregion, iso3, reported)

  f4.12_txt_rep_glob <- f4.12_txt_reps %>%
    summarise(rep = sum(reported),
              all =  n())

  f4.12_txt_rep_reg <- f4.12_txt_reps %>%
    group_by(g_whoregion) %>%
    summarise(rep = sum(reported),
              all =  n()) %>%
    ungroup()

  # Clear up
  rm(f4.12_reg, f4.12_glob)

}


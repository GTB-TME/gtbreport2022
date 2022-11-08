# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data preparation script for ch3-1.rmd
# Hazim Timimi, July 2022
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


# Load chapter 3 and 4 packages, settings and data
source(here::here('report/ch3_ch4_load_data.r'))



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: table 3.1.1 ----
# (Notifications of TB, HIV-positive TB, and DR-TB cases, globally and for WHO regions)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

t3.1.1_data <- filter(notification, year == (report_year - 1)) %>%
  select(g_whoregion,
         c_notified,
         c_newinc,
         new_labconf, new_clindx, new_ep,
         ret_rel_labconf, ret_rel_clindx, ret_rel_ep,
         newrel_hivpos,
         conf_rr_nfqr,
         conf_rr_fqr) %>%

  # calculate regional aggregates
  group_by(g_whoregion) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%

  # merge with regional names
  inner_join(who_region_shortnames, by = "g_whoregion") %>%
  arrange(entity) %>%
  select(-g_whoregion)


# Add global summary to the regional summary
t3.1.1_global <- t3.1.1_data %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  mutate(entity="Global")


t3.1.1_data <- rbind(t3.1.1_data, t3.1.1_global)


# Calculate total pulmonary and %ages that are bac confirmed and that are extrapulmonary
t3.1.1_data <- t3.1.1_data %>%
  mutate( newrel_pulm = new_labconf + new_clindx + ret_rel_labconf + ret_rel_clindx,
          newrel_pulm_conf_pct = (new_labconf + ret_rel_labconf) * 100
          /
            (new_labconf + new_clindx + ret_rel_labconf + ret_rel_clindx),
          newrel_ep_pct = (new_ep + ret_rel_ep) * 100
          /
            (c_newinc)
  ) %>%
  # Restrict to variables needed in the final output
  select(entity,
         c_notified,
         c_newinc,
         newrel_pulm,
         newrel_pulm_conf_pct,
         newrel_ep_pct,
         newrel_hivpos,
         conf_rr_nfqr,
         conf_rr_fqr)

# summary dataset for the text
t3.1.1_region <- filter(t3.1.1_data, entity!="Global") %>%
  arrange(desc(c_notified)) %>%
  mutate(c_total_p = c_notified/t3.1.1_global$c_notified*100) %>%
  mutate(c_newinc_cum_p = cumsum(c_total_p))

t3.1.1_txt <- filter(t3.1.1_data, entity=="Global") %>%
  mutate(c_pulm_p = newrel_pulm/c_newinc *100)

t3.1.1_txt <- filter(t3.1.1_region, entity=="Western Pacific Region") %>%
  select(c_newinc_cum_p = c_newinc_cum_p) %>%
  cbind(t3.1.1_txt)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.1.1 ----
# (Panel plot of incidence estimates compared to notifications by WHO region and globally since 2000)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if(show_estimates){

  inc_global <- filter(est_global, year>=2000) %>%
    mutate(entity = 'Global') %>%
    select(year,
           entity,
           e_inc_num = inc.num,
           e_inc_num_lo = inc.lo.num,
           e_inc_num_hi = inc.hi.num,
           c_newinc = c.newinc)

  inc_regional <- filter(est_regional, year>=2000) %>%
    select(year,
           g_whoregion = g.whoregion,
           e_inc_num = inc.num,
           e_inc_num_lo = inc.lo.num,
           e_inc_num_hi = inc.hi.num,
           c_newinc = c.newinc) %>%

    # merge with regional names
    inner_join(who_region_shortnames, by = "g_whoregion") %>%
    select(-g_whoregion)

  # Combine the two sets
  f3.1.1_data <- rbind(inc_regional, inc_global) %>%

    # Set the entity order for plotting
    mutate(entity = factor(entity,
                           levels = c("Global", "African Region", "Region of the Americas", "South-East Asia Region",
                                      "European Region", "Eastern Mediterranean Region", "Western Pacific Region")))

  # Summary dataset for simple quoting of numbers in the text of section 3.1
  f3.1.1_past3y <- filter(inc_global, year>=report_year-3) %>%
    mutate(c_newinc_p   = lag(c_newinc)) %>%
    mutate(c_newinc_p2 = lag(lag(c_newinc))) %>%
    mutate(pct_dif      = (c_newinc - c_newinc_p)*100/c_newinc_p) %>%
    mutate(pct_dif2    = (c_newinc - c_newinc_p2)*100/c_newinc_p2)

  f3.1.1_txt <- filter(f3.1.1_past3y, year==report_year-3) %>%
    select(c_newinc_y3 = c_newinc)

  f3.1.1_txt <- filter(f3.1.1_past3y, year==report_year-2) %>%
    select(c_newinc_y2 = c_newinc,
           pct_dif) %>%
    cbind(f3.1.1_txt)

  f3.1.1_txt <- filter(f3.1.1_past3y, year==report_year-1) %>%
    select(c_newinc_y1 = c_newinc,
           pct_dif2) %>%
    cbind(f3.1.1_txt)

}



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.1.2 ----
# (Panel plot of incidence estimates compared to notifications for 30 countries since 2000)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if(show_estimates){

  f3.1.2_data <- filter(est_country, iso3 %in% hbc30$iso3)  %>%
    select(year,
           iso3,
           e_inc_num = inc.num,
           e_inc_num_lo = inc.lo.num,
           e_inc_num_hi = inc.hi.num) %>%

    inner_join(select(notification, iso3, year, country, c_newinc), by = c("iso3", "year")) %>%
    arrange(country)


  # Summary dataset for simple quoting of numbers in the text of section 3.1
  f3.1.2_txt <- filter(f3.1.2_data, year>=report_year-2 & iso3 %in% c("IND", "IDN", "PHL")) %>%
    mutate(c_newinc_p = lag(c_newinc)) %>%
    mutate(pct_dif = (c_newinc - c_newinc_p)*100/c_newinc_p) %>%
    filter(year==report_year-1) %>%
    select(iso3, pct_dif) %>%
    pivot_wider(names_from = iso3,
                values_from = pct_dif)

  # Calculate % change for India and Indonesia between 2013 and 2019
  # and add it to the summary
  f3.1.2_txt <- filter(f3.1.2_data, year %in% c(2013, 2019) & iso3 %in% c("IND", "IDN")) %>%
    select(year, c_newinc) %>%
    mutate(c_newinc_p = lag(c_newinc)) %>%
    group_by(year) %>%
    summarise(across(starts_with("c_"), sum, na.rm=TRUE)) %>%
    mutate(diff = c_newinc - c_newinc_p) %>%
    filter(year==2019) %>%
    select(IND_IDN_2013 = diff) %>%
    cbind(f3.1.2_txt)

}




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.1.3 ----
# (Irwin's doughnuts -- % completion of UNHLM targets)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Summary dataset for simple quoting of numbers in the text of section 3.1
f3.1.3_txt <- filter(notification, year >= 2018) %>%
  summarise(across(c(c_newinc,
                     c_new_014,
                     conf_rrmdr_tx,
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
         -conf_rr_fqr_tx)

f3.1.3_data <- f3.1.3_txt %>%

  # Calculate percentage complete for each UNHLM 2018-2022 target
  mutate(c_newinc_pct = c_newinc * 100/ 40e6,  # target: 40 million notified
         c_new_014_pct = c_new_014 * 100 / 3.5e6, # target: 3.5 million children notified
         rr_treated_pct  = rr_treated * 100 / 1.5e6,  # target: 1.5 million treated for drug-resistant TB
         rrmdr_014_tx_pct = rrmdr_014_tx * 100 / 115e3 # target: 115 thousand children treated for drug-resistant TB
         ) %>%

  select(contains("_pct"))  %>%

  pivot_longer(cols = contains("_pct"),
               names_to = "target_completion")

# Supplementary data for quoting in the text
f3.1.3_kids_txt <- filter(notification, year >= 2019) %>%
  group_by(year) %>%
  summarise(sum_kids = sum(c_new_014, na.rm = TRUE)) %>%
  ungroup %>%
  pivot_wider(names_from = year,
              names_prefix = "kids_",
              values_from = sum_kids) %>%
  mutate(kids_change_pct = abs(kids_2020 - kids_2019) * 100 / kids_2019)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.1.4 ----
# (Bar chart showing numbers of adults and children notified with TB each year since 2015)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

f3.1.4_data <- filter(notification, year >= 2015) %>%

  select(iso2, year, c_newinc, c_new_014) %>%

  group_by(year) %>%
  summarise(across(starts_with("c_new"), sum, na.rm = TRUE)) %>%

  # calculate the "adult" fraction
  mutate(c_new_15plus = c_newinc - c_new_014) %>%

  # switch to long format for plotting
  pivot_longer(cols = starts_with("c_new_"),
               names_to = "age_group",
               values_to = "how_many")

f3.1.4_txt <- filter(notification, year == report_year-1) %>%

  select(iso2, year, c_newinc, c_new_014, newrel_f15plus, newrel_m15plus) %>%

  group_by(year) %>%
  summarise(across(contains("new"), sum, na.rm = TRUE)) %>%

  # calculate the "adult" fraction
  mutate(c_new_15plus = c_newinc - c_new_014) %>%

  #calculate pct
  mutate(pct_m = newrel_m15plus/c_newinc * 100,
         pct_f = newrel_f15plus/c_newinc * 100,
         pct_c = c_new_014/c_newinc * 100)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.1.5 ----
# (Panel plot of incidence rate estimates and notification rates by age group and sex as population pyramid-style bar charts)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if(show_estimates){

  agesex_notifs <- filter(notification, year == (report_year - 1)) %>%

    select(iso2,
           g_whoregion,
           c_newinc,
           newrel_m04, newrel_m514, newrel_m1524, newrel_m2534, newrel_m3544, newrel_m4554, newrel_m5564, newrel_m65,
           newrel_f04, newrel_f514, newrel_f1524, newrel_f2534, newrel_f3544, newrel_f4554, newrel_f5564, newrel_f65) %>%

    # Calculate total notifications in the 10-year age groups for ages 15 and over for each country
    rowwise() %>%
    mutate(tot_adults = sum(c_across(newrel_m1524:newrel_m65), na.rm = TRUE) +
             sum(c_across(newrel_f1524:newrel_f65), na.rm = TRUE)
    ) %>%
    mutate(tot_children = sum(c_across(newrel_m04:newrel_m514), na.rm = TRUE) + sum(c_across(newrel_f04:newrel_f514), na.rm = TRUE)) %>%
    ungroup()

  # remove countries that had not reported cases in the 10-year age groups for ages 15 and over
  # But keep a record of them for a footnote

  f3.1.5_excluded <- filter(agesex_notifs, tot_adults == 0) %>%
    summarise(countries=n())

  f3.1.5_excluded <- filter(agesex_notifs, tot_adults > 0) %>%

    summarise(across(c_newinc, sum, na.rm=TRUE)) %>%

    mutate(agesex_pct = c_newinc * 100 /
             sum(notification[notification$year==report_year -1, "c_newinc"], na.rm=TRUE)) %>%
    select(agesex_pct) %>%
    cbind(f3.1.5_excluded)


  agesex_notifs <- filter(agesex_notifs, tot_adults > 0)

  # Get population by age and sex
  agesex_pop    <- filter(estimates_population, year == (report_year - 1)) %>%
    select(iso2,
           g_whoregion,
           e_pop_m04, e_pop_m514, e_pop_m1524, e_pop_m2534, e_pop_m3544, e_pop_m4554, e_pop_m5564, e_pop_m65,
           e_pop_f04, e_pop_f514, e_pop_f1524, e_pop_f2534, e_pop_f3544, e_pop_f4554, e_pop_f5564, e_pop_f65)

  # Join the two tables
  agesex <- inner_join(agesex_notifs, agesex_pop, by = c("iso2", "g_whoregion"))

  # Calculate global and regional aggregates and combine them
  agesex_global <- agesex %>%
    summarise(across(newrel_m04:e_pop_f65, sum, na.rm = TRUE)) %>%
    mutate(group_name = "global")

  agesex_regional <- agesex %>%
    group_by(g_whoregion) %>%
    summarise(across(newrel_m04:e_pop_f65, sum, na.rm = TRUE)) %>%
    rename(group_name = g_whoregion)

  # Put the aggregates together
  agesex_agg_notifs_wide <- rbind(agesex_global, agesex_regional)
  rm(agesex_global, agesex_regional, agesex_notifs, agesex)

  # Calculate case notification rates
  agesex_agg_notifs_wide <- agesex_agg_notifs_wide %>%
    mutate(m04 = newrel_m04 * 1e5 / e_pop_m04,
           m514 = newrel_m514 * 1e5 / e_pop_m514,
           m1524 = newrel_m1524 * 1e5 / e_pop_m1524,
           m2534 = newrel_m2534 * 1e5 / e_pop_m2534,
           m3544 = newrel_m3544 * 1e5 / e_pop_m3544,
           m4554 = newrel_m4554 * 1e5 / e_pop_m4554,
           m5564 = newrel_m5564 * 1e5 / e_pop_m5564,
           m65 = newrel_m65 * 1e5 / e_pop_m65,

           f04 = newrel_f04 * 1e5 / e_pop_f04,
           f514 = newrel_f514 * 1e5 / e_pop_f514,
           f1524 = newrel_f1524 * 1e5 / e_pop_f1524,
           f2534 = newrel_f2534 * 1e5 / e_pop_f2534,
           f3544 = newrel_f3544 * 1e5 / e_pop_f3544,
           f4554 = newrel_f4554 * 1e5 / e_pop_f4554,
           f5564 = newrel_f5564 * 1e5 / e_pop_f5564,
           f65 = newrel_f65 * 1e5 / e_pop_f65)

  # Flip into long format
  agesex_agg_notifs <- agesex_agg_notifs_wide %>%

    select(group_name:f65) %>%

    pivot_longer(cols = m04:f65,
                 names_to = c("sex", "age_group"),
                 names_pattern = "(.)(.*)",
                 values_to = "cnr")

  # Get the aggregated incidence estimates by age group and sex (already in long format)
  agesex_agg_inc <- filter(db_estimates_group,
                           year == (report_year - 1) &
                             measure == "inc" &
                             unit == "num" &
                             sex != "a" &
                             age_group %in% c("0-4","5-14", "15-24", "25-34", "35-44", "45-54", "55-64", "65plus") ) %>%

    select(group_name,
           age_group,
           sex,
           inc_num = best) %>%

    # Match format of age_group with that for notifications
    mutate(age_group = str_remove(age_group,"-")) %>%
    mutate(age_group = ifelse(age_group=="65plus", "65", age_group))

  # Calculate aggregate populations so can convert incidence numbers to rates
  agesex_agg_pop <- agesex_pop %>%
    summarise(across(e_pop_m04:e_pop_f65, sum, na.rm = TRUE)) %>%
    mutate(group_name = "global")

  agesex_agg_pop <- agesex_pop %>%
    group_by(g_whoregion) %>%
    summarise(across(e_pop_m04:e_pop_f65, sum, na.rm = TRUE)) %>%
    rename(group_name = g_whoregion) %>%
    rbind(agesex_agg_pop) %>%

    # flip to long format
    pivot_longer(cols = e_pop_m04:e_pop_f65,
                 names_to = c("sex", "age_group"),
                 names_pattern = "e_pop_(.)(.*)",
                 values_to = "pop")


  # merge aggregate incidence and population and convert numbers to rates
  agesex_agg_inc <- agesex_agg_inc %>%
    inner_join(agesex_agg_pop, by = c("group_name", "age_group", "sex")) %>%

    # calculate the incidence rate for each agesex combination
    mutate(inc_100k = inc_num * 1e5 / pop) %>%

    # keep the incidence rate only
    select(-inc_num, -pop)

  # Merge incidence and notification records
  f3.1.5_data <- agesex_agg_notifs %>%
    inner_join(agesex_agg_inc, by=c("group_name", "sex", "age_group")) %>%

    # Define factors, labels and sorting ready for plotting
    mutate(age_group = factor(age_group,
                              levels=c("04","514", "1524", "2534", "3544", "4554", "5564", "65"),
                              labels=c("0-4","5-14", "15-24", "25-34", "35-44", "45-54", "55-64", "\u226565"))) %>%

    mutate(sex = factor(sex,
                        levels = c("f", "m"),
                        labels=c("Female", "Male"))) %>%

    # merge with regional names
    left_join(who_region_shortnames, by = c("group_name" = "g_whoregion")) %>%
    mutate(entity = ifelse(group_name=="global", "Global", as.character(entity))) %>%

    # get rid of unneeded variables
    select(-group_name) %>%

    # Change the order of the entities
    mutate(entity = factor(entity,
                           levels = c("Global", "African Region", "Region of the Americas", "South-East Asia Region",
                                      "European Region", "Eastern Mediterranean Region", "Western Pacific Region")))



  # Create summary stats for the text
  f3.1.5_txt <- agesex_agg_notifs_wide %>%
    rowwise() %>%
    mutate(men = sum(across(starts_with("newrel_m"))) - newrel_m04 - newrel_m514,
           women = sum(across(starts_with("newrel_f"))) - newrel_f04 - newrel_f514,
           all_disagg = sum(across(starts_with("newrel_"))),
           mf_ratio = sum(across(starts_with("newrel_m")))/sum(across(starts_with("newrel_f")))) %>%

    mutate(pct_men = round(men*100/all_disagg),
           pct_women = round(women*100/all_disagg)) %>%

    # do the usual trick to avoid % not adding up to 100
    mutate(pct_kids = 100 - pct_men - pct_women) %>%

    select(entity = group_name,
           pct_men,
           pct_women,
           pct_kids,
           mf_ratio) %>%
    filter(entity == "global")

}





# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.1.6 ----
# (Map showing percentage of new and relapse TB cases that were children)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

kids_data <- filter(notification, year>=report_year-2) %>%

  select(iso3,
         country,
         year,
         c_new_014,
         newrel_m15plus,
         newrel_mu,
         newrel_f15plus,
         newrel_sexunk15plus,
         newrel_fu) %>%

  # calculate % of children in the age/sex data
  rowwise() %>%
  mutate(agesex_tot = sum(c_across(c_new_014:newrel_fu), na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(kids_pct = ifelse(agesex_tot > 0,
                           c_new_014 * 100 / agesex_tot,
                           NA)) %>%

  # Assign the categories for the map
  mutate(var = cut(kids_pct,
                   c(0, 5.0, 10.0, 15.0, Inf),
                   c('0-4.9', '5-9.9', '10-14.9', '\u226515'),
                   right=FALSE))

# Find the countries with empty data for latest year and see if there are data for the previous year
kids_prev_year_data <- kids_data %>%
  filter(year == report_year - 1 & is.na(kids_pct)) %>%
  select(iso3) %>%
  inner_join(filter(kids_data, year == report_year - 2), by = "iso3") %>%
  filter(!is.na(kids_pct))

# Now combine into one dataframe, with previous data used if latest year's data are not available
f3.1.6_data <- kids_data %>%
  filter(year == report_year - 1) %>%
  anti_join(kids_prev_year_data, by= "iso3") %>%
  rbind(kids_prev_year_data)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.1.7 ----
# (Map showing percentage of extrapulmonary cases among new and relapse TB cases)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

ep_data <- notification %>%
  filter(year  >= report_year - 2) %>%
  select(iso3,
         country,
         year,
         new_labconf, new_clindx, new_ep,
         ret_rel_labconf, ret_rel_clindx, ret_rel_ep) %>%

  # calculate % of extrapulmonary cases
  rowwise() %>%
  mutate(newrel_tot = sum(c_across(new_labconf:ret_rel_ep), na.rm = TRUE)) %>%
  mutate(ep_tot = sum(c_across(contains("_ep")), na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(ep_pct = ifelse(newrel_tot > 0,
                         ep_tot * 100 / newrel_tot,
                         NA)) %>%

  # Assign the categories for the map
  mutate(var = cut(ep_pct,
                   c(0, 10, 20, 30, Inf),
                   c('0-9.9', '10-19', '20-29', '\u226530'),
                   right=FALSE))

# Find the countries with empty data for latest year and see if there are data for the previous year
ep_prev_year_data <- ep_data %>%
  filter(year == report_year - 1 & is.na(ep_pct)) %>%
  select(iso3) %>%
  inner_join(filter(ep_data, year == report_year - 2), by = "iso3") %>%
  filter(!is.na(ep_pct))

# Now combine into one dataframe, with previous data used if latest year's data are not available
f3.1.7_data <- ep_data %>%
  filter(year == report_year - 1) %>%
  anti_join(ep_prev_year_data, by= "iso3") %>%
  rbind(ep_prev_year_data)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.1.8 ----
# (Panel plot of line charts showing percentage contribution of the private sector to TB notifications by year since 2010)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Get PPM data for 7 priority countries
f3.1.8_data <- filter(strategy, year >= 2010 & iso2 %in% c('BD', 'IN', 'ID', 'MM', 'NG', 'PK', 'PH')) %>%
  select(iso2, year, country, priv_new_dx) %>%

  # Merge with notifications
  inner_join(select(notification, iso2, year, c_notified), by = c("iso2", "year")) %>%

  # Calculate percent contributions
  mutate(private_pcnt = ifelse(c_notified > 0,
                               priv_new_dx * 100 / c_notified,
                               NA))





# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.1.9 ----
# (Map showing percentage of management units with community contributions to case finding and/or treatment support)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Get the variable to identify countries requested to report on community indicators
comm_datarequest <- data_collection %>%
  filter(datcol_year==report_year) %>%
  select(iso3,
         dc_engage_community_display)

f3.1.9_data <- strategy %>%
  filter(year==report_year - 1) %>%
  select(iso3,
         country,
         bmu,
         bmu_community_impl,
         community_data_available)%>%
  mutate(comm_pct = ifelse(bmu > 0,
                           bmu_community_impl * 100 / bmu,
                           NA))  %>%

  # Assign the categories for the map
  mutate(var = cut(comm_pct,
                   c(0, 25, 50, 75, Inf),
                   c('0-24', '25-49', '50-74', '\u226575'),
                   right=FALSE)) %>%

  # Link to data request status (used in the footnote for the figure)
  inner_join(comm_datarequest, by = "iso3")





# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.1.10 ----
# (Map showing which countries have case-based TB surveillance systems)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

cb_data <- strategy %>%
  filter(year >= report_year - 2) %>%
  select(iso3,
         country,
         year,
         caseb_err_nat)

# Make sure all ECDC countries are marked as having case-based surveillance for all TB cases
cb_ecdc <- data_collection %>%
  filter(datcol_year == report_year & dc_ecdc == 1) %>%
  select(iso3)

cb_data <- cb_data %>%
  mutate(caseb_err_nat = ifelse(iso3 %in% cb_ecdc$iso3, 42, caseb_err_nat)) %>%

  # UK hasn't responded, but we know it has a case-based system, so fudge it
  mutate(caseb_err_nat = ifelse(iso3=="GBR", 42, caseb_err_nat)) %>%

  # Assign the categories for the map
  mutate(var = factor(caseb_err_nat,
                      levels = c(0, 44, 43, 42),
                      labels = c("None", "Partially (in transition)", "MDR-TB patients only", "All TB patients")))

# Find the countries with empty data for latest year and see if there are data for the previous year
cb_prev_year_data <- cb_data %>%
  filter(year == report_year - 1 & is.na(caseb_err_nat)) %>%
  select(iso3) %>%
  inner_join(filter(cb_data, year == report_year - 2), by = "iso3") %>%
  filter(!is.na(caseb_err_nat))

# Now combine into one dataframe, with previous data used if latest year's data are not available
f3.1.10_data <- cb_data %>%
  filter(year == report_year - 1) %>%
  anti_join(cb_prev_year_data, by= "iso3") %>%
  rbind(cb_prev_year_data)


# Simple summary for the section text
f3.1.10_txt <- filter(f3.1.10_data, caseb_err_nat == 42) %>%
  select(iso3) %>%
  inner_join(filter(notification, year==report_year-1), by = "iso3") %>%
  summarise(across(c_newinc, sum, na.rm=TRUE))

f3.1.10_txt <- filter(notification, year==report_year-1) %>%
  summarise(across(c_newinc, sum, na.rm=TRUE)) %>%
  select(c_newinc_glob = c_newinc) %>%
  cbind(f3.1.10_txt) %>%
  mutate(pct_notif_caseb = c_newinc*100/c_newinc_glob) %>%
  select(pct_notif_caseb)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.1.additional ----
# (Map showing % of foreign borne cases)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

f3.1.add_data <- notification %>%
  select(iso3,country,year,c_notified,notif_foreign) %>%
  filter(year==report_year-1|year==report_year-2) %>%
  mutate(pct_foreign = notif_foreign/c_notified * 100) %>%
  pivot_wider(names_from = year, values_from = c_notified:pct_foreign) %>%
  mutate(pct_foreign = ifelse(is.na(pct_foreign_2021),pct_foreign_2020,pct_foreign_2021)) %>%
# Assign the categories for the map
  mutate(var = cut(pct_foreign,
                   c(0, 5, 25, 50, 75, Inf),
                   c('0\u20134','5\u201324','25\u201349','50\u201374','\u226575'),
                   right=FALSE)) 
  
#   whomap(colours = brewer.pal(5, "YlOrRd"),
#          legend.title = "Percentage (%)",
#          na.col = "white",
#          water.col = "white")

## title: Percentage of people notified with TB (new, relapse and retreatment) who were born outside the country of report, 2021
## texts: At times the age-group and sex distribution of people with TB varies between countries because of different demographic characteristics between native and non-native populations. In most industrialised countries, the majority of TB notifications occur in people born outside the country of report (Fig. 3.1.7)
# ggsave(f3.1.add_plot,file=paste0("./report/ch3_data/f3.1.additional_",as.character(Sys.Date()),".pdf"),width=14, height = 6)
# f3.1.add_data %>% write_csv(here("./report/ch3_data/fig3.1.additional.csv"))

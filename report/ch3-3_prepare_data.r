# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data preparation script for ch3-3.rmd
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
# Data: figure 3.3.1 ----
# (Forest plot of TB treatment coverage in 30 countries, regionally and globally)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if(show_estimates){

  # A. Countries
  # - - - - - - - -
  coverage_inc_country <- filter(est_country, year==report_year-1) %>%
    select(year,
           iso3,
           e_inc_num = inc.num,
           e_inc_num_lo = inc.lo.num,
           e_inc_num_hi = inc.hi.num) %>%

    # restrict to high burden countries
    inner_join(hbc30, by = "iso3")

  coverage_country <- filter(notification, year==report_year-1) %>%
    select(entity = country,
           iso3,
           c_newinc)  %>%
    inner_join(coverage_inc_country, by = "iso3") %>%
    select(-iso3) %>%
    mutate(c_cdr = c_newinc * 100 / e_inc_num,
           c_cdr_lo = c_newinc * 100  / e_inc_num_hi,
           c_cdr_hi = c_newinc * 100  / e_inc_num_lo,
           # highlight countries with no data
           entity = ifelse(is.na(c_newinc), paste0(entity, "*"), entity )) %>%
    select(entity,
           c_cdr,
           c_cdr_lo,
           c_cdr_hi) %>%
    arrange(desc(c_cdr))


  # B. Regions
  # - - - - - - - -
  coverage_inc_region <- filter(est_regional, year==report_year-1) %>%
    select(g_whoregion = g.whoregion,
           e_inc_num = inc.num,
           e_inc_num_lo = inc.lo.num,
           e_inc_num_hi = inc.hi.num)

  coverage_region <- filter(notification, year==report_year-1) %>%
    select(g_whoregion,
           c_newinc) %>%
    # calculate regional aggregates
    group_by(g_whoregion) %>%
    summarise(across(c_newinc, sum, na.rm=TRUE)) %>%
    ungroup() %>%

    # merge with incidence estimates
    inner_join(coverage_inc_region, by = "g_whoregion") %>%

    # Calculate coverage
    mutate(c_cdr = c_newinc * 100 / e_inc_num,
           c_cdr_lo = c_newinc * 100  / e_inc_num_hi,
           c_cdr_hi = c_newinc * 100  / e_inc_num_lo) %>%

    # merge with regional names and simplify to match structure of country table
    inner_join(who_region_shortnames, by = "g_whoregion") %>%
    select(entity,
           c_cdr,
           c_cdr_lo,
           c_cdr_hi) %>%
    arrange(desc(c_cdr))


  # C. Global (Calculate for two years as CDR for the earlier year is needed for the text)
  # - - - - - - - -
  coverage_inc_global <- filter(est_global, year>=report_year-3) %>%
    select(year,
           e_inc_num = inc.num,
           e_inc_num_lo = inc.lo.num,
           e_inc_num_hi = inc.hi.num) %>%
    mutate(entity="Global")

  coverage_global <- filter(notification, year>=report_year-3) %>%
    select(year,
           c_newinc) %>%
    # calculate global aggregate
    group_by(year) %>%
    summarise(across(c_newinc, sum, na.rm=TRUE)) %>%
    mutate(entity="Global") %>%
    ungroup() %>%
    inner_join(coverage_inc_global, by=c("entity", "year")) %>%

    # Calculate coverage
    mutate(c_cdr = c_newinc * 100 / e_inc_num,
           c_cdr_lo = c_newinc * 100  / e_inc_num_hi,
           c_cdr_hi = c_newinc * 100  / e_inc_num_lo)
  # %>%
  #   select(entity,
  #          year,
  #          c_cdr,
  #          c_cdr_lo,
  #          c_cdr_hi)

  coverage_global_latest <- filter(coverage_global, year==report_year-1) %>%
    select(entity,
           c_cdr,
           c_cdr_lo,
           c_cdr_hi)

  # D. Bring them all together
  # - - - - - - - - - - - - -

  # Create dummy records so can see a horizontal line in the output to separate countries, regions and global parts
  coverage_dummy1 <- data.frame(entity = " ", c_cdr = NA, c_cdr_lo = 0, c_cdr_hi = 100)
  coverage_dummy2 <- data.frame(entity = "  ", c_cdr = NA, c_cdr_lo = 0, c_cdr_hi = 100)


  # Create combined dataframe in order of countries then regional and global estimates
  f3.3.1_data <- rbind(coverage_country,
                       coverage_dummy1,
                       coverage_region,
                       coverage_dummy2,
                       coverage_global_latest) %>%

    # The dataframe is in the order I want, so make entity an ordered factor based on
    # what I already have. That way ggplot will not reorder by entity name
    # But I need to reverse order for plotting

    mutate(entity = factor(entity,
                           levels = rev(entity)))


  # Create summary for the text
  f3.3.1_txt <- coverage_global %>%
    select(-entity, -e_inc_num_lo, -e_inc_num_hi) %>%
    # Calculate the gap between incidence and notifications
    mutate(global_gap = e_inc_num - c_newinc) %>%
    pivot_wider(names_from = year,
                values_from = c(c_cdr, c_cdr_lo, c_cdr_hi, e_inc_num, c_newinc, global_gap)) %>%
    # Calculate the % change in notifications
    mutate(pct_change_c_newinc = abs(c_newinc_2021 - c_newinc_2020) * 100 / c_newinc_2020)

  # Add EUR and EMR
  f3.3.1_txt <- filter(coverage_region, entity %in% c("European Region", "Eastern Mediterranean Region", "Region of the Americas")) %>%
    select(entity, c_cdr) %>%
    pivot_wider(names_from = entity,
                values_from = c_cdr)  %>%
    # handle spaces in column names
    select(c_cdr_EUR = `European Region`,
           c_cdr_EMR = `Eastern Mediterranean Region`,
           c_cdr_AMR = `Region of the Americas`
    ) %>%
    cbind(f3.3.1_txt)

  # Add Mozambique pct pulmonary bac confirmed
  f3.3.1_txt <- filter(notification, year==report_year-1 & iso3 == "MOZ") %>%
    select(new_labconf, new_clindx,
           ret_rel_labconf, ret_rel_clindx) %>%
    mutate(bacconf_pct_MOZ = (new_labconf + ret_rel_labconf) * 100 / (new_labconf + new_clindx + ret_rel_labconf + ret_rel_clindx)) %>%
    select(bacconf_pct_MOZ) %>%
    cbind(f3.3.1_txt)


  # Create ist of countries with cdr < 50% for the text
  f3.3.1_txt_list_hi <- filter(coverage_country, c_cdr > 74, entity != "Mozambique") %>% # hi TC in MOZ may be due to over reliance on clinical dx
    select(entity) #%>%
    # arrange(entity)

  f3.3.1_txt_list_lo <- filter(coverage_country, c_cdr < 50) %>%
    select(entity) %>%
    arrange(entity)

  # remove the temporary dataframes
  rm(list=ls(pattern = "^coverage"))

}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.3.2 ----
# (Forest plot of TB treatment coverage in 30 countries, regionally and globally)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## for 0-14
# A. Countries
# - - - - - - - -
coverage_inc_country_014 <- filter(db_estimates_country, year==report_year-1, age_group=="0-14", sex=="a") %>%
  select(year,
         iso3,
         e_inc_num_014 = best,
         e_inc_num_014_lo = lo,
         e_inc_num_014_hi = hi) %>%

  # restrict to high burden countries
  inner_join(hbc30, by = "iso3")

coverage_country_014 <- filter(notification, year==report_year-1) %>%
  select(entity = country,
         iso3,
         c_new_014)  %>%
  inner_join(coverage_inc_country_014, by = "iso3") %>%
  select(-iso3) %>%
  mutate(c_cdr = c_new_014 * 100 / e_inc_num_014,
         c_cdr_lo = c_new_014 * 100  / e_inc_num_014_hi,
         c_cdr_hi = c_new_014 * 100  / e_inc_num_014_lo,
         # highlight countries with no data
         entity = ifelse(is.na(c_new_014), paste0(entity, "*"), entity )) %>%
  select(entity,
         c_cdr,
         c_cdr_lo,
         c_cdr_hi) %>%
  arrange(desc(c_cdr))


# B. Regions
# - - - - - - - -
coverage_inc_region_014 <- filter(db_estimates_group, year==report_year-1, group_type=="g_whoregion", age_group=="0-14", sex=="a") %>%
  select(g_whoregion = group_name ,
         e_inc_num_014 = best,
         e_inc_num_014_lo = lo,
         e_inc_num_014_hi =hi)

coverage_region_014 <- filter(notification, year==report_year-1) %>%
  select(g_whoregion,
         c_new_014) %>%
  # calculate regional aggregates
  group_by(g_whoregion) %>%
  summarise(across(c_new_014, sum, na.rm=TRUE)) %>%
  ungroup() %>%

  # merge with incidence estimates
  inner_join(coverage_inc_region_014, by = "g_whoregion") %>%

  # Calculate coverage
  mutate(c_cdr = c_new_014 * 100 / e_inc_num_014,
         c_cdr_lo = c_new_014 * 100  / e_inc_num_014_hi,
         c_cdr_hi = c_new_014 * 100  / e_inc_num_014_lo) %>%

  # merge with regional names and simplify to match structure of country table
  inner_join(who_region_shortnames, by = "g_whoregion") %>%
  select(entity,
         c_cdr,
         c_cdr_lo,
         c_cdr_hi) %>%
  arrange(desc(c_cdr))


# C. Global (Calculate for two years as CDR for the earlier year is needed for the text)
# - - - - - - - -
coverage_inc_global_014 <- filter(db_estimates_group, year==report_year-1, group_type=="global", age_group=="0-14", sex=="a") %>%
  select(year,
         e_inc_num_014 = best,
         e_inc_num_014_lo = lo,
         e_inc_num_014_hi =hi) %>%
  mutate(entity="Global")

coverage_global_014 <- filter(notification, year>=report_year-2) %>%
  select(year,
         c_new_014) %>%
  # calculate global aggregate
  group_by(year) %>%
  summarise(across(c_new_014, sum, na.rm=TRUE)) %>%
  mutate(entity="Global") %>%
  ungroup() %>%
  inner_join(coverage_inc_global_014, by=c("entity", "year")) %>%

  # Calculate coverage
  mutate(c_cdr = c_new_014 * 100 / e_inc_num_014,
         c_cdr_lo = c_new_014 * 100  / e_inc_num_014_hi,
         c_cdr_hi = c_new_014 * 100  / e_inc_num_014_lo) %>%
  select(entity,
         c_cdr,
         c_cdr_lo,
         c_cdr_hi) %>%
  arrange(desc(c_cdr))


# D. Bring them all together
# - - - - - - - - - - - - -

# Create dummy records so can see a horizontal line in the output to separate countries, regions and global parts
coverage_014_dummy1 <- data.frame(entity = " ", c_cdr = NA, c_cdr_lo = 0, c_cdr_hi = 100)
coverage_014_dummy2 <- data.frame(entity = "  ", c_cdr = NA, c_cdr_lo = 0, c_cdr_hi = 100)


# Create combined dataframe in order of countries then regional and global estimates
f3.3.2a_data <- rbind(coverage_country_014,
                      coverage_014_dummy1,
                      coverage_region_014,
                      coverage_014_dummy2,
                      coverage_global_014) %>%

  # The dataframe is in the order I want, so make entity an ordered factor based on
  # what I already have. That way ggplot will not reorder by entity name
  # But I need to reverse order for plotting

  mutate(entity = factor(entity,
                         levels = rev(entity)))

# texts
f3.3.2a_txt <- f3.3.2a_data %>%
  filter(entity=="Global")


## for 15plus
# A. Countries
# - - - - - - - -
coverage_inc_country_15plus <- filter(db_estimates_country, year==report_year-1, age_group=="15plus", sex=="a") %>%
  select(year,
         iso3,
         e_inc_num_15plus = best,
         e_inc_num_15plus_lo = lo,
         e_inc_num_15plus_hi = hi) %>%

  # restrict to high burden countries
  inner_join(hbc30, by = "iso3")

coverage_country_15plus <- filter(notification, year==report_year-1) %>%
  mutate(c_new_15plus = NZ(newrel_f15plus) + NZ(newrel_m15plus)) %>%
  select(entity = country,
         iso3,
         c_new_15plus)  %>%
  inner_join(coverage_inc_country_15plus, by = "iso3") %>%
  select(-iso3) %>%
  mutate(c_cdr = c_new_15plus * 100 / e_inc_num_15plus,
         c_cdr_lo = c_new_15plus * 100  / e_inc_num_15plus_hi,
         c_cdr_hi = c_new_15plus * 100  / e_inc_num_15plus_lo,
         # highlight countries with no data
         entity = ifelse(is.na(c_new_15plus), paste0(entity, "*"), entity )) %>%
  select(entity,
         c_cdr,
         c_cdr_lo,
         c_cdr_hi) %>%
  arrange(desc(c_cdr))


# B. Regions
# - - - - - - - -
coverage_inc_region_15plus <- filter(db_estimates_group, year==report_year-1, group_type=="g_whoregion", age_group=="15plus", sex=="a") %>%
  select(g_whoregion = group_name ,
         e_inc_num_15plus = best,
         e_inc_num_15plus_lo = lo,
         e_inc_num_15plus_hi =hi)

coverage_region_15plus <- filter(notification, year==report_year-1) %>%
  mutate(c_new_15plus = NZ(newrel_f15plus) + NZ(newrel_m15plus)) %>%
  select(g_whoregion,
         c_new_15plus) %>%
  # calculate regional aggregates
  group_by(g_whoregion) %>%
  summarise(across(c_new_15plus, sum, na.rm=TRUE)) %>%
  ungroup() %>%

  # merge with incidence estimates
  inner_join(coverage_inc_region_15plus, by = "g_whoregion") %>%

  # Calculate coverage
  mutate(c_cdr = c_new_15plus * 100 / e_inc_num_15plus,
         c_cdr_lo = c_new_15plus * 100  / e_inc_num_15plus_hi,
         c_cdr_hi = c_new_15plus * 100  / e_inc_num_15plus_lo) %>%

  # merge with regional names and simplify to match structure of country table
  inner_join(who_region_shortnames, by = "g_whoregion") %>%
  select(entity,
         c_cdr,
         c_cdr_lo,
         c_cdr_hi) %>%
  arrange(desc(c_cdr))


# C. Global (Calculate for two years as CDR for the earlier year is needed for the text)
# - - - - - - - -
coverage_inc_global_15plus <- filter(db_estimates_group, year==report_year-1, group_type=="global", age_group=="15plus", sex=="a") %>%
  select(year,
         e_inc_num_15plus = best,
         e_inc_num_15plus_lo = lo,
         e_inc_num_15plus_hi =hi) %>%
  mutate(entity="Global")

coverage_global_15plus <- filter(notification, year>=report_year-2) %>%
  mutate(c_new_15plus = NZ(newrel_f15plus) + NZ(newrel_m15plus)) %>%
  select(year,
         c_new_15plus) %>%
  # calculate global aggregate
  group_by(year) %>%
  summarise(across(c_new_15plus, sum, na.rm=TRUE)) %>%
  mutate(entity="Global") %>%
  ungroup() %>%
  inner_join(coverage_inc_global_15plus, by=c("entity", "year")) %>%

  # Calculate coverage
  mutate(c_cdr = c_new_15plus * 100 / e_inc_num_15plus,
         c_cdr_lo = c_new_15plus * 100  / e_inc_num_15plus_hi,
         c_cdr_hi = c_new_15plus * 100  / e_inc_num_15plus_lo) %>%
  select(entity,
         c_cdr,
         c_cdr_lo,
         c_cdr_hi) %>%
  arrange(desc(c_cdr))


# D. Bring them all together
# - - - - - - - - - - - - -

# Create dummy records so can see a horizontal line in the output to separate countries, regions and global parts
coverage_15plus_dummy1 <- data.frame(entity = " ", c_cdr = NA, c_cdr_lo = 0, c_cdr_hi = 100)
coverage_15plus_dummy2 <- data.frame(entity = "  ", c_cdr = NA, c_cdr_lo = 0, c_cdr_hi = 100)


# Create combined dataframe in order of countries then regional and global estimates
f3.3.2b_data <- rbind(coverage_country_15plus,
                      coverage_15plus_dummy1,
                      coverage_region_15plus,
                      coverage_15plus_dummy2,
                      coverage_global_15plus) %>%

  # The dataframe is in the order I want, so make entity an ordered factor based on
  # what I already have. That way ggplot will not reorder by entity name
  # But I need to reverse order for plotting

  mutate(entity = factor(entity,
                         levels = rev(entity)))


# texts
f3.3.2b_txt <- f3.3.2b_data %>%
  filter(entity=="Global")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.3.3 ----
# (Bubble map of difference between notifications and estimated incidence for 10 countries)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if(show_estimates){

  f3.3.3_data  <- filter(est_country, year == report_year - 1) %>%
    select(iso3,
           year,
           e_inc_num = inc.num
    ) %>%

    # Link to notifications
    inner_join(notification, by =c("year", "iso3")) %>%
    select(iso3,
           country,
           e_inc_num,
           c_newinc) %>%

    # Calculate the gap and use that for the bubble sizes
    mutate(size = e_inc_num - c_newinc) %>%

    # limit to the top 10 by size of gap
    top_n(10, size) %>%

    # sort in descending order so can list the country names in the figure footnote
    arrange(desc(size))


  # Summary number of gaps for the section text
  # Get global incidence
  f3.3.3_txt <- filter(est_global, year == report_year-1) %>%
    select(e_inc_num = inc.num)

  # Add global notification
  f3.3.3_txt <- filter(notification, year == report_year-1) %>%
    select(year,
           c_newinc) %>%
    # calculate global aggregate
    group_by(year) %>%
    summarise(across(c_newinc, sum, na.rm=TRUE)) %>%
    ungroup() %>%
    cbind(f3.3.3_txt) %>%

    # Calculate the global gap and drop the other variables
    mutate(gap = e_inc_num - c_newinc) %>%
    select(gap) %>%

    # Calculate % of global gap contributed by the top 10 countries
    cbind(f3.3.3_data) %>%
    mutate(pct_gap = size * 100 / gap) %>%

    # flip wider for easy quoting in text
    select(iso3, pct_gap) %>%
    pivot_wider(names_from = iso3,
                names_prefix = "pct_gap_",
                values_from = pct_gap) %>%

    # Add total % accounted for by all ten countries
    mutate(pct_gap_top_ten = rowSums(.))

}




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.3.4 ----
# (Line and ribbon chart of estimated TB/HIV incidence, number notified and number on antiretroviral therapy)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if(show_estimates){

  inctbhiv_data <- filter(est_global, year >= 2004) %>%
    select(year,
           e_inc_tbhiv_num = inc.h.num,
           e_inc_tbhiv_num_lo = inc.h.lo.num,
           e_inc_tbhiv_num_hi = inc.h.hi.num)



  f3.3.4_data <- filter(TBHIV_for_aggregates, year>=2004) %>%
    select(year,
           hivtest_pos,
           hiv_art,
           # new variables from 2015 onwards
           newrel_hivpos,
           newrel_art) %>%
    group_by(year) %>%
    summarise(across(hivtest_pos:newrel_art, sum, na.rm=TRUE)) %>%

    # Merge pre/post 2014 variables
    mutate(hivpos = ifelse(year < 2015,
                           hivtest_pos,
                           newrel_hivpos),
           art = ifelse(year < 2015,
                        hiv_art,
                        newrel_art)) %>%
    select(year,
           hivpos,
           art) %>%

    inner_join(inctbhiv_data, by = "year")

  # Summary of numbers for the section text
  f3.3.4_txt <- filter(f3.3.4_data, year>=report_year-2) %>%
    mutate(c_art_notified = art * 100 / hivpos,
           c_art_estimated = art * 100 / e_inc_tbhiv_num,
           gap_hivpos_estimated = 100 - hivpos * 100 / e_inc_tbhiv_num) %>%
    select(year,
           c_art_notified,
           c_art_estimated,
           gap_hivpos_estimated) %>%
    pivot_wider(names_from = year,
                values_from = c(c_art_notified, c_art_estimated, gap_hivpos_estimated))

  f3.3.4_num_txt <- filter(f3.3.4_data, year>=report_year-2) %>%
    mutate(c_art_notified =  hivpos - art,
           c_art_estimated = e_inc_tbhiv_num - art) %>%
    select(year,
           c_art_notified,
           c_art_estimated) %>%
    pivot_wider(names_from = year,
                values_from = c(c_art_notified, c_art_estimated))
  

  #-------------------
  ### regional data 
  inctbhiv_region_data <- filter(est_regional, year >= 2004) %>%
    select(year,
           g_whoregion = g.whoregion,
           e_inc_tbhiv_num = inc.h.num,
           e_inc_tbhiv_num_lo = inc.h.lo.num,
           e_inc_tbhiv_num_hi = inc.h.hi.num)
  
  
  
  f3.3.4_region_data <- filter(TBHIV_for_aggregates, year>=2004) %>%
    select(year,
           g_whoregion,
           hivtest_pos,
           hiv_art,
           # new variables from 2015 onwards
           newrel_hivpos,
           newrel_art) %>%
    group_by(year,g_whoregion) %>%
    summarise(across(hivtest_pos:newrel_art, sum, na.rm=TRUE)) %>%
    
    # Merge pre/post 2014 variables
    mutate(hivpos = ifelse(year < 2015,
                           hivtest_pos,
                           newrel_hivpos),
           art = ifelse(year < 2015,
                        hiv_art,
                        newrel_art)) %>%
    select(year,
           g_whoregion,
           hivpos,
           art) %>%
    
    inner_join(inctbhiv_region_data, by = c("year","g_whoregion"))
  
  # Summary of numbers for the section text
  f3.3.4_region_txt <- filter(f3.3.4_region_data, year>=report_year-2) %>%
    mutate(c_art_notified = art * 100 / hivpos,
           c_art_estimated = art * 100 / e_inc_tbhiv_num,
           gap_hivpos_estimated = 100 - hivpos * 100 / e_inc_tbhiv_num) %>%
    select(year,
           g_whoregion,
           c_art_notified,
           c_art_estimated,
           gap_hivpos_estimated) %>%
    pivot_wider(names_from = year,
                values_from = c(c_art_notified, c_art_estimated, gap_hivpos_estimated))
  
  f3.3.4_region_num_txt <- filter(f3.3.4_region_data, year>=report_year-2) %>%
    mutate(c_art_notified =  hivpos - art,
           c_art_estimated = e_inc_tbhiv_num - art) %>%
    select(year,
           g_whoregion,
           c_art_notified,
           c_art_estimated) %>%
    pivot_wider(names_from = year,
                values_from = c(c_art_notified, c_art_estimated))
  
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.3.5 ----
# (Forest plot of TB/HIV treatment coverage in 30 countries, regionally and globally)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if(show_estimates){

  # A. Countries
  # - - - - - - - -
  coverage_inc_country <- filter(est_country, year==report_year-1) %>%
    select(year,
           iso3,
           e_inc_tbhiv_num = inc.h.num,
           e_inc_tbhiv_num_lo = inc.h.lo.num,
           e_inc_tbhiv_num_hi = inc.h.hi.num) %>%

    # restrict to high burden countries
    inner_join(hbtbhiv30, by = "iso3")

  coverage_country <- filter(TBHIV_for_aggregates, year==report_year-1) %>%
    select(entity = country,
           iso3,
           newrel_art)  %>%
    inner_join(coverage_inc_country, by = "iso3") %>%
    select(-iso3) %>%
    mutate(c_art = newrel_art * 100 / e_inc_tbhiv_num,
           c_art_lo = newrel_art * 100  / e_inc_tbhiv_num_hi,
           c_art_hi = newrel_art * 100  / e_inc_tbhiv_num_lo,
           # highlight countries with no data
           entity = ifelse(is.na(newrel_art), paste0(entity, "*"), entity )) %>%
    select(entity,
           c_art,
           c_art_lo,
           c_art_hi) %>%
    arrange(desc(c_art))


  # B. Regions
  # - - - - - - - -
  coverage_inc_region <- filter(est_regional, year==report_year-1) %>%
    select(g_whoregion = g.whoregion,
           e_inc_tbhiv_num = inc.h.num,
           e_inc_tbhiv_num_lo = inc.h.lo.num,
           e_inc_tbhiv_num_hi = inc.h.hi.num)

  coverage_region <- filter(TBHIV_for_aggregates, year==report_year-1) %>%
    select(g_whoregion,
           newrel_art) %>%
    # calculate regional aggregates
    group_by(g_whoregion) %>%
    summarise(across(newrel_art, sum, na.rm=TRUE)) %>%
    ungroup() %>%

    # merge with incidence estimates
    inner_join(coverage_inc_region, by = "g_whoregion") %>%

    # Calculate coverage
    mutate(c_art = newrel_art * 100 / e_inc_tbhiv_num,
           c_art_lo = newrel_art * 100  / e_inc_tbhiv_num_hi,
           c_art_hi = newrel_art * 100  / e_inc_tbhiv_num_lo) %>%

    # merge with regional names and simplify to match structure of country table
    inner_join(who_region_shortnames, by = "g_whoregion") %>%
    select(entity,
           c_art,
           c_art_lo,
           c_art_hi) %>%
    arrange(desc(c_art))

  # C. Global
  # - - - - - - - -
  coverage_inc_global <- filter(est_global, year==report_year-1) %>%
    select(e_inc_tbhiv_num = inc.h.num,
           e_inc_tbhiv_num_lo = inc.h.lo.num,
           e_inc_tbhiv_num_hi = inc.h.hi.num) %>%
    mutate(entity="Global")

  coverage_global <- filter(TBHIV_for_aggregates, year==report_year-1) %>%
    select(newrel_art) %>%
    # calculate global aggregate
    summarise(across(newrel_art, sum, na.rm=TRUE)) %>%
    mutate(entity="Global") %>%
    inner_join(coverage_inc_global, by="entity") %>%

    # Calculate coverage
    mutate(c_art = newrel_art * 100 / e_inc_tbhiv_num,
           c_art_lo = newrel_art * 100  / e_inc_tbhiv_num_hi,
           c_art_hi = newrel_art * 100  / e_inc_tbhiv_num_lo) %>%
    select(entity,
           c_art,
           c_art_lo,
           c_art_hi)

  # D. Bring them all together
  # - - - - - - - - - - - - -

  # Create dummy records so can see a horizontal line in the output to separate countries, regions and global parts
  coverage_dummy1 <- data.frame(entity = " ", c_art = NA, c_art_lo = 0, c_art_hi = 100)
  coverage_dummy2 <- data.frame(entity = "  ", c_art = NA, c_art_lo = 0, c_art_hi = 100)


  # Create combined dataframe in order of countries then regional and global estimates
  f3.3.5_data <- rbind(coverage_country,
                       coverage_dummy1,
                       coverage_region,
                       coverage_dummy2,
                       coverage_global) %>%

    # The dataframe is in the order I want, so make entity an ordered factor based on
    # what I already have. That way ggplot will not reorder by entity name
    # But I need to reverse order for plotting

    mutate(entity = factor(entity,
                           levels = rev(entity)))


  # Summary of numbers for the section text: Find top and bottom country
  # (Use c(1,30) as selector since data are already in descending order)
  f3.3.5_txt <- coverage_country[c(1,30),] %>%
    mutate(pos = c("top", "bottom")) %>%
    select(pos, entity, c_art) %>%
    pivot_wider(names_from = pos,
                values_from = c(entity, c_art))

  # Add number of countries with coverage >= 50%
  f3.3.5_txt <- filter(coverage_country, c_art >= 50) %>%
    summarise(over_50 = n()) %>%
    cbind(f3.3.5_txt)

  # remove the temporary dataframes
  rm(list=ls(pattern = "^coverage"))

}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.3.6 ----
# (Horizontal bar chart showing TB treatment outcomes for WHO regions and globally)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# A. Regional aggregates
# - - - - - - - - - - - -
txout_region  <- filter(outcomes, year==report_year - 2) %>%

  select(iso2,
         g_whoregion,
         newrel_coh,
         newrel_succ,
         newrel_fail,
         newrel_died,
         newrel_lost,
         c_newrel_neval) %>%

  group_by(g_whoregion) %>%
  summarise(across(newrel_coh:c_newrel_neval, sum, na.rm=TRUE)) %>%
  ungroup() %>%

  # merge with regional names
  inner_join(who_region_shortnames, by = "g_whoregion") %>%
  select(-g_whoregion) %>%

  # Calculate outcome proportions for plotting as stacked bars
  calculate_outcomes_pct("newrel_") %>%

  # Sort regions in descending order of success rate
  arrange(desc(`Treatment success`))


# B. Global aggregates
# - - - - - - - - - - - -
txout_global  <- filter(outcomes, year==report_year - 2) %>%

  select(iso2,
         newrel_coh,
         newrel_succ,
         newrel_fail,
         newrel_died,
         newrel_lost,
         c_newrel_neval) %>%

  summarise(across(newrel_coh:c_newrel_neval, sum, na.rm=TRUE)) %>%
  ungroup()  %>%
  mutate(entity="Global")  %>%

  # Calculate outcome proportions for plotting as stacked bars
  calculate_outcomes_pct("newrel_")


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
f3.3.6_data <- rbind(txout_region, txout_dummy, txout_global) %>%

  # Keep record of current order (in reverse) so plot comes out as we want it
  mutate(entity = factor(entity, levels=rev(entity))) %>%

  # Drop the actual numbers and keep percentages
  select(-coh, -succ, -fail, -died, -lost, -c_neval) %>%

  # Flip into long mode for stacked bar plotting
  pivot_longer(cols = `Treatment success`:`Not evaluated`,
               names_to = "outcome")

# Create summary for section text
f3.3.6_txt <- filter(f3.3.6_data, outcome == "Treatment success" &
                       entity %in% c("Global", "Eastern Mediterranean Region", "European Region")) %>%
  pivot_wider(names_from = entity,
              names_prefix = "tsr_",
              values_from = value) %>%
  #simplify variable names
  select(tsr_EMR = `tsr_Eastern Mediterranean Region`,
         tsr_EUR = `tsr_European Region`,
         tsr_Global)

# tidy up
rm(list = ls(pattern = "^txout_"))




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.3.7 ----
# (Panel of 3 horizontal bar charts showing TB treatment outcomes globally by year since 2012 for TB, TB/HIV and MDR/RR-TB)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

txout_tb_data <- outcomes %>%
  filter(between(year, 2012, report_year - 2)) %>%
  select(iso2,
         year,
         newrel_coh,
         newrel_succ,
         newrel_fail,
         newrel_died,
         newrel_lost,
         c_newrel_neval) %>%

  # calculate global aggregates
  group_by(year) %>%
  summarise(across(newrel_coh:c_newrel_neval, sum, na.rm=TRUE)) %>%
  ungroup()%>%

  # Calculate outcome proportions for plotting as stacked bars
  calculate_outcomes_pct("newrel_") %>%

  # Drop the actual numbers and keep percentages
  select(-coh, -succ, -fail, -died, -lost, -c_neval) %>%

  # Add tx group type
  mutate(subgroup = "New and relapse TB cases")

txout_hiv_data <- outcomes %>%
  filter(between(year, 2012, report_year - 2)) %>%
  select(iso2,
         year,
         country,
         contains("tbhiv_")) %>%

  # calculate global aggregates
  group_by(year) %>%
  summarise(across(tbhiv_coh:c_tbhiv_neval, sum, na.rm=TRUE)) %>%
  ungroup() %>%

  # Calculate outcome proportions for plotting as stacked bars
  calculate_outcomes_pct("tbhiv_") %>%

  # Drop the actual numbers and keep percentages
  select(-coh, -succ, -fail, -died, -lost, -c_neval) %>%

  # Add tx group type
  mutate(subgroup = "New and relapse HIV-positive TB cases")

# Combine the two data frames
f3.3.7_data <- rbind(txout_tb_data, txout_hiv_data) %>%

  # flip into long format
  pivot_longer(cols = `Treatment success`:`Not evaluated`,
               names_to = "outcome") %>%

  # Determine the order of subgroup for plotting
  mutate(subgroup = factor(subgroup,
                           levels = c("New and relapse TB cases",
                                      "New and relapse HIV-positive TB cases")))

# tidy up
rm(list = ls(pattern = "^txout_"))

# text
f3.3.7_txt <- f3.3.7_data %>%
  filter(year==2020, outcome=="Treatment success", subgroup == "New and relapse HIV-positive TB cases")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.3.8 ----
# (Horizontal bar chart showing TB treatment outcomes in HIV-positive cases for WHO regions and globally)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# A. Regional aggregates
# - - - - - - - - - - - -
txout_region  <- filter(outcomes, year==report_year - 2) %>%

  select(iso2,
         g_whoregion,
         contains("tbhiv_")) %>%

  group_by(g_whoregion) %>%
  summarise(across(tbhiv_coh:c_tbhiv_neval, sum, na.rm=TRUE)) %>%
  ungroup() %>%

  # merge with regional names
  inner_join(who_region_shortnames, by = "g_whoregion") %>%
  select(-g_whoregion) %>%

  # Calculate outcome proportions for plotting as stacked bars
  calculate_outcomes_pct("tbhiv_") %>%

  # Sort regions in descending order of success rate
  arrange(desc(`Treatment success`))


# B. Global aggregates
# - - - - - - - - - - - -
txout_global  <- filter(outcomes, year==report_year - 2) %>%

  select(iso2,
         contains("tbhiv_")) %>%

  summarise(across(tbhiv_coh:c_tbhiv_neval, sum, na.rm=TRUE)) %>%
  ungroup()  %>%
  mutate(entity="Global")  %>%

  # Calculate outcome proportions for plotting as stacked bars
  calculate_outcomes_pct("tbhiv_")


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
f3.3.8_data <- rbind(txout_region, txout_dummy, txout_global) %>%

  # Keep record of current order (in reverse) so plot comes out as we want it
  mutate(entity = factor(entity, levels=rev(entity))) %>%

  # Drop the actual numbers and keep percentages
  select(-coh, -succ, -fail, -died, -lost, -c_neval) %>%

  # Flip into long mode for stacked bar plotting
  pivot_longer(cols = `Treatment success`:`Not evaluated`,
               names_to = "outcome")

# Create summary for section text
f3.3.8_txt <- filter(f3.3.8_data, outcome == "Treatment success" & entity == "Global") %>%
  #simplify variable name
  select(tsr_tbhiv_Global = value)

# tidy up
rm(list = ls(pattern = "^txout_"))




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.3.9 ----
# (Horizontal bar chart showing TB treatment success rates in children for WHO regions and globally)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# A. Regional aggregates
# - - - - - - - - - - - -
txout_region  <- filter(outcomes, year==report_year - 2) %>%

  select(iso2,
         g_whoregion,
         newrel_014_coh,
         newrel_014_succ) %>%

  group_by(g_whoregion) %>%
  summarise(across(newrel_014_coh:newrel_014_succ, sum, na.rm=TRUE)) %>%
  ungroup() %>%

  # merge with regional names
  inner_join(who_region_shortnames, by = "g_whoregion") %>%
  select(-g_whoregion) %>%

  # Calculate treatment success rate
  mutate(c_tsr_014 = newrel_014_succ * 100/newrel_014_coh) %>%

  # Sort regions in descending order of success rate
  arrange(desc(c_tsr_014))

# B. Global aggregates
# - - - - - - - - - - - -
txout_global  <- filter(outcomes, year==report_year - 2) %>%

  select(iso2,
         newrel_014_coh,
         newrel_014_succ) %>%

  summarise(across(newrel_014_coh:newrel_014_succ, sum, na.rm=TRUE)) %>%
  ungroup()  %>%
  mutate(entity="Global")  %>%

  # Calculate treatment success rate
  mutate(c_tsr_014 = newrel_014_succ * 100/newrel_014_coh)

# Create a dummy record to make a gap in the output to separate regional and global parts
txout_dummy <- data.frame(entity = " ", newrel_014_coh = NA, newrel_014_succ = NA, c_tsr_014 = NA)

# Create combined table in order of countries then regional and global estimates
f3.3.9_data <- rbind(txout_region, txout_dummy, txout_global) %>%

  # Keep record of current order (in reverse) so plot comes out as we want it
  mutate(entity = factor(entity, levels=rev(entity))) %>%

  # Drop the actual numbers and keep percentages
  select(-newrel_014_coh,
         -newrel_014_succ)

# Create summary for section text and footnote
f3.3.9_txt <- filter(f3.3.9_data, entity == "Global") %>%
  #simplify variable name
  select(tsr_014_Global = c_tsr_014)

# Add number of countries that reported and total cohort
f3.3.9_txt <- filter(outcomes, year==report_year - 2 &
                       !is.na(newrel_014_coh) & !is.na(newrel_014_succ)) %>%
  summarise(countries = n(),
            kids_coh = sum(newrel_014_coh, na.rm=TRUE)) %>%
  cbind(f3.3.9_txt)

# Calculate percent of notified that had an outcome reported
f3.3.9_txt <- filter(notification, year==report_year-2) %>%
  summarise(kids_notified = sum(c_new_014, na.rm=TRUE)) %>%
  cbind(f3.3.9_txt) %>%
  mutate(kids_outcome_pct = kids_coh * 100 / kids_notified)

# tidy up
rm(list = ls(pattern = "^txout_"))



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 3.3.9 ----
# (Panel of bar charts showing treatment outcomes in absolute numbers by year since 2000 globally and for WHO regions)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# txoutnum_regional <- outcomes %>%
#   filter(year >= 2000 & year <= report_year - 2) %>%
#   select(year,
#          g_whoregion,
#          new_sp_coh,
#          new_sp_cur,
#          new_sp_cmplt,
#          c_new_sp_neval,
#          new_snep_coh,
#          new_snep_cmplt,
#          c_new_snep_neval,
#          newrel_coh,
#          newrel_succ,
#          c_newrel_neval)%>%
#   group_by(year, g_whoregion) %>%
#   summarise(across(new_sp_coh:c_newrel_neval, sum, na.rm=TRUE)) %>%
#   ungroup() %>%
#
#   # merge with regional names
#   inner_join(who_region_shortnames, by = "g_whoregion") %>%
#   select(-g_whoregion)
#
# txoutnum_global <- outcomes %>%
#   filter(year >= 2000 & year <= report_year - 2) %>%
#   select(year,
#          g_whoregion,
#          new_sp_coh,
#          new_sp_cur,
#          new_sp_cmplt,
#          c_new_sp_neval,
#          new_snep_coh,
#          new_snep_cmplt,
#          c_new_snep_neval,
#          newrel_coh,
#          newrel_succ,
#          c_newrel_neval)%>%
#   group_by(year) %>%
#   summarise(across(new_sp_coh:c_newrel_neval, sum, na.rm=TRUE)) %>%
#   ungroup() %>%
#
#   mutate(entity = "Global")
#
# #Combine regional and global data and reorganise
# f3.3.9_data <- rbind(txoutnum_regional, txoutnum_global) %>%
#
#   mutate(entity = factor(entity,
#                          levels = c("African Region", "Region of the Americas", "South-East Asia Region",
#                                     "European Region", "Eastern Mediterranean Region", "Western Pacific Region",
#                                     "Global") )) %>%
#
#   # Simplify the data for plotting
#   mutate(success = (new_sp_cur + new_sp_cmplt + new_snep_cmplt + newrel_succ) / 1e6,
#          neval = (c_new_sp_neval + c_new_snep_neval + c_newrel_neval) / 1e6,
#          coh = (new_sp_coh + new_snep_coh + newrel_coh) / 1e6) %>%
#   mutate(fail_other = coh - success - neval) %>%
#   select(entity,
#          year,
#          `Treatment success` = success,
#          `Failure/Died/Lost to follow-up` = fail_other,
#          `Not evaluated` = neval) %>%
#
#   # Flip into long mode for stacked bar plotting
#   pivot_longer(cols = `Treatment success`:`Not evaluated`,
#                names_to = "outcome")
#
# # tidy up
# rm(list = ls(pattern = "^txoutnum_"))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: Table 3.3.1 ----
# (Table showing cumulative number of lives saved by WHO region and globally)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if(show_estimates){
  library("splitstackshape")

  ftb2 <- Vectorize(function(x) {
    # formatter according to GTB rounding rules
    # https://docs.google.com/document/d/1cu_syknBiF3scAX9d7hEfN0LZwdG40i8ttN6yua2xTQ/edit
    #' @param x vector of values
    #' @export
    stopifnot(!is.character(x))

    fstring = "-"

    if (!is.na(x) & is.numeric(x) & x < 2e9) {

      fstring <-  ifelse(x==0, "0",
                         ifelse(signif(x, 1) < 0.1, formatC(signif(x,3), format="f", digits=3),
                                ifelse(signif(x, 2) < 1, formatC(signif(x,2), format="f", digits=2),
                                       ifelse(signif(x, 2) < 10, formatC(signif(x,2), format="f", digits=1),
                                              ifelse(signif(x, 3) < 100, formatC(signif(x, 2), big.mark=" ", format="d"),
                                                     formatC(signif(x, 3), big.mark=" ", format="d"))))))




    }

    return(fstring)
  }, 'x')

  beg <- "("
  endash <- "–"
  end <- ")"

  # Uses the lives_saved table which Hazim tweaked from Philippe's original
  # to standardise numbers to 2 significant figures to match what Irwin is showing
  # in table E2
  t3.3.1_data <- lives_saved %>%
    # Add long region names
    left_join(who_region_shortnames, by = c("Region" = "g_whoregion")) %>%
    mutate(entity = ifelse(Region=="Global", "Global", as.character(entity))) %>%

    # Put entity at the beginning and drop Region
    select(entity,
           saved.hivneg,
           saved.hivneg.ui,
           saved.hivpos,
           saved.hivpos.ui,
           saved,
           saved.ui,
           -Region)  %>%

    # Change the order based on the bizarre method chosen by WHO Press ...
    mutate(entity = factor(entity,
                           levels = c("African Region", "Region of the Americas", "South-East Asia Region",
                                      "European Region", "Eastern Mediterranean Region", "Western Pacific Region",
                                      "Global"))) %>%
    arrange(entity)

  t3.3.1_data <- t3.3.1_data %>%
    cSplit(c('saved.hivneg.ui','saved.hivpos.ui','saved.ui'), sep = "-", type.convert=F) %>%
    rename(saved.hivneg.lo=5, saved.hivneg.hi=6, saved.hivpos.lo=7, saved.hivpos.hi=8, saved.lo=9, saved.hi=10) %>%
    cSplit(c('saved.hivneg.lo','saved.hivneg.hi','saved.hivpos.lo','saved.hivpos.hi','saved.lo','saved.hi'), sep = ")", type.convert=F) %>%
    rename(saved.hivneg.lo=5, saved.hivneg.hi=6, saved.hivpos.lo=7, saved.hivpos.hi=8, saved.lo=9, saved.hi=10) %>%
    cSplit(c('saved.hivneg.lo','saved.hivpos.lo','saved.lo'), sep = "(", type.convert=F) %>%
    rename(saved.hivneg.lo=9, saved.hivpos.lo=11,  saved.lo=13) %>%
    select(entity:saved.hi,saved.hivneg.lo,saved.hivpos.lo,saved.lo) %>%
    mutate_if(is.character, as.numeric) %>%
    mutate_if(is.numeric, ftb2) %>%
    mutate(saved.hivneg.ui=paste0(beg,saved.hivneg.lo,endash,saved.hivneg.hi,end),
           saved.hivpos.ui=paste0(beg,saved.hivpos.lo,endash,saved.hivpos.hi,end),
           saved.ui=paste0(beg,saved.lo,endash,saved.hi,end)
    ) %>%
    select(entity,saved.hivneg,saved.hivneg.ui,saved.hivpos,saved.hivpos.ui,saved,saved.ui)


  # Get summary for easy quoting in the text
  t3.3.1_txt <- filter(t3.3.1_data, entity=="Global") %>%
    select(saved) %>%
    as.numeric()

}


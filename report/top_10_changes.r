# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Table for Katherine: the top 10 countries with
# (a) the largest increases in the estimated number of TB deaths between 2020 and 2021 and
# (b) the largest increases in the estimated number of incident cases between 2020 and 2021
#
# Hazim Timimi, October 2022
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

library(here)

# Load the data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
source(here('report/ch3_ch4_load_data.r'))


mort <- est_country %>%
  filter(year >= 2020) %>%
  select(iso3, year, mort.nh.num) %>%
  arrange(iso3, year) %>%
  mutate(mort.nh.num.prev = ifelse(year == 2021, lag(mort.nh.num), NA),
         delta = ifelse( year == 2021, mort.nh.num - lag(mort.nh.num), NA)) %>%
  filter(year == 2021) %>%
  arrange(desc(delta)) %>%
  head(10) %>%
  inner_join(list_iso3_country, by = "iso3") %>%
  select(country,
         year,
         mortality_hiv_neg = mort.nh.num,
         mortality_hiv_neg_2020 = mort.nh.num.prev,
         change = delta)

print(mort)




inc <- est_country %>%
  filter(year >= 2020) %>%
  select(iso3, year, inc.num) %>%
  arrange(iso3, year) %>%
  mutate(inc.num.prev = ifelse(year == 2021, lag(inc.num), NA),
         delta = ifelse( year == 2021, inc.num - lag(inc.num), NA)) %>%
  filter(year == 2021) %>%
  arrange(desc(delta)) %>%
  head(10) %>%
  inner_join(list_iso3_country, by = "iso3") %>%
  select(country,
         year,
         incidence = inc.num,
         incidence_2020 = inc.num.prev,
         change = delta)


print(inc)


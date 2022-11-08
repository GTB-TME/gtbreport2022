# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data preparation script for ch1.rmd
# Takuya Yamanaka, Augusut 2022
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Load chapter 3 and 4 packages, settings and data
source(here::here('report/ch3_ch4_load_data.r'))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 1.1 and 1.2 ----
# (Global trend in case notifications of people newly diagnosed with TB, 2016–2021)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
f1.2_data <- filter(notification, year >= 2015) %>%
  select(year,
         g_whoregion,
         c_newinc) %>%

  # calculate regional aggregates
  group_by(g_whoregion,year) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%

  # merge with regional names
  inner_join(who_region_shortnames, by = "g_whoregion") %>%
  arrange(entity) %>%
  select(-g_whoregion)

# to adjust yaxis range in facet_wrap
f1.2_data <- data.table(f1.2_data)
f1.2_data[g_whoregion == "AFR",y_min := 1.25*1e6]
f1.2_data[g_whoregion == "AFR",y_max := 1.5*1e6]
f1.2_data[g_whoregion == "AMR",y_min := 0.19*1e6]
f1.2_data[g_whoregion == "AMR",y_max := 0.24*1e6]
f1.2_data[g_whoregion == "SEA",y_min := 2.5*1e6]
f1.2_data[g_whoregion == "SEA",y_max := 3.5*1e6]
f1.2_data[g_whoregion == "EUR",y_min := 0.15*1e6]
f1.2_data[g_whoregion == "EUR",y_max := 0.25*1e6]
f1.2_data[g_whoregion == "EMR",y_min := 0.39*1e6]
f1.2_data[g_whoregion == "EMR",y_max := 0.55*1e6]
f1.2_data[g_whoregion == "WPR",y_min := 1*1e6]
f1.2_data[g_whoregion == "WPR",y_max := 1.5*1e6]

# texts!
f1.2_text <- f1.2_data %>%
  filter(entity=="African Region", year==2019|year==2020) %>%
  pivot_wider(names_from = year,
              values_from = c_newinc) %>%
  mutate(pct_decline=(1-`2020`/`2019`)*100)

f1.2_txt <- filter(f1.2_data, year>=2019) %>%
  mutate(c_newinc_p = lag(c_newinc)) %>%
  mutate(pct_dif = (c_newinc - c_newinc_p)*100/c_newinc_p) %>%
  filter(year==2020, entity=="African Region")


# Add global summary to the regional summary
f1.1_data <- f1.2_data %>%
  group_by(year) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  mutate(entity="Global")

f1.1_txt <- filter(f1.1_data, year>=2019) %>%
  mutate(c_newinc_p = lag(c_newinc)) %>%
  mutate(pct_dif = (c_newinc - c_newinc_p)*100/c_newinc_p) %>%
  filter(year==2020)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 1.3 ----
# (Countries with the largest contributions to global shortfalls in TB notifications in 2020 and 2021, compared with 2019)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# grouping countries by changes in 2019-2020 and 2020-2021
segmentation <- notification %>%
  select(iso3,year,country,g_whoregion,c_newinc) %>%
  pivot_wider(names_from = year, values_from = c_newinc) %>%
  mutate(pct1920 = `2020`/`2019`-1,
         pct1921 = `2021`/`2019`-1,
         pct2021 = `2021`/`2020`-1,
         num1920 = `2019`-`2020`,
         num1921 = `2019`-`2021`,
         num2021 = `2020`-`2021`) %>%
  select(iso3,country,g_whoregion,`2019`:`2021`,pct1920:num2021) %>%
  mutate(col=ifelse(pct1920<(-0.05)&pct2021>=0&num1921>0&pct1920,"a", #impact in 2020 and some recovery in 2021
                    ifelse(pct1920<(-0.05)&pct2021>=0&num1921<=0,"b", #impact in 2020 and recovery in 2021 exceeding 2019
                           ifelse(pct1920>=(-0.05)&pct2021<0,"d", #impact only in 2021
                                  ifelse(pct1920<(-0.05)&pct2021<0,"c", #impact both in 2020 and 2021
                                         ifelse(pct1920>0&pct2021>0,"e", #increases in 2020 and 2021
                                                "f")))))) %>% # the rest
  left_join(list_hbcs_plus_wl,by="iso3") %>%
  rename(country=country.x, country_hbc=country.y) %>%
  mutate(hbc=ifelse(is.na(country_hbc),0,1))

## Percentage of shortfall
f1.3_data <- segmentation %>%
  mutate(tot2020=sum(num1920,na.rm=T),
         tot2021=sum(num1921,na.rm=T))

## for 2020
f1.3a_data <- f1.3_data %>%
  arrange(desc(num1920)) %>%
  slice(1:10) %>%
  mutate(pct_contribute2020=num1920/tot2020) %>%
  mutate(cumsum=cumsum(pct_contribute2020)) %>%
  # mutate(hit90=ifelse(cumsum>=0.9&cumsum<=0.902,"yes","no"))
  mutate(hit90=ifelse(cumsum<=0.902,"yes","no")) %>%
  mutate(country = ifelse(iso3=="CHN"|iso3=="ZAF",paste0(country,"\u1D43"),country)) %>%
  select(iso3, country, pct_contribute2020, cumsum, hit90)

f1.3a_sel_order <-
  f1.3a_data %>%
  arrange(pct_contribute2020) %>%
  mutate(country = factor(country))

f1.3a_txt <- f1.3a_data %>%
  select(iso3,country,pct_contribute2020) %>%
  filter(country=="India"|country=="Indonesia"|country=="Philippines") %>%
  mutate(cumsum=cumsum(pct_contribute2020)*100) %>%
  select(cumsum) %>%
  slice(3)

## for 2021
f1.3b_data <- f1.3_data %>%
  arrange(desc(num1921)) %>%
  slice(1:10) %>%
  mutate(pct_contribute2021=num1921/tot2021) %>%
  mutate(cumsum=cumsum(pct_contribute2021)) %>%
  # mutate(hit90=ifelse(cumsum>=0.9&cumsum<=0.93,"yes","no"))
  mutate(hit90=ifelse(cumsum<=0.93,"yes","no")) %>%
  mutate(country = ifelse(iso3=="CHN"|iso3=="ZAF",paste0(country,"\u1D43"),country)) %>%
  select(iso3, country, pct_contribute2021, cumsum, hit90)

f1.3b_sel_order <-
  f1.3b_data %>%
  arrange(pct_contribute2021) %>%
  mutate(country = factor(country))

f1.3b_txt <- f1.3b_data %>%
  select(iso3,country,pct_contribute2021) %>%
  filter(country=="India"|country=="Indonesia"|country=="Philippines") %>%
  mutate(cumsum=cumsum(pct_contribute2021)*100) %>%
  select(cumsum) %>%
  slice(3)

f1.3_region_text <- notification %>% # this is for main PDF report
  select(iso3,year,country,g_whoregion,c_newinc) %>%
  filter(year>2018) %>%
  pivot_wider(names_from = year, values_from = c_newinc) %>%
  group_by(g_whoregion) %>%
  summarize(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  mutate(pct1920 = `2020`/`2019`-1,
         pct1921 = `2021`/`2019`-1,
         pct2021 = `2021`/`2020`-1,
         num1920 = `2019`-`2020`,
         num1921 = `2019`-`2021`,
         num2021 = `2020`-`2021`) %>%
  mutate(tot2020=sum(num1920,na.rm=T),
         tot2021=sum(num1921,na.rm=T)) %>%
  mutate(pct_contribute2020=num1920/tot2020) %>%
  mutate(pct_contribute2021=num1921/tot2021) %>%
  arrange(desc(pct_contribute2020)) %>%
  ungroup() %>%
  slice(1:2) %>%
  mutate(cumsum2020=cumsum(pct_contribute2020)*100) %>%
  mutate(cumsum2021=cumsum(pct_contribute2021)*100) %>%
  select(cumsum2020,cumsum2021) %>%
  slice(2)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  Data: figure 1.4 ----
#  Changes in national case notificationsa of people newly diagnosed with TB (%), 2019–2020 and 2020–2021
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

f1.4_data <- segmentation %>%
  filter(`2019`>100) %>%
  select(iso3, country, hbc, country_hbc,  pct1920, pct1921, pct2021, cn2021=`2021`)

data_breaks <- data.frame(xstart = c(min(f1.4_data$pct1920)*100-5,min(f1.4_data$pct1920)*100-5,0,0),  # Create data with breaks
                          xend   = c(0,0,max(f1.4_data$pct1920)*100+5,max(f1.4_data$pct1920)*100+5),
                          ystart = c(min(f1.4_data$pct2021,na.rm = T)*100-5,0,min(f1.4_data$pct2021,na.rm = T)*100-5,0),
                          yend   = c(0,max(f1.4_data$pct2021,na.rm = T)*100+5,0,max(f1.4_data$pct2021,na.rm = T)*100+5),
                          colors = factor(1:4))

data_breaks2 <- data.frame(xstart = c(min(f1.4_data$pct1920)*100-5,min(f1.4_data$pct1920)*100-5,0,0),  # Create data with breaks
                          xend   = c(0,0,max(f1.4_data$pct1920)*100+5,max(f1.4_data$pct1920)*100+5),
                          ystart = c(min(f1.4_data$pct1921,na.rm = T)*100-5,0,min(f1.4_data$pct1921,na.rm = T)*100-5,0),
                          yend   = c(0,max(f1.4_data$pct1921,na.rm = T)*100+5,0,max(f1.4_data$pct1921,na.rm = T)*100+5),
                          colors = factor(1:4))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 1.5 ----
# (Trends in case notifications of people newly diagnosed with TB, 30 high TB burden countries, 2016-2021)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
f1.5_data <-
  filter(segmentation, iso3 %in% list_hbcs_plus_wl$iso3)  %>%
  select(iso3,country,g_whoregion,pct1920,pct1921) %>%
  pivot_longer(cols = 'pct1920':'pct1921', names_to = "shortfall", values_to = "value") %>%
  mutate(shortfall=ifelse(shortfall=='pct1920',"pct2020","pct2021"),
         value=value+1)

f1.5_sel_order <-
  f1.5_data %>%
  filter(shortfall == "pct2020") %>%
  arrange(value) %>%
  mutate(country = factor(country))

f1.5_txt <- filter(f1.5_sel_order, value<0.8)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 1.6 ----
# (Trends in case notifications of people newly diagnosed with TB, 30 high TB burden countries, 2016-2021)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# pct_change <- notification %>%
#   select(iso3,year,country,g_whoregion,c_newinc) %>%
#   pivot_wider(names_from = year, values_from = c_newinc) %>%
#   mutate(`2019_2` = round(`2019`/`2018`*100-100,2),
#          `2020` = round(`2020`/`2019`*100-100,2),
#          `2021` = round(`2021`/`2019`*100-100,2)) %>%
#   select(iso3,country,g_whoregion,`2019_2`,`2020`,`2021`) %>%
#   rename(`2019`=`2019_2`) %>%
#   pivot_longer(names_to = "year", cols = `2019`:`2021`, values_to = "pct_change") %>%
#   mutate(year=as.numeric(year))

pct_change <- notification %>%
  select(iso3,year,country,g_whoregion,c_newinc) %>%
  pivot_wider(names_from = year, values_from = c_newinc) %>%
  mutate(
         `2016_2` = round(`2016`/`2015`*100-100,2),
         `2017_2` = round(`2017`/`2016`*100-100,2),
         `2018_2` = round(`2018`/`2017`*100-100,2),
         `2019_2` = round(`2019`/`2018`*100-100,2),
         `2020_2` = round(`2020`/`2019`*100-100,2),
         `2021`   = round(`2021`/`2020`*100-100,2)) %>%
  select(iso3,country,g_whoregion,`2016_2`,`2017_2`,`2018_2`,`2019_2`,`2020_2`,`2021`) %>%
  rename(`2016`=`2016_2`,`2017`=`2017_2`,`2018`=`2018_2`,`2019`=`2019_2`,`2020`=`2020_2`) %>%
  pivot_longer(names_to = "year", cols = `2016`:`2021`, values_to = "pct_change") %>%
  mutate(year=as.numeric(year))


f1.6_data <- notification %>%
  select(iso3,year,country,g_whoregion,c_newinc) %>%
  subset(year>=2015) %>%
  filter(iso3 %in% list_hbcs_plus_wl$iso3) %>%
  left_join(segmentation %>% select(iso3,col)) %>%
  mutate(col=ifelse(iso3=="ETH"|iso3=="NAM"|iso3=="PRK"|iso3=="ZAF"|iso3=="CHN"#iso3=="RUS"||iso3=="ZWE"|iso3=="PNG"|iso3=="GAB"
                      ,"f",col)) %>%
  mutate(col=ifelse(iso3=="SLE","b",col)) %>%
  left_join(pct_change %>% select(iso3,year,pct_change)) %>%
  mutate(country = ifelse(country=="Russian Federation","Russian Federation\u1D47",country)) %>%
  mutate(country = ifelse(country=="China","China\u1D43",country)) %>%
  mutate(country = ifelse(country=="Namibia","Namibia\u1D43",country))

f1.6a_sel_order <-
  f1.6_data %>%
  filter(col=="a") %>%
  filter(year == 2020) %>%
  arrange(pct_change) %>%
  mutate(country = factor(country))

f1.6a_txt <- f1.6a_sel_order %>%
  select(pct_change)

f1.6b_sel_order <-
  f1.6_data %>%
  filter(col=="b") %>%
  filter(year == 2020) %>%
  arrange(pct_change) %>%
  mutate(country = factor(country))

f1.6b_txt <- f1.6b_sel_order %>%
  select(pct_change)

f1.6c_sel_order <-
  f1.6_data %>%
  filter(col=="c") %>%
  filter(year == 2020) %>%
  arrange(pct_change) %>%
  mutate(country = factor(country))

f1.6c_txt <- f1.6c_sel_order %>%
  select(pct_change)

f1.6c_txt_RUS_2020 <- f1.6_data %>%
  filter(iso3=="RUS", year==2020)

f1.6c_txt_RUS <- f1.6_data %>%
  filter(iso3=="RUS", year<2020)

f1.6d_sel_order <-
  f1.6_data %>%
  filter(col=="d") %>%
  filter(year == 2021) %>%
  arrange(pct_change) %>%
  mutate(country = factor(country))

f1.6d_txt <-
  f1.6_data %>%
  filter(col=="d") %>%
  filter(year == 2021) %>%
  mutate(pct_change=pct_change*(-1))

f1.6f_txt_CHN_2020 <- f1.6_data %>%
  filter(iso3=="CHN", year==2020)

f1.6f_txt_CHN_2019 <- f1.6_data %>%
  filter(iso3=="CHN", year==2019)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: figure 1.7 ----
# (Trends in monthly case notification in India 2019-2022)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# # Download the provisional notifications data using this code:
# # -------- Start of provisional notifications download code -----
#
# url <- "https://extranet.who.int/tme/generateJSON.asp?ds=c_newinc&iso2=IN"
# 
# json_IND <- fromJSON(readLines(url, warn = FALSE, encoding = 'UTF-8'))
# 
# json_IND_prov <- json_IND$c_newinc_prov %>% rowwise() %>% mutate(sum= sum(c_across(m_01:m_12), na.rm = T))
# 
# # Save to the CSV folder.
# write.csv(json_IND_prov, here(paste0("csv/prov_notifs_IND_", Sys.Date(), '.csv')), row.names = FALSE)
#
# # -------- End of provisional notifications download code -----

# Get India's total notifications for 2019
IND_2019 <- filter(notification, iso3=='IND' & year ==2019) %>% select(c_newinc)

# Read the previously downloaded provisional notifications for India
f1.7_data <- read.csv(here('csv/prov_notifs_IND_2022-10-06.csv')) %>%
  select(year, m_01:m_12) %>%
  pivot_longer(!year, names_to = "month", values_to = "value") %>%
  mutate(ave_2019 = IND_2019$c_newinc/12) %>%
  mutate(month = month.abb[as.numeric(as.factor(month))]) %>%
  mutate(month = factor(month, levels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct",
                                                   "Nov","Dec")))



# # Download the provisional notifications data using this code:
# # -------- Start of provisional notifications download code -----
# #
# url <- "https://extranet.who.int/tme/generateJSON.asp?ds=c_newinc&iso2=ZM"
# 
# json_ZMB <- fromJSON(readLines(url, warn = FALSE, encoding = 'UTF-8'))
# 
# json_ZMB_prov <- json_ZMB$c_newinc_prov %>% rowwise() %>% mutate(sum= sum(c_across(m_01:m_12), na.rm = T))
# 
# # Save to the CSV folder.
# write.csv(json_ZMB_prov, here(paste0("csv/prov_notifs_ZMB_", Sys.Date(), '.csv')), row.names = FALSE)
#
# # -------- End of provisional notifications download code -----

# Get India's total notifications for 2019
ZMB_2019 <- filter(notification, iso3=='ZMB' & year ==2019) %>% select(c_newinc)

# Read the previously downloaded provisional notifications for India
ft_inn_data <- read.csv(here('csv/prov_notifs_ZMB_2022-10-24.csv')) %>%
  select(year, m_01:m_12) %>%
  pivot_longer(!year, names_to = "month", values_to = "value") %>%
  mutate(ave_2019 = ZMB_2019$c_newinc/12) %>%
  mutate(month = month.abb[as.numeric(as.factor(month))]) %>%
  mutate(month = factor(month, levels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct",
                                        "Nov","Dec")))



report_year <- 2022
# And the latest year for which the available data are being displayed in graphics
latest_year <- 2021
# Kill any attempt at using factors, unless we explicitly want them!
options(stringsAsFactors=FALSE)

csv_folder_name = "Outputs_R/CSV"

# Load needed data objects (light summary data)

Fig5.1 <- readRDS("./ch5_data/Fig5_1.rds")

Fig5.2 <- readRDS("./ch5_data/Fig5_2.rds")
countries <- Fig5.2  %>% filter( g_income!= "HIC") %>% dim() 
countries <- countries[1]

burden_tot <- Fig5.2  %>% group_by(var) %>%  
  select(c_notified) %>%  summarise_all(sum , na.rm = T)
# Percentage of all notified, that are accounted for in the current gtb report
burden_tot$prop <- NA
burden_tot$prop  <- round(100*prop.table(burden_tot[,2]), digits = 0) %>% unlist()
burden_inc <- burden_tot$prop[1]

Fig5.3 <- readRDS("./ch5_data/Fig5_3.rds")

Fig5.4 <- readRDS("./ch5_data/Fig5_4.rds")
Fig5.4 <- Fig5.4 %>%
  filter(year<report_year)
Fig5.5 <- readRDS("./ch5_data/Fig5_5.rds")

Fig5.6 <- readRDS("./ch5_data/Fig5_6.rds")
Fig5.7 <- readRDS("./ch5_data/Fig5_7.rds")
Fig5.8 <- readRDS("./ch5_data/Fig5_8.rds")

Fig5.9 <- readRDS("./ch5_data/Fig5_9.rds")
# Fig5.9_inc_grp <- readRDS(here::here(outputs_folder_name,"Fig5_9_incgrp.rds")) 

# Total gap in financing in latest year referenced
latest_year_budgetgap_bn = Fig5.9 %>% ungroup() %>% 
  filter(gap_tot > 0 & !is.na(gap_tot)) %>% 
  filter(year == latest_year) %>% 
  summarise(sum(gap_tot)/1000) %>% 
  unlist() %>% 
  round(1)

# how many countries report having gaps/ shortfalls in financing?
# If gap (or any value) is less than zero, it is zeroed
countries_with_gaps <- Fig5.9 %>% ungroup() %>% 
  filter(gap_tot > 0 & !is.na(gap_tot)) %>% 
  filter(year == latest_year) %>% 
  summarise_at(vars(gap_tot), length) %>% 
  unlist()

# What are the top 5 country shortfalls in latest_year?
table_country_gaps <- Fig5.9 %>%
  filter(gap_tot > 0 & !is.na(gap_tot)) %>% 
  filter(year == latest_year ) %>% 
  group_by(year) %>% 
  arrange(desc(gap_tot)) %>% 
  select(country, gap_tot) %>% 
  slice(1:5) %>% 
  mutate(gap_tot  = round(gap_tot,0))

# Among LICs: How many countries? How much gap?
latest_year_budgetgap_lic <- Fig5.9 %>% ungroup() %>% 
  filter(gap_tot > 0 & !is.na(gap_tot)) %>% 
  filter(g_income == "LIC" & year == latest_year) %>% 
  summarise(sum(gap_tot)) %>% 
  unlist() %>% 
  round(0)

lics_with_gaps <- Fig5.9 %>% ungroup() %>% 
  filter(gap_tot > 0 & !is.na(gap_tot)) %>% 
  filter(g_income == "LIC" & year == latest_year) %>% 
  summarise_at(vars(gap_tot), length) %>% 
  unlist()

Fig5_10_1 <- readRDS("./ch5_data/Fig5_10_1.rds")
Fig5_10_1_data <- Fig5_10_1$data
Fig5_10_1_data$country <- trimws(Fig5_10_1_data$country)
Fig5_10_2 <- readRDS("./ch5_data/Fig5_10_2.rds")
Fig5_10_2_data <- Fig5_10_2$data
Fig5_10_2_data$country <- trimws(Fig5_10_2_data$country)
Fig5_10_3 <- readRDS("./ch5_data/Fig5_10_3.rds")
Fig5_10_3_data <- Fig5_10_3$data
Fig5_10_3_data$country <- trimws(Fig5_10_3_data$country)

Fig5.11p <- readRDS("./ch5_data/Fig5.11p.rds")
Fig5.12p <- readRDS("./ch5_data/Fig5.12p.rds")

dstb_cpp_no <- read.csv("./ch5_data/CSV/fig5.11_data.csv")
mdr_cpp_no <- read.csv("./ch5_data/CSV/fig5.12_data.csv") %>% dim()
dstb_cpp_no <- dstb_cpp_no[1]
mdr_cpp_no <- mdr_cpp_no[1]

# Median provider cost per case notified (DSTB)
c_pp_dots <- Fig5.11p$data %>% ungroup() %>% summarise(median(c_pp_dots)) %>% round(0) %>%  unlist()
# Median provider cost per case notified (DSTB)
c_pp_mdr <- Fig5.12p$data %>% ungroup() %>% summarise(median(c_pp_mdr)) %>% round(0) %>% unlist()



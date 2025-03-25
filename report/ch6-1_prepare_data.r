# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data preparation script for ch6-1.rmd
# Takuya Yamanaka, August 2022
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Load chapter 6 packages, settings and data
source(here::here('report/ch6_load_data.r'))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: fig 6.1.1 ----
# (Trends in the UHC service coverage index in WHO regions and World Bank income groups, 2000–2019)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

palatte_fig6.1.1a = c("#FC1C1E","#8FD314","#37F546","#FF8C1F","#4ABAFC","#CC0ED5","#124CFC") 
# AFRO 0, 89, 88, 1   #FC1C1E 
# AMRO 37, 7, 91, 11 #8FD314 
# EURO 0, 45, 88, 0 #FF8C1F 
# EMRO 78, 2, 72, 2 #37F546 
# SEARO 71, 27, 1, 0 #4ABAFC 
# WPRO 13, 94, 9, 8 #CC0ED5 
# Global 93, 70, 0, 1 #124CFC 

f6.1.1a_data <-
  read.csv(here::here("./csv/sci.csv")) %>%
  dplyr::rename(who_reg=7,year=10,SCI=FactValueNumeric) %>%
  dplyr::select(who_reg,year,SCI) %>%
  mutate(who_reg=factor(recode(who_reg,"AFR"="African Region",
                               "AMR"="Region of the Americas",
                               "EMR"="Eastern Mediterranean Region",
                               "EUR"="European Region",
                               "SEAR"="South-East Asia Region",
                               "WPR"="Western Pacific Region",
                               "GLOBAL"="GLOBAL"),
                        levels=c("African Region",
                                 "Region of the Americas",
                                 "Eastern Mediterranean Region",
                                 "European Region",
                                 "South-East Asia Region",
                                 "Western Pacific Region","GLOBAL")))
# (a) WHO regions
f6.1.1a_data <- 
  f6.1.1a_data %>%
  mutate(year=as.numeric(year)) %>% 
  arrange(who_reg) 

## text!
f6.1.1a_txt <- f6.1.1a_data %>%
  filter(who_reg == "GLOBAL", year == 2000|year == report_year-3) %>%
  pivot_wider(names_from = "year", values_from = "SCI") %>%
  rename(sci_2019 = `2019`, sci_2000 = `2000`)

# (b) WB income classifications
palatte_fig6.1.1b = c("#66D6FF","#E63E13","#8FD314","#814550","#124CFC") 
#                
# High income 60, 16, 0, 0 #66D6FF 
# UMI 3, 74, 92, 7 #E63E13 
# LMI 37, 7, 91, 11 #8FD314 
# LI 42, 69, 64, 13 #814550 
# Global 93, 70, 0, 1 #124CFC 

f6.1.1b_data <- 
  read.csv(here::here("./csv/sci_wb.csv")) %>%
  dplyr::rename(income_group=8,year=10,SCI=FactValueNumeric) %>%
  dplyr::select(income_group,year,SCI) 

f6.1.1b_data <-
  subset(f6.1.1a_data,who_reg=="GLOBAL") %>%
  dplyr::rename(income_group=1) %>% rbind.data.frame(.,f6.1.1b_data) %>%
  mutate(income_group=factor(
    income_group,
    levels=c("High-income","Upper-middle-income","Lower-middle-income",
             "Low-income","GLOBAL"))) 

f6.1.1b_data <- 
  f6.1.1b_data %>%
  mutate(year=as.numeric(year)) %>% 
  arrange(income_group) 

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: fig 6.1.2 ----
# (UHC service coverage index by country, 2019)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

palatte_fig6.1.2 = c("#B11219","#ED1D24","#D2745B","#D98D93","#FDCEB8","#F8E5CF") 

f6.1.2_data <- 
  sdg %>% filter(indicator_id == "UHC_INDEX_REPORTED", sex=="a") %>% 
  slice_max(year) %>% 
  select(country,iso3,value,year) %>%
  mutate(var=cut(value,breaks=c(0,40,50,60,70,80,100),right = FALSE,
                 labels=c("<40","40–49","50–59","60–69","70\u201379","\u226580")))  

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: fig 6.1.3 ----
# (Percentage of the general population facing catastrophic health expenditure)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

palatte_fig6.1.3 = c("#D2EEF6","#8BB2CF","#5187B3","#13ABD7","#026FBA","#033270")

f6.1.3_data <- 
  sdg %>% filter(indicator_id == "FINPROTECTION_CATA_TOT_10_POP", sex=="a") %>% 
  group_by(iso3) %>%
  slice(which.max(year)) %>%  ## select latest available year
  ungroup() %>% 
  select(country,iso3,year,value)  %>% 
  mutate(var=cut(value,breaks=c(0,3,6,9,12,15,100),right = FALSE,
                 labels=c("<3%","3–5%","6–8%","9–11%","12–14%","\u226515%")))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: fig 6.1.4 ----
# (UHC service coverage index (SDG 3.8.1))
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

label_income <- c(LIC="Low-income",LMC="Lower-middle-income",UMC="Upper-middle-income")

uhc_index <- 
  sdg %>% filter(indicator_id == "UHC_INDEX_REPORTED", sex=="a") %>% 
  group_by(iso3) %>%
  slice(which.max(year)) %>%  ## select latest available year
  ungroup() %>% 
  mutate(value=value) %>% 
  select(iso3,value,year) %>% 
  rename(uhc_index=2,uhc_index_yr=3)

uhc_ce10 <- 
  sdg %>% filter(indicator_id == "FINPROTECTION_CATA_TOT_10_POP", sex=="a") %>% 
  group_by(iso3) %>%
  slice(which.max(year)) %>%  ## select latest available year
  ungroup() %>% 
  select(iso3,value,year) %>% 
  rename(uhc_ce10=2,uhc_ce10_yr=3)

hgf <- hgf %>%
  select(!SDG.3.8.2) %>%
  rename(iso3=1,uhc_index_yr=2,country=3,uhc_index=4,g_income=5) %>%
  filter(!is.na(uhc_index)) %>%
  mutate(g_income=ifelse(g_income=="Lower-middle income","LMC",
                ifelse(g_income=="Upper-middle income","UMC",
                       ifelse(g_income=="Low income","LIC","HIC"))))

f6.1.4_data <- 
  # list_hbcs_plus_wl_income %>% 
  # left_join(uhc_index,by='iso3') %>% 
  hgf %>%
  left_join(uhc_ce10,by='iso3')

f6.1.4_data %>% 
  filter(!is.na(uhc_ce10)) %>% 
  arrange(g_income) 


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: fig 6.1.5 ----
# (Current health expenditure per capita, 30 high TB burden countries, 2000–2019)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

f6.1.5_data <- 
  sdg %>% filter(iso3 %in% iso3_hbc) %>% 
  filter(indicator_id=="GHED_CHE_pc_PPP_SHA2011") %>% 
  arrange(country) %>% 
  left_join(list_hbcs_income) %>% 
  select(iso3,year,value,country,g_income) 


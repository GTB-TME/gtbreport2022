# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data preparation script for ch6-3.rmd
# Takuya Yamanaka, August 2022
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Load chapter 6 packages, settings and data
source(here::here('report/ch6_load_data.r'))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: fig 6.3.1 ----
# (The relationship between GDP per capita and the prevalence of undernutrition, and TB incidence per 100 000 population)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

f6.3.1a_data <-
  sdg %>% arrange(-year) %>%
  group_by(iso3,indicator_id) %>% slice(1) %>%
  filter(indicator_id == 'NY.GDP.PCAP.PP.KD') %>%
  left_join(est %>% filter(year==report_year-1) %>% select(iso3,inc)) %>%
  mutate(gdp=value/1e3) %>%
  ungroup() %>%
  select(country,iso3,gdp,inc) 

# f6.3.1a2_data <-
#   read_csv(here::here("report/ch6_data/gdp.csv")) %>%
#   pivot_longer(cols = `1960`:`2021`, names_to = "year", values_to = "value") %>%
#   rename(country = 1, iso3 = 2, indicator_id = 4) %>%
#   select(country, iso3, indicator_id, year, value) %>%
#   filter(!is.na(value)) %>%
#   arrange(-as.numeric(year)) %>%
#   group_by(iso3,indicator_id) %>% slice(1) %>%
#   left_join(est %>% filter(year==report_year-1) %>% select(iso3,inc)) %>%
#   mutate(gdp=value/1e3) %>%
#   ungroup() %>%
#   select(country,iso3,gdp,inc) %T>%
#   write_csv(paste0("./ch6_data/out/f6.3.1a_data_",as.character(Sys.Date()),".csv"))
  

f6.3.1b_data <- 
  sdg %>% arrange(-year) %>% 
  group_by(iso3,indicator_id) %>% slice(1) %>% 
  filter(indicator_id == 'SN.ITK.DEFC.ZS') %>%
  left_join(est %>% filter(year==report_year-1) %>% select(iso3,inc)) %>% 
  mutate(nut=value) %>% 
  ungroup() %>% 
  select(country,iso3,nut,inc) 

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: fig 6.3.2 ----
# (Global estimates of the number of TB cases attributable to selected risk factors)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

f6.3.2_data <- rf_global %>%
  filter(group_type=="global",sex=="a",risk_factor!="all") %>%
  mutate(risk_factor=factor(risk_factor, labels=c("Alcohol use disorders","Diabetes","HIV infection","Smoking","Undernourishment")))%>%
  mutate(risk_factor=fct_rev(risk_factor)) 


f6.3.2_data$risk_factor <- forcats::fct_reorder(f6.3.2_data$risk_factor,f6.3.2_data$best,min) 
  
f6.3.2_txt <- f6.3.2_data %>%
  select(group_type, risk_factor, best) %>%
  pivot_wider(names_from = risk_factor, values_from = best) %>%
  rename(alcohol = `Alcohol use disorders`, hiv = `HIV infection`)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: fig 6.3.3 ----
# (Estimated number of TB cases attributable to five risk factors per 100 000 population)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# merge pop to att

att2 <- rf_country %>%
  subset((age=="15+"|age=="18+"|age=="a")&sex=="a")

# pop %>%
#   subset(year==report_year-1) %>%
#   select(iso3,g_whoregion,e_pop_15plus,e_pop_num) %>% right_join(att2,by="iso3") -> att2

db_estimates_country %>%
  filter((age_group=="15plus"|age_group=="a")&sex=="a") %>%
  select(iso3, age_group, best) %>%
  pivot_wider(names_from = age_group, values_from = best) %>%
  rename(adult=2,all=3) %>%
  right_join(att2,by="iso3") -> att2

att2 <- att2 %>%
  mutate(paf_pct=ifelse(age=="a",best/all*100,best/adult*100)) %>%
  right_join(list_iso3_country,by="iso3")

f6.3.3_data <- 
  att2 %>% #subset(sex=="a"&risk_factor!="all") %>% 
  # slice_max(year) %>% 
  select(country,iso3,year,age,risk_factor,best,paf_pct) 

#### subset for malnutrition
palatte_fig6.3.3a = c(#"#FFEDA0",
  "#FED976","#FEB24C","#FD8D3C",#"#FC4E2A",
  "#E31A1C") 

f6.3.3a_data <- 
  f6.3.3_data %>% subset(risk_factor=="undernutrition") %>%
  mutate(var=cut(paf_pct,breaks=c(0,10,15,20,Inf),right = FALSE,labels=c("<10","10\u201314","15\u201319","\u226520"))
  )


#### subset for alcohol use
palatte_fig6.3.3b = c(#"#ECE7F2",
  "#D0D1E6","#A6BDDB","#74A9CF",#"#3690C0",
  "#0570B0") 

f6.3.3b_data <- 
  f6.3.3_data %>% subset(risk_factor=="alcohol") %>%
  mutate(var=cut(paf_pct,breaks=c(0,5,10,15,Inf),right = FALSE,labels=c("<5","5\u20139","10\u201314","\u226515"))
)

### subset for HIV
palatte_fig6.3.3c = c(#"#F7FCB9",
  "#D9F0A3","#ADDD8E","#78C679",#"#41AB5D",
  "#238443") 

f6.3.3c_data <- 
  f6.3.3_data %>% subset(risk_factor=="hiv") %>%
  mutate(var=cut(paf_pct,breaks=c(0,1,5,10,Inf),right = FALSE,labels=c("<1","1\u20134","5\u20139","\u226510"))
)

### subset for Diabetes
palatte_fig6.3.3d = c(#"#FDE0DD",
  "#FCC5C0","#FA9FB5","#F768A1",#"#DD3497",
  "#AE017E") 
f6.3.3d_data <- 
  f6.3.3_data %>% subset(risk_factor=="diabetes") %>%
  mutate(var=cut(paf_pct,breaks=c(0,3,4,5,Inf),right = FALSE,labels=c("<3","3\u20134","4\u20135","\u22655"))
  )

# subset for smoking
palatte_fig6.3.3e = c(#"#EFEDF5",
  "#DADAEB","#BCBDDC","#9E9AC8",#"#807DBA",
  "#6A51A3") 
f6.3.3e_data <- 
  f6.3.3_data %>% subset(risk_factor=="smoking") %>%
  mutate(var=cut(paf_pct,breaks=c(0,5,10,15,Inf),right = FALSE,labels=c("<5","5\u20139","10\u201314","\u226515"))
  )


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: fig 6.3.4 ----
# (Estimated number of TB cases attributable to five risk factors, 30 high TB burden countries and three global TB watchlist countries)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

palatte_f6.3.4 = c("#084EA2","#0491D1","#ED1D24","#B92270","#91A93E") 
# Blue (alcohol), Skyblue (diabetes), Red (HIV), purple (smoking), green (undernourishment)

#' # Risk Ratios (mean, sd)
#'
# rr.alc <-
#   c(3.33, (2.14 - 5.19) / 3.94) # alcohol, Eur Resp Dis 2017
# rr.dia <-
#   c(1.5, (1.76 - 1.28) / 3.92) # diabetes, Trop Med Int Health 2018
# rr.und <- c(3.2, 0.2 / 3.92) # under-nourishment, GTB 2018
# rr.smk <- c(1.57, (2.1 - 1.18) / 3.92) # smoking, GTB 2019
# #HIV globally: mean 19.2; SD 1.5
# rr.hiv <- c(19.2,1.5)

f6.3.4_data <- 
  rf_country %>% subset(sex=="a"&risk_factor!="all") %>% 
  right_join(list_iso3_country,by="iso3") %>%
  filter(iso3 %in% iso3_hbc_plus_wl, sex=='a') %>%
  select(country,iso3,year,age,risk_factor,best,lo,hi) %>%
  # mutate(across(6:8, ~./1000)) %>% 
  mutate(riskgrp=factor(risk_factor,levels=c('alcohol','diabetes','hiv','smoking','undernutrition'))) %>% 
  left_join(list_hbcs_plus_wl) %>% 
  mutate(risk_factor=fct_recode(risk_factor,'Alcohol use disorders'='alcohol',Smoking='smoking',Diabetes='diabetes',
                                HIV='hiv',Undernourishment='undernutrition')) %>%
  arrange(country,risk_factor)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: fig 6.3.5 ----
# (Status of selected SDG indicators beyond SDG 3 by country, latest available year)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Maps for SDG indicators 
sdgind_radar <- c("EG.CFT.ACCS.ZS","SN.ITK.DEFC.ZS","EN_LND_SLUM","per_allsp.cov_pop_tot",
                  "SI_POV_DAY1","SI.POV.GINI")

f6.3.5_data <- 
  sdg %>% 
  arrange(-year) %>% 
  group_by(iso3,indicator_id) %>% slice(1) %>% 
  filter(indicator_id %in% sdgind_radar) %>%
  pivot_wider(id_cols = c(iso3,country), names_from = indicator_id, values_from = value) %>% 
  mutate("Clean fuels"=EG.CFT.ACCS.ZS, "Undernourishment"=SN.ITK.DEFC.ZS, 
         "Living in slums"=EN_LND_SLUM, "Social protection"=per_allsp.cov_pop_tot, 
         "Living in poverty"=SI_POV_DAY1, "Income inequality"=SI.POV.GINI) %>% 
  select(iso3,country,'Clean fuels':'Income inequality') %>% 
  arrange(country)

f6.3.5_data <- f6.3.5_data %>% 
  mutate(country = ifelse(iso3=='PRK',"Democratic People's\nRepublic of Korea",country)) %>% 
  mutate(country = ifelse(iso3=='COD',"Democratic Republic\nof the Congo",country))

#### subset for clean fuels
f6.3.5a_data <- 
  f6.3.5_data %>% 
  mutate(var=cut(`Clean fuels`,breaks=c(0,30,60,80,Inf),right = FALSE,
                 labels=c("<30","30–59","60–79","\u226580")))  

#### subset for income equality
f6.3.5b_data <- 
  f6.3.5_data %>% 
  mutate(var=cut(`Income inequality`,breaks=c(0,32,35,42,Inf),right = FALSE,
                 labels=c("<32","32–34","35–41","\u226542")))  

### subset for living in poverty
f6.3.5c_data <- 
  f6.3.5_data %>% 
  mutate(var=cut(`Living in poverty`,breaks=c(0,1,2,20,Inf),right = FALSE,
                 labels=c("<1","1–1.9","2–19","\u226520"))) 

# subset for social protection
f6.3.5d_data <- 
  f6.3.5_data %>% 
  mutate(var=cut(`Social protection`,breaks=c(0,20,40,60,Inf),right = FALSE,
                 labels=c("<20","20–39","40–59","\u226560"))) 

# subset for living in slums
f6.3.5e_data <- 
  f6.3.5_data %>% 
  mutate(var=cut(`Living in slums`,breaks=c(0,15,30,50,Inf),right = FALSE,
                 labels=c("<15","15-29","30–49","\u226550")))  

# subset for malnutrition
f6.3.5f_data <- 
  f6.3.5_data %>% 
  mutate(var=cut(`Undernourishment`,breaks=c(0,3,5,15,Inf),right = FALSE,
                 labels=c("<3","3-4.9","5–14","\u226515")))  


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: fig 6.3.6 ----
# (Status of selected SDG indicators beyond SDG 3 in 30 high TB burden and 3 global TB watchlist countries, latest available year)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

f6.3.6_data <- 
  sdg %>% filter(iso3 %in% iso3_hbc_plus_wl) %>% 
  arrange(-year) %>% 
  group_by(iso3,indicator_id) %>% slice(1) %>% 
  filter(indicator_id %in% sdgind_radar) %>%
  pivot_wider(id_cols = c(iso3,country), names_from = indicator_id, values_from = value) %>%
  mutate("Clean fuels"=EG.CFT.ACCS.ZS, "Nutrition"=100-SN.ITK.DEFC.ZS, 
         "Not in slums"=100-EN_LND_SLUM, "Social protection"=per_allsp.cov_pop_tot, 
         "Not in poverty"=100-SI_POV_DAY1, "Income equality"=100-SI.POV.GINI) %>% 
  select(iso3,country,'Clean fuels':'Income equality') %>% 
  arrange(country) %>%
  pivot_longer(cols = 'Clean fuels':'Income equality', names_to = "sdg", values_to = "value") 

f6.3.6_data <- f6.3.6_data %>% 
  mutate(country = ifelse(iso3=='PRK',"Democratic People's\nRepublic of Korea",country)) %>% 
  mutate(country = ifelse(iso3=='COD',"Democratic Republic\nof the Congo",country))

f6.3.6_data <- f6.3.6_data %>%
  mutate(value=as.integer(value))



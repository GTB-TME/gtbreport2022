# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data preparation script for ch6-2.rmd
# Takuya Yamanaka, August 2022
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Load chapter 6 packages, settings and data
source(here::here('report/ch6_load_data.r'))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: fig 6.2.1 ----
# (National surveys of costs faced by TB patients and their households since 2015: progress and plans)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

palatte_fig6.2.1 = c("#8FC63F","#0491D1","#8FD8F8")

#! run by TY to generate .csv !#

# f6.2.1_data <- 
#   read_excel(here::here('report/ch6_data_raw/Quarterly update table PCS July 2022.xlsx'),sheet="for map") %>%
#   filter(var!=0) %>% 
#   mutate(var = factor(var, levels=c("Completed","Ongoing","Planned"))) %>%
#   #  left_join(iso3_country) %>% 
#   #  select(country, iso3, var) %T>%
#   select(iso3, var) %T>%
#   write_csv(here::here('csv/f6.2.1_data.csv'))

f6.2.1_data <- read.csv(here::here('csv/f6.2.1_data.csv'))

num <- table(f6.2.1_data$var)
labs <- c(
  paste0("Completed (n=",num[1],")"),
  paste0("Ongoing (n=",num[2],")"),
  paste0("Planned (n=",num[3],")"))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: fig 6.2.3 ----
# (Selected baseline results from national surveys^a^ of costs faced by TB patients and their households, latest year   )
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

palatte_f6.2.3 = c("#084EA2","#00225C","#ED1D24","#8C304C") # blue, dark blue, red, dark red

# PCS book
# f6.2.3_data <- read_excel(here::here('report/ch6_data_raw/TBPCS_ALL_database_v23.xlsx')) %>% #from PCS book
#   subset((variable=="cc1"|variable=="cc2"|variable=="cc3")&
#            quintile=="all"&sex=="all") %>%
#   left_join(list_iso3_country) %>%
#   dplyr::select(iso3,country,year,dr,value,lci,uci) %>%
#   dplyr::rename(grp=dr,catast_pct=value,catast_pct_lo=lci,catast_pct_hi=uci) %>%
#   mutate_at(c("catast_pct","catast_pct_lo","catast_pct_hi"), ~.*100)
# 
# # country's reports
# f6.2.3b_data <- read_excel(here::here('report/ch6_data_raw/catast indicator summary 2022_v1.xlsx')) %>% # reported data from countries
#   select(-country) %>% left_join(list_iso3_country) %>%
#   select(iso3, country, year,grp,catast_pct,catast_pct_lo,catast_pct_hi) %>%
#   subset(iso3=="BEN"|iso3=="LSO"|iso3=="MDA"|iso3=="SLV"|iso3=="TLS"|iso3=="COL"|iso3=="NER"|iso3=="ZAF") %>% #Colombia/Niger/South Africa added in 2022
#   mutate(grp=factor(grp,levels=c('all','ds','dr'))) %>%
#   mutate(grp=factor(grp,labels=c('overall','TB (first-line treatment)','Drug-resistant TB'))) %>%
#   filter(iso3!="MDA") %>%
#   arrange(country)
# 
# f6.2.3_data <- f6.2.3_data %>%
#   rbind.data.frame(.,f6.2.3b_data) %>%
#   mutate(catast_pct_lo = ifelse(iso3=="MLI"&grp=="Drug-resistant TB", 80.73638, catast_pct_lo),
#          catast_pct_lo = ifelse(iso3=="UGA"&grp=="Drug-resistant TB", 92.01535, catast_pct_lo),
#          catast_pct_hi = ifelse(iso3=="UGA"&grp=="Drug-resistant TB", 100.0, catast_pct_hi)) %T>%
#   write_csv(paste0(here::here('csv/f6.2.3_data.csv')))  # write csv
# 
# write csv for GTB Database (pass on to Hazim)
# f6.2.3_data %>%
#   mutate(across(5:7, round, 0)) %>%
#   write_csv(paste0(here::here('report/ch6_data_raw/'),"/_GTBD_tbpcs_catacost_gtr",report_year,".csv"))   # write csv

f6.2.3_data <- read.csv(here::here('csv/f6.2.3_data.csv'))

## Subset data for All TB to estimate pooled average
tb <- tb %>%
  mutate(c.notified=c_newinc) %>%
  mutate(c_ds=ifelse(year<2020,c.notified-conf_rrmdr,c.notified-(conf_rr_nfqr + conf_rr_fqr))) %>%
  mutate(conf_rrmdr=ifelse(year<2020,conf_rrmdr,conf_rr_nfqr + conf_rr_fqr))

f6.2.3a_data <- f6.2.3_data %>% 
  filter(grp=='overall') %>% 
  mutate(year=ifelse(year>report_year-1,report_year-1,year)) %>% # to be changed to 2021
  mutate(country = ifelse(country=="Niger","Niger\u1D9C",country))

tb %>% select(iso3,year,c.notified) %>% right_join(f6.2.3a_data,by=c('iso3','year')) -> f6.2.3a_data

fit_all <-
  rma(
    yi = catast_pct,
    sei = (catast_pct_hi - catast_pct_lo)/3.92,
    data = f6.2.3a_data, 
    weights = c.notified
  )

f6.2.3a_data <- f6.2.3a_data %>% 
  add_row(iso3="AVE",country="Pooled average", 
          grp="ave",
          catast_pct    = as.numeric(fit_all$b),
          catast_pct_lo = fit_all$ci.lb,
          catast_pct_hi = fit_all$ci.ub) #%>% 
  # mutate(country=factor(country,levels=rev(country))) #%T>%  # factorize in the order of rows 
  # write_csv(paste0("./ch6_data/out/f6.2.3_data_all-TB_",as.character(Sys.Date()),".csv"))   

f6.2.3_sel_order <- 
  f6.2.3a_data %>% 
  arrange(catast_pct) %>% 
  arrange(grp) %>%
  mutate(country_a = ifelse(iso3=="NER","Niger",country)) %>% 
  mutate(country = factor(country),
         country_a = factor(country_a)) 

f6.2.3a_data <- f6.2.3a_data %>% 
  mutate(country = ifelse(iso3=="NER","Niger",country))

## Subset data for DS-TB 
f6.2.3b_data <- f6.2.3_data %>% 
  filter(grp=='TB (first-line treatment)'|iso3=="TLS"|iso3=="SLV"|iso3=="FJI"|iso3=="SLB") %>% 
  mutate(year=ifelse(year>2021,2021,year))  %>%
  mutate(country = ifelse(country=="Niger","Niger\u1D9C",country))

tb %>% select(iso3,year,c_ds) %>% right_join(f6.2.3b_data,by=c('iso3','year')) -> f6.2.3b_data

fit_all <-
  rma(
    yi = catast_pct,
    sei = (catast_pct_hi - catast_pct_lo)/3.92,
    data = f6.2.3b_data, 
    weights = c_ds
  )

f6.2.3b_data <- f6.2.3b_data %>% 
  add_row(iso3="AVE",country="Pooled average", 
          grp="ave",
          catast_pct    = as.numeric(fit_all$b),
          catast_pct_lo = fit_all$ci.lb,
          catast_pct_hi = fit_all$ci.ub) %>% 
  mutate(country=factor(country,levels=rev(country))) 

f6.2.3b_data <- f6.2.3b_data %>% 
  mutate(grp=ifelse(grp=="overall","TB (first-line treatment)",grp))


## Subset data for DR-TB 
f6.2.3c_data <- f6.2.3_data %>% 
  filter(grp=='Drug-resistant TB') %>% 
  mutate(year=ifelse(year>report_year-1,report_year-1,year)) %>%
  mutate(country = ifelse(country=="Niger","Niger\u1D9C",country))

tb %>% select(iso3,year,conf_rrmdr) %>% right_join(f6.2.3c_data,by=c('iso3','year')) -> f6.2.3c_data

fit_all <-
  rma(
    yi = catast_pct,
    sei = (catast_pct_hi - catast_pct_lo)/3.92,
    data = f6.2.3c_data, 
    weights = conf_rrmdr   
  )

f6.2.3c_data <- f6.2.3c_data %>% 
  add_row(iso3="AVE",country="Pooled average", 
          grp="ave",
          catast_pct    = as.numeric(fit_all$b),
          catast_pct_lo = fit_all$ci.lb,
          catast_pct_hi = fit_all$ci.ub) %>% 
  mutate(country=factor(country,levels=rev(country))) %>%  # factorize in the order of rows 
  add_row(iso3="SLB",grp="Drug-resistant TB",country="Solomon Islands") %>%
  add_row(iso3="FJI",grp="Drug-resistant TB",country="Fiji") %>%
  add_row(iso3="TLS",grp="Drug-resistant TB",country="Timor-Leste") %>%
  add_row(iso3="SLV",grp="Drug-resistant TB",country="El Salvador") 


# extract pooled averages for texts
f6.2.3a_txt <- f6.2.3a_data %>%
  subset(iso3=="AVE") %>%
  mutate(grp="Overall\n(End TB Strategy indicator)") 

f6.2.3a_txt_lo <- f6.2.3a_data %>%
  arrange(catast_pct) %>%
  slice(1)

f6.2.3a_txt_hi <- f6.2.3a_data %>%
  arrange(desc(catast_pct)) %>%
  slice(1)

f6.2.3a_txt_num <- f6.2.3a_data %>%
  filter(iso3!="AVE")

f6.2.3c_txt <- f6.2.3c_data %>%
  subset(iso3=="AVE") %>%
  mutate(grp="Drug-resistant TB") 

f6.2.3c_txt_num <- f6.2.3c_data %>%
  filter(iso3!="AVE", !is.na(catast_pct))



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data: fig 6.2.4 ----
# (Distribution of costs faced by TB patients and their households in 25 national surveys completed 2016–2021)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
palatte_f6.2.4 = c("goldenrod2","dodgerblue1","darkblue") # green, blue, skyblue

# # from PCS book
# f6.2.4_data <- read_excel(here::here('report/ch6_data/Fig3.3.1.cost drivers_v6_2022-06-14.xlsx')) %>%
#   rename(cat=3) %>%
#   mutate(value=value*100)
# 
# # from country's reports
# f6.2.4_data <- read_excel(here::here('report/ch6_data/catast cost breakdown summary 2022_v1.xlsx'))  %>% 
#   select(-country) %>% left_join(list_iso3_country) %>% 
#   filter(grp=="all") %>% select(-source,-grp) %>% 
#   #  mutate(lab=case_when(
#   #    method_indirect=="hca" ~ paste0(country,"^a"),
#   #    TRUE ~ as.character(country))) %>% 
#   #  mutate(lab=str_replace_all(lab," ","~")) %>% 
#   #  mutate(lab=fct_reorder(lab,pcost_med_pct)) %>% 
#   arrange(pcost_med_pct) %>%
#   mutate(country=fct_reorder(country,pcost_med_pct)) %>%
#   select(iso3,country,pcost_med_pct,pcost_nmed_pct,pcost_indirect_pct) %>%
#   subset(iso3=="BEN"|iso3=="LSO"|iso3=="MDA"|iso3=="SLV"|iso3=="TLS"|iso3=="COL"|iso3=="NER"|iso3=="ZAF") %>%
#   rename(p_med=3,p_nmed=4,p_indirect=5) %>%
#   pivot_longer(p_med:p_indirect,names_to="cat") %>% rbind.data.frame(.,f6.2.4_data) %T>% 
#   write_csv(paste0(here::here('csv/f6.2.4_data.csv')))

f6.2.4_data <- read.csv(here::here('csv/f6.2.4_data.csv'))

f6.2.4_sel_order <- 
  f6.2.4_data %>% 
  filter(cat == "p_med") %>% 
  arrange(value) %>% 
  mutate(country = factor(country))

f6.2.4_txt_num <- f6.2.4_data %>%
  filter(cat == "p_med")

f6.2.4_txt_med <- f6.2.4_data %>%
  filter(cat == "p_med", value > 25) %>%
  arrange(desc(value))

f6.2.4_txt_nmed <- f6.2.4_data %>%
  filter(cat == "p_nmed", value > 50) %>%
  arrange(as.character(country))

f6.2.4_txt_indirect <- f6.2.4_data %>%
  filter(cat == "p_indirect", value > 43.6) %>% # find cutoff indirect > nmed
  arrange(as.character(country))

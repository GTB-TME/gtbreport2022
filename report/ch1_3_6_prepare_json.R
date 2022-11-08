# Load output packages ----
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
library(dplyr)
library(here)
library(jsonlite)

source(here('report/ch1_prepare_data.r'))

# chapter 1: COVID-19 and TB
## fig 1.1
f1.1_data %>% mutate(c_newinc = c_newinc/1e6) %>% select(year, c_newinc) %>% toJSON("rows")

## fig 1.2
f1.2_data %>% mutate(c_newinc = c_newinc/1e6) %>% select(year, entity, c_newinc) %>% toJSON("rows")

## fig 1.3 (a)
f1.3a_data %>% mutate(value = pct_contribute2020*100) %>% select(country, value) %>% toJSON("rows")

## fig 1.3 (b)
f1.3b_data %>% mutate(value = pct_contribute2021*100) %>% select(country, value) %>% toJSON("rows")

## fig 1.4
f1.4_data %>% mutate(x = pct1920*100, y = pct1921*100, size = cn2021, 
                     color = ifelse(hbc==1, "#E41B17", "dodgerblue")) %>% select(country, x, y, size, color) %>% toJSON("rows")

## fig 1.5
f1.5_data %>% mutate(value = value*100) %>% select(country, year=shortfall, value) %>%
  pivot_wider(names_from = year, values_from = value) %>% 
  arrange(pct2020) %>% mutate(pct2019=100) %>% toJSON("rows")

## fig 1.6 
f1.6_data %>% filter(col=="a") %>% select(country, year, c_newinc) %>% toJSON("rows")

f1.6_data %>% filter(col=="b") %>% select(country, year, c_newinc) %>% toJSON("rows")

f1.6_data %>% filter(col=="c") %>% select(country, year, c_newinc) %>% toJSON("rows")

f1.6_data %>% filter(col=="d") %>% select(country, year, c_newinc) %>% toJSON("rows")

f1.6_data %>% filter(col=="e") %>% select(country, year, c_newinc) %>% toJSON("rows")

f1.6_data %>% filter(col=="f") %>% select(country, year, c_newinc) %>% toJSON("rows")

## fig 1.7
f1.7_data  %>% pivot_wider(names_from = year, values_from = value) %>% rename(cn2019=2, cn2020=3, cn2021=4, cn2022=5) %>% mutate(month = month.name[c(1:12)], month_abb = month.abb[c(1:12)]) %>% toJSON("rows")

# chapter 3: case notifications

## section 3.1
source(here('report/ch3-1_prepare_data.r'))

## fig 3.1.1 global
f3.1.1_data %>% filter(entity=="Global") %>% 
  mutate(e_inc_num = e_inc_num/1e6,
         e_inc_num_lo = e_inc_num_lo/1e6,e_inc_num_hi = e_inc_num_hi/1e6,
         c_newinc = c_newinc/1e6) %>% toJSON("rows")

## fig 3.1.1 region
f3.1.1_data %>% filter(entity!="Global") %>%
  mutate(e_inc_num = e_inc_num/1e6,
         e_inc_num_lo = e_inc_num_lo/1e6,e_inc_num_hi = e_inc_num_hi/1e6,
         c_newinc = c_newinc/1e6) %>% toJSON("rows")

## fig 3.1.2 
f3.1.2_data %>%  mutate(e_inc_num = e_inc_num,e_inc_num_lo = e_inc_num_lo,e_inc_num_hi = e_inc_num_hi,c_newinc = c_newinc) %>% toJSON("rows")

## fig 3.1.4
f3.1.4_data %>% pivot_wider(names_from = age_group, values_from = how_many) %>% select(year, c_new_014, c_new_15plus) %>% toJSON("rows")

## fig 3.1.4
f3.1.5_data %>% pivot_wider(names_from = sex, values_from = cnr:inc_100k) %>% mutate(cnr_Female = cnr_Female*(-1), inc_100k_Female = inc_100k_Female*(-1)) %>% mutate(diff_Male = inc_100k_Male-cnr_Male, diff_Female=inc_100k_Female-cnr_Female) %>% arrange(rev(age_group)) %>% toJSON("rows")

## fig 3.1.8
f3.1.8_data %>% select(year, country, private_pcnt) %>% toJSON("rows")

## section 3.2
source(here('report/ch3-2_prepare_data.r'))

## fig 3.2.1 
f3.2.1_data %>% select(year,entity,value=bacconf_pct) %>% toJSON("rows")

## fig 3.2.3 **need to find summarized data to draw boxplot with Kendo UI!
### function to find outliers
is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

f3.2.3_data_outliers <- f3.2.3_data %>% filter(year == report_year-1 & !is.na(income) & bacconf_pct_denominator>=100) %>% mutate(bacconf_pct = bacconf_pct*100) %>% select(year,country,income,bacconf_pct) %>%
  group_by(income) %>% mutate(outlier = ifelse(is_outlier(bacconf_pct), bacconf_pct, as.numeric(NA))) %>% filter(!is.na(outlier))

f3.2.3_data_in <- f3.2.3_data %>% filter(year == report_year-1 & !is.na(income) & bacconf_pct_denominator>=100) %>% mutate(bacconf_pct = bacconf_pct*100) %>% select(year,country,income,bacconf_pct) %>%
  group_by(income) %>% mutate(outlier = ifelse(is_outlier(bacconf_pct), bacconf_pct, as.numeric(NA))) %>% filter(is.na(outlier))

### produce json data without outliers
f3.2.3_data_in %>%
  summarise(lower = min (bacconf_pct), q1 = quantile(bacconf_pct, 0.25), median = quantile(bacconf_pct, 0.5), q3 = quantile(bacconf_pct, 0.75), upper = max(bacconf_pct)) %>%
  mutate(color = c("#BDD7E7","#6BAED6","#3182BD","#08519C")) %>% toJSON("rows") 

### write json data for outliers manually
f3.2.3_data_outliers %>% select(income,outlier) %>% arrange(income)%>% ungroup() %>% select(outlier) %>% toJSON("rows") 

# [{"income":"Low-income","lower":46.079,"q1":65.1627,"median":70.5893,"q3":81.0107,"upper":95.843,"color":"#BDD7E7", "outliers": [38.3]},{"income":"Lower-middle-income","lower":36.3136,"q1":62.6421,"median":77.092,"q3":85.8785,"upper":100,"color":"#6BAED6"},{"income":"Upper-middle-income","lower":55.8426,"q1":69.1873,"median":76.7474,"q3":79.0273,"upper":88.0893,"color":"#2171B5", "outliers": [96.5,98.7,51.1,94.8,50.8]},{"income":"High-income","lower":72.3739,"q1":84.1652,"median":88.2377,"q3":92.656,"upper":96.6102,"color":"#08519C", "outliers": [67.2,58.2,56.1,67.1]}] 

## fig 3.2.4 
f3.2.4_data  %>% toJSON("rows")

## fig 3.2.7
f3.2.7b_data  %>% select(entity,median) %>% toJSON("rows")

## fig 3.2.8
f3.2.8_data %>% toJSON("rows")

## fig 3.2.9
f3.2.9b_data  %>% select(entity,median) %>% toJSON("rows")

## fig 3.2.10
f3.2.10_data  %>% select(year, entity, value=hivstatus_pct) %>% toJSON("rows")

## fig 3.2.12
f3.2.12_data  %>% select(year, entity, value=dst_pcnt) %>% toJSON("rows")

## fig 3.2.14
f3.2.14_data  %>% select(year, entity, value=fqdst_pct) %>% toJSON("rows")


## section 3.3
source(here('report/ch3-3_prepare_data.r'))

## fig 3.3.1 
f3.3.1_data %>% mutate(value = c_cdr, lo = c_cdr_lo, hi = c_cdr_hi) %>% toJSON("rows")

## fig 3.3.2
f3.3.2a_data %>% toJSON("rows")
f3.3.2b_data %>% toJSON("rows")

## fig 3.3.4
f3.3.4_data %>% toJSON("rows")

## fig 3.3.5
f3.3.5_data %>% mutate(value = c_art, lo = c_art_lo, hi = c_art_hi) %>% toJSON("rows")

## fig 3.3.6
f3.3.6_data %>% pivot_wider(names_from = outcome, values_from = value) %>% rename(succ=2, fail=3, died=4, ltfu=5, neval=6) %>% toJSON("rows")

## fig 3.3.7
f3.3.7_data %>% pivot_wider(names_from = outcome, values_from = value) %>% rename(succ=3, fail=4, died=5, ltfu=6, neval=7) %>% arrange(rev(year)) %>% toJSON("rows")

## fig 3.3.8
f3.3.8_data %>% pivot_wider(names_from = outcome, values_from = value) %>% rename(succ=2, fail=3, died=4, ltfu=5, neval=6) %>% toJSON("rows")

## fig 3.3.9
f3.3.9_data  %>% toJSON("rows")


## section 3.4
source(here('report/ch3-4_prepare_data.r'))

## fig 3.4.1 
f3.4.1_data %>% pivot_wider(names_from = age_group, values_from = how_many) %>% rename(age_014=3, age_15plus=4) %>% toJSON("rows")

## fig 3.4.2 
f3.4.2_data %>% select(country, year, rr_detected, rr_treated) %>% toJSON("rows")

## fig 3.4.4 
f3.4.4_data %>% toJSON("rows")

## fig 3.4.5
f3.4.5_data %>% rename(value = c_rr_coverage , lo = c_rr_coverage_lo, hi = c_rr_coverage_hi) %>% toJSON("rows")

## fig 3.4.7
f3.4.7_data %>% pivot_wider(names_from = outcome, values_from = value) %>% arrange(rev(year)) %>% rename(entity = 1, succ=3, fail=4, died=5, ltfu=6, neval=7) %>% toJSON("rows")


# chapter 6: TB determinants

## section 6.1
source(here('report/ch6-1_prepare_data.r'))

## fig 6.1.1 
f6.1.1a_data %>% pivot_wider(names_from = who_reg, values_from = SCI) %>% rename(afro=2, amro=3, emro=4, euro=5, searo=6, wpro=7) %>% arrange(year) %>% toJSON("rows")

f6.1.1b_data %>% pivot_wider(names_from = income_group, values_from = SCI) %>% rename(hic=2, umic=3, lmic=4, lic=5) %>% arrange(year) %>% toJSON("rows")

## fig 6.1.4
f6.1.4_data %>% mutate(size=5) %>% toJSON("rows")

## fig 6.1.5
f6.1.5_data %>% select (year,country,value) %>% add_row(country = "Democratic People's Republic of Korea", year = rep(2000:2019, 1), value = NA) %>% toJSON("rows")


## section 6.2
source(here('report/ch6-2_prepare_data.r'))

## fig 6.2.3
f6.2.3a_data %>% filter(iso3!="AVE") %>% arrange(desc(catast_pct)) %>% add_row() %>% rbind.data.frame(filter(f6.2.3a_data, iso3=="AVE"))  %>% toJSON("rows")

f6.2.3b_data %>% arrange(desc(factor(country, levels = f6.2.3_sel_order$country, ordered = TRUE)))  %>% filter(iso3!="AVE") %>% add_row() %>% rbind.data.frame(filter(f6.2.3b_data, iso3=="AVE"))  %>% toJSON("rows")

## fig 6.2.4
f6.2.4_data %>% pivot_wider(names_from = cat, values_from = value) %>% arrange(desc(p_med)) %>% toJSON("rows")

## section 6.3
source(here('report/ch6-3_prepare_data.r'))

## fig 6.3.2
f6.3.2_data %>% select(risk_factor,best,lo,hi) %>% mutate(color = c("dodgerblue","deeppink","goldenrod","firebrick","green")) %>% toJSON("rows")

## fig 6.3.4
f6.3.4_data %>% select(country,risk_factor,riskgrp,best,lo,hi) %>% mutate(color = ifelse(riskgrp=="alc","#084EA2",ifelse(riskgrp=="dia","#0491D1",ifelse(riskgrp=="hiv","#ED1D24",ifelse(riskgrp=="smk","#B92270","#91A93E"))))) %>% add_row(country = c("Angola","Central African Republic","Gabon","Uganda"), risk_factor = "Smoking") %>% add_row(country = c("Uganda","Zambia","Zimbabwe"), risk_factor = "Undernourishment")  %>% add_row(country = c("Uganda"), risk_factor = "Alcohol use disorders")  %>% add_row(country = c("Uganda"), risk_factor = "Diabetes") %>% arrange(factor(risk_factor, levels = c("Alcohol use disorders","Diabetes","HIV","Smoking","Undernourishment"), ordered = TRUE)) %>% toJSON("rows")

## fig 6.3.6
f6.3.6_data %>% mutate(country = ifelse(iso3=='PRK',"Democratic People's Republic of Korea",country)) %>%  mutate(country = ifelse(iso3=='COD',"Democratic Republic of the Congo",country)) %>% add_column(color = rep(brewer.pal(6, "Paired"), 33))  %>% toJSON("rows")



## section 2.1
source(here('report/ch1_2_dataprep.R'))

## fig 2.1.1 
names(global) <- gsub('[.]', '_', names(global))
global %>% select(year,inc_num,inc_lo_num,inc_hi_num,c_newinc,inc_h_num,inc_h_hi_num,inc_h_lo_num) %>% toJSON("rows")

global %>% select(year,inc,inc_lo,inc_hi,newinc,inc_h,inc_h_hi,inc_h_lo) %>% mutate(milestone = as.numeric(select(filter(global,year==2015),inc))*0.8) %>% toJSON("rows")

## fig 2.1.5 
globsplt %>% pivot_wider(names_from = sex, values_from = inc:pop) %>% mutate(inc_F = inc_F*(-1), newrel_F = newrel_F*(-1)) %>% mutate(diff_M = inc_M-newrel_M, diff_F=inc_F-newrel_F) %>% arrange(rev(age))%>% mutate(age=ifelse(age=="65plus","\u226565",as.character(age))) %>% toJSON("rows")
regsplt %>% select(!pop) %>% pivot_wider(names_from = sex, values_from = inc:newrel) %>% mutate(inc_F = inc_F*(-1), newrel_F = newrel_F*(-1)) %>% mutate(diff_M = inc_M-newrel_M, diff_F=inc_F-newrel_F) %>% arrange(rev(age)) %>% mutate(age=ifelse(age=="65plus","\u226565",as.character(age))) %>% toJSON("rows")

## fig 2.1.7
names(regional) <- gsub('[.]', '_', names(regional))

regional %>% select(region,year,inc,inc_lo,inc_hi,newinc,inc_h,inc_h_hi,inc_h_lo,inc_milestone) %>% toJSON("rows")

## fig 2.1.8
names(hest) <- gsub('[.]', '_', names(hest))

hest %>% select(country,iso3,year,inc,inc_lo,inc_hi,newinc,inc_h,inc_h_hi,inc_h_lo,inc_milestone) %>% toJSON("rows")

## fig 2.1.9
names(dta) <- gsub('[.]', '_', names(dta))
dta %>% select(country,iso3,year,inc,inc_lo,inc_hi,newinc,inc_h,inc_h_hi,inc_h_lo,inc_milestone) %>% toJSON("rows")



## section 2.2
source(here('report/ch1_2_dataprep.R'))

## fig 2.2.1 
names(global) <- gsub('[.]', '_', names(global))
global %>% select(year,mort_num,mort_lo_num,mort_hi_num,mort_h_num,mort_h_hi_num,mort_h_lo_num,mort_nh_num,mort_nh_hi_num,mort_nh_lo_num) %>% mutate(milestone = as.numeric(select(filter(global,year==2015),mort_num))*0.65) %>% toJSON("rows")

global %>% select(year,mort,mort_lo,mort_hi,mort_h,mort_h_hi,mort_h_lo,mort_nh,mort_nh_hi,mort_nh_lo)  %>% toJSON("rows")

## fig 2.2.3
top10  %>% select(cause,deaths,tbhiv) %>% toJSON("rows")

## fig 2.2.4
cod  %>% arrange(desc(n)) %>% toJSON("rows")

## fig 2.2.5
dta %>% select(year,mort_h_num,mort_h_hi_num,mort_h_lo_num,mort_nh_num,mort_nh_hi_num,mort_nh_lo_num,r_num,r_lo_num,r_hi_num) %>% toJSON("rows")

## fig 2.2.6
regional %>% select(region,year,mort_h,mort_h_hi,mort_h_lo,mort_nh,mort_nh_hi,mort_nh_lo)  %>% toJSON("rows")

## fig 2.2.7
regional %>% select(region,year,mort_num,mort_hi_num,mort_lo_num,mort_milestone)  %>% toJSON("rows")

## fig 2.2.9
hest %>% select(country,year,mort_num,mort_hi_num,mort_lo_num,mort_milestone)  %>% toJSON("rows")

## fig 2.2.10
dta %>% select(country,year,mort_num,mort_hi_num,mort_lo_num,mort_milestone)  %>% toJSON("rows")



## section 2.3
source(here('report/ch2-3_prepare_data.r'))

## fig 2.3.1
f2_3_01_data %>% toJSON("rows")

## fig 2.3.1
f2_3_02_data %>% toJSON("rows")

## fig 2.3.3
f2_3_03_data %>% toJSON("rows")

## fig 2.3.4
f2_3_04_data %>% rename(entity=country) %>% toJSON("rows")

## section 5
source(here::here('report/ch5_prepare_data_kendo.r'))

## fig 5.1
Fig5.1 %>%  toJSON("rows")

## fig 5.3 
Fig5.3 %>% pivot_wider(names_from = name, values_from = value) %>% mutate(year=as.numeric(as.character(year))) %>% rename(domestic=2,international=3) %>% add_row(year=c(2014,2022)) %>% mutate(target=13) %>% arrange(year)  %>% toJSON("rows")

## fig 5.4
Fig5.4 %>%  toJSON("rows")

## fig 5.5
Fig5.5 %>%   filter(year <= latest_year) %>% group_by(year, g_brics) %>% summarise_at(vars(DSTB,MDR),sum,na.rm = T ) %>% mutate(g_brics=ifelse(as.character(g_brics)=="High TB burden and global TB watchlist countries outside BRICSᵃ (n=28)","High TB burden and global TB watchlist countries\noutside BRICS\u1D43 (n=28)",as.character(g_brics))) %>% toJSON("rows")

## fig 5.6
Fig5.6 %>% filter(year <= latest_year) %>% group_by(year) %>% summarise_at(vars(int, ext, tot), function(x) sum(x, na.rm = T)) %>% toJSON("rows")

## fig 5.7
Fig5.7 %>%   filter(year <= latest_year) %>% group_by(year,grp) %>% summarise_all(sum, na.rm = T) %>% toJSON("rows")

## fig 5.8
Fig5.8 %>%   filter(name != "gap_tot")  %>% pivot_wider(names_from = name, values_from = value) %>% rename(int=3,ext=4) %>%  toJSON("rows")

## fig 5.9
Fig5.9 %>% filter( year <= latest_year) %>%   group_by(year, grp = g_income) %>%   mutate_at(vars(gap_tot), function(x) ifelse(x < 0, 0 , x)) %>%   summarise_at(vars(gap_tot), sum, na.rm = T) %>% pivot_wider(names_from = grp, values_from = gap_tot) %>%  toJSON("rows")

Fig5.9 %>% filter( year <= latest_year) %>% group_by(year, grp = g_whoregion) %>% mutate_at(vars(gap_tot), function(x) ifelse(x < 0, 0 , x)) %>% summarise_at(vars(gap_tot), sum, na.rm = T) %>% pivot_wider(names_from = grp, values_from = gap_tot) %>%  toJSON("rows")


## fig 5.10
Fig5_10_1_data %>% pivot_wider(names_from = variable,values_from = value) %>% arrange(desc(int_pct)) %>% mutate(country=str_wrap(country, width = 23,indent = 2)) %>% toJSON("rows")
Fig5_10_2_data %>% pivot_wider(names_from = variable,values_from = value) %>% arrange(desc(int_pct)) %>% mutate(country=str_wrap(country, width = 23,indent = 2)) %>% toJSON("rows")
Fig5_10_3_data %>% pivot_wider(names_from = variable,values_from = value) %>% arrange(desc(int_pct)) %>% toJSON("rows")


## MAF
maf1 %>% mutate(year=factor(year, labels=c(2020,2021,2022))) %>% select(!cat2) %>% pivot_wider(names_from = year,values_from = value) %>% rename(pct2020=2,pct2021=3,pct2022=4) %>% mutate(cat=str_wrap(cat, width = 24,indent = 2)) %>% toJSON("rows") 


## CSE
cse2 %>% pivot_wider(names_from = cat2, values_from = as.numeric(value)) %>% mutate(cat1=str_wrap(cat1, width = 34,indent = 2)) %>% arrange(desc(value)) %>% toJSON("rows") 

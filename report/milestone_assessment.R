#Some code to determine which countries and WHO regions
#have met the 2020 End TB milestones (for the current year) compared to 2015.

# Clear the decks ----
rm(list=ls())

library(tidyverse)
library(reshape2)
library(plyr)
library(dplyr)
library(data.table)
library(here)
library(whomap)
library(RColorBrewer)

source(here('report/ch1_2_dataprep.R'))

#by country
dta.inc.2015 <- est[year==2015, .(iso3, country, g.hbc,  inc, year, mort.num, mort.milestone, inc.milestone, g.whoregion)]
dta.inc.2021 <- est[year==2021, .(iso3, country, g.hbc,  inc, year, mort.num, mort.milestone, inc.milestone, g.whoregion)]

dta1 <- rbind(dta.inc.2015, dta.inc.2021)
dta1_wide<-dcast(melt(dta1, id.vars=c("country", "year")), country~variable+year)
dta1_wide <- dta1_wide[, c(1,2,4,6,7,8,9,10,12,14)]
dta1_wide[, 4:9] <- lapply(dta1_wide[, 4:9], as.numeric)

dta1_wide$inc.percent.change<-((dta1_wide$inc_2015-dta1_wide$inc_2021)/dta1_wide$inc_2015)*100
dta1_wide$mort.percent.change<-((dta1_wide$mort.num_2015-dta1_wide$mort.num_2021)/dta1_wide$mort.num_2015)*100

names(dta1_wide)[names(dta1_wide)=="mort.milestone_2015"] <- "mort.milestone"
names(dta1_wide)[names(dta1_wide)=="inc.milestone_2015"] <- "inc.milestone"

names(dta1_wide)[names(dta1_wide)=="iso3_2015"] <- "iso3"
names(dta1_wide)[names(dta1_wide)=="g.hbc_2015"] <- "g.hbc"

dta1_wide$mort.target<-dta1_wide$mort.milestone-dta1_wide$mort.num_2021
dta1_wide$inc.target<-dta1_wide$inc.milestone-dta1_wide$inc_2021

#If inc.met or mort.met = 1, then the 2020 milestone has been met

dta1_wide <-dta1_wide %>%
  mutate(inc.met = case_when(inc.target <0~ 0, #condition 1
                             inc.target >0~ 1, #condition 2
  ))%>%
  mutate(mort.met = case_when(mort.target <0~ 0, #condition 1
                              mort.target >0~ 1, #condition 2
  ))

#amend the path to your choosing
dta1_wide %>% write.csv(here("report/html_drafts/country.csv"))   # write csv

#make map of countries which have met the 2020 milestone

#incidence milestone

dta <- dta1_wide[ .(iso3, inc.met)]
dta$var<-dta$V2
dta$iso3<-dta$country

dta$var<-factor(dta$var, levels = c(0, 1),
                labels = c("Not met", "Met"))

p1<-whomap(dta, colours= c("white","#66CC99"),
          na.col="#999999",
          na.label = "Not applicable",
          water.col = 'white',
          hidef=F)
p1

pdf(file = here("./report/html_drafts/inc_targets.pdf"), width = 10, height = 6)
p1
dev.off()


#mortality milestone

dta <- dta1_wide[ .(iso3, mort.met)]
dta$var<-dta$V2
dta$iso3<-dta$country

dta$var<-factor(dta$var, levels = c(0, 1),
       labels = c("Not met", "Met"))

p2<-whomap(dta, colours= c("white","#9999CC"),
           na.col="#999999",
           na.label = "Not applicable",
           water.col = 'white',
           hidef=F)
p2

pdf(file = here("./report/html_drafts/mort_targets.pdf"), width = 10, height = 6)
p2
dev.off()


#by WHO region

dta2015 <- regional[year==2015, .(year, g.whoregion, inc, mort.num, mort.milestone, inc.milestone)]
dta2021 <- regional[year==2021, .(year, g.whoregion, inc, mort.num, mort.milestone, inc.milestone)]

dta3 <- rbind(dta2015, dta2021)

dta3_wide<-dcast(melt(dta3, id.vars=c("g.whoregion", "year")), g.whoregion~variable+year)
dta3_wide <- dta3_wide[, c(1,2,3,4,5,6,8)]

dta3_wide$inc.percent.change<-((dta3_wide$inc_2015-dta3_wide$inc_2021)/dta3_wide$inc_2015)*100
dta3_wide$mort.percent.change<-((dta3_wide$mort.num_2015-dta3_wide$mort.num_2021)/dta3_wide$mort.num_2015)*100

names(dta3_wide)[names(dta3_wide)=="mort.milestone_2015"] <- "mort.milestone"
names(dta3_wide)[names(dta3_wide)=="inc.milestone_2015"] <- "inc.milestone"

dta3_wide$mort.target<-dta3_wide$mort.milestone-dta3_wide$mort.num_2021
dta3_wide$inc.target<-dta3_wide$inc.milestone-dta3_wide$inc_2021

dta3_wide <- dta3_wide[, c(1,2,3,4,5,7,6,8,9,10,11)]

#If inc.met or mort.met = 1, then the 2020 milestone has been met
dta3_wide <-dta3_wide %>%
  mutate(inc.met = case_when(inc.target <0~ 0, #condition 1
                             inc.target >0~ 1, #condition 2
  ))%>%
  mutate(mort.met = case_when(mort.target <0~ 0, #condition 1
                              mort.target >0~ 1, #condition 2
  ))

#amend the path to your choosing
dta3_wide %>% write.csv(here("report/html_drafts/regional.csv"))   # write csv







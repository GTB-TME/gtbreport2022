## example plots of RRTB outputs

## testing - uncomment & run para to load data+libraries:
rm(list=ls())
library(ggplot2); library(gtbreport); library(scales)
library(data.table); library(here); library(glue)


## incidence data from dboutput folder
load(here('drtb/dboutput/db_dr_country.rda'))
load(here('drtb/dboutput/db_dr_group.rda'))
load(here('drtb/dboutput/drhbc.rda'))
load(here('data/est.rda'))

gh <- function(x) glue(here(x))

## utility for formating scales
absspace <- function(x,...) {             #works
  format(abs(x), ..., big.mark=" ",scientific = FALSE, trim = TRUE)
}

regkey <- data.table(
  name = c('Africa','The Americas','Eastern Mediterranean',
           'Europe','South-East Asia','Western Pacific'),
  group_name=c('AFR','AMR','EMR','EUR','SEA','WPR'))
ckey <- unique(est[,.(iso3,name=country)])


## ======================= SUMMARY PLOTS ============

## ------ global total
GP <- ggplot(db_dr_group[group_name=='global'],
             aes(x=year,y=e_inc_rr_num,
                 ymin=e_inc_rr_num_lo,
                 ymax=e_inc_rr_num_hi
                 ))+
  geom_line()+geom_ribbon(col=NA,alpha=0.25) +
  scale_y_continuous(label=absspace,limits = c(0,NA)) +
  xlab('Year') + ylab('RR-TB incidence')+
  ggtitle('')+
  theme_gtb()

ggsave(GP,file=gh('drtb/reportoutput/global.inc.pdf'),w=5,h=5)
## ggsave(GP,file=gh('drtb/reportoutput/Fig_16_inMasterlist.pdf'),w=5,h=5)



## ------ global proportion
GP <- ggplot(db_dr_group[group_name=='global'],
             aes(x=year,
                 y=prop.rr,
                 ymin=prop.rr.lo,
                 ymax=prop.rr.hi
                 ))+
  geom_line()+geom_ribbon(col=NA,alpha=0.25) +
  scale_y_continuous(label=percent,limits = c(0,NA)) +
  xlab('Year') + ylab('Proportion of all incidence that is RR-TB')+
  ggtitle('')+
  theme_gtb()

ggsave(GP,file=gh('drtb/reportoutput/global.prop.pdf'),w=5,h=5)
## ggsave(GP,file=gh('drtb/reportoutput/Fig_15_inMasterlist.pdf'),w=5,h=5)


## ------ regional total
tmp <- db_dr_group[group_name!='global']
tmp <- merge(tmp,regkey,by='group_name',all.x = TRUE)

GP <- ggplot(tmp,
             aes(x=year,y=e_inc_rr_num,
                 ymin=e_inc_rr_num_lo,
                 ymax=e_inc_rr_num_hi
                 ))+
  geom_line()+geom_ribbon(col=NA,alpha=0.25) +
  scale_y_continuous(label=absspace,limits = c(0,NA)) +
  xlab('Year') + ylab('RR-TB incidence')+
  ggtitle('')+
  facet_wrap(~name)+
  theme_gtb()

ggsave(GP,file=gh('drtb/reportoutput/regional.inc.pdf'),w=5*3,h=5*2)


## ------ regional proportion
tmp <- db_dr_group[group_name!='global']
tmp <- merge(tmp,regkey,by='group_name',all.x = TRUE)

GP <- ggplot(tmp,
             aes(x=year,
                 y=prop.rr,
                 ymin=prop.rr.lo,
                 ymax=prop.rr.hi
                 ))+
  geom_line()+geom_ribbon(col=NA,alpha=0.25) +
  scale_y_continuous(label=percent,limits = c(0,NA)) +
  xlab('Year') + ylab('Proportion of all incidence that is RR-TB')+
  ggtitle('')+
  facet_wrap(~name,scales = 'free_y')+
  theme_gtb()

ggsave(GP,file=gh('drtb/reportoutput/regional.prop.pdf'),w=5*3,h=5*2)


## ------ 30HBC total
tmp <- merge(drhbc,ckey,by='iso3',all.x=TRUE,all.y=FALSE)

GP <- ggplot(tmp,
             aes(x=year,y=e_inc_rr_num,
                 ymin=e_inc_rr_num_lo,
                 ymax=e_inc_rr_num_hi
                 ))+
  geom_line()+geom_ribbon(col=NA,alpha=0.25) +
  scale_y_continuous(label=absspace,limits = c(0,NA)) +
  xlab('Year') + ylab('RR-TB incidence')+
  ggtitle('')+
  facet_wrap(~name,ncol = 5,nrow=6,scales='free_y')+
  theme_gtb()

ggsave(GP,file=gh('drtb/reportoutput/hbc30.inc.pdf'),w=5*5,h=5*6)

## ------ 30HBC proportions

## new
tmp <- merge(drhbc,ckey,by='iso3',all.x=TRUE,all.y=FALSE)

GP <- ggplot(tmp,
             aes(x=year,
                 y=e_rr_prop_new,
                 ymin=e_rr_prop_new_lo,
                 ymax=e_rr_prop_new_hi
                 ))+
  geom_line()+geom_ribbon(col=NA,alpha=0.25) +
  scale_y_continuous(label=percent,limits = c(0,NA)) +
  xlab('Year') + ylab('Proportion of all incidence that is RR-TB')+
  ggtitle('New & relapse patients')+
  facet_wrap(~name,ncol = 5,nrow=6)+ #,scales = 'free_y'
  theme_gtb()

ggsave(GP,file=gh('drtb/reportoutput/hbc30.prop.new.pdf'),w=5*5,h=5*6)


## ret
tmp <- merge(drhbc,ckey,by='iso3',all.x=TRUE,all.y=FALSE)

GP <- ggplot(tmp,
             aes(x=year,
                 y=e_rr_prop_ret,
                 ymin=e_rr_prop_ret_lo,
                 ymax=e_rr_prop_ret_hi
                 ))+
  geom_line()+geom_ribbon(col=NA,alpha=0.25) +
  scale_y_continuous(label=percent,limits = c(0,NA)) +
  xlab('Year') + ylab('Proportion of all incidence that is RR-TB')+
  ggtitle('Retreatment patients')+
  facet_wrap(~name,ncol = 5,nrow=6)+ #,scales = 'free_y'
  theme_gtb()

ggsave(GP,file=gh('drtb/reportoutput/hbc30.prop.ret.pdf'),w=5*5,h=5*6)


## overall
tmp <- merge(drhbc,ckey,by='iso3',all.x=TRUE,all.y=FALSE)

GP <- ggplot(tmp,
             aes(x=year,
                 y=rr.prop,
                 ymin=rr.prop.lo,
                 ymax=rr.prop.hi
                 ))+
  geom_line()+geom_ribbon(col=NA,alpha=0.25) +
  scale_y_continuous(label=percent,limits = c(0,NA)) +
  xlab('Year') + ylab('Proportion of all incidence that is RR-TB')+
  ggtitle('')+
  facet_wrap(~name,ncol = 5,nrow=6,scales = 'free_y')+
  theme_gtb()

ggsave(GP,file=gh('drtb/reportoutput/hbc30.prop.pdf'),w=5*5,h=5*6)




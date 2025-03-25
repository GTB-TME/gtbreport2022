## generating incidence plots & tables
## libraries
rm(list=ls())
library(here)
library(data.table)
library(ggplot2)
library(scales)
library(glue)

## main data
load(file=here('drtb/data/rhofnr.Rdata'))
load(file=here('data/est.rda'))

## utilities
source(here('drtb/R/utils/incdr.R'))
source(here('drtb/R/utils/ftb.R'))
gh <- function(x) glue(here(x))
ssum <- function(x,...) sqrt(sum(x^2,...)) #sqrt sum sqrs
brkt <- function(x,y) paste0(ftb(x),' (',
                             ftb(pmax(0,x-1.96*y)),' - ',
                             ftb(x+1.96*y),')')

## list of models to graph

## dir(path=here('drtb/outdata'),pattern='KO')
rootlst <- c(
  "separate_global_2",
  "separate_global_4",
  "separate_globalH_2",
  "separate_globalH_4",
  "separate_LerouxIntercept_2",
  "separate_LerouxIntercept_4",
  "separate_LerouxInterceptLerouxSlope_2",
  "separate_LerouxInterceptLerouxSlope_4",
  "together_LerouxIntercept_2",
  "together_LerouxIntercept_4",
  "together_LerouxInterceptLerouxSlope_2",
  "together_LerouxInterceptLerouxSlope_4"
)

## ============================================
## ======= each model graphs and tables =======
## ============================================

## NOTE
NARM <- TRUE #dangerous if TRUE as sweeps over input problems

## loop over models
ylimz <- c(2015,est[,max(year)])
TLG <- TLR <- list()
for(modrt in rootlst){
  cat('working on model ',modrt,'................\n')

  ## example proportions estimates
  load(gh('drtb/outdata/KO.{modrt}.Rdata'))

  ## create incidence
  tout <- addIncidence(KO[year>1999])

  ## === global summaries
  ## -- totals
  glob <- tout[,.(inc=sum(inc.num,na.rm = NARM),
                  inc.sd=ssum(inc.num.sd,na.rm = NARM),
                  inc.rr=sum(inc.rr,na.rm = NARM),
                  inc.rr.sd=ssum(inc.rr.sd,na.rm = NARM)),
               by=year]
  glob[,c('inc.rr.lo','inc.rr.hi'):=.(pmax(0,inc.rr-1.96*inc.rr.sd),
                                      inc.rr+1.96*inc.rr.sd)]

  GP <- ggplot(glob,aes(x=year,y=inc.rr,
                        ymin=inc.rr.lo,ymax=inc.rr.hi))+
    geom_line()+geom_ribbon(col=NA,alpha=0.25) +
    scale_y_continuous(label=comma,limits = c(0,NA)) +
    theme_light() +
    xlab('Year') + ylab('RR-TB incidence')+
    ggtitle(modrt)
  ## GP

  ggsave(GP,file=gh('drtb/plots/tot.glob.long.{modrt}.png'),w=5,h=5)

  mxlim <- glob[year>=2015,max(inc.rr.hi)]
  GP <- GP + xlim(ylimz) +
    scale_y_continuous(label=comma,
                       limits = c(0,mxlim))

  ggsave(GP,file=gh('drtb/plots/tot.glob.short.{modrt}.png'),w=5,h=5)


  ## --- proportions
  glob[,prop.rr:=inc.rr/inc]
  glob[,prop.rr.sd:=prop.rr * sqrt((inc.rr.sd/inc.rr)^2+(inc.sd/inc)^2)]
  glob[,c('prop.rr.lo','prop.rr.hi'):=.(pmax(0,prop.rr-1.96*prop.rr.sd),
                                      prop.rr+1.96*prop.rr.sd)]

  GP <- ggplot(glob,aes(x=year,y=prop.rr,
                        ymin=prop.rr.lo,ymax=prop.rr.hi))+
    geom_line()+geom_ribbon(col=NA,alpha=0.25) +
    scale_y_continuous(label=percent,limits = c(0,NA)) +
    theme_light() +
    xlab('Year') + ylab('RR-TB as proportion of all incidence')+
    ggtitle(modrt)
  ## GP

  ggsave(GP,file=gh('drtb/plots/prop.glob.long.{modrt}.png'),w=5,h=5)

  mxlim <- glob[year>=2015,max(prop.rr.hi)]
  GP <- GP + xlim(ylimz) +
    scale_y_continuous(label=percent,
                       limits = c(0,mxlim))

  ggsave(GP,file=gh('drtb/plots/prop.glob.short.{modrt}.png'),w=5,h=5)


  ## === regional  summaries

  ## -- totals
  reg <- tout[,.(inc=sum(inc.num,na.rm = NARM),
                 inc.sd=ssum(inc.num.sd,na.rm = NARM),
                 inc.rr=sum(inc.rr,na.rm = NARM),
                 inc.rr.sd=ssum(inc.rr.sd,na.rm = NARM)),
              by=.(year,g_whoregion)]
  reg[,c('inc.rr.lo','inc.rr.hi'):=.(pmax(0,inc.rr-1.96*inc.rr.sd),
                                      inc.rr+1.96*inc.rr.sd)]

  GP <- ggplot(reg,aes(x=year,y=inc.rr,
                       ymin=inc.rr.lo,
                       ymax=inc.rr.hi))+
    geom_line()+geom_ribbon(col=NA,alpha=0.25) +
    scale_y_continuous(label=comma,limits = c(0,NA)) +
    theme_light() +
    xlab('Year') + ylab('RR-TB incidence')+
    facet_wrap(~g_whoregion)+
    ggtitle(modrt)
  ## GP


  ggsave(GP,file=gh('drtb/plots/tot.reg.long.{modrt}.png'),w=5*3,h=5*2)

  mxlim <- reg[year>=2015,max(inc.rr.hi)]
  GP <- GP + xlim(ylimz) +
    scale_y_continuous(label=comma,
                       limits = c(0,mxlim))

  ggsave(GP,file=gh('drtb/plots/tot.reg.short.{modrt}.png'),w=5*3,h=5*2)

  ## --- proportions
  reg[,prop.rr:=inc.rr/inc]
  reg[,prop.rr.sd:=prop.rr * sqrt((inc.rr.sd/inc.rr)^2+(inc.sd/inc)^2)]
  reg[,c('prop.rr.lo','prop.rr.hi'):=.(pmax(0,prop.rr-1.96*prop.rr.sd),
                                        prop.rr+1.96*prop.rr.sd)]

  GP <- ggplot(reg,aes(x=year,y=prop.rr,
                       ymin=prop.rr.lo,
                       ymax=prop.rr.hi))+
    geom_line()+geom_ribbon(col=NA,alpha=0.25) +
    scale_y_continuous(label=percent,limits = c(0,NA)) +
    theme_light() +
    xlab('Year') + ylab('RR-TB as proportion of all incidence')+
    facet_wrap(~g_whoregion)+
    ggtitle(modrt)
  ## GP

  ggsave(GP,file=gh('drtb/plots/prop.reg.long.{modrt}.png'),w=5*3,h=5*2)

  mxlim <- reg[year>=2015,max(prop.rr.hi)]
  GP <- GP + xlim(ylimz) +
    scale_y_continuous(label=percent,
                       limits = c(0,mxlim))

  ggsave(GP,file=gh('drtb/plots/prop.reg.short.{modrt}.png'),w=5*3,h=5*2)

  ## === tables

  ## --- global
  glob[,inc.rr.prty:=brkt(inc.rr,inc.rr.sd)] #pretty incidence
  glob[,prop.rr.prty:=brkt(1e2*prop.rr,1e2*prop.rr.sd)] #pretty props
  glob[year>=2015]

  tmp <- glob[year>=2015,.(year,
                           `RR-TB incidence`=inc.rr.prty,
                           `RR as % TB incidence`=prop.rr.prty)]
  tmp[,model:=modrt]
  TLG[[modrt]] <- tmp
  fwrite(tmp,file=gh('drtb/plots/tabs.glob.short.{modrt}.csv'))


  ## --- regional
  reg[,inc.rr.prty:=brkt(inc.rr,inc.rr.sd)] #pretty incidence
  reg[,prop.rr.prty:=brkt(1e2*prop.rr,1e2*prop.rr.sd)] #pretty props

  tmp <- reg[year==2020,
             .(year,g_whoregion,
               `RR-TB incidence`=inc.rr.prty,
               `RR as % TB incidence`=prop.rr.prty)][order(g_whoregion)]
  tmp[,model:=modrt]
  TLR[[modrt]] <- tmp
  fwrite(tmp,file=gh('drtb/plots/tabs.reg.recent.{modrt}.csv'))

} #end of model loop

## ============================================
## ======= comparison graphs and tables =======
## ============================================

TLG <- rbindlist(TLG)
TLR <- rbindlist(TLR)


TLG <- TLG[order(year)]
TLR <- TLR[order(g_whoregion)]

fwrite(TLG,file=gh('drtb/outdata/tabs.glob.compare.csv'))
fwrite(TLR,file=gh('drtb/outdata/tabs.reg.compare.csv'))

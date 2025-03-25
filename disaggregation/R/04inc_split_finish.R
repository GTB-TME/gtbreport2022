#' ---
#' title: Generating age/sex incidence, completion stages
#' author: Pete Dodd
#' date: 18 July, 2022
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---

#' #Pre-amble
#' (Last updated: `r Sys.Date()`)
#'
#' This file is for generating age/sex disaggregation of incidence: making final adjustments and outputting in correct format.
#'
#' This is therefore to be run after the analysis of 03inc_split.R
#'
#' N.B. This file should render to html with `rmarkdown::render('04inc_split_finish.R')` or from the command line with ` R -q -e "rmarkdown::render(\"04inc_split_finish.R\",output_dir=\"../html\")"`
#'
#' #Read in necessary data
#'
#' Relevant libraries for loading and manipulating data:
#'
rm(list=ls())
library(here)                           # file locations
source(here('disaggregation/R/utilities.R'))

#' Various pre-requisite data
load(here('data/est.rda'))         # current burden estimates
load(here('data/pop.rda'))         # populations
load(here('disaggregation/output/incsplits/data/tbimp.Rdata')) # notifications with child splits
load(here('disaggregation/R/02children/data/K.Rdata')) #child estimates
load(here('disaggregation/output/incsplits/data/incsplit.Rdata')) #incidnce splits

Sys.setlocale(locale="C")               #to avoid report writing UTF8 error

## check 3 age-cat countries
incsplit[,age:=gsub('_','-',age_group)]



#' Introduces `inc.num` if not there
if('inc.num' %ni% names(est)){
  cat('Warning: making own inc.num!\n')
  est[,c('inc.num',
         'inc.lo.num',
         'inc.hi.num'):=.(pop*inc*1e-5,pop*inc.lo*1e-5,pop*inc.hi*1e-5)]
}

#' If Philippe's `regional.Rdata` file is ready, use it; otherwise, make one.
fn <- here('data/regional.rda')
if(!file.exists(fn)){
  cat('Warning: making own regional estimate data!\n')
  regional <- est[,{
    inc.num=sum(inc.num);
    ss = sqrt(sum((inc.hi.num-inc.lo.num)^2))/3.92;
    inc.lo.num=inc.num - 1.96*ss; inc.hi.num=inc.num + 1.96*ss;
    list(inc.num=inc.num,inc.lo.num=inc.lo.num,inc.hi.num=inc.hi.num)
  },by=.(g.whoregion,year)]
} else {load(fn);}

#' Similarly, if the `global.Rdata` file is ready, use it; otherwise, make one.
fn <- here('data/global.rda')
if(!file.exists(fn)){
  cat('Warning: making own global estimate data!\n')
  global <- est[,{
    inc.num=sum(inc.num);
    ss = sqrt(sum((inc.hi.num-inc.lo.num)^2))/3.92;
    inc.lo.num=inc.num - ss*1.96; inc.hi.num=inc.num + ss*1.96;
    list(inc.num=inc.num,inc.lo.num=inc.lo.num,inc.hi.num=inc.hi.num)
  },by=.(year)]
} else {load(fn);}

#'
#' #Recalculation of uncertainty
#'
#' Suppose
#'
#' $$\sum_i X_{i} = X$$
#'
#' assuming the $$X_i$$ are independent,
#'
#' $$ var(X) = \sum_i var(X_i) $$
#'
#' Let $$p_i = E(X_i)/E(X)$$, $$F_i = sd(X_i)/E(X_i)$$ and $$F = sd(X)/E(X)$$. Then
#'
#' $$F^2 = \sum_i p_i^2 F_i^2$$
#'
#' To disaggregate correctly, we need to choose $$F_i$$ to satisfy this equation.
#'
## ' One solution is to take $$F_i = F / sqrt(\sum_i p_i^2) $$.
## '
## ' Another approach is to take $$F_i^2 \propto p_i$$, in which case $$F_i = sqrt(p_i)F / sqrt(\sum_i p_i^3)$$
## ' 

#'
#' A new section for creating incidence splits from Hazim's unified
#' database format: https://www.dropbox.com/s/0yxmkz6vnmywpj6/Estimates_files_spec.docx?dl=0
#'
#'
#'
#' First, some basic consistency checks:
#'
#' 
#' Check consistency with above 

TC <- est[year==estyr,.(iso3,g.whoregion,inc.num,inc.lo.num,inc.hi.num)]
TR <- regional[year==estyr][,.(g.whoregion,inc.num,
                               inc.lo.num,inc.hi.num)]
TG <- global[year==estyr][,.(inc.num,inc.lo.num,inc.hi.num)]

## -- midpoint checks
## regional
TC[,sum(inc.num),by='g.whoregion'];
TR[,.(g.whoregion,inc.num)]

## global
TC[,sum(inc.num)]
TR[,sum(inc.num)]
TG[,sum(inc.num)]

## -- uncertainty checks: some slight differences
## regional
TC[,sum((inc.lo.num-inc.hi.num)^2),by='g.whoregion'];
TR[,.(g.whoregion,(inc.lo.num-inc.hi.num)^2)] #

## global
TC[,sum((inc.lo.num-inc.hi.num)^2)]
TR[,sum((inc.lo.num-inc.hi.num)^2)]
TG[,((inc.lo.num-inc.hi.num)^2)]


#' regional outputs for database (done first to allow non-incsplit countries below)
## g.whoregion
db2 <- merge(incsplit[,.(iso3,sex=tolower(sex),age_group=age,prop)],
             est[year==estyr,.(iso3,year,inc.num,
                               inc.lo.num,inc.hi.num,g.whoregion)],
             by='iso3',all.y = FALSE)
db2[,inc:=inc.num*prop]            #incidence in this country/cat

## ## start temp
## (bad <- db2[!is.finite(inc),unique(iso3)])
## db2[bad,inc:=1e-6]

## ## end temp
db2[,prop:=inc/(sum(inc)+1e-15)] #renormalize to whole world
db2a <- db2[,.(prop=sum(prop)),by=.(g.whoregion,year,sex,age_group)]
db2a[,pattern:='a']; db2a[,sum(prop)]          #check

## other patterns
db2[,kid:=ifelse(age_group %in% c('0-4','5-14'),'0-14','15plus')] #kid or not
db2b <- db2[,.(prop=sum(prop)),by=.(g.whoregion,year,sex,age_group=kid)]
db2[,kid:=NULL]                         #whoregion/sex/<>15
db2b[,pattern:='b']; db2b[,sum(prop)]                        #check
db2c <- db2b[,.(prop=sum(prop)),by=.(g.whoregion,year,age_group)];
db2c[,pattern:='c']; db2c[,sex:='a']; db2c[,sum(prop)]       #whoregion/<>15
db2d <- db2b[,.(prop=sum(prop)),by=.(g.whoregion,year,sex)];
db2d[,pattern:='d']; db2d[,age_group:='a']; db2d[,sum(prop)] #whoregion/sex
db2e <- db2d[,.(prop=sum(prop)),by=.(g.whoregion,year)]
db2e[,c('age_group','sex','pattern'):=list('a','a','e')]; db2e[,sum(prop)] #whoregion


## join together
dbf <- rbindlist(list(db2a,db2b,
                      db2c[,names(db2a),with=FALSE],
                      db2d[,names(db2a),with=FALSE],
                      db2e[,names(db2a),with=FALSE]))
dbf[,sum(prop),by=pattern]              #check

dbf <- merge(dbf,regional[year==estyr,.(g.whoregion,inc.num,
                                        inc.lo.num,inc.hi.num)],
             by='g.whoregion',
             all.x=TRUE)
dbf[,prop:=prop/sum(prop),by=.(g.whoregion,pattern)] #normalize by pattern x whoregion

## make the relevant estimates
dbf[,bsd:=0.5*(inc.hi.num-inc.lo.num)/(1e-15+inc.num*sqrt(sum(prop^2))),
    by=.(pattern,g.whoregion)]; # disaggregate so as to sum
dbf[,best:=inc.num*prop]; dbf[,lo:=best*pmax(0,1-bsd)]; dbf[,hi:=best*(1+bsd)];
dbf[,c('measure','unit'):=.('inc','num')]

#' Checks:
#' 
## sums by region
dbf[,.SD[,.(test=sum(best),aggr=inc.num[1])],by=.(g.whoregion,pattern)];

## uncertainty aggregation:
dbf[,F2:=((inc.hi.num - inc.lo.num)/inc.num)^2]
(dbf[,.(Fsq=sum(prop^2*(hi-lo)^2/best^2),F2=F2[1]),by=.(g.whoregion,pattern)])
## dbf[,.(sum((hi-lo)^2),((inc.hi.num - inc.lo.num)^2)[1]),by=.(g.whoregion,pattern)]

## to keep
db2o <- dbf[,.(group_type='g_whoregion',group_name=g.whoregion,
               year,measure,
               unit,age_group,sex,best,lo,hi)]

## global decomposion
dbf[,sum(prop),by=.(g.whoregion,pattern)]
dbf[,incR:=inc.num[1],by=.(g.whoregion,pattern)] #regional incidence
dbf[,propR:=incR/dbf[sex=='a' & age_group=='a',sum(incR)]] #proportion in region
db3f <- dbf[,.(prop=sum(propR*prop)),by=.(year,sex,age_group,pattern)]
db3f <- cbind(db3f,global[year==estyr,.(inc.num,inc.lo.num,inc.hi.num)]) #global data
db3f[,bsd:=0.5*(inc.hi.num-inc.lo.num)/(1e-15+inc.num*sqrt(sum(prop^2))),
    by=.(pattern)]; # disaggregate so as to sum

db3f[,best:=inc.num*prop]; db3f[,lo:=best*pmax(0,1-bsd)]; db3f[,hi:=best*(1+bsd)];

db3f[,c('measure','unit'):=.('inc','num')]
db3f[,c('group_type','group_name'):=list('global','global')]

##checks
db3f[,.SD[,.(test=sum(best),inc.num[1])],by=.(pattern)] #check
db3f[,F2:=((inc.hi.num - inc.lo.num)/inc.num)^2]
db3f[,.(Fsq=sum(prop^2*(hi-lo)^2/best^2),F2=F2[1]),by=.(pattern)]

## join
db2o <- rbind(db2o,db3f[,names(db2o),with=FALSE])

## checks on joint output cross - consistency
tmp1 <- db2o[group_name!='global',
             .(suminc=sum(best)),by=.(age_group,sex)]
tmp2 <- db2o[group_name=='global',
             .(glob=sum(best)),by=.(age_group,sex)]
(tmp <- merge(tmp1,tmp2))

## save out
db_estimates_group <- copy(db2o)
attr(db_estimates_group, "timestamp") <- Sys.Date() #set date
fn2 <- here('disaggregation/dboutput/db_estimates_group.Rdata')
save(db_estimates_group,file=fn2)


#' Country incsplits file for database
db1a <- incsplit[,.(iso3,sex=tolower(sex),age_group=age,prop)]
db1a[,pattern:='a']
db1a[,kid:=ifelse(age_group %in% c('0-4','5-14'),'0-14','15plus')] #kid or not
db1b <- db1a[,.(prop=sum(prop)),by=.(iso3,sex,age_group=kid)]
db1a[,kid:=NULL]; db1b[,pattern:='b']
db1c <- db1b[,.(prop=sum(prop)),by=.(iso3,age_group)]
db1c[,sex:='a']; db1c[,pattern:='c']
db1d <- db1b[,.(prop=sum(prop)),by=.(iso3,sex)]
db1d[,age_group:='a']; db1d[,pattern:='d']
db1e <- db1d[,.(prop=sum(prop)),by=.(iso3)]
db1e[,c('age_group','sex','pattern'):=list('a','a','e')]

db1f <- rbindlist(list(db1a,db1b,
                      db1c[,names(db1a),with=FALSE],
                      db1d[,names(db1a),with=FALSE],
                      db1e[,names(db1a),with=FALSE]))
db1f[,analysed:=TRUE]                    #those with data
dbreg <- merge(db1f,est[year==estyr][,.(iso3,g.whoregion,inc.num)],
               by='iso3')
## dbreg <- dbreg[,.(prop=mean(prop)),by=.(g.whoregion,sex,age_group,pattern)] #reg averages
dbreg <- dbreg[,.(prop=weighted.mean(prop,inc.num)),
               by=.(g.whoregion,sex,age_group,pattern)] #reg averages
## countries with data
db1f <- merge(db1f,est[year==estyr][,.(iso3,year,inc.num,
                                       inc.lo.num,inc.hi.num)],
              by='iso3',all.x=TRUE)

## countries using regional averages
db10 <- merge(est[year==estyr][,.(iso3,year,inc.num,
                                  inc.lo.num,inc.hi.num,g.whoregion)],
              dbreg,
              by='g.whoregion',all.x=TRUE,all.y=TRUE,allow.cartesian = TRUE)
db10[,c('g.whoregion','analysed'):=list(NULL,FALSE)]
db10 <- db10[! iso3 %in% db1f[,unique(iso3)]]
db1f <- rbind(db1f,db10[,names(db1f),with=FALSE]) #join

## uncertainty
db1f[,bsd:=0.5*(inc.hi.num-inc.lo.num)/(1e-15+inc.num*sqrt(sum(prop^2))),
     by=.(iso3,pattern)];
db1f[,best:=inc.num*prop]; db1f[,lo:=best*pmax(0,1-bsd)]; db1f[,hi:=best*(1+bsd)];
db1f[,c('measure','unit'):=.('inc','num')]

##checks
db1f[,.(nm=sum(prop)),by=.(pattern,iso3)][,summary(nm)]
db1f[,.SD[,.(A=sum(best),B=inc.num[1])],by=.(iso3,pattern)][,any(abs(A-B)>.5)] #NOTE
db1f[,.(sum((hi-lo)^2),((inc.hi.num - inc.lo.num)^2)[1]),by=.(iso3,pattern)]
db1f[,F2:=((inc.hi.num - inc.lo.num)/inc.num)^2] #uncertainty
db1f[,.(Fsq=sum(prop^2*(hi-lo)^2/(1e-12+best)^2),F2=F2[1]),by=.(iso3,pattern)]
## db1f[,.(Fsq=sum(prop^2*(2*bsd)^2),F2=F2[1]),by=.(iso3,pattern)]

## checks on joint output cross - consistency
db1f[,length(unique(iso3))]; est[,length(unique(iso3))]
setdiff(est[,(unique(iso3))],db1f[,(unique(iso3))]) #
db1f <- merge(db1f,est[year==estyr][,.(iso3,g.whoregion)],by='iso3')
tmp1 <- db1f[,.(suminc=sum(best)),by=.(age_group,sex,g.whoregion)]
tmp2 <- db2o[group_name!='global',.(reg=(best)),
             by=.(age_group,sex,g.whoregion=group_name)]
(tmp <- merge(tmp1,tmp2))
tmp[,ratio:=suminc/reg];tmp[,ratidf:=suminc-reg]
summary(tmp)

## save out
db_estimates_country <- db1f[order(iso3),.(iso3,year,measure,
                                           unit,age_group,sex,
                                           best,lo,hi,analysed)]
attr(db_estimates_country, "timestamp") <- Sys.Date() #set date
fn2 <- here('disaggregation/dboutput/db_estimates_country.Rdata')
save(db_estimates_country,file=fn2)

#' # Clean up
#' 
#' Save out in number format
tmpD <- dcast(db1f[pattern=='d',.(iso3,age_group,sex,best,lo,hi)],
              iso3~age_group+sex,
              value.var = c('best','lo','hi'))
tmpC <- dcast(db1f[pattern=='c',.(iso3,age_group,sex,best,lo,hi)],
              iso3~age_group+sex,
              value.var = c('best','lo','hi'))
tmpB <- dcast(db1f[pattern=='b',.(iso3,age_group,sex,best,lo,hi)],
              iso3~age_group+sex,
              value.var = c('best','lo','hi'))
tmpA <- merge(tmpB,tmpC)
tmpA <- merge(tmpA,tmpD)
tmpA <- merge(tmpA,est[year==estyr][,.(iso3,inc.num,
                                       inc.lo.num,inc.hi.num)])
tmpA <- merge(tmpA,pop[year==estyr][,.(iso3,year,country,g.whoregion,
                                      e.pop.m014,e.pop.f014,
                                      e.pop.014=e.pop.m014+e.pop.f014,
                                      e.pop.m15plus,e.pop.f15plus,
                                      e.pop.15plus = e.pop.m15plus+e.pop.f15plus,
                                      e.pop.f=e.pop.f014+e.pop.f1524+e.pop.f2534+
                                        e.pop.f3544+e.pop.f4554+e.pop.f5564+e.pop.f65,
                                      e.pop.m=e.pop.m014+e.pop.m1524+e.pop.m2534+
                                        e.pop.m3544+e.pop.m4554+e.pop.m5564+e.pop.m65,
                                      e.pop.num)])
names(tmpA) <- gsub('-','_',names(tmpA))

## join whoregion, populations
PN <- tmpA[,.( iso3,g.whoregion,
             year, country, inc.num,
             inc.num.lo=inc.lo.num,  inc.num.hi=inc.hi.num,
             inc.num.f=best_a_f,
             inc.num.f.lo=lo_a_f,inc.num.f.hi=hi_a_f,inc.num.m=best_a_m,
             inc.num.m.lo=lo_a_m, inc.num.m.hi=hi_a_m, inc.num.15plus=best_15plus_a,
             inc.num.15plus.lo=lo_15plus_a,
             inc.num.15plus.hi=hi_15plus_a,inc.num.15plus.f=best_15plus_f,
             inc.num.15plus.f.lo=lo_15plus_f,inc.num.15plus.f.hi=hi_15plus_f,
             inc.num.15plus.m=best_15plus_m,
             inc.num.15plus.m.lo=lo_15plus_m,inc.num.15plus.m.hi=hi_15plus_m,
             inc.num.014=best_0_14_a,
             inc.num.014.lo=lo_0_14_a,  inc.num.014.hi=hi_0_14_a,
             inc.num.014.f=best_0_14_f,
             inc.num.014.f.lo=lo_0_14_f,  inc.num.014.f.hi=hi_0_14_f,
             inc.num.014.m=best_0_14_m,
             inc.num.014.m.lo=lo_0_14_m,inc.num.014.m.hi=hi_0_14_m,
             e.pop.m014,
             e.pop.f014,  e.pop.014, e.pop.m15plus,
             e.pop.f15plus, e.pop.15plus, e.pop.f,
             e.pop.m,  e.pop.num
             )]

PNn <- copy(PN)                         #numeric version
popz <- grep('pop',names(PN))
PN[,c(popz):=lapply(.SD,function(x)x*1e-3),.SDcols=popz] #pops in thousands
## format
PN[,5:40 := lapply(.SD, ftb), .SDcols=5:40]
PN[,5:40 := lapply(.SD, zapdotlt), .SDcols=5:40]

## save
A_country_incidence_disaggregated_age_sex_num <- PN
attr(A_country_incidence_disaggregated_age_sex_num,"timestamp") <- Sys.Date()

fn <- here('disaggregation/dboutput/A_country_incidence_disaggregated_age_sex_num.Rdata')
save(A_country_incidence_disaggregated_age_sex_num,file=fn)

#' save out in rate format
PR <- PNn[,.( iso3,g.whoregion,
             year, country, inc=inc.num/e.pop.num,  
             inc.sd=(inc.num.hi - inc.num.lo)/(3.92*e.pop.num),
             inc.f=inc.num.f/e.pop.f,
             inc.f.sd=(inc.num.f.hi - inc.num.f.lo)/(3.92*e.pop.f),
             inc.m=inc.num.m/e.pop.m,
             inc.m.sd=(inc.num.m.hi - inc.num.m.lo)/(3.92*e.pop.m),
             inc.15plus=inc.num.15plus/e.pop.15plus,
             inc.15plus.sd=(inc.num.15plus.hi-inc.num.15plus.lo)/(3.92*e.pop.15plus),
             inc.15plus.f=inc.num.15plus.f/e.pop.f15plus,
             inc.15plus.f.sd=(inc.num.15plus.f.hi-inc.num.15plus.f.lo)/
               (3.92*e.pop.f15plus),
             inc.15plus.m=inc.num.15plus.m/e.pop.m15plus,
             inc.15plus.m.sd=(inc.num.15plus.m.hi-inc.num.15plus.m.lo)/
               (3.92*e.pop.m15plus),
             inc.014=inc.num.014/e.pop.014,
             inc.014.sd=(inc.num.014.hi-inc.num.014.lo)/(3.92*e.pop.014),
             inc.014.f=inc.num.014.f/e.pop.f014,
             inc.014.f.sd=(inc.num.014.f.hi-inc.num.014.f.lo)/(3.92*e.pop.f014),
             inc.014.m=inc.num.014.m/e.pop.m014,
             inc.014.m.sd=(inc.num.014.m.hi-inc.num.014.m.lo)/(3.92*e.pop.m014),
             e.pop.m014,
             e.pop.f014,  e.pop.014, e.pop.m15plus,
             e.pop.f15plus, e.pop.15plus ,e.pop.f,
             e.pop.m,  e.pop.num
             )]


PRn <- copy(PR)
popz <- grep('pop',names(PR))
nonpopz <- setdiff(names(PR),c(names(PR)[popz],names(PR)[1:4]))
PR[,c(nonpopz):=lapply(.SD,function(x)x*1e5),.SDcols=c(nonpopz)] #per 100k
PR[,c(popz):=lapply(.SD,function(x)x*1e-3),.SDcols=popz] #pops in thousands
PR[,5:ncol(PR) := lapply(.SD, ftb), .SDcols=5:ncol(PR)]
PR[,5:ncol(PR) := lapply(.SD, zapdotlt), .SDcols=5:ncol(PR)]

## save
A_country_incidence_disaggregated_age_sex_rate <- PR
attr(A_country_incidence_disaggregated_age_sex_rate,'timestamp') <- Sys.Date()
fn <- here('disaggregation/dboutput/A_country_incidence_disaggregated_age_sex_rate.Rdata')
save(A_country_incidence_disaggregated_age_sex_rate,file=fn)

#' New version of regional aggregate save out, first formatting data:
## col names for different patterns
dnmz <- c(t(outer(c('inc.num.f','inc.num.m'),
                  c("",".lo",".hi"),paste0)))
cnmz <- c(t(outer(c('inc.num.15plus','inc.num.014'),
                  c("",".lo",".hi"),paste0)))
bnmz <- c(t(outer(c('inc.num.15plus.f','inc.num.15plus.m',
                    'inc.num.014.f','inc.num.014.m'),
                  c("",".lo",".hi"),paste0)))

## M/F
tmpD <- dcast(dbf[pattern=='d',.(g.whoregion,age_group,sex,best,lo,hi)],
              g.whoregion~age_group+sex,
              value.var = c('best','lo','hi'))
tmpD <- tmpD[,.(g.whoregion,best_a_f,lo_a_f,hi_a_f,best_a_m,lo_a_m,hi_a_m)]
names(tmpD)[2:ncol(tmpD)] <- dnmz

## kid/adult
tmpC <- dcast(dbf[pattern=='c',.(g.whoregion,age_group,sex,best,lo,hi)],
              g.whoregion~age_group+sex,
              value.var = c('best','lo','hi'))
tmpC <- tmpC[,.(g.whoregion,best_15plus_a,lo_15plus_a,hi_15plus_a,
                `best_0-14_a`,`lo_0-14_a`,`hi_0-14_a`)]
names(tmpC)[2:ncol(tmpC)] <- cnmz

## age x sex
tmpB <- dcast(dbf[pattern=='b',.(g.whoregion,age_group,sex,best,lo,hi)],
              g.whoregion~age_group+sex,
              value.var = c('best','lo','hi'))
tmpB <- tmpB[,.(g.whoregion,
                best_15plus_f,lo_15plus_f,hi_15plus_f,
                best_15plus_m,lo_15plus_m,hi_15plus_m,
                `best_0-14_f`,`lo_0-14_f`,`hi_0-14_f`,
                `best_0-14_m`,`lo_0-14_m`,`hi_0-14_m`)]
names(tmpB)[2:ncol(tmpB)] <- bnmz


## join
tmpA <- merge(tmpD,tmpC)
tmpA <- merge(tmpA,tmpB)

tmpA <- merge(regional[year==estyr][,.(g.whoregion,inc.num,
                                       inc.lo.num,inc.hi.num)],
              tmpA,
              by='g.whoregion')
tmpA <- cbind(year=estyr,tmpA)
tmpA[,3:ncol(tmpA) := lapply(.SD, ftb), .SDcols=3:ncol(tmpA)] #format


## save
A_regional_incidence_disaggregated_age_sex <- tmpA
attr(A_regional_incidence_disaggregated_age_sex,'timestamp') <- Sys.Date()

fn <- here('disaggregation/dboutput/A_regional_incidence_disaggregated_age_sex.Rdata')
save(A_regional_incidence_disaggregated_age_sex,file=fn)



#' Save out aggregated format (global)

## M/F
Tmpd <- dcast(db3f[pattern=='d',.(group_name,age_group,sex,best,lo,hi)],
              group_name~age_group+sex,
              value.var = c('best','lo','hi'))
Tmpd <- Tmpd[,.(group_name,best_a_f,lo_a_f,hi_a_f,best_a_m,lo_a_m,hi_a_m)]
names(Tmpd)[2:ncol(Tmpd)] <- dnmz

## kid/adult
Tmpc <- dcast(db3f[pattern=='c',.(group_name,age_group,sex,best,lo,hi)],
              group_name~age_group+sex,
              value.var = c('best','lo','hi'))
Tmpc <- Tmpc[,.(group_name,best_15plus_a,lo_15plus_a,hi_15plus_a,
                `best_0-14_a`,`lo_0-14_a`,`hi_0-14_a`)]
names(Tmpc)[2:ncol(Tmpc)] <- cnmz

## age x sex
Tmpb <- dcast(db3f[pattern=='b',.(group_name,age_group,sex,best,lo,hi)],
              group_name~age_group+sex,
              value.var = c('best','lo','hi'))
Tmpb <- Tmpb[,.(group_name,
                best_15plus_f,lo_15plus_f,hi_15plus_f,
                best_15plus_m,lo_15plus_m,hi_15plus_m,
                `best_0-14_f`,`lo_0-14_f`,`hi_0-14_f`,
                `best_0-14_m`,`lo_0-14_m`,`hi_0-14_m`)]
names(Tmpb)[2:ncol(Tmpb)] <- bnmz


## join
Tmpa <- merge(Tmpd,Tmpc)
Tmpa <- merge(Tmpa,Tmpb)
Tmpa[,group_name:=NULL]
Tmpa <- cbind(global[year==estyr][,.(inc.num,inc.lo.num,inc.hi.num)],Tmpa)
Tmpa <- cbind(year=estyr,g.whoregion='Global',Tmpa)
Tmpa[,3:ncol(Tmpa) := lapply(.SD, ftb), .SDcols=3:ncol(Tmpa)]


## save
A_global_incidence_disaggregated_age_sex <- Tmpa
attr(A_global_incidence_disaggregated_age_sex,'timestamp') <- Sys.Date()
fn <- here('disaggregation/dboutput/global_incidence_disaggregated_age_sex.Rdata')
save(A_global_incidence_disaggregated_age_sex,file = fn)

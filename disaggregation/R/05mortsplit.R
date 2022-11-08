#' ---
#' title: Generating age/sex mortality disaggregation
#' author: Pete Dodd
#' date: 20 July, 2022
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---

#' # Pre-amble
#' (Last updated: `r Sys.Date()`)
#'
#' This file is for generating age/sex disaggregation of mortality.
#'
#' This is therefore to be run after the analysis of prevalence data and after the child TB data has been generated (incidence and mortality), and after the incidence age/sex splits have been run.
#'
#' N.B. This file should render to html with `rmarkdown::render('05mortsplit.R')` or from the command line with ` R -q -e "rmarkdown::render(\"05mortsplit.R\",output_dir=\"../html\")"`
#'
#' 
#' # Read in necessary data
#'
#' Relevant libraries for loading and manipulating data:
#'
rm(list=ls())
library(here)                           # file locations
source(here('disaggregation/R/utilities.R'))


#' encoding etc
gypt <- TRUE

#' Various pre-requisite data 
load(here('data/est.rda'))         # current burden estimates NOTE will be est
load(here('data/pop.rda'))         # populations for country name
load(here('data/unaids.rda'))          # includes notifications
load(here('disaggregation/R/02children/data/KM.Rdata'))   # Child FR estimates
## load(here('disaggregation/prevsplits/prevsplits.Rdata'))   # from prevalence split work
load(here('disaggregation/output/incsplits/data/incsplit.Rdata')) #from prevalence split work
## names(prevsplits)[2] <- 'age_group'

#' Create output directory if missing
fn <- here('disaggregation/output/mortsplits')
if(!file.exists(fn)) dir.create(fn)
fn <- here('disaggregation/output/mortsplits/plots')
if(!file.exists(fn)) dir.create(fn)

#' Check which countries to get IHME data for:
fn <- here('disaggregation/output/mortsplits/vr_ihme_iso3.txt')
if(!file.exists(fn)){
  ihme <- est[source.mort=='IHME',unique(iso3)]
  cat(ihme,file=fn)
} else {
  ihme <- scan(fn,'character')
}
ihme

#' If Philippe's `regional.rda` file is ready, use it; otherwise, make one.
fn <- here('data/regional.rda')
if(!file.exists(fn)){
  cat('Warning: making own regional estimate data!\n')
  regional <- est[,{
    ## incidence
    inc.num=sum(inc.num);
    ss = sqrt(sum((inc.hi.num-inc.lo.num)^2))/3.92;
    inc.lo.num=inc.num - 1.96*ss; inc.hi.num=inc.num + 1.96*ss;
    ## mort HIV-ve
    mort.nh.num=smn(mort.nh.num);
    ss = sqrt(smn((mort.nh.hi.num-mort.nh.lo.num)^2))/3.92;
    mort.nh.lo.num=mort.nh.num - 1.96*ss;
    mort.nh.hi.num=mort.nh.num + 1.96*ss;
    ## mort HIV+ve
    mort.h.num=smn(mort.h.num);
    ss = sqrt(smn((mort.h.hi.num-mort.h.lo.num)^2))/3.92;
    mort.h.lo.num=mort.h.num - 1.96*ss;
    mort.h.hi.num=mort.h.num + 1.96*ss;
    list(inc.num=inc.num,inc.lo.num=inc.lo.num,inc.hi.num=inc.hi.num,
       mort.nh.num=mort.nh.num,mort.nh.lo.num=mort.nh.lo.num,
       mort.nh.hi.num=mort.nh.hi.num,
       mort.h.num=mort.h.num,mort.h.lo.num=mort.h.lo.num,
       mort.h.hi.num=mort.h.hi.num)
  },by=.(g.whoregion,year)]
} else {load(fn);}

#' Similarly, if the `global.Rdata` file is ready, use it; otherwise, make one.
fn <- here('data/global.rda')
if(!file.exists(fn)){
  cat('Warning: making own global estimate data!\n')
  global <- est[,{
    ## incidence
    inc.num=sum(inc.num);
    ss = sqrt(sum((inc.hi.num-inc.lo.num)^2))/3.92;
    inc.lo.num=inc.num - ss*1.96; inc.hi.num=inc.num + ss*1.96;
    ## mort HIV-ve
    mort.nh.num=smn(mort.nh.num);
    ss = sqrt(smn((mort.nh.hi.num-mort.nh.lo.num)^2))/3.92;
    mort.nh.lo.num=mort.nh.num - 1.96*ss;
    mort.nh.hi.num=mort.nh.num + 1.96*ss;
    ## mort HIV+ve
    mort.h.num=smn(mort.h.num);
    ss = sqrt(smn((mort.h.hi.num-mort.h.lo.num)^2))/3.92;
    mort.h.lo.num=mort.h.num - 1.96*ss;
    mort.h.hi.num=mort.h.num + 1.96*ss;
    list(inc.num=inc.num,inc.lo.num=inc.lo.num,inc.hi.num=inc.hi.num,
         mort.nh.num=mort.nh.num,mort.nh.lo.num=mort.nh.lo.num,
         mort.nh.hi.num=mort.nh.hi.num,
         mort.h.num=mort.h.num,mort.h.lo.num=mort.h.lo.num,
         mort.h.hi.num=mort.h.hi.num)
  },by=.(year)]
} else {load(fn);}


#' Changes to data:
#' 
#' Introduces `inc.num` if not there & adds country column
if('inc.num' %ni% names(est)){
  cat('Warning: calculating own inc.num!\n')
  est[,c('inc.num',
         'inc.lo.num',
         'inc.hi.num'):=.(pop*inc*1e-5,pop*inc.lo*1e-5,pop*inc.hi*1e-5)]
}
if('country' %ni% names(est)){
  cat('Adding country name to est!\n')
  est <- merge(est,unique(pop[,.(iso3,country)]),all.x=TRUE)
}

## a function to sum without NAs
smn <- function(x) sum(x,na.rm=TRUE)

#' Restrict to current year
regional <- regional[year==estyr]
global <- global[year==estyr]

#' Namings
hbcsh <- unique(pop[,.(iso3,name=country)])
hbcsh <- merge(hbcsh,unique(est[,.(iso3,g.hbc)]))
hbcsh <- hbcsh[g.hbc==TRUE]
## shorter names for graphs
hbcsh[iso3=='COD',name:="DR Congo"]
hbcsh[iso3=='PRK',name:="DPR Korea"]
hbcsh[iso3=='TZA',name:="UR Tanzania"]

#' Some preparation with `incsplit`
## incsplit$age <- NULL#factor(incsplit$age,levels=agz3,ordered=TRUE)
## incsplit$sex <- factor(incsplit$sex,levels=c('M','F'),ordered=TRUE)
incsplit <- incsplit[order(iso3,sex,age_group),]
setkey(incsplit,iso3)

#' And `KM`
KM[,h:=HIVinDeaths/1e3]

#'
#' # VR data
#'
#' First a function for reformatting
refrm <- function(indat){
  indat <- indat[,.SD[Year==max(Year)],by=Country]             #most recent year
  ## rename & aggregate
  indat <- indat[,.(Country,name,Year,icd,cause1,Sex,
                    `0-4`=Deaths2+Deaths3+Deaths4+Deaths5+Deaths6,
                    `5-14`=Deaths7+Deaths8,
                    `15-24`=Deaths9+Deaths10,
                    `25-34`=Deaths11+Deaths12,
                    `35-44`=Deaths13+Deaths14,
                    `45-54`=Deaths15+Deaths16,
                    `55-64`=Deaths17+Deaths18,
                    `65plus`=Deaths19+Deaths20+
                      Deaths21+Deaths22+Deaths23+Deaths24+Deaths25)]
  ## reshape
  MM <- melt(indat,id = c("Country","name","Year","icd","cause1","Sex"))
  ## separate for TB/ill/tot
  M3a <- MM[cause1=='TB',.(Country,name,Year,icd,Sex,variable,value)]
  M3b <- MM[cause1=='ill',.(Country,name,Year,icd,Sex,variable,ill=value)]
  M3c <- MM[cause1=='tot',.(Country,name,Year,icd,Sex,variable,tot=value)]
  ## M3h <- M3[cause1=='HIV',.(Country,name,Year,icd,Sex,variable,HIV=value)] #2agg 
  ## merge & calculate
  mv <- c('Country','name','Year','icd','Sex','variable')
  M3a <- merge(M3a,M3b,all.x=TRUE)
  M3a <- merge(M3a,M3c,all.x=TRUE)
  ## print(M3a);print(M3h)
  ## M3a <- merge(M3a,M3h,all.x=TRUE,all.y = FALSE,by=mv)
  M3a[,g:=ill/tot]
  M3a[,value:=value/(1-g)]
  M3a[,keep:=!is.na(sum(value)),by=name]
  M3a <- M3a[keep==TRUE,]
  MM <- M3a[Sex %in% c(1,2),]
  MM <- MM[order(name),]
  MM <- MM[order(name),]
  MM$sex <- c('M','F')[as.numeric(MM$Sex)]
  MM$sex <- factor(MM$sex)
  MM$age <- factor(MM$variable,levels=agz3,ordered=TRUE)
  MM[,age_group:=gsub('-','_',age)]
  MM[,age:=NULL]
  MM
}



#' Extract and reformat each of the sheets in the spreadsheet in turn:
#'
typz <- c('text','text','text','text',
          rep('text',5),rep('numeric',26)) #specify column types to avoid warning
xlfn <- here('input/VR/info_TB_june2022.xls')


#' big loop to generate cleaned/reshaped file if not there:
mfn <- here('disaggregation/output/mortsplits/VR.Rdata')
if(!file.exists(mfn)){
  ## reading in excel data

  ## -------- ICD>10.1
  M1 <- as.data.table(read_excel(xlfn,sheet = 3,
                                 skip=0,col_types=typz)) #ICD > 10.1
  M1 <- refrm(M1)
  ## -------- ICD 9
  M2 <- as.data.table(read_excel(xlfn,sheet = 5,
                                 skip=0,col_types=typz)) #ICD = 9
  M2 <- refrm(M2)
  ## -------- ICD 8
  M3 <- as.data.table(read_excel(xlfn,sheet = 6,
                                 skip=0,col_types=typz)) #ICD = 8
  M3 <- refrm(M3)
  ## -------- ICD 10.1
  M4 <- as.data.table(read_excel(xlfn,sheet = 4,
                                 skip=0,col_types=typz)) #ICD = 10.1
  M4 <- refrm(M4)

  ## --- join ---
  VR <- rbind(M1,M2,M3,M4)

  ## known fixes
  tky <- pop[iso3=='TUR',country][1] ## "T<U+00FC>rkiye"
  VR[name=='Turkey',name:=tky] #new name
  VR[name=='Turkey'] #old name
  VR[name==tky] #new name
  pop[country==tky]

  ## Differences in names to be done by hand:
  (vrbad <- setdiff(VR[,as.character(unique(name))],
                   pop[,as.character(unique(country))]))

  ## direct renamings:
  done <- c()
  for(i in seq_along(vrbad)){
    newnm <- grep(vrbad[i],pop[,unique(country)],value=TRUE)
    if(length(newnm)>0){
      print(newnm)
      VR[name==vrbad[i],name:=newnm]      #rename
      done <- c(done,i)
    }
  }
  vrbad <- vrbad[-done]

  ## others for renaming
  vrbad
  (newnm <- grep('Czech',pop[,unique(country)],value=TRUE))
  VR[name==grep('Czech',vrbad,value=TRUE),name:=newnm]
  (newnm <- grep('Macedonia',pop[,unique(country)],value=TRUE)[1])
  VR[name==grep('Mace',vrbad,value=TRUE),name:=newnm]

  ## those still bad
  vrbad <- vrbad[!str_detect(vrbad,'Cze|Serb|Mace')]

  ## sub-countries
  VR[name %in% c("French Guiana","Martinique","Reunion","Mayotte","Guadeloupe"),
     name:="France"]
  ## VR[name %in% c("Rodrigues"),name:="Mauritius"]

  (fulln <- pop[grepl('Vincent',country),country][1])
  VR[grepl('Vincent',name),name:=fulln]

  ## check
  setdiff(VR[,unique(name)],pop[,unique(country)]) #should be none missing

  ## aggregate
  VR <- VR[,.(value=sum(value*(1-g)),ill=sum(ill),tot=sum(tot)),
           by=.(name,Sex,variable)]
  VR[,g:=ill/tot]                         #reintroduce
  VR[,value:=value/(1-g)]

  ## Add iso3, and tidy up
  VR <- merge(VR,pop[year==estyr][,.(iso3,country)],
              by.x = 'name',by.y="country",all.x=TRUE)
  (nmz <- VR[is.na(iso3),unique(name)])   #should be none
  ccnt <- rowSums(VR[,table(iso3,Sex)]) #some more than once or not at all
  ccnt <- data.table(iso3=names(ccnt),N=ccnt)
  VR <- VR[iso3 %in% ccnt[N>0,as.character(unique(iso3))],]
  ## VR <- VR[,.(value=sum(value)),by=.(iso3,name,Country,Sex,variable,sex,age_group)]
  (rs <- rowSums(VR[,table(iso3,Sex)]))
  (bad <- names(rs)[rs<16])
  VR[iso3 %in% bad,.(name,Sex,variable,value)]
  ## probably no great loss in dropping these:
  VR <- VR[iso3 %ni% bad]
  ## tidy:
  VR$sex <- factor(c('M','F')[as.numeric(as.character(VR$Sex))],
                   levels=c('M','F'),ordered=TRUE)
  VR[,Sex:=NULL]
  VR[,age_group:=variable]
  VR$age_group <- factor(gsub("-","_",VR$age_group),levels=agz2,ordered=TRUE)
  VR[,variable:=NULL]
  VR$iso3 <- factor(VR$iso3)
  setkey(VR,iso3)                         #80 countries!
  rm(M1,M2,M3,M4)

  ## load IHME country age data
  load(here('disaggregation/R/02children/mortinput/IHME.Rdata'))

  ## are these in VR?
  ihme[which(ihme %in% VR$iso3)]
  VR[iso3=='ZAF',sum(value)]
  GP <- ggplot(VR[iso3=='ZAF'],aes(x=age_group,y=value,col=sex)) +
    geom_point() +
    geom_point(data=IHME[iso3=='ZAF'],shape=2)
  ggsave(GP,
         file=here('disaggregation/output/mortsplits/plots/ihme_ZAF.pdf'))

  ## replace values in VR
  VR <- VR[iso3 %ni% ihme]
  VR <- rbind(VR,IHME)

  ## save out!
  save(VR,file=mfn)

} else {
    load(file=mfn)
}


#'
#' # Plotting
#'
#' A plotting sub-function:
muplot <- function(indat,cnm){
    if('age' %ni% names(indat)) indat[,age:=gsub("_","-",age_group)]
    indat[age=='65plus',age:=rev(agz4)[1]]
    indat$age <- factor(indat$age,levels=agz4,ordered=TRUE)
    indat$sex <- factor(indat$sex,levels=c('F','M'),ordered=TRUE)
    ## plot
    mp <- prodplot(data=indat,mort ~ sex + age,divider=mosaic()) +
        aes(fill=sex)
    mp <- mp + theme(axis.line=element_blank(),
                     axis.text.x=element_text(angle=90),
                     axis.text.y=element_text(),
                     axis.ticks=element_blank(),
                     axis.title.x=element_blank(),
                     axis.title.y=element_blank(),
                     legend.title=element_blank(),
                     legend.position="none",
                     panel.background=element_blank(),
                     panel.border=element_blank(),
                     panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank(),
                     plot.background=element_blank()) +
        ggtitle(cnm)
    mp
}


#'
#' # Process for VR countries
#'
#' There is an over-ride that can use VR data in preference to `KM` for children.
#'
#' This will need regional averages of VR-based mortality splits:
VRS <- est[year==estyr][,.(iso3,source.mort)]
VR <- merge(VR,est[year==estyr][,.(iso3,g.whoregion)],by='iso3')
VRreg <- VR[,.(mort=sum(value)),by=.(sex,age_group,g.whoregion)]
VRreg[,mort:=mort/sum(mort),by=g.whoregion]
setkey(VRreg,g.whoregion)

#' Have a look at these distributions (NB only countries with VR data):
#' 

plts <- list()
for(reg in VRreg[,as.character(unique(g.whoregion))]){
  plts[[reg]] <- muplot(VRreg[reg],reg)
}
ggarrange(plotlist=plts,ncol=3,nrow=2)



#'
#' A function for constructing splits from VR data
#'
#'

getVRMortnh <- function(cn,KMoverride=TRUE,allplt=FALSE){
    el <- est[cn][year==estyr][,.(inc.num,inc.nh.num,
                                  mort.nh.num,country)] #incidence & mortality ests
    tmp <- VR[cn]           # split of mortality in VR data
    tmpi <- incsplit[cn]    # split of incidence
    km <- KM[cn]            # child mortality data
    tmp[,mort:=-9.9]        # for output - initialized as dble
    tmp[,src.mort:='VR']    # source for mortality disaggregation
    if(tmp[,is.na(tot[1])]) tmp[,src.mort:='IHME VR']
    tmp <- merge(tmp,tmpi[,.(age_group,sex,inc)],by=c('sex','age_group'))
    ## kids: mf split young & older from incidence
    mfy <- tmpi[age_group=='0_4' & sex=='M',inc]/
        tmpi[age_group=='0_4',sum(inc)+1e-15]
    mfo <- tmpi[age_group=='5_14' & sex=='M',inc]/
        tmpi[age_group=='5_14',sum(inc)+1e-15]
    tmp[age_group=='0_4' & sex=='M',
        mort := mfy * km[,deathsY*(1-h)]]   # young,m
    tmp[age_group=='0_4' & sex=='F',
        mort := (1-mfy) * km[,deathsY*(1-h)]]   # young,f
    tmp[age_group=='5_14' & sex=='M',
        mort := mfo * km[,(deaths-deathsY)*(1-h)]]  # old,m
    tmp[age_group=='5_14' & sex=='F',
        mort := (1-mfo) * km[,(deaths-deathsY)*(1-h)]] # old,f
    ## adults
    amort <- el$mort.nh.num - km[,deaths*(1-h)]      # adult deaths
    ## add to data
    M <- tmp[age_group %ni% kds,sum(value)]
    tmp[age_group %ni% kds,mort:= amort * value/M]
    if(KMoverride){                       # uses VR data directly
        tmp[,mort:=value]
        ## if(el$mort.nh.num<1e2){# regional average for low totals
        ##   cat('Fewer than 100 TB deaths in',cn,'- using regional VR average!\n')
        ##   reg <- as.character(tmp$g.whoregion[1])
        ##   tmp[,mort:=VRreg[reg,mort]]
        ## }
    }
    tmp[,mort:=mort/(sum(mort)+1e-15)]
    tmp[,CFR:=mort*el[,mort.nh.num/(inc.nh.num+1e-15)]] #CFRs by age gender
    if(any(tmp$CFR>.5)) warning(paste0('a CFR>50% in ',cn,'!')) #NB approximate
    tmp <- tmp[,.(iso3,sex,age_group,mort,src.mort,CFR)]
    tmp[,mort:=mort/(sum(mort)+1e-15)]
    ## plot
    mp <- muplot(tmp,el[,country])
    mph <- mpnh <- NA
    if(allplt) mpnh <- mp
    list(dat=tmp,
         HN=tmp[,.(iso3,sex,age_group,mort,src.mort)],
         plt=mp,plth=mph,pltnh=mpnh)
}


#' Testing...
## getVRMortnh('PHL')
## ## getVRMort('SWE')
## getVRMortnh('AUS')


#' Looking at the IHME places
#' 
## getVRMortnh('ZAF',KMoverride = FALSE)                      #pretty suspicious of this.
## getVRMortnh('ZAF',KMoverride = TRUE)
## getVRMortnh('BOL')
## getVRMortnh('HTI')

## ## India one of these
## getVRMortnh('IND',KMoverride = TRUE)
## getVRMortnh('IND',KMoverride = FALSE)


#'
#' # Adding CFR-based splits 
#'
#' 
#' Where VR data have not been used for estimating mortality, we use the incidence age/sex split and apply CFRs to this
#'
#' ## HIV by age
#'
#' There is a gender split in Philippe's data, which it is natiural to preserve.
#' However, we still need to split the HIV prevalence by age group (and make simple assumptions about non-differential CDR and ART coverage and CFR).
#'
#' This is problematic as UNAIDS public age-disaggregations from http://aidsinfo.unaids.org are different age groups (as well as annoyingly formatted).
#'
#' 
## ------ HIV data 
## total

HA <- unaids[year==estyr][,.(iso3,hiv=hiv.num)]
HA <- merge(HA,est[year==estyr][,.(iso3,Country=country)],by='iso3')
## kids
HK <- unaids[year==estyr][,.(iso3,hivk=h014)]
HK <- merge(HK,est[year==estyr][,.(iso3,Country=country)],by='iso3')
HA <- merge(HA,HK,by=c('iso3','Country'))

## young people 15-24
bn <- here('disaggregation/R/02children/mortinput')
fn <- here(bn,'People living with HIV_People living with HIV - Young people (15-24)_Population_ All young people (15-24).csv')
HYP <- fread(fn,header=TRUE)
HYP <- HYP[,.(Country,`2020`)]
HYP <- HYP[,hivyp:=as.numeric(gsub("[[:space:]+|<+|\\.+]","",`2020`))]
## older folk 50+
bn <- here('disaggregation/R/02children/mortinput')
fn <- here(bn,'People living with HIV_People living with HIV - People aged 50 and over_Population_ All people aged 50 and over.csv')
HO <- fread(fn,header=TRUE)
HO <- HO[,.(Country,`2020`)]
HO <- HO[,hivo:=as.numeric(gsub("[[:space:]+|<+|\\.+]","",`2020`))]
HO[,`2020`:=NULL]; HYP[,`2020`:=NULL]
HA <- merge(HA,HYP,by='Country')
HA <- merge(HA,HO,by='Country')

HA <- HA[,.(Country,hivyp=hivyp/(hiv-hivk),
            hivmid=(hiv-hivk-hivo-hivyp)/(hiv-hivk),
            hivo=hivo/(hiv-hivk))]

#' Here yp = 15-24, mid = 25-50, and old=50+
#' Need some ad hoc rules to arrive at our age splits of:
#' 15-24 (trivial), 25-34, 35-44, 44-55, 55-64, 65plus
#' Need to parametrize 25-34, 35-44, 44-55, 55-64, 65plus in terms of two variables
#' as only two inputs. Assume:
#'
#' * 25-34 and 35-44 are the same (x)
#' * 44-55 is at x/2
#' * 55-64 is y
#' * 65plus is at y/2
#'
#' Then (assuming even ages across 44-55):
#'
#' $$mid=(2+1/4)\times x$$
#' $$old=x/4 + (1+1/2)y$$
#'
#' allowing:
#'
#' $$x=4.mid/9$$
#' $$y=2(old-mid/9)/3$$
#'
#' Calculate:

HA[,x:=4*hivmid/9]
HA[,y:=2*(hivo-hivmid/9)/3]
HA[y<0,y:=0]                            #couple
HA[,`15_24`:=hivyp]
HA[,`25_34`:=x]
HA[,`35_44`:=x]
HA[,`45_54`:=x/2]
HA[,`55_65`:=y]
HA[,`65plus`:=y/2]
## ditch unnecessary
HA[,hivyp:=NULL];HA[,hivmid:=NULL];HA[,hivo:=NULL];HA[,x:=NULL];HA[,y:=NULL];
HA <- melt(HA,id='Country')             #reshape
HA <- HA[order(Country,variable),]

#'
#' Merge by name (danger!): 
#'
(bad <- setdiff(HA[,unique(Country)],est[,unique(country)])) #none


HA <- merge(HA,est[year==estyr][,.(iso3,country)],
            by.x='Country',by.y='country',all.x = TRUE)
## drop NAs
HA[,keep:=!is.na(sum(value)),by=Country]
HA <- HA[keep==TRUE,]                   #drop NAs
setkey(HA,iso3)
HA                                      #inspect

#' Create an average for countries without data:
hav <- HA[,.(value=mean(value,na.rm = FALSE)),by=variable]
hav

#'
#' ## CFR with HIV-weighting
#'
#' 
#' This function just fetches the age-split and sex-split of HIV in adults
#' 

getHIVsplit <- function(cn){
  sr <- 1                               #default sex ratio
  if(cn %in% unaids$iso3) sr <- unaids[cn][year==estyr][,h15m/h15f]
  sr <- sr/(1+sr)
  pa <- hav$value                       #default age split
  if(cn %in% HA$iso3) pa <- HA[cn]$value
  list(pa=pa,sr=sr)
}

#' This function modifies the CFR approach above to apply CFRs that differ by age/sex in order to appropriately average between the HIV+ CFR and the HIV- CFR in this country

getCFRMortH <- function(cn,allplt=FALSE){
    ##incidence & mortality ests
    elr <- est[year==estyr][, .(inc.num,country,inc.h.num, inc.nh.num,
                                mort.h.num = mort.h.num,mort.nh.num,
                                iso3,g.whoregion)]  #countries in same region
    el <- elr[iso3==cn]
    elr <- elr[g.whoregion==el$g.whoregion & iso3%ni%VR$iso3 & iso3!=cn] #regional data
    tmp <- incsplit[cn]                   #split of incidence
    km <- KM[cn]                          #child mortality data
    tmp[,mort:=-9.9]                      #to hold output - initialized as dble
    tmp[,src.mort:='CFR']                 #source for mortality disaggregation
  ## kids
    mfy <- tmp[age_group=='0_4' & sex=='M',inc]/
        tmp[age_group=='0_4',(sum(inc)+1e-15)] #mf split young
    mfo <- tmp[age_group=='5_14' & sex=='M',inc]/
        tmp[age_group=='5_14',(sum(inc)+1e-15)] #mf split old
    tmp[age_group=='0_4' & sex=='M',
        mort := mfy * km$deathsY ]  #young,m
    tmp[age_group=='0_4' & sex=='F',
        mort := (1-mfy) * km$deathsY ]#young,f
    tmp[age_group=='5_14' & sex=='M',
        mort := mfo * km[,deaths-deathsY] ]          #old,m
    tmp[age_group=='5_14' & sex=='F',
        mort := (1-mfo) * km[,deaths-deathsY] ]      #old,f
    ## adults
    amort <- el$mort.nh.num+el$mort.h.num - km$deaths      #adult deaths
    ainc <- el$inc.num*tmp[age_group %ni% kds,sum(prop)] #adult inc
    CFR <- amort / (ainc+1e-15)             #CFR in adults
    ## slight bodge as not HIV split for children (but usually small)
    h <- el[,inc.h.num/(inc.num+1e-15)]           #HIV in TB
    CFR.h <- (el$mort.h.num - km$deathsH) /(ainc*h+1e-15)
    CFR.nh <- (el$mort.nh.num - km$deaths + km$deathsH) /
        (ainc*(1-h)+1e-15)
    if(CFR.nh<=0){                        #use regional average CFR
        cat(cn,':CFR.nh < 0: using regional non-VR average!\n')
        CFR.nh <- elr[,mean(mort.nh.num/(inc.nh.num+1e-15),na.rm = TRUE)]
    }
    hs <- getHIVsplit(cn)
    pa <- hs$pa; sr <- hs$sr;
    hm <- h * sr * pa;
    hf <- h * (1-sr) * pa; #HIV split among M/F TB cases
    CFRm <- hm*CFR.h + (1-hm)*CFR.nh
    CFRf <- hf*CFR.h + (1-hf)*CFR.nh
    ## add to data
    tmp[age_group %ni% kds & sex=='F',mort:= CFRf * inc]
    tmp[age_group %ni% kds & sex=='M',mort:= CFRm * inc]
    tmp[,mort:=mort/(sum(mort)+1e-15)]
    tmp[,CFR:=mort*(el$mort.nh.num+el$mort.h.num)/(inc+1e-15)] #CFRs by age gender (total)
    tmph <- copy(tmp);
    tmpnh <- copy(tmp)        #separate copies for HIV+/- below
    tmp <- tmp[,.(iso3,sex,age_group,mort,src.mort,CFR)]
    ## --- HIV + / - splits separately (for tables -- see below)
    ## -- HIV +ve
    tmph[,mort:=0.0]
    tmph[age_group=='0_4' & sex=='M',
         mort := mfy * km[,deathsY*h] ]   #young,m
    tmph[age_group=='0_4' & sex=='F',
         mort := (1-mfy) * km[,deathsY*h] ]     #young,f
    tmph[age_group=='5_14' & sex=='M',
         mort := mfo * km[,(deaths-deathsY)*h] ] #old,m
    tmph[age_group=='5_14' & sex=='F',
         mort := (1-mfo) * km[,(deaths-deathsY)*h] ]#old,f
    tmph[,mort:=mort/(sum(mort)+1e-15)]           #safety
    PK <- km$deathsH/(el$mort.h.num+1e-15) #prop deaths in kids 
    tmph[,mort:=PK*mort]
    tmph[age_group %ni% kds & sex=='F',mort:= hf * inc]
    tmph[age_group %ni% kds & sex=='M',mort:= hm * inc]
    fac <- tmph[age_group %ni% kds,sum(mort)]
    tmph[age_group %ni% kds,mort:= (1-PK)*mort/fac]
    ## -- HIV -ve
    tmpnh[,mort:=0.0]
    tmpnh[age_group=='0_4' & sex=='M',
          mort := mfy * km[,deathsY*(1-h)] ]   #young,m
    tmpnh[age_group=='0_4' & sex=='F',
          mort := (1-mfy) * km[,deathsY*(1-h)] ]     #young,f
    tmpnh[age_group=='5_14' & sex=='M',
          mort := mfo * km[,(deaths-deathsY)*(1-h)] ] #old,m
    tmpnh[age_group=='5_14' & sex=='F',
          mort := (1-mfo)*km[,(deaths-deathsY)*(1-h)] ]#old,f
    tmpnh[,mort:=mort/(sum(mort)+1e-15)]   #safety
    PK <- km[,deaths-deathsH]/(el$mort.nh.num+1e-15)                #prop deaths in kids
    if(PK>=1){ 
        cat(cn,':PK>=1: using regional non-VR average!\n')
        elr <- merge(elr,KM,by='iso3',all.y=FALSE)
        PK <- elr[,mean((deaths-deathsH)/(mort.nh.num+1e-15),na.rm=TRUE)]
    }
    tmpnh[,mort:=PK*mort]
    tmpnh[age_group %ni% kds & sex=='F',mort:= (1-hf) * inc]
    tmpnh[age_group %ni% kds & sex=='M',mort:= (1-hm) * inc]
    fac <- tmpnh[age_group %ni% kds,sum(mort)]
    tmpnh[age_group %ni% kds,mort:= (1-PK)*mort/fac]
    ## plot
    mp <- mpnh <- mph <- NA
    mp <- muplot(tmp,cn)
    if(allplt){
        mpnh <- muplot(tmpnh,cn) #el[,country]
        mph <- muplot(tmph,cn)
    }
    ans <- list(
        dat=tmp,
        HN=tmpnh[,.(iso3,sex,age_group,mort,src.mort)],
        HP=tmph[,.(iso3,sex,age_group,mort,src.mort)],
        plt=mp)
    if(allplt){
        ans[['pltnh']] <- muplot(ans$HN,cn)
        ans[['plth']] <- muplot(ans$HP,cn)
    }
    ans
}


## testing
## getCFRMortH('AGO')
## tmp <- getCFRMortH('MOZ',allplt = TRUE)
## tmp$pltnh
## tmp$plt
## tmp$HN[age_group %in% kds,sum(mort)]
## muplot(tmp$HN,'ZAF, HIV -ve')
## muplot(tmp$HP,'ZAF, HIV +ve')
#'
#' # Unified function and application
#'
#' This just wraps the above depending on data availability:
#' 

getMort <- function(cn,verbose=FALSE,allplt=FALSE){
    VRagogo <- cn %in% VR[!is.na(g),unique(iso3)]  #
    ans <- getCFRMortH(cn,allplt = allplt)
    if(VRagogo){
        if(verbose) cat("Using VR...\n")
        ansVR <- getVRMortnh(cn,allplt = allplt) 
        ansVR$HP <- ans$HP
        ansVR$plth <- ans$plth
        return(ansVR)
    } else {
        if(verbose) cat("Using CFR...\n")
        return(ans)
    }
}

#' Testing:
#' 
## getMort('AGO')
## tmp <- getCFRMortH('MOZ',allplt = TRUE)

tmp <- getMort('MOZ',verbose = TRUE,allplt = TRUE)
print(tmp$plt)
print(tmp$pltnh)
tmp$HN
muplot(tmp$HN,'MOZ test')
tmp$HN[age_group %in% kds,sum(mort)]

## getMort('ZAF',verbose=TRUE)             #CFR w/HIV NB VR from IHME
## getCFRMortH('ZAF')                       #CFR w/o HIV
## getVRMortnh('ZAF')                        #VR HIV
## getVRMortnh('GBR')  

#' The countries to work on

CNZ <- unique(incsplit$iso3)
notinK <- setdiff(CNZ,KM$iso3)          #these have missing notification data
notinK                                  #
CNZ <- CNZ[CNZ %ni% notinK]

#' Countries with missing estimates to be excluded
#' 

missco <- est[year==estyr][,.(iso3,inc.num,inc.h.num,inc.nh.num,
                              mort.h.num,mort.nh.num)]
missco[,tot:=inc.num+inc.h.num+inc.nh.num+mort.h.num+mort.nh.num]
missco <- missco[is.na(tot),iso3]
missco
CNZ <- setdiff(CNZ,missco)

#' Looping through countries
## all

HPsplit <- HNsplit <- mortsplit <- allcns <-
    allcnsNH <- allcnsH <- list()
length(CNZ)        #
for(cn in CNZ){
    print(cn)
    tmp <- getMort(cn,allplt=TRUE)
    mortsplit[[cn]] <- tmp$dat
    HNsplit[[cn]] <- tmp$HN
    HPsplit[[cn]] <- tmp$HP
    allcns[[cn]] <- tmp$plt
    allcnsH[[cn]] <- tmp$plth
    allcnsNH[[cn]] <- tmp$pltnh
}
allcns[["CIV"]] <- allcns[["CIV"]] + ggtitle("Cote d'Ivoire") #avoid encoding error
allcnsH[["CIV"]] <- allcnsH[["CIV"]] + ggtitle("Cote d'Ivoire") #avoid encoding error
allcnsNH[["CIV"]] <- allcnsNH[["CIV"]] + ggtitle("Cote d'Ivoire") #avoid encoding error

length(CNZ)
length(HNsplit)
length(HPsplit)
## names(HPsplit)

## save(allcns,file=here('disaggregation/output/mortsplits/allcns.Rdata'))
save(allcnsH,
     file=here('disaggregation/output/mortsplits/plots/allcnsH.Rdata'))
save(allcnsNH,
     file=here('disaggregation/output/mortsplits/plots/allcnsNH.Rdata'))
save(HNsplit,
     file=here('disaggregation/output/mortsplits/plots/HNsplitL.Rdata'))
save(HPsplit,
     file=here('disaggregation/output/mortsplits/plots/HPsplitL.Rdata'))

#' Country plots

## HIV-ve
fn <- glue(here('disaggregation/output/mortsplits/plots/nhMort_'))
j <- k <- 0
for(i in seq_along(names(HNsplit))){
  if(k==0) plots <- list()  # new empty list
  k <- k+1
  plots[[k]] <- allcnsNH[[i]]
  if(k==30 | i==length(CNZ)){
    W <- H <- 15
    ## if(k!=30){ W <- 15; H <- 10}
    j <- j+1
    print(j)
    GP <- ggarrange(plotlist=plots, ncol = 5,nrow=6)
    fnl <- fn + j + '.pdf'
    if(gypt) ggsave(filename=fnl,GP,height=H,width=W,device=cairo_pdf)
    k <- 0
  }
}

## HIV+ve
fn <- glue(here('disaggregation/output/mortsplits/plots/hMort_'))
j <- k <- 0
for(i in seq_along(names(HPsplit))){
  if(k==0) plots <- list()  # new empty list
  k <- k+1
  plots[[k]] <- allcnsH[[i]]
  if(k==30 | i==length(CNZ)){
    W <- H <- 15
    ## if(k!=30){ W <- 15; H <- 10}
    j <- j+1
    print(j)
    GP <- ggarrange(plotlist=plots, ncol = 5,nrow=6)
    fnl <- fn + j + '.pdf'
    if(gypt) ggsave(filename=fnl,GP,height=H,width=W,device=cairo_pdf)
    k <- 0
  }
}


#'
#' 
#' save data:
## HNsplit <- do.call('rbind',HNsplit)

HNsplit <-  rbindlist(HNsplit,fill=TRUE)
HNsplit[is.na(age),age:=rep(agz4,nrow(HNsplit[is.na(age)])/length(agz4))]
unique(HNsplit[,.(age_group,age)])

#' Check
#' 

HNsplit[is.na(mort)]
HNsplit[is.na(mort),mort:=0]


HNsplit <- merge(HNsplit,est[year==estyr][,.(iso3,g.whoregion)],
                 by='iso3',all.y = FALSE)
save(HNsplit,file=here('disaggregation/output/mortsplits/HNsplit.Rdata'))


HPsplit <- do.call('rbind',HPsplit)
HPsplit <- merge(HPsplit,est[year==estyr][,.(iso3,g.whoregion)],
                 by='iso3',all.y = FALSE)
HPsplit[is.na(mort)]
HPsplit[is.na(mort),mort:=0]
HPsplit[mort<0]
HPsplit[mort<0,mort:=0]
save(HPsplit,file=here('disaggregation/output/mortsplits/HPsplit.Rdata'))


#' Checking:

hbc <- as.character(hbcsh[order(as.character(name)),iso3])
hbcn <- as.character(hbcsh[order(as.character(name)),name])
hbc <- c(t(matrix(hbc,byrow = TRUE,ncol=5)))         #re-order for plot
hbcn <- c(t(matrix(hbcn,byrow = TRUE,ncol=5)))         #re-order for plot



#' Special plots for 30 HBCs - straight to HIV-ve
## 30 HBC

pltlst <- list()
for(i in 1:30){
  plt <- allcnsNH[[hbc[i]]] + theme(legend.position="none")
  if(!(i%%5==1)) plt <- plt + theme(axis.text.y=element_blank())
  pltlst[[i]] <- plt + ggtitle(hbcn[i])
}


## Order required
## AGO, BGD, BRA, KHM, CAR
## CHN, Congo, PRK, DRC, ETH
## IND, IDN, KEN, LSO, LIB
## MOZ, MMR, NMB, NIG, PAK
## PNG, PHL, RUS, SLE, ZAF,
## THA, TZA, VNM, ZMB, ZWE

H <- 15
W <- H*.75
GP <- ggarrange(plotlist = pltlst, ncol = 5,nrow=6)
fn <- here('disaggregation/output/mortsplits/plots/mu30_NH.pdf')
if(gypt) ggsave(filename=fn,GP,height=H,width=W,device=cairo_pdf)




## resizing
wz <- rep(1,5);wz[1] <- 1.05
fn <- here('disaggregation/output/mortsplits/plots/mu30_NHb.pdf')
GP <- grid.arrange(pltlst[[1]],pltlst[[2]],pltlst[[3]],pltlst[[4]],
                   pltlst[[5]],
                   pltlst[[6]],pltlst[[7]],pltlst[[8]],pltlst[[9]],
                   pltlst[[10]],
                   pltlst[[11]],pltlst[[12]],pltlst[[13]],pltlst[[14]],
                   pltlst[[15]],
                   pltlst[[16]],pltlst[[17]],pltlst[[18]],pltlst[[19]],
                   pltlst[[20]],
                   pltlst[[21]],pltlst[[22]],pltlst[[23]],pltlst[[24]],
                   pltlst[[25]],
                   pltlst[[26]],pltlst[[27]],pltlst[[28]],pltlst[[29]],
                   pltlst[[30]],
             ncol=5,widths=wz)
if(gypt) ggsave(filename=fn,GP,height=H,width=W,device=cairo_pdf)

## rotated version
Pltlst <- pltlst
for(i in 1:length(Pltlst)){
  Pltlst[[i]] <- Pltlst[[i]]+ coord_flip() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
}

fn <- here('disaggregation/output/mortsplits/plots/mu30_NHc.pdf')
GP <- grid.arrange(Pltlst[[1]],Pltlst[[2]],Pltlst[[3]],
                   Pltlst[[4]],Pltlst[[5]],
                   Pltlst[[6]],Pltlst[[7]],Pltlst[[8]],
                   Pltlst[[9]],Pltlst[[10]],
                   Pltlst[[11]],Pltlst[[12]],Pltlst[[13]],
                   Pltlst[[14]],Pltlst[[15]],
                   Pltlst[[16]],Pltlst[[17]],Pltlst[[18]],
                   Pltlst[[19]],Pltlst[[20]],
                   Pltlst[[21]],Pltlst[[22]],Pltlst[[23]],
                   Pltlst[[24]],Pltlst[[25]],
                   Pltlst[[26]],Pltlst[[27]],Pltlst[[28]],
                   Pltlst[[29]],Pltlst[[30]],
             ncol=5,widths=wz)
if(gypt) ggsave(filename=fn,GP,height=H,width=W,device=cairo_pdf)



#' NB: countries like this (LSO), that appear to have better child mortality proportions,
#' despite poor CDRs are because the adult CFR is so much worse. TODO: consider a CFR in children that reflects the overall health system performance and the adult CFR
## getMort("LSO")                          #

#'
#' # Regional aggregations
#'
#' Regional weighted by mortality
#'

## HIV-ve only
HNsplit[,age:=gsub('_','-',age_group)]
HNsplit <- merge(HNsplit,
                 est[year==estyr][,.(iso3,mort.nh.num)],by='iso3')
mortregHN <- HNsplit[,.(mort=weighted.mean(mort,w=mort.nh.num)),
                     by=.(sex,age,g.whoregion)]
mortgloHN <- HNsplit[,.(mort=weighted.mean(mort,w=mort.nh.num)),
                     by=.(sex,age)]

mortgloHN[,mort:=mort/sum(mort)]
mortgloHN[age=='65plus',age:=rev(agz4)[1]]
mortregHN[,mort:=mort/sum(mort)]
mortregHN[age=='65plus',age:=rev(agz4)[1]]

#'
#' ## Plotting
#'
#' Global:

f <- 1.5
hh <- f*7.5/2; ww <- f*10/3;
fn <- here('disaggregation/output/mortsplits/plots/muGlobal_NH.pdf')
GP <- muplot(mortgloHN,'Global')
if(gypt) ggsave(filename=fn,GP,height=hh,width=ww,device=cairo_pdf)

## rotated
fn <- here('disaggregation/output/mortsplits/plots/muGlobal_NHc.pdf')
GP <- muplot(mortgloHN,'Global') +
  coord_flip() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
if(gypt) ggsave(filename=fn,GP,height=hh,width=ww,device=cairo_pdf)

## save out data for report plots
Mglobsplt <- mortgloHN
fn <- here('disaggregation/reportoutput/Mglobsplt.rda')
save(Mglobsplt,file=fn)

#' Regional:
#' 

regs <- sort(as.character(est[,unique(g.whoregion)]))
whozt <- c('Africa','The Americas','Eastern Mediterranean','Europe',
           'South-East Asia',
           'Western Pacific')

## save out data for report plots
Mregsplt <- mortregHN
fn <- here('disaggregation/reportoutput/Mregsplt.rda')
save(Mregsplt,file=fn)


#' Regional HN:
#' 

pltlst <- list()
for(i in 1:6){
  plt <- muplot(mortregHN[g.whoregion==regs[i]],whozt[i]) +
    theme(legend.position="none")
  if(!(i%%3==1)) plt <- plt + theme(axis.text.y=element_blank())
  pltlst[[i]] <- plt ## + ggtitle(regs[i])
}

## save out
fn <- here('disaggregation/output/mortsplits/plots/muRegional_HN.pdf')
GP <- ggarrange(plotlist = pltlst, ncol = 3,nrow=2)
ggsave(GP,filename=fn,h=7.5,w=10,device=cairo_pdf)


## rotated
Pltlst <- pltlst
for(i in 1:length(Pltlst)){
  Pltlst[[i]] <- Pltlst[[i]]+ coord_flip() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
}

## save out
fn <- here('disaggregation/output/mortsplits/plots/muRegional_HNc.pdf')
GP <- ggarrange(plotlist = Pltlst, ncol = 3,nrow=2)
ggsave(GP,filename=fn,h=7.5,w=10,device=cairo_pdf)


#'
#' ## Tables
#' 
#' What are the sources of mortality for countries with HIV?

isz <- muz <- as.character(HA[,unique(iso3)])
for(i in 1:length(muz)) muz[i] <- HNsplit[isz[i],src.mort][1]
srcs <- data.table(iso3=isz,src.mort=muz)
srcs <- srcs[!is.na(src.mort),]
kable(srcs[src.mort=='VR'])


#' The rest are CFR, so not outrageous to base HIV mortality on these (see above).
#'
#' 
#' HIV +ve:

HPsplit[,age:=age_group]
HPsplit[age_group %in% kds,age:='0-14']
HPsplit[age_group %ni% kds,age:='15plus']
HPsplit <- HPsplit[,.(mort=sum(mort)),by=.(iso3,sex,age,g.whoregion)]

## odd countries:
HPsplit[,keep:=is.finite(sum(mort)),by=iso3]
HPsplit[age=='0-14' & sex=='M' & keep==FALSE,iso3] #oddities
HPsplit <- HPsplit[keep==TRUE,]

## aggregate
HPsplit <- merge(HPsplit,est[year==estyr][,.(iso3,mort.h.num=mort.h.num)],by='iso3')
HPreg <- HPsplit[,.(mort=weighted.mean(mort,w=mort.h.num,na.rm=TRUE)),
                 by=.(sex,age,g.whoregion)]
HPglo <- HPsplit[,.(mort=weighted.mean(mort,w=mort.h.num,na.rm=TRUE)),
                 by=.(sex,age)]


#' HIV -ve:

HNsplit[,age:=age_group]
HNsplit[age_group %in% kds,age:='0-14']
HNsplit[age_group %ni% kds,age:='15plus']
HNsplit <- HNsplit[,.(mort=sum(mort)),by=.(iso3,sex,age,g.whoregion)]

## odd countries:
HNsplit[,keep:=is.finite(sum(mort)),by=iso3]
HNsplit[age=='0-14' & sex=='M' & keep==FALSE,iso3] #oddities
HNsplit <- HNsplit[keep==TRUE,]
## aggregate
HNsplit <- merge(HNsplit,est[year==estyr][,.(iso3,mort.nh.num)],
                 by='iso3')
HNreg <- HNsplit[,.(mort=weighted.mean(mort,w=mort.nh.num,na.rm=TRUE)),
                 by=.(sex,age,g.whoregion)]
HNglo <- HNsplit[,.(mort=weighted.mean(mort,w=mort.nh.num,na.rm=TRUE)),
                 by=.(sex,age)]


#' Now make tables that are consistent with Philippe's envelopes and uncertainty and
#' write out
#'
#'

## === HIV-ve
HNregb <- merge(HNreg,
                regional[year==estyr,
                         .(g.whoregion,mort.nh.num,
                           mort.nh.lo.num,mort.nh.hi.num)],
                all.y = TRUE)

HNregb[,mid:=mort.nh.num*mort]
## HNregb[,hi:=(mort.nh.hi.num)*mort]      #actually too narrow
## HNregb[,lo:=(mort.nh.lo.num)*mort]

HNregb[,hi:=mid+0.5*(mort.nh.hi.num-mort.nh.lo.num)*
            mort/sqrt(sum(mort^2)),
       by=g.whoregion]#w_i=p_i W/sqrt(sum(p^2))
HNregb[,lo:=mid-0.5*(mort.nh.hi.num-mort.nh.lo.num)*mort/
            sqrt(sum(mort^2)),
       by=g.whoregion]

## tests
tmp <- regional[year==estyr,
                .(g.whoregion,mort.nh.num,mort.nh.lo.num,mort.nh.hi.num)]
tmp2 <- HNregb[,.(mid=sum(mid),W=sqrt(sum((hi-lo)^2))),by=g.whoregion]
test <- merge(tmp,tmp2,by='g.whoregion')
test[,mid-mort.nh.num]
test[,(mort.nh.hi.num-mort.nh.lo.num)^2 - W^2]

## global
HNglo[,c("mort.nh.num","mort.nh.lo.num","mort.nh.hi.num"):=
           global[year==estyr,.(mort.nh.num,
                                mort.nh.lo.num,
                                mort.nh.hi.num)]]

HNglo[,mid:=mort.nh.num*mort]
HNglo[,hi:=mid+0.5*(mort.nh.hi.num-mort.nh.lo.num)*mort/
           sqrt(sum(mort^2))]
HNglo[,lo:=mid-0.5*(mort.nh.hi.num-mort.nh.lo.num)*mort/
           sqrt(sum(mort^2))]

## tests
HNglo[,sum(mid)]
HNregb[,sum(mid)]
HNglo[,sum((hi-lo)^2)]
HNglo[1,(mort.nh.hi.num-mort.nh.lo.num)^2]

HNrego <- dcast(HNregb,g.whoregion ~ sex+age,
                value.var = c('mid','lo','hi'))

HNregout <- cbind(
    regional[year==estyr][order(g.whoregion),
                          .(g.whoregion,
                            mid=mort.nh.num,
                            lo=mort.nh.lo.num,
                            hi=mort.nh.hi.num)],
    HNrego[order(g.whoregion),
         .(M014.mid=`mid_M_0-14`,M014.lo=`lo_M_0-14`,M014.hi=`hi_M_0-14`,
           F014.mid=`mid_F_0-14`,F014.lo=`lo_F_0-14`,F014.hi=`hi_F_0-14`,
           M15plus.mid=mid_M_15plus,M15plus.lo=lo_M_15plus,
           M15plus.hi=hi_M_15plus,
           F15plus.mid=mid_F_15plus,F15plus.lo=lo_F_15plus,
           F15plus.hi=hi_F_15plus)])

cbind(HNregout[,mid],HNregout[,M014.mid+F014.mid+
                               M15plus.mid+F15plus.mid])
HNregout[,.(sum(M014.mid),sum(F014.mid),sum(M15plus.mid),sum(F15plus.mid))]
HNUR <- copy(HNregout)

## apply rounding
HNregout[,2:ncol(HNregout):=lapply(.SD, ftb),
         .SDcols=2:ncol(HNregout)]

## write out
fn <- here('disaggregation/output/mortsplits/HNregional.csv')
fwrite(HNregout,file=fn)

## === HIV+ve
HPregb <- merge(HPreg,
                regional[year==estyr,
                         .(g.whoregion,mort.h.num,
                           mort.h.lo.num,
                           mort.h.hi.num)],
                all.y = TRUE)


HPregb[,mid:=mort.h.num*mort]
HPregb[,hi:=mid+0.5*(mort.h.hi.num-mort.h.lo.num)*mort/sqrt(sum(mort^2)),
       by=g.whoregion]#w_i=p_i W/sqrt(sum(p^2))
HPregb[,lo:=mid-0.5*(mort.h.hi.num-mort.h.lo.num)*mort/sqrt(sum(mort^2)),
       by=g.whoregion]

## tests
tmp <- regional[year==estyr,
                .(g.whoregion,mort.h.num,mort.h.lo.num,mort.h.hi.num)]
tmp2 <- HPregb[,.(mid=sum(mid),W=sqrt(sum((hi-lo)^2))),by=g.whoregion]
test <- merge(tmp,tmp2,by='g.whoregion')
test[,mid-mort.h.num]
test[,(mort.h.hi.num-mort.h.lo.num)^2 - W^2]

## global
HPglo[,c("mort.h.num","mort.h.lo.num","mort.h.hi.num"):=
         global[year==estyr,.(mort.h.num,mort.h.lo.num,mort.h.hi.num)]]

HPglo[,mid:=mort.h.num*mort]
HPglo[,hi:=mid+0.5*(mort.h.hi.num-mort.h.lo.num)*mort/sqrt(sum(mort^2))]
HPglo[,lo:=mid-0.5*(mort.h.hi.num-mort.h.lo.num)*mort/sqrt(sum(mort^2))]

## tests
HPglo[,sum(mid)]
HPregb[,sum(mid)]
HPglo[,sum((hi-lo)^2)]
HPglo[1,(mort.h.hi.num-mort.h.lo.num)^2]

HPrego <- dcast(HPregb,g.whoregion ~ sex+age,
                value.var = c('mid','lo','hi'))

HPregout <- cbind(
  regional[year==estyr][order(g.whoregion),.(g.whoregion,
                                             mid=mort.h.num,
                                             lo=mort.h.lo.num,
                                             hi=mort.h.hi.num)],
  HPrego[order(g.whoregion),
         .(M014.mid=`mid_M_0-14`,M014.lo=`lo_M_0-14`,M014.hi=`hi_M_0-14`,
           F014.mid=`mid_F_0-14`,F014.lo=`lo_F_0-14`,F014.hi=`hi_F_0-14`,
           M15plus.mid=mid_M_15plus,M15plus.lo=lo_M_15plus,
           M15plus.hi=hi_M_15plus,
           F15plus.mid=mid_F_15plus,F15plus.lo=lo_F_15plus,
           F15plus.hi=hi_F_15plus)])

cbind(HPregout[,mid],
      HPregout[,M014.mid+F014.mid+M15plus.mid+F15plus.mid])
HPregout[,.(sum(M014.mid),sum(F014.mid),
            sum(M15plus.mid),sum(F15plus.mid))]
HPUR <- copy(HPregout)

## round
HPregout[,2:ncol(HPregout):=lapply(.SD, ftb),
         .SDcols=2:ncol(HPregout)]

## write out
fn <- here('disaggregation/output/mortsplits/HPregional.csv')
fwrite(HPregout,file=fn)

## globals
HNglo[,g.whoregion:='Global']
HPglo[,g.whoregion:='Global']
HPgloo <- dcast(HPglo,g.whoregion ~ sex+age,
                value.var = c('mid','lo','hi'))
HNgloo <- dcast(HNglo,g.whoregion ~ sex+age,
                value.var = c('mid','lo','hi'))

HPglout <- cbind(
    global[year==estyr][,.(mid=mort.h.num,
                           lo=mort.h.lo.num,
                           hi=mort.h.hi.num)],
  HPgloo[,
         .(M014.mid=`mid_M_0-14`,M014.lo=`lo_M_0-14`,M014.hi=`hi_M_0-14`,
           F014.mid=`mid_F_0-14`,F014.lo=`lo_F_0-14`,F014.hi=`hi_F_0-14`,
           M15plus.mid=mid_M_15plus,M15plus.lo=lo_M_15plus,
           M15plus.hi=hi_M_15plus,
           F15plus.mid=mid_F_15plus,F15plus.lo=lo_F_15plus,
           F15plus.hi=hi_F_15plus)])

HNglout <- cbind(
    global[year==estyr][,.(mid=mort.nh.num,
                           lo=mort.nh.lo.num,
                           hi=mort.nh.hi.num)],
  HNgloo[,
         .(M014.mid=`mid_M_0-14`,M014.lo=`lo_M_0-14`,M014.hi=`hi_M_0-14`,
           F014.mid=`mid_F_0-14`,F014.lo=`lo_F_0-14`,F014.hi=`hi_F_0-14`,
           M15plus.mid=mid_M_15plus,M15plus.lo=lo_M_15plus,
           M15plus.hi=hi_M_15plus,
           F15plus.mid=mid_F_15plus,F15plus.lo=lo_F_15plus,
           F15plus.hi=hi_F_15plus)])


## round
HPglout[,1:ncol(HPglout):=lapply(.SD, ftb), .SDcols=1:ncol(HPglout)]
HNglout[,1:ncol(HNglout):=lapply(.SD, ftb), .SDcols=1:ncol(HNglout)]

HNglout <- cbind(g.whoregion='Global',HNglout)
HPglout <- cbind(g.whoregion='Global',HPglout)

## test sums
HPUR[,.(mid=sum(mid),M014.mid=sum(M014.mid),F014.mid=sum(F014.mid),
        M15plus.mid=sum(M15plus.mid),F15plus.mid=sum(F15plus.mid))]
HPglout[,.(mid,M014.mid,F014.mid,M15plus.mid,F15plus.mid)]

## write out
fwrite(HNglout,
       file=here('disaggregation/output/mortsplits/HNglobal.csv'))
fwrite(HPglout,
       file=here('disaggregation/output/mortsplits/HPglobal.csv'))

#'
#' # Ad hoc additions and corrections
#'
#' This section adds in the missing countries as regional averages
#'

load(here('disaggregation/output/mortsplits/HNsplit.Rdata')) #to be overwritten
HNsplit <- merge(HNsplit,est[year==estyr][,.(iso3,mort.nh.num)],
                 by='iso3')
HNR <- HNsplit[,.(mort=weighted.mean(mort,w=mort.nh.num,na.rm=TRUE)),
               by=.(sex,age_group,g.whoregion)] #regional averages
extras <- est[year==estyr][iso3 %ni%HNsplit[,iso3],
                           .(iso3,g.whoregion,mort.nh.num)]
extras[,src.mort:='regional average']
extras <- merge(x=extras,y=HNR,by='g.whoregion',allow.cartesian = TRUE)
extras <- extras[,.(iso3,sex,age_group,mort,src.mort,g.whoregion)]

HNsplit <- rbind(HNsplit[,.(iso3,sex,
                            age_group,mort,
                            src.mort,g.whoregion)],extras)
HNsplit <- HNsplit[order(as.character(iso3)),]
save(HNsplit,
     file=here('disaggregation/output/mortsplits/HNsplit.Rdata'))





## ------- add uncertainty and output
HNsplit[,unique(src.mort)]              #check
HNsplit <- merge(HNsplit,
                 est[year==estyr][,.(iso3,mort.nh.num,
                                     mort.nh.lo.num,mort.nh.hi.num)],
                 by='iso3')

## add uncertainty
HNsplit[,mid:=mort*mort.nh.num]
HNsplit[,F:=sqrt(mort*(mort.nh.hi.num-mort.nh.lo.num)^2) ]
## HNsplit[,lo:=pmax(mid-0.5*F,0)]
HNsplit[,hi:=mid+0.5*F]
HNsplit[,lo:=mid-0.5*F]
HNsplit[lo<0,hi:=hi-lo]
HNsplit[lo<0,lo:=0]


## check:
## midpoint
test <- HNsplit[,.(tt=(mort.nh.num[1] - sum(mid))^2),by=iso3]
test[,summary(tt)]
(nac <- test[is.na(tt)])

## uncertainty
test <- HNsplit[,.(tt=(hi-lo)^2,tt2=(mort.nh.hi.num-mort.nh.lo.num)^2),
                by=iso3]
test[,summary(tt)]
test[is.na(tt),unique(iso3)]
test2 <- test[,.(v1=sum(tt),v2=tt2[1]),by=iso3]
test2[,summary(v1-v2)]

HNsplit[,best:=mid]
HNsplit[,year:=estyr]
HNsplit[,measure:='mortality']
HNsplit[,unit:='number']
HNsplit[,age_group:=gsub('_','-',age_group)]
HNsplit[,sex:=tolower(sex)]
HNsplit[,c('F','mort.nh.hi.num','mort.nh.lo.num','mort.nh.num',
           'mort','mid','src.mort'):=NULL]

attr(HNsplit,'timestamp') <- Sys.Date()

db_hn_mortality_disaggregated <- HNsplit
attr(db_hn_mortality_disaggregated,'timestamp') <- Sys.Date()

fn <- here('disaggregation/output/mortsplits/db_hn_mortality_disaggregated.Rdata')
save(db_hn_mortality_disaggregated,file=fn)
fn <- here('disaggregation/dboutput/db_hn_mortality_disaggregated.Rdata')
save(db_hn_mortality_disaggregated,file=fn)

## also save out parts of HPsplit
MPglobsplit <- HPsplit
MPglobsplit[,acat:='15+']
## MPglobsplit <- merge(MPglobsplit,
##                      est[year==estyr,.(iso3,mort.h.num)],
##                      by='iso3')
## MPglobsplit[age %in% c('0-4','5-14'),acat:='0-14']
MPglobsplit[age %in% c('0-14'),acat:='0-14']
MPglobsplit[,deaths:=mort * mort.h.num]
MPglobsplit <- MPglobsplit[is.finite(deaths) & deaths>0,
                           .(deaths=sum(deaths)),
                           by=.(sex,acat)]
MPglobsplit[,deaths.pc:=1e2*deaths/sum(deaths)]
MPglobsplit[,deaths:=NULL]
attr(MPglobsplit,'timestamp') <- Sys.Date()
fn <- here('disaggregation/reportoutput/MPglobsplit.rda')
save(MPglobsplit,file=fn)


#' CFR checks:

load(here('disaggregation/output/mortsplits/HNsplit.Rdata'))
load(here('disaggregation/output/incsplits/data/incsplit.Rdata'))

HNsplit <- merge(HNsplit,est[year==estyr][,.(iso3,
                                             mort.nh.num,
                                             mort.nh.lo.num,
                                             mort.nh.hi.num)],
                 by='iso3')
HNsplit[,mid:=mort*mort.nh.num]


HNC <- merge(HNsplit[,.(iso3,age_group,sex,mid)],
             incsplit[,.(iso3,age_group,sex,inc)],
             by=c('iso3','age_group','sex'))

HNC[,CFR:=mid/inc]
HNC[inc>0,summary(CFR)]

kable(HNC[inc>0 & (CFR>1)])
(odd <- HNC[inc>0 & (CFR>1),iso3])
## HNsplit[iso3 %in% odd]

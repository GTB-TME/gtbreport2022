#' ---
#' title: VR data
#' author: Philippe Glaziou, Anna Dean
#' date: 2022/06/07
#' output:
#'    html_document:
#'      mode: selfcontained
#'      toc: true
#'      toc_depth: 3
#'      toc_float: true
#'      number_sections: true
#'      theme: flatly
#'      highlight: zenburn
#'      df_print: paged
#'      code_folding: hide
#' ---

#' Last updated: `r Sys.Date()`
#'

# Load libraries and data
library(data.table)
library(imputeTS)
library(haven) # read_dta
library(readxl)
library(stringr)
library(here)



# load data
load(here('data/old.rda'))   
load(here('data/cty.rda'))
load(here('data/pop.rda'))
load(here('data/unaids.rda'))  #!!! introduced dep on 04/08/2021 
load(here('data/ovr.rda'))  # vr from last year (this year, data from 2000-2004 are missing from the excel file)

source(here('R/fun.R'))

m <- 1e5
yr <- 2021

# convert VR data in excel
# adapting code from Pete (2019)
#
# function for reformatting
refrm <- function(indat) {
  # rename & aggregate
  indat <- indat[, .(
    Country,
    name,
    Year,
    icd,
    Cause,
    cause1,
    Sex,
    `0-4` = rowSums(cbind(
      Deaths2, Deaths3, Deaths4, Deaths5, Deaths6
    ), na.rm = T),
    `5-14` = rowSums(cbind(Deaths7, Deaths8), na.rm = T),
    `15-24` = rowSums(cbind(Deaths9, Deaths10), na.rm = T),
    `25-34` = rowSums(cbind(Deaths11, Deaths12), na.rm = T),
    `35-44` = rowSums(cbind(Deaths13, Deaths14), na.rm = T),
    `45-54` = rowSums(cbind(Deaths15, Deaths16), na.rm = T),
    `55-64` = rowSums(cbind(Deaths17, Deaths18), na.rm = T),
    `65plus` = rowSums(
      cbind(
        Deaths19,
        Deaths20,
        Deaths21,
        Deaths22,
        Deaths23,
        Deaths24,
        Deaths25
      ),na.rm = T))
    ]

  # Sequelae
  seq <- c("B90", "B900", "B901", "B902", "B908", "B909", "B077")
  indat[Cause %in% seq, cause1 := "tbseq"]
  indat[, Cause := NULL]

  ## reshape
  MM <-
    melt(indat, id = c("Country", "name", "Year", "icd", "cause1", "Sex"))
  MM$sex <- c('M', 'F')[as.numeric(MM$Sex)]
  MM[is.na(sex), sex := 'U']
  MM$sex <- factor(MM$sex)
  MM$age <- factor(MM$variable, levels = agz3, ordered = TRUE)
  MM[, age_group := gsub('-', '_', age)]
  MM[, age := NULL]
  MM
}

# Some useful age range vectors:
agz <-
  c('04', '514', '1524', '2534', '3544', '4554', '5564', '65') #for extract
agz2 <-
  c('0_4',
    '5_14',
    '15_24',
    '25_34',
    '35_44',
    '45_54',
    '55_64',
    '65plus') #for labels
agz3 <- gsub('_', '-', agz2)
agz4 <- c(rev(rev(agz3)[-1]), "\u2265 65")
kds <- c('0_4', '5_14')
kdz <- c('04', '514')
AA <-
  data.table(
    a1 = agz,
    age_group = agz2,
    age = agz3,
    agegt = agz4
  ) #for conversion


typz <- c('text',
          'text',
          'text',
          'text',
          rep('text', 5),
          rep('numeric', 26)) #specify column types to avoid warning
xlfn <- 'input/VR/info_TB_june2022.xls'



## -------- ICD>10.1
M1 <-
  as.data.table(read_excel(
    xlfn,
    sheet = 3,
    #skip = 1,
    col_names = TRUE,
    col_types = typz
  ))
M1 <- refrm(M1)


## -------- ICD 10.1
M2 <-
  as.data.table(read_excel(
    xlfn,
    sheet = 4,
    #skip = 1,
    col_names = TRUE,
    col_types = typz
  ))
M2 <- refrm(M2)


## -------- ICD 9
M3 <-
  as.data.table(read_excel(
    xlfn,
    sheet = 5,
    #skip = 1,
    col_types = typz
  ))
# remove double counted B02
M3[, dble := sums(str_count(Cause, 'B02[0-9]')), by = .(name, Year, Sex)]
M3[Cause == 'B02' & dble > 0, drop := TRUE]
(M3[, table(drop)])
M3 <- M3[is.na(drop)]
M3[, drop := NULL]
M3[, dble := NULL]
M3 <- refrm(M3)


## -------- ICD 8
M4 <-
  as.data.table(read_excel(
    xlfn,
    sheet = 6,
    #skip = 1,
    col_types = typz
  ))
M4 <- refrm(M4)


## --- join ---
VR <- rbind(M1, M2, M3, M4)


## Differences in names:
(vrbad <- setdiff(VR[, unique(name)],
                  pop[, unique(country)]))

## direct renamings:
done <- c()
for (i in seq_along(vrbad)) {
  newnm <- grep(vrbad[i], pop[, unique(country)], value = TRUE)
  if (length(newnm) > 0) {
    print(newnm)
    VR[name == vrbad[i], name := newnm]
    done <- c(done, i)
  }
}
vrbad <- vrbad[-done]

## others for renaming
vrbad
(newnm <- grep('Czech', pop[, unique(country)], value = TRUE))
VR[name == grep('Czech', vrbad, value = TRUE), name := newnm]
(newnm <- grep('Serbia', pop[, unique(country)], value = TRUE)[1])
#VR[name == grep('Serb', vrbad, value = TRUE), name := newnm]
(newnm <-
    grep('Macedonia', pop[, unique(country)], value = TRUE)[1])
VR[name == grep('Mace', vrbad, value = TRUE), name := newnm]
(newnm <- grep('Vincent', pop[, unique(country)], value = TRUE)[1])
VR[name == grep('Vincent', vrbad, value = TRUE), name := newnm]
(newnm <- grep('Libya', pop[, unique(country)], value = TRUE))
VR[name == grep('Libya', vrbad, value = TRUE), name := newnm]
(newnm <- grep('Palestinian', pop[, unique(country)], value = TRUE))
VR[name == grep('Palestinian', vrbad, value = TRUE), name := newnm]

## sub-countries
VR[name %in% c(
  "French Guiana",
  "Martinique",
  "Reunion",
  "Mayotte",
  "Guadeloupe",
  "Saint Pierre and Miquelon"
), name := "France"]
VR[name %in% c("Rodrigues"), name := "Mauritius"]
VR[name %in% c("Virgin Islands (USA)"), name := "United States of America"]

## check
(dropname <- setdiff(VR[, unique(name)], cty[, unique(country)]))
VR <- VR[name %ni% dropname]

VR[, year := as.integer(Year)]
VR[, Year := NULL]

## Add iso3, and tidy up
VR <- merge(VR,
            cty[, .(iso3, country)],
            by.x = 'name',
            by.y = "country",
            all.x = TRUE)
(VR[is.na(iso3), unique(name)])

VR[, Sex := NULL]
VR[, age_group := variable]
VR$age_group <-
  factor(gsub("-", "_", VR$age_group),
         levels = agz2,
         ordered = TRUE)
VR[, variable := NULL]

setkey(VR, iso3)
rm(M1, M2, M3, M4)

# save(VR, file = 'data/VR.rda')



# aggregate and reshape long to wide
#
vr <- VR[, .(value = sums(value)),
         by = .(iso3, year, cause1)]

vr <- dcast(vr, ... ~ cause1)
setkey(vr, iso3)
setnames(vr, c('iso3','year','hiv','tb','ill_def','tb_seq','total'))





# # Country data additions to WHO database
#

# * Azerbaijan additions
#
# source: epi review, the Hague 2017
#

vr['AZE']
NewYEARS <- 2010:2015
addAZE <- as.data.frame(matrix(nrow=length(NewYEARS),ncol=ncol(vr)))
names(addAZE) <- names(vr)
addAZE$iso3 <- "AZE"
addAZE$year[1:nrow(addAZE)] <- NewYEARS

dim(vr)
vr2 <- rbind(vr, addAZE, use.names = TRUE)
setkey(vr2, iso3)
sel <- vr2$iso3 == 'AZE' & vr2$year %in% 2010:2015

# vr2$vr.coverage[sel] <- rep(1, 6)
vr2$total[sel] <- c(53580, 53726, 55017, 54383, 55648, 54697)
vr2$tb[sel] <- c(709, 577, 373, 378, 372, 485)
vr2$ill_def[sel] <- c(1343, 1771, 1836, 1892, 2440, 1864)
vr2['AZE']

# checking that each combination iso3-year appears only once in the database
tCheck <- as.data.frame(table(paste(vr2$iso3,vr2$year)))
nrow(tCheck[ tCheck[,2]>1 , ])==0

#!!!
# check whether the following 2 blocks do not overwrite newly reported VR data
#!!!

# RUS 2020
#
# Total deaths = 2138586
# Ill-defined = 142370
# TB deaths = 6841
addRUS <- vr2['RUS'][year==2019]
addRUS$year <- 2020
addRUS$total <- 2138586
addRUS$ill_def <- 142370
addRUS$tb <- 6841
addRUS$hiv <- unaids['RUS'][year==2020, mort.hiv.num]
vr3 <- rbind(vr2, addRUS, use.names = TRUE)
setkey(vr3, iso3)

vr <- copy(vr3)
(dim(vr))
# 
# 
# 
# RUS 2021
#
# Total deaths =
# Ill-defined =
# TB deaths =
addRUS <- vr['RUS'][year==2020]
addRUS$year <- 2021
addRUS$total <- 2441594
addRUS$ill_def <- 137231
addRUS$tb <- 6313
addRUS$hiv <- unaids['RUS'][year==2021, mort.hiv.num]
vr4 <- rbind(vr, addRUS, use.names = TRUE)
setkey(vr4, iso3)

vr <- copy(vr4)
(dim(vr))





# check that the additions do not end-up duplicating country-year entries
#
sum(duplicated(vr[,.(iso3,year)]))==0



# # process raw VR data
#
(dim(vr))

# missing tb deaths
miss.tb <- unique(vr$iso3[is.na(vr$tb)])

vr[iso3 %in% miss.tb, tb.nm := sum(!is.na(tb)), by = iso3]
vr[, table(tb.nm)]
vr[tb.nm == 0]

(dim(vr))
vr <- vr[tb.nm > 0 | iso3 %ni% miss.tb]
(dim(vr))

vr[tb.nm == 1] # not yet imputable

vr[, prop.tb := tb / total]
vr[tb.nm >= 2, tb.imp := na_interpolation(round(prop.tb * total)), by =
     iso3]
(summary(vr$tb.imp))
vr[!is.na(tb) & is.na(tb.imp), tb.imp := tb]
(vr['ZAF'])
(vr['RUS'])

# clean up
#
vr[, tb.nm := NULL]
rm(miss.tb)



# missing ill-defined
miss.ill <- unique(vr$iso3[is.na(vr$ill_def)])

vr[iso3 %in% miss.ill, ill.nm := sum(!is.na(ill_def)), by = iso3]
(vr[, table(ill.nm)])

vr[, garbage := ill_def / total]

vr[ill.nm >= 2, garbage.imp := na_interpolation(garbage), by = iso3]
(summary(vr$garbage.imp))
vr[!is.na(garbage) & is.na(garbage.imp), garbage.imp := garbage]
(vr['ZAF'])
(vr['RUS'])



# clean up
#
vr[, ill.nm := NULL]
rm(miss.ill)


# missing tb_seq
#
vr[, sum(is.na(tb_seq))]


# total TB deaths
#
vr[, totaltb := rowSums(cbind(vr$tb.imp, vr$tb_seq), na.rm = TRUE)]

# proportion sequelae
#
vr[, seq.prop := tb_seq / totaltb]


setkey(vr, iso3, year)



# proportion of TB deaths out of well documented deaths
#
vr[, tb.prop := totaltb / (total - ill_def)]



# check
#
vr[!is.na(tb.prop), test.isbinom(tb.prop)]





# # VR quality
vrqual <- fread('input/VR/vrqual.csv')
vrqual <- vrqual[country != ""]
(vrqual)


## Differences in names:
(vrqbad <- setdiff(vrqual[, unique(country)],
                  cty[, country]))

(newnm <- grep('Macedonia', cty[, country], value = TRUE))
vrqual[country == grep('Macedonia', vrqbad, value = TRUE), country := newnm]
(newnm <- grep('Britain', cty[, country], value = TRUE))
vrqual[country == grep('United Kingdom', vrqbad, value = TRUE), country := newnm]


vrqual <- merge(vrqual, cty[,.(country, iso3)], by='country', all.x =TRUE)
sum(is.na(vrqual$iso3))==0 # check

vr <-
  merge(vr, vrqual[, .(iso3,
                       codqual = quality,
                       min.coverage,
                       max.coverage,
                       min.usability,
                       max.usability,
                       use.ghe)], by = 'iso3')


vr[codqual %in% c("high", "medium"), keep.vr := T]
vr[codqual %ni% c("high", "medium"), keep.vr := F]
vr[, summary(keep.vr)]

remove.pc <- function(x) as.numeric(gsub("%", "", x))

vr[, min.coverage := remove.pc(min.coverage)]
vr[, max.coverage := remove.pc(max.coverage)]
vr[, min.usability := remove.pc(min.usability)]
vr[, max.usability := remove.pc(max.usability)]

# vr[, vr.coverage := min.coverage]







# GHE
# dths1 = males
# dths2 = females
#
ghe <- fread('input/GHE/dths_total.csv')
ghe[causename=='All Causes', ghe.env := rowSums(cbind(dths1, dths2), na.rm=TRUE)]

vr2 <- merge(vr, ghe[!is.na(ghe.env),.(iso3, year, ghe.env)], by=c('iso3','year'), all.x=TRUE, all.y=FALSE)
(dim(vr))
(dim(vr2))


vr2[, env := ghe.env]                    # default envelope
vr2[iso3 %in% c('DMA','PRI','SMR','KNA'), env := total]

vr2[, env.nm := sum(!is.na(env)), by=iso3]
vr2[, table(env.nm)]
vr2[env.nm >=2, env := na_interpolation(env), by=iso3]
vr2[, sum(is.na(env))]


vr2[, vr.coverage := pmin(1, total / env)] # VR coverage
vr2[is.na(vr.coverage), vr.coverage := total] # non-GHE estimate where missing GHE


# update RUS
vr2[iso3=='RUS' & year==yr-1, ghe.env := total]
vr2[iso3=='RUS' & year==yr-1, env := total]






# proportion of TB deaths out of well documented deaths
#
vr2[, tb.prop := totaltb / (total - ill_def)]


#
# adjusted TB deaths
#
vr2[, tb.adj := env * tb.prop]

vr <- copy(vr2)

# check
#
vr[!is.na(tb.prop), test.isbinom(tb.prop)]




# # SDs
#
# assume TB deaths between 0.5 and 1.5 times observed $t$ rate among
# garbage $g$ and non covered $c$
#
# $t_{adj} = \frac{t}{c(1-g)}$
#
# $\text{SD}(t_{adj}) = \frac{t}{4} \left(\frac{1}{c(1-g)} - 1\right)$
#
vr[, tb.adj.sd := (tb.adj / 4) * (1 / (vr.coverage * (1 - garbage)) - 1)]

vr[!is.na(tb.adj.sd), test.ispos(tb.adj.sd)]
vr[keep.vr == T & is.na(tb.adj.sd), tb.adj.sd := tb.adj * .2]
vr[keep.vr == T, summary(tb.adj.sd / tb.adj)]


vr2  <-
  merge(vr[iso3 %ni% 'VIR'], pop[, .(iso3, year, pop = e.pop.num)], by =
          c('iso3', 'year'), all.x = TRUE)
dim(vr[iso3 %ni% 'VIR'])
dim(vr2)
vr2[is.na(pop)]

# exclude MNE and SRB 2000:2004
vr <- copy(vr2[!is.na(pop)])







# Codecorrect from IHME (GBD 2020, accessed June 2021)
#
load('data/ihmetb.rda')   # see ihme.R
load('data/ihmeall.rda')


# envelope ratios GHO/IHME
#
ihme.all <-
  ihmeall[year >= 2000 &
         cause == 'All causes' & measure == 'Deaths']
ihme.all <- ihme.all[!is.na(iso3)]
gbd <-
  merge(ihme.all[, .(iso3,
                     year,
                     ihme = val,
                     ihme.lo = lower,
                     ihme.hi = upper)],
        ghe,
        by = c('iso3', 'year'),
        all.y = T)


# missing IHME envelopes?
#
gbd[, sapply(.SD, function(x)
  sum(is.na(x)))]


# WHO/IHME env ratio
#
gbd[, env.ratio := ghe.env / ihme]

vr2 <- merge(vr, gbd[!is.na(env.ratio),.(iso3,year,ihme.env=ihme,env.ratio)], by=c('iso3','year'), all.x=TRUE)
(dim(vr))
(dim(vr2))

vr <- copy(vr2)

vr[, prop.tb := NULL]





#!!! 
# add missing 2000-2004 data 
# (next year, check for completeness of excel file)
dim(vr); dim(ovr[year<2005])
vr2 <- rbind(ovr[year<2005], vr)
dim(vr2)

vr <- copy(vr2)

sum(duplicated(vr[,.(iso3,year)]))==0 
setkey(vr, iso3, year)
#
#!!!



# simple checks
ovr[, table(year)]
vr[, table(year)]

ovr[, sum(tb, na.rm=TRUE), by=year]
vr[, sum(tb, na.rm=TRUE), by=year]




#  save vr dataset
#
save(vr, file = 'data/vr.rda')
fwrite(vr, file = paste0('csv/vr_', Sys.Date(), '.csv'))






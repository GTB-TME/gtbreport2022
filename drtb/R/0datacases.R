#' ---
#' title: Organizing DRTB data
#' author: Pete Dodd
#' date: 4 July, 2022
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' 
#' #Pre-amble
#' (Last updated: `r Sys.Date()`)
#' 
#' This is clearning and reshaping survey and surveillance data on drug-resistant TB for onward analysis
#'
#' N.B. This file should render to html with `rmarkdown::render('0datacases.R',output_dir='../html')` or from the command line with `R -q -e "rmarkdown::render(\"0datacases.R\",output_dir=\"../html\")"`
#'
#'
#' # Relevant libraries
#'
#' 

rm(list=ls())
library(here)
library(data.table)
library(ggplot2)
library(ggrepel)
library(glue)
gh <- function(x) glue(here(x))
rot45 <- theme(axis.text.x = element_text(angle = 45, hjust = 1))


#' Include subnational data from: PNG, BRA, CAF, PRK, but
#' we will exclude other subnational data
subnationalkeepers <- c("PNG", "BRA", "CAF", "PRK")

#' # INH-rated data
#'
#' Read in
## new
load(here('data/drhnew.rda'))
Hn <- drhnew
names(Hn) <- gsub('\\.','_',names(drhnew))

## ret
load(here('data/drhret.rda'))
Hr <- drhret
names(Hr) <- gsub('\\.','_',names(drhret))


#' ## ZAF INH data exclusion
#' 
#' Exclude INH data for South Africa for 2021 from surveillance
#'  because.
#' The data report that bacteriologically confirmed cases also have a test result for isoniazid. This is
#' not correct – many cases are not actually tested for isoniazid (discussion with NICD). Additionally, the reported
#' percentages with isoniazid resistance are quite low (new: 1.6%).
Hn <- Hn[!(iso3=='ZAF' & year_new==2021)]
Hr <- Hr[!(iso3=='ZAF' & year_ret==2021)]


#' # Keys & joining
#' 
#' Make data keys so can distinguish if already used comparing vs RR data
Hn[all_areas_covered_new==1 | iso3 %in% subnationalkeepers,
   id:=paste0('N_',ifelse(source_new=='Survey','Y','L'),
              '_',iso3,'_',year_new)]
Hn[,.(iso3,source_new,year_new,all_areas_covered_new,id)]
Hr[all_areas_covered_ret==1 | iso3 %in% subnationalkeepers,
   id:=paste0('R_',ifelse(source_ret=='Survey','Y','L'),
              '_',iso3,'_',year_ret)]
Hr[,.(iso3,source_ret,year_ret,all_areas_covered_ret,id)]

#' drop subnational (barring exceptions)
Hn <- Hn[!is.na(id)]
Hr <- Hr[!is.na(id)]

#' remove new/ret labels from names and combine
names(Hn) <- gsub("_new","",names(Hn))
names(Hr) <- gsub("_ret","",names(Hr))
Hn[,patients:='new']
Hr[,patients:='ret']
H <- rbind(Hn,Hr)

#' check have included relevant subnations
H[iso3 %in% subnationalkeepers]
subnationalkeepers

#' check for duplicate ids
H[,.(length(id),length(unique(id)))]

save(H,file=here('drtb/data/H.Rdata'))


#' # Rif-rated data
#'
#' Read in
## new
load(here('data/drnew.rda'))
Rn <- drnew
names(Rn) <- gsub('\\.','_',names(drnew))

## ret
load(here('data/drret.rda'))
Rr <- drret
names(Rr) <- gsub('\\.','_',names(drret))

#' make data keys so can distinguish if already used comparing vs RR data
Rn[all_areas_covered_new==1 | iso3 %in% subnationalkeepers,
   id:=paste0('N_',ifelse(source_new=='Survey','Y','L'),
              '_',iso3,'_',year_new)]
Rn[,.(iso3,source_new,year_new,all_areas_covered_new,id)]
Rr[all_areas_covered_ret==1 | iso3 %in% subnationalkeepers,
   id:=paste0('R_',ifelse(source_ret=='Survey','Y','L'),
              '_',iso3,'_',year_ret)]
Rr[,.(iso3,source_ret,year_ret,all_areas_covered_ret,id)]

#' drop subnational
Rn <- Rn[!is.na(id)]
Rr <- Rr[!is.na(id)]

#' remove new/ret labels from names and combine
names(Rn) <- gsub("_new","",names(Rn))
names(Rr) <- gsub("_ret","",names(Rr))
Rn[,patients:='new']
Rr[,patients:='ret']
R <- rbind(Rn,Rr)

#' check have included relevant subnations
R[iso3 %in% subnationalkeepers]
subnationalkeepers

#' check for duplicate ids
R[,.(length(id),length(unique(id)))]

save(R,file=here('drtb/data/R.Rdata'))

#' Fewer data points are rated for use to estimate HR than to estimate RR:
(HvR.rated <- c(nrow(R), nrow(H)))
cat(HvR.rated,file=here('drtb/data/HvR.rated.txt'))

rm(Rr,Rn,R,Hr,Hn,H)

#'
#' # Checking & labeling cases
#' 
#' Load
load(file=here('drtb/data/H.Rdata'))
load(file=here('drtb/data/R.Rdata'))

#' NOTE data edits
#' denominator but all answers NA - assume zero or NA?
H[source!='Survey' & is.na(r_rlt) & !is.na(dst_rlt) &
           is.na(dr_h_nr) & is.na(dr_r_nh) & is.na(mdr) &
           is.na(dst_rlt_hr),
  .(id,r_rlt,dst_rlt,dr_h_nr,dr_r_nh,mdr,dst_rlt_hr)]

#' Assume that for this entry, the NA numerators are in fact 0:
H[id=='R_L_MUS_2015',c('dr_h_nr','dr_r_nh','mdr'):=0]


#'
#' There are 3 different denominators, with potential numerators:
#'
#'
#'  dst_rlt:
#'    dr_h_nr, dr_r_nh, mdr,
#'    dst_rlt_hr, dst_rlt_rr, 
#'
#' r_rlt:
#'    rr
#'
#' (non-overlapping group so can always use)
#' xpert:
#'    xpert_dr_r
#'
#' Check their patterns of co-occurence:
#' 

H[source!='Survey',.N,
  by=.(dabs=is.na(dst_rlt),rabs=is.na(r_rlt),xabs=is.na(xpert))]

#' hr there only when dst & rrlt is:
H[source!='Survey',.N,
  by=.(dabs=is.na(dst_rlt),rabs=is.na(r_rlt),hrabs=is.na(dst_rlt_hr))]

H[source!='Survey' & !is.na(dst_rlt_hr),.N,
  by=.(is.na(dr_h_nr), is.na(dr_r_nh), is.na(mdr))]

H[source!='Survey' & !is.na(dst_rlt_hr) &
  !is.na(dr_h_nr)& !is.na(dr_r_nh)& !is.na(mdr)] #none

H[source!='Survey',.N,
  by=.(dabs=is.na(dst_rlt),rabs=is.na(r_rlt),xabs=is.na(xpert))]

#' Remove overlap between H-rated and R-rated data (ie since all H-rated is R-rated, exclude these records from R-rated data):
H[,Hrated:='yes']
R[,Hrated:='no']
R <- R[!id %in% H$id] #exclude overlap: don't need in R data if H-rated

#' Variable names to keep:
nmz <- c('iso3','year','patients','Hrated','id',
         'dr_h_nr', 'dr_r_nh', 'mdr',
         'dst_rlt',
         'dst_rlt_hr', 'dst_rlt_rr',
         'rr','r_rlt',
         'xpert_dr_r','xpert'
         )

#' Join into single surveillance dataset:
B <- rbind(H[source!='Survey',..nmz],R[source!='Survey',..nmz])


B[!is.na(dst_rlt_hr),.N,
  by=.(is.na(dr_r_nh), is.na(mdr))]


#' NOTE exclude where all denoms are only NA or 0
totden <- rowSums(as.matrix(B[,.(dst_rlt,r_rlt,xpert)]),na.rm=TRUE)
length(totden)
totden0 <- B[totden==0,.(iso3,year,patients,id,dst_rlt,r_rlt,xpert,
                         dr_h_nr, dr_r_nh, mdr,
                         dst_rlt_hr, dst_rlt_rr,
                         rr,xpert_dr_r)]

totden0[,length(unique(iso3))] #30
totden0[,length((iso3))] #108
fwrite(totden0,file=here('drtb/data/totden0.csv'))

#' Print for html output
totden0[,unique(paste0(iso3,':',year))]

#' NOTE boils down to dst_rlt
B <- B[!id %in% totden0$id] #NOTE drop!
nrow(B)

#' hr cases
B[Hrated=='yes',.N,by=.(is.na(r_rlt),is.na(dst_rlt_hr),
                        is.na(dr_h_nr),is.na(dr_r_nh),is.na(mdr))]

B[Hrated=='yes' & !is.na(r_rlt),.N,
  by=.(is.na(dst_rlt_hr),is.na(dr_h_nr),is.na(dr_r_nh),is.na(mdr))]

B[Hrated=='yes' & !is.na(r_rlt) & is.na(dst_rlt_hr),.N,
  by=.(is.na(dr_h_nr),is.na(dr_r_nh),is.na(mdr))]

B[Hrated=='yes' & !is.na(r_rlt) & is.na(dst_rlt_hr) & !is.na(mdr)]


#' r.x?
B[,.N,by=.(is.na(r_rlt),is.na(xpert))] #xpert always missing with rrlt



#'
#' ## Label cases
#' 
#' Set cases:
cases <- B[,{
  if(!is.na(r_rlt)){
    if(Hrated=='no'){
      cz <- 'r'
    } else { #H-rated
      if(!is.na(dst_rlt_hr)){
        cz <- 'r.h'
      } else {
        if(!is.na(mdr)){
          cz <- 'r.m'
        }
      }
    }
  } else {
    if(!is.na(dst_rlt)){
      if(Hrated=='yes'){
        ## dh0
        if(!is.na(dr_h_nr) & !is.na(dr_r_nh) & !is.na(mdr)){
          if(is.na(xpert)){
            cz <- 'dh0'
          } else {
            cz <- 'dh0.x'
          }
        }
        ## dh1
        if(!is.na(dr_h_nr) & is.na(dr_r_nh) & !is.na(mdr)){
          if(is.na(xpert)){
            cz <- 'dh1r'
          } else {
            cz <- 'dh1r.x'
          }
        }
        ## dh2
        if(is.na(dr_h_nr) & is.na(dr_r_nh) & !is.na(mdr)){
          if(is.na(xpert)){
            cz <- 'dh2'
          } else {
            cz <- 'dh2.x'
          }
        }
      } ## else { #not H-rated
      ##   if(is.na(dr_r_nh) & is.na(mdr)){
      ##     if(is.na(xpert)){
      ##       cz <- 'dr'
      ##     } else {
      ##       cz <- 'dr.x'
      ##     }
      ##   }
      ## } # NOTE doesn't happen
    } #no r_rlt, but dst_rlt
  }
  list(case=cz)
},by=id]


cases

#' merge back in
B <- merge(B,cases,by='id')

(CY <- B[,.(countryyears=.N),by=.(patients,case)])
CY[,sum(countryyears)]==nrow(B) #NOTE all there?
nrow(B)

(CC <- B[,.(countries=length(unique(iso3))),by=.(patients,case)])
CC <- merge(CC,CY,by=c('patients','case'))

cord <- c('r','r.h','r.m',
          'dh0','dh0.x','dh1r','dh1r.x','dh2','dh2.x')

CC$case <- factor(CC$case,levels=cord,ordered = TRUE)

setkey(CC,patients,case)
CC

setdiff(cord,CC$case)


fwrite(CC,file=here('drtb/data/L.CC.csv'))


Lcases <- cases
save(Lcases,file=here('drtb/data/Lcases.Rdata'))
LB <- B
save(LB,file=here('drtb/data/LB.Rdata'))


LB[,.(length(unique(iso3)),.N),by=patients]

#' ----- surveys

nmzy <- c('iso3','year','patients','Hrated','id',
         'dst_rlt','r_rlt',
         grep('pcnt',names(H),value=TRUE)
         )

#' Make YB
YB <- rbind(H[source=='Survey',..nmzy],R[source=='Survey',..nmzy])
nrow(YB)


#' CHECK
lonm <- grep('lo',nmzy,value=TRUE)
hinm <- gsub('lo','hi',lonm)
mdnm <- gsub('_lo','',lonm)

YB[!is.na(mdr_pcnt),.N,by=.(is.na(mdr_pcnt_lo),
                            is.na(mdr_pcnt_hi))] #OK
YB[!is.na(xpert_dr_r_pcnt),.N,by=.(is.na(xpert_dr_r_pcnt_lo),
                                   is.na(xpert_dr_r_pcnt_hi))] #OK
YB[!is.na(rr_pcnt),.N,by=.(is.na(rr_pcnt_lo),is.na(rr_pcnt_hi))] #OK
YB[!is.na(dst_rlt_hr_pcnt),.N,by=.(is.na(dst_rlt_hr_pcnt_lo),
                                   is.na(dst_rlt_hr_pcnt_hi))] #OK
YB[!is.na(dst_rlt_rr_pcnt),.N,by=.(is.na(dst_rlt_rr_pcnt_lo),
                                   is.na(dst_rlt_rr_pcnt_hi))] #OK

#' all H-rated
YB[!is.na(dr_h_nr_pcnt) & Hrated=='yes',.N,
   by=.(is.na(dr_h_nr_pcnt_lo),
        is.na(dr_h_nr_pcnt_hi))] #BUG


YB[ !is.na(dr_h_nr_pcnt) &
     Hrated=='yes' &
             is.na(dr_h_nr_pcnt_lo) & is.na(dr_h_nr_pcnt_hi)] #BUG

#' DROP YB (above)
YB[iso3=='TUN' & year==2012,dr_h_nr_pcnt:=NA] #TODO change


#' VIM::aggr(YB[,.(mdr_pcnt, dr_h_nr_pcnt, dr_r_nh_pcnt,
#'                 rr_pcnt, dst_rlt_hr_pcnt, dst_rlt_rr_pcnt)],
#'           combined=TRUE,prop=TRUE,cex.axis=.7)

YB[,.N,by=.(is.na(xpert_dr_r_pcnt),is.na(r_rlt))] #use rr_pcnt for case
YB[,.N,by=.(is.na(rr_pcnt),is.na(r_rlt))] #use rr_pcnt for case
YB[,.N,by=.(is.na(rr_pcnt),is.na(dst_rlt_rr_pcnt))] #use rr_pcnt for case


mdnm
YB[!is.na(rr_pcnt),.N,by=.(is.na(xpert_dr_r_pcnt),
                           is.na(dst_rlt_hr_pcnt),
                           is.na(dst_rlt_rr_pcnt),
                           is.na(dr_h_nr_pcnt),
                           is.na(dr_r_nh_pcnt),
                           is.na(mdr_pcnt))]

#' from r_rlt:
#'   rr (R-res regardless of H)
#' 
#' from dst_rlt
#'   dst_rlt_hr (H-res regardless of R),
#'   dst_rlt_rr (R-res regardless of H),
#'   dr_h_nr (H-res not R-res),
#'   dr_r_nh (R-res not H-res),
#'   mdr

YB[!is.na(rr_pcnt) & Hrated=='yes',.N,
   by=.(is.na(xpert_dr_r_pcnt),
        is.na(dst_rlt_hr_pcnt),
        is.na(dst_rlt_rr_pcnt),
        is.na(dr_h_nr_pcnt),
        is.na(dr_r_nh_pcnt),
        is.na(mdr_pcnt))]

YB[is.na(rr_pcnt),.N,
   by=.(
     is.na(dst_rlt_hr_pcnt),
     is.na(dst_rlt_rr_pcnt))]


YB[is.na(rr_pcnt),.N,
   by=.(
        ## is.na(dst_rlt_hr_pcnt),
        ## is.na(dst_rlt_rr_pcnt),
        is.na(dr_h_nr_pcnt),
        is.na(dr_r_nh_pcnt),
        is.na(mdr_pcnt))]


#' dst_rlt_hr_pcnt
#' dst_rlt_rr_pcnt
#' dr_h_nr_pcnt
#' dr_r_nh_pcnt
#' mdr_pcnt

YB[,.(length(unique(id)),length(id))]

 
Ycases <- YB[,{
  cz <- ''
  if(!is.na(rr_pcnt)){
    if(Hrated=='no'){
      cz <- 'r'
    } else { #H-rated
      if(!is.na(dr_h_nr_pcnt)){
        cz <- 'r.hnr'
      } else if(!is.na(dst_rlt_hr_pcnt)){
        cz <- 'r.hr'
      } else{
        cz <- 'r'
      }
    }
  } else{
    if(Hrated=='yes'){
      if(!is.na(dr_r_nh_pcnt) & !is.na(dr_h_nr_pcnt) & !is.na(mdr_pcnt)){
        cz <- 'dh0'
      }
      if(is.na(dr_r_nh_pcnt) & !is.na(dr_h_nr_pcnt) & !is.na(mdr_pcnt)){
        cz <- 'dh1.r'
      }
      if(!is.na(dr_r_nh_pcnt) & is.na(dr_h_nr_pcnt) & !is.na(mdr_pcnt)){
        cz <- 'dh1.h'
      }
      if(is.na(dr_r_nh_pcnt) & is.na(dr_h_nr_pcnt) & !is.na(mdr_pcnt)){
        cz <- 'dh2'
      }
    }
  }
  ## Xpert cases! 
  if(!is.na(xpert_dr_r_pcnt)) cz <- paste0(cz,".x")
  if(cz=='.x') cz <- 'x'
  list(case=cz)
},by=id]

#' merge back in
YB <- merge(YB,Ycases,by='id')

(YCY <- YB[,.(countryyears=.N),by=.(patients,case)])
YCY[,sum(countryyears)]==nrow(YB) #NOTE all there?
nrow(YB)


(YCC <- YB[,.(countries=length(unique(iso3))),by=.(patients,case)])
YCC <- merge(YCC,YCY,by=c('patients','case'),all=TRUE)
YCC <- YCC[case!='']

cordy <- c('r','r.hnr','r.hr',
           'dh0','dh1.r','dh1.h','dh2')
cordyx <- paste0(cordy,".x")
cordyb <- c(cordy,'x',cordyx)

unique(YCC$case)

YCC$case <- factor(YCC$case,levels=cordyb,ordered = TRUE)

setkey(YCC,patients,case)
YCC

setdiff(cordyb,YCC$case) #NOTE not occuring


fwrite(YCC,file=here('drtb/data/Y.CC.csv'))
save(Ycases,file=here('drtb/data/Ycases.Rdata'))
save(YB,file=here('drtb/data/YB.Rdata'))


YB[,.(length(unique(iso3)),.N),by=patients]



YB[Hrated=='no' & is.na(rr_pcnt)] #all have Xpert data


library(MCMCpack)                       # for inverse wishart etc
library(lattice)                        # MCMC chain plots
library(readxl)                         # for reading excel
library(data.table)                     # data handling
library(ggplot2)                        # plotting
library(scales)                         # plotting
library(productplots)                   # for mosaic plots
library(gridExtra)                      # grouping plots
library(ggpubr)                         # plot arrangents
library(glue)                           # strings
library(knitr)                          # for table formatting
library(stringr)                        # string searching
library(metafor)                        # for MA of sex splits in kids
library(gtbreport)                      # formatting glaziou/gtbreport

estyr <- 2021

Sys.setlocale(locale="C")               #to avoid report UTF8 error
rot45 <- theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
knitr::opts_chunk$set(fig.width=12, fig.height=8) #make the figure size larger

#' A function for making absolute labels with space separators
absspace <- function(x,...) {             #works
   format(abs(x), ..., big.mark=" ",scientific = FALSE, trim = TRUE)
}

## not int
`%ni%` <- Negate(`%in%`)

## function for time stamps in right format:
getm <- function()glue(format(Sys.time(), "%Y-%m-%d-%H%M"))
getm()

## Function for date stamps in right format:
getdate <- function()glue(gsub('-','-',Sys.Date()))
getdate()


whoz <- c('AFR','AMR','EMR','EUR','SEA','WPR')
whozt <- c('Africa','The Americas',
           'Eastern Mediterranean','Europe','South-East Asia',
           'Western Pacific')


#' Some useful age range vectors:
agz <- c('04','514','1524','2534','3544','4554','5564','65') #for extract
agz2 <- c('0_4','5_14','15_24','25_34','35_44','45_54','55_64','65plus') #for labels
agz3 <- gsub('_','-',agz2)
agz4 <- c(rev(rev(agz3)[-1]),"\u2265 65")
kds <- c('0_4','5_14')
kdz <- c('04','514')
kdznew <- c('04','59','1014')
kdzboth <- unique(c(kdz,kdznew))
AA <- data.table(a1=agz,age_group=agz2,age=agz3,agegt=agz4) #for conversion

## MCEE age groups
agzmc <- c('04','59','1014','1519','2024','2534','3544','4554','5564','65') #for MCEE
agzmc2 <- c('0_4','5_9','10_14','15_19','20_24','25_34','35_44','45_54','55_64','65plus') #for labels
agzmc3 <- gsub('_','-',agzmc2)
agzmc4 <- c(rev(rev(agzmc3)[-1]),"\u2265 65")

#' A function get notification data (new & relapse) by age for a given country
#'
#' Needs tb.Rdata to be loaded and est.Rdata and pop.Rdata
#' 
#'
getNotes <- function(cn,m=TRUE,mostrecent=TRUE,silent=FALSE){
  if(!silent) cat('Country=',cn,'\n')
  notes <- matrix(nrow=2,ncol=8)
  yr <- myr <- max(tb$year)
  if(mostrecent){
    while(all(is.na(tb[iso3==cn&year==yr,paste0('newrel.m',agz),with=FALSE]))
          & yr>1980){
      yr <- yr-1
      if(!silent) cat('Check notification data! Using a previous year due to NAs!\n')
    }
  }
  here <- tb[iso3==cn]
  tmp <- unlist(here[year==yr,paste0('newrel.m',agz),with=FALSE])
  notes[1,] <- tmp
  tmp <- unlist(here[year==yr,paste0('newrel.f',agz),with=FALSE])
  notes[2,] <- tmp
  colnames(notes) <- agz2
  rownames(notes) <- c('M','F')
  if(sum(notes[,-c(1:2)],na.rm=TRUE)==0 &
     sum(notes[,c(1:2)],na.rm=TRUE)>0 &
     here[year==yr,!is.na(newrel.m15plus)] &
     here[year==yr,!is.na(newrel.f15plus)]){
    print('3 age groups!')
    notes <- notes[,1:3]                #case where only >15 for 3 age cats
    colnames(notes)[3] <- '15plus'
    notes[1,3] <- here[year==yr,newrel.m15plus] 
    notes[2,3] <- here[year==yr,newrel.f15plus] 
  }
  if(m){                              #melt first
    hagz <- colnames(notes)
    notes <- data.table(newrel=c(notes),
                        sex=rep(c('M','F'),ncol(notes)),
                        age_group=rep(hagz,each=2))
    notes$age_group <- factor(notes$age_group,levels=hagz,ordered=TRUE)
    notes$sex <- factor(notes$sex,levels=c('M','F'),ordered=TRUE)
  }
  notes[,c('iso3','year'):=.(cn,yr)]
  setkeyv(notes,c('age_group','sex'))
  notes
}

## ## test
## print(getNotes('ZAF'))
## print(getNotes('AUT'))
## print(getNotes('AUT',mostrecent=FALSE))
## print(getNotes('BES'))
## tt <- getNotes('NGA'); print(tt)
## tt <- getNotes('MOZ'); print(tt)


niplot <- function(cndt,label.problem=FALSE){
  cn <- cndt[,iso3][1]
  cndt[,age:=gsub("_","-",age_group)];
  hagz3 <- agz3; hagz4 <- agz4;         #local age cats
  if(length(unique(cndt$age))==3){
    hagz3[3] <- '15plus'; hagz4[3] <- "\u2265 15"
  } 
  cndt$age <- factor(cndt$age,levels=hagz3,ordered=TRUE)
  cndt$sex <- factor(cndt$sex,levels=c('M','F'),ordered=FALSE)
  ans <- ggplot(data=cndt,aes(x=age,y=newrel,fill=sex)) +
    coord_flip() + grids() +
    geom_bar(data=cndt[sex=='M'],stat='identity',aes(y=newrel))+
    geom_bar(data=cndt[sex=='F'],stat='identity',aes(y=-newrel))+
    ylab('New & relapse cases') + xlab('Age group')+
    scale_y_continuous(labels = absspace)
  if('inc' %in% names(cndt)){
    ans <- ans +
      geom_bar(data=cndt[sex=='M'],stat='identity',
               aes(x=age,y=inc),fill='transparent',col=1)+
      geom_bar(data=cndt[sex=='F'],stat='identity',
               aes(x=age,y=-inc),fill='transparent',col=1)
  }
  ans <- ans + scale_x_discrete(labels=hagz4)
  nm <- tb[iso3==cn,country[1]]
  if(cn=='CIV') nm <- "Cote d'Ivoire"   #otherwise encoding pbm
  tlt <- nm
  if('problem' %in% names(cndt)){
    if(cndt[,any(problem)] & label.problem)
      tlt <- paste0(tlt,': -- PROBLEM!')
    if('note' %in% names(cndt)){
      if(cndt$note[1]!='')
        tlt <- paste0(tlt,' (',cndt$note[1],')')
    }
  }
  ans <- ans + ggtitle(tlt)
  ans
}

## tmpn <- getNotes('NGA'); print(tmpn)
## tmp <- getNotes('MOZ'); print(tmp)
## niplot(tt)
## niplot(getNotes('ZAF'))


## function for expanding 3 age cats to full
## NOTE may change from integer to real
ni3cat <- function(dat){
  ## c(t(cbind(1:6,7:12))) alternate MF
  nwz <- c(t(cbind(rep(dat[age_group=='15plus' & sex=='M',newrel]/6,6),
                   rep(dat[age_group=='15plus' & sex=='F',newrel]/6,6))
             ))
  nwz <- c(dat[1:4,newrel],nwz)
  datnew <- data.table(newrel=nwz,
                       sex=rep(c('M','F'),length(agz2)),
                       age_group=rep(agz2,each=2),
                       iso3=rep(dat[1,iso3],2*length(agz2)),
                       year=rep(dat[1,year],2*length(agz2))
                       )
  datnew
}

## tmp3 <- ni3cat(tmp)
## tmp


##' A plotting sub-function:
muplot <- function(indat,cnm){
  if('age' %ni% names(indat)) indat[,age:=gsub("_","-",age_group)]
  ## indat[age=='65plus',age:=rev(agz4)[1]]
  indat$age <- factor(indat$age,levels=mcagl,ordered=TRUE)
  indat$sex <- factor(indat$sex,levels=c('F','M'),ordered=TRUE)
  ## plot
  mp <- prodplot(data=indat,mort ~ sex + age,divider=mosaic()) + aes(fill=sex)
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


rowSD <- function(M) sqrt(rowMeans((M-rowMeans(M))^2))
colSD <- function(M) sqrt(colMeans((M-colMeans(M))^2))

getLNparms <- function(m,v){
  lvv <- log(1+v/m^2)
  rbind(mu=log(m)-lvv/2,sig2=lvv)
}

## ## test
## getLNparms(1,.2)
## getLNparms(rep(1,3),rep(.2,3))
## tt <- getLNparms(1,1)
## mean(rlnorm(1e4,tt[1,1],sqrt(tt[2,1])))
## sd(rlnorm(1e4,tt[1,1],sqrt(tt[2,1])))


## Function for rounding anything left with a decimal point and removing '<'
## (for number formatting for tables)
zapdotlt <- function(x){
    x <- gsub("<","",x)
    got <- grep("\\.",x)
    x[got] <- as.character(round(as.numeric(x[got])))
    x
}

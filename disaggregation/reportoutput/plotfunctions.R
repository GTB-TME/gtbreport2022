## functions to make plots of age/sex disaggregation
##
## NOTES:
##
## assumes following libraries are loaded:
## ggplot2, data.table, gtbreport, here (in repo for examples)
##
## assumes the following libraries are installed:
## ggpubr, gridExtra
##
## this is functions only; input data path by each function

## ## testing - uncomment & run para to load data+libraries:
## rm(list=ls())
## library(ggplot2); library(data.table); library(here);
## library(productplots)                   # for mosaic plots
## library(gtbreport)
## ## incidence data
## load(here('disaggregation/reportoutput/h30splt.rda'))
## load(here('disaggregation/reportoutput/globsplt.rda'))
## load(here('disaggregation/reportoutput/regsplt.rda'))
## ## mortality data
## load(here('disaggregation/reportoutput/Mglobsplt.rda'))
## load(here('disaggregation/reportoutput/Mregsplt.rda'))

## utility for formating scales
absspace <- function(x,...) {             #works
    format(abs(x), ..., big.mark=" ",scientific = FALSE, trim = TRUE)
}


## NOTE the following functions have examples commented out above them
## which can be run having loaded the above data & illustrate use

## ===== incidence =========

## disaggregated age/sex incidence at HBC30 level

## ## expects D = h30split
## ## use as:
## ## absolute values, variable scale
## (plt <- disag.plot.inc.h30(h30splt)) 
## ggsave(filename=here('disaggregation/reportoutput/plotexx/h30.absv.pdf'),
##        plot=plt,height=15,width=15*.75,device = cairo_pdf)
## ## absolute values, fixed scale
## (plt <- disag.plot.inc.h30(h30splt,fixed = TRUE))
## ggsave(filename=here('disaggregation/reportoutput/plotexx/h30.absf.pdf'),
##        plot=plt,height=15,width=15*.75,device = cairo_pdf)
## ## per capita
## (plt <- disag.plot.inc.h30(h30splt[,.(iso3,sex,age,name,
##                                       inc=1e5*inc/pop,
##                                       newrel=1e5*newrel/pop)],
##                            fixed = TRUE))
## ggsave(filename=here('disaggregation/reportoutput/plotexx/h30.perc.pdf'),
##        plot=plt,height=15,width=15*.75,device = cairo_pdf)

disag.plot.inc.h30 <- function(D,fixed=FALSE){
    ## ages and colors etc
    agz3 <- c('0-4','5-14','15-24','25-34','35-44','45-54',
              '55-64','65plus') #for labels
    agz4 <- c(rev(rev(agz3)[-1]),"\u2265 65")
    clz <- c('F'=palette_gtb('female'),
             'M'=palette_gtb('male'))
    wz <- rep(1,5);wz[1] <- 1.15; wz <- wz/1.15 #widths
    whoz <- sort(D[,unique(iso3)])
    whosh <- unique(D[,.(iso3,name)])
    max.dat <- 1.05*max(D$inc)
    ## plot construction
    nplst <- list()
    for(i in seq_along(whoz)){
        reg <- D[iso3==whoz[i]]
        plt <-
            ggplot(data=reg,aes(x=age,y=newrel,fill=sex)) +
            coord_flip() +
            geom_bar(stat='identity',aes(y=ifelse(sex=='M',
                                                  newrel,
                                                  -newrel)))+

            scale_fill_manual(values=clz)+
            geom_bar(stat='identity',fill='transparent',col=1,
                     aes(y=ifelse(sex=='M',
                                  inc,
                                  -inc)))+
            ggtitle(whosh[iso3==whoz[i],name]) +
            theme_gtb()+
            theme(legend.position="none") + xlab('') + ylab('') +
            scale_x_discrete(labels=agz4)
        if(fixed){
            plt <- plt + scale_y_continuous(labels = absspace,
                                            limits=c(-max.dat,max.dat))
        } else {
            plt <- plt + scale_y_continuous(labels = absspace)
        }
        if(!(i%%5==1)) plt <- plt + theme(axis.text.y=element_blank())
        nplst[[i]] <- plt
    }
    gridExtra::grid.arrange(nplst[[1]],nplst[[2]],nplst[[3]],
                            nplst[[4]],nplst[[5]],nplst[[6]],
                            nplst[[7]],nplst[[8]],nplst[[9]],
                            nplst[[10]],nplst[[11]],nplst[[12]],
                            nplst[[13]],nplst[[14]],nplst[[15]],
                            nplst[[16]],nplst[[17]],nplst[[18]],
                            nplst[[19]],nplst[[20]],nplst[[21]],
                            nplst[[22]],nplst[[23]],nplst[[24]],
                            nplst[[25]],nplst[[26]],nplst[[27]],
                            nplst[[28]],nplst[[29]],nplst[[30]],
                        ncol=5,widths=wz)
}




## ## disaggregated age/sex incidence at WHO regional level
## ## expects D = regsplt
## ## here('disaggregation/reportoutput/regsplt.rda')
## ## use as:
## ## variable scales, absolute
## (plt <- disag.plot.inc.regional(regsplt))
## ggsave(filename=here('disaggregation/reportoutput/plotexx/reg.absv.pdf'),
##        plot=plt,height=10,width=15,device = cairo_pdf)
## ## same scales, absolute
## (plt <- disag.plot.inc.regional(regsplt,fixed=TRUE)) #same scales
## ggsave(filename=here('disaggregation/reportoutput/plotexx/reg.absf.pdf'),
##        plot=plt,height=10,width=15,device = cairo_pdf)
## ## per capita
## (plt <- disag.plot.inc.regional(regsplt[,.(g.whoregion,sex,age,name,
##                                            inc=1e5*inc/pop,
##                                            newrel=1e5*newrel/pop)],
##                                 fixed=TRUE))
## ggsave(filename=here('disaggregation/reportoutput/plotexx/reg.perc.pdf'),
##        plot=plt,height=10,width=15,device = cairo_pdf)

disag.plot.inc.regional <- function(D,fixed=FALSE ){
    ## ages and colors etc
    agz3 <- c('0-4','5-14','15-24','25-34','35-44','45-54',
              '55-64','65plus') #for labels
    agz4 <- c(rev(rev(agz3)[-1]),"\u2265 65")
    clz <- c('F'=palette_gtb('female'),
             'M'=palette_gtb('male'))
    whoz <- sort(D[,unique(g.whoregion)])
    whosh <- unique(D[,.(g.whoregion,name)])
    max.dat <- 1.05*max(D$inc)
    ## plot construction
    nplst <- list()
    for(i in seq_along(whoz)){
        reg <- D[g.whoregion==whoz[i]]
        plt <-
            ggplot(data=reg,aes(x=age,y=newrel,fill=sex)) +
            coord_flip() +
            geom_bar(stat='identity',aes(y=ifelse(sex=='M',
                                                  newrel,
                                                  -newrel)))+
                     scale_fill_manual(values=clz)+
            geom_bar(stat='identity',fill='transparent',col=1,
                     aes(y=ifelse(sex=='M',
                                  inc,
                                  -inc)))+
            ggtitle(whosh[g.whoregion==whoz[i],name]) +
            theme_gtb()+
            theme(legend.position="none") + xlab('') + ylab('') +
            scale_x_discrete(labels=agz4)
        if(fixed){
            plt <- plt + scale_y_continuous(labels = absspace,
                                            limits=c(-max.dat,max.dat))
        } else {
            plt <- plt + scale_y_continuous(labels = absspace)
        }
        nplst[[i]] <- plt
    }
    ggpubr::ggarrange(plotlist = nplst,ncol=3,nrow=2)
}


## ## disaggregated age/sex incidence at global level
## ## expects D = globsplt
## ## here('disaggregation/reportoutput/globsplt.rda')
## ## use as:
## ## absolute
## (plt <- disag.plot.inc.global(globsplt))
## ggsave(filename=here('disaggregation/reportoutput/plotexx/glob.abs.pdf'),
##        plot=plt,height=7,width=7,device = cairo_pdf)
## ##per capita
## (plt <- disag.plot.inc.global(globsplt[,.(sex,age,
##                                    inc=1e5*inc/pop,
##                                    newrel=1e5*newrel/pop)]))
## ggsave(filename=here('disaggregation/reportoutput/plotexx/glob.perc.pdf'),
##        plot=plt,height=7,width=7,device = cairo_pdf)

disag.plot.inc.global <- function(D){
    ## ages and colors etc
    agz3 <- c('0-4','5-14','15-24','25-34','35-44','45-54',
              '55-64','65plus') #for labels
    agz4 <- c(rev(rev(agz3)[-1]),"\u2265 65")
    clz <- c('F'=palette_gtb('female'),
             'M'=palette_gtb('male'))
    ## plot construction
    plt <- ggplot(data=D,aes(x=age,y=newrel,fill=sex)) +
        coord_flip() +
        geom_bar(stat='identity',aes(y=ifelse(sex=='M',
                                              newrel,
                                              -newrel)))+
        scale_fill_manual(values=clz)+
        geom_bar(stat='identity',fill='transparent',col=1,
                 aes(y=ifelse(sex=='M',
                              inc,
                              -inc)))+
        theme_gtb() +
        theme(legend.position="none") + xlab('') + ylab('') +
        scale_x_discrete(labels=agz4)+
        scale_y_continuous(labels = absspace)
    plt
}



## ===== mortality =========

## utility plot function used in other plotters
muplot <- function(indat,cnm){
    agz2 <- c('0_4','5_14','15_24','25_34','35_44',
              '45_54','55_64','65plus') #for labels
    agz3 <- gsub('_','-',agz2)
    agz4 <- c(rev(rev(agz3)[-1]),"\u2265 65")
    if(!'age' %in% names(indat)) indat[,age:=gsub("_","-",age_group)]
    indat[age=='65plus',age:=rev(agz4)[1]]
    indat$age <- factor(indat$age,levels=agz4,ordered=TRUE)
    indat$sex <- factor(indat$sex,levels=c('F','M'),ordered=TRUE)
    ## plot
    mp <- prodplot(data=indat,
                   mort ~ sex + age,
                   divider=mosaic()) +
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


## ## disaggregated age/sex HIV-ve mortality at global level
## ## expects D = Mglobsplt
## ## here('disaggregation/reportoutput/Mglobsplt.rda')
## ## use as:
## ## absolute
## (plt <- disag.plot.mort.global(Mglobsplt,'Global'))
## ggsave(filename=here('disaggregation/reportoutput/plotexx/glob.mort.pdf'),
##        plot=plt,height=7,width=7,device = cairo_pdf)

disag.plot.mort.global <- function(D,title){
    ## ages and colors etc
    agz3 <- c('0-4','5-14','15-24','25-34','35-44','45-54',
              '55-64','65plus') #for labels
    agz4 <- c(rev(rev(agz3)[-1]),"\u2265 65")
    clz <- c('F'=palette_gtb('female'),
             'M'=palette_gtb('male'))
    ## plot construction
    plt <- muplot(D,title) +
        coord_flip() +
        scale_fill_manual(values=clz)+
        theme_gtb()+
        theme(axis.title.x=element_blank(), #remove
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              legend.position = 'none')
    plt
}



## ## disaggregated age/sex HIV-ve mortality at regional level
## ## expects D = Mglobsplt
## ## here('disaggregation/reportoutput/Mregsplt.rda')
## ## use as:
## (plt <- disag.plot.mort.regional(Mregsplt))
## ggsave(filename=here('disaggregation/reportoutput/plotexx/reg.mort.pdf'),
##        plot=plt,height=7.5,width=10,device=cairo_pdf)

disag.plot.mort.regional <- function(D){
    ## ages and colors etc
    agz3 <- c('0-4','5-14','15-24','25-34','35-44','45-54',
              '55-64','65plus') #for labels
    agz4 <- c(rev(rev(agz3)[-1]),"\u2265 65")
    whozt <- c('Africa','The Americas','Eastern Mediterranean','Europe',
               'South-East Asia',
               'Western Pacific')
    regs <- c('AFR','AMR','EMR','EUR','SEA','WPR')
    clz <- c('F'=palette_gtb('female'),
             'M'=palette_gtb('male'))
    ## plot construction
    Pltlst <- list()
    for(i in 1:6){
        plt <- muplot(D[g.whoregion==regs[i]],whozt[i]) +
            theme(legend.position="none")
        if(!(i%%3==1)) plt <- plt + theme(axis.text.y=element_blank())
        plt <- plt + coord_flip() +
            scale_fill_manual(values=clz)+
            theme_gtb()+
            theme(axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank(),
                  legend.position = 'none')
        Pltlst[[i]] <- plt
    }
    ggpubr::ggarrange(plotlist = Pltlst, ncol = 3,nrow=2)
}

## this uses the stuff in Xdatanfunction.R
rm(list=ls())
library(here)
source(here('disaggregation/R/02children/Xpaeddatanfunctions.R'))

## running...
LHsample <- randomLHS(1e4,LHN)                #latin hypercube sample
LHparms <- makeParameters2(LHsample,metaparm) #turn the raw LHS into parameters

## qbeta warnings...
system.time({                           #around 100 sec for 1e4
    LAcl <- loopAns2(LHparms,metaparm,method=c('comm','lat'))
})

Dcl <- make2cat(LAcl)                      #aggregate output
testof <- as.data.table(Dcl)

if(!file.exists(here('disaggregation/R/02children/plots'))) dir.create(here('disaggregation/R/02children/plots'))
## save(testof,file=here('disaggregation/R/02children/plots/testof.Rdata'))
##27M for 1e4, which seems a bit pointless

## Make BB age split data
## regional average where missing
fn <- here('disaggregation/R/02children/data/BB.Rdata')
BB <- testof[,prop:=value/sum(value),by=.(iso3,replicate)]
BB[replicate==1 & iso3=='AFG']
BB <- BB[age=='0-4',.(iso3,replicate,prop)]
BB <- BB[,.(propu5=mean(prop),propu5.sd=sd(prop)),by=iso3]
BB <- merge(BB,unique(est[,.(iso3,g.whoregion)]),by='iso3',
            all.x=TRUE,all.y=FALSE)
BBR <- BB[,.(propu5=mean(propu5),propu5.sd=mean(propu5.sd)),
          by=g.whoregion]
BBR <- merge(unique(est[,.(iso3,g.whoregion)]),BBR,
             by='g.whoregion',
             all.x=TRUE,all.y = TRUE,allow.cartesian = TRUE)
setkey(BBR,iso3)
(missing <- setdiff(BBR$iso3,BB$iso3))
BB <- rbind(BB,BBR[missing])  #add missing with regional averages
BB[,g.whoregion:=NULL]
## add beta parms
BB[,tot:=propu5*(1-propu5)/propu5.sd^2-1]
BB[,summary(tot)]
BB[,c('aa','ab'):=.(tot*propu5,tot*(1-propu5))]
BB[,tot:=NULL]
## save
setkey(BB,iso3)
save(BB,file=fn)


## reporting function
## NB have got rid of ARI/LAT
reporting <- function(testdf){          #takes data.table
    ## least disaggregated
    report4 <- testdf[,list(median=round(quantile(value,0.5)),
                            LQ=round(quantile(value,0.25)),
                            UQ=round(quantile(value,0.75))),
                      by=list(country,iso3,age)]
    report4 <- report4[order(report4$iso3,report4$age),] #order by iso
    cat('1\n')
    fn <- here('disaggregation/R/02children/plots/Report_leastdisag.csv')
    write.table(file=fn,
                report4,row.names = FALSE,col.names = TRUE,sep=',')
    ##join yng and old
    tmp <- testdf[,.(value=sum(value)),
                  by=.(iso3,country,replicate)] #testob
    report <- tmp[,.(median=round(quantile(value,0.5)),
                     LQ=round(quantile(value,0.25)),
                      UQ=round(quantile(value,0.75))),
                by=.(country,iso3)]
    report <- report[order(report$iso3),]          #order by country
    cat('2\n')
    fn <- here('disaggregation/R/02children/plots/Report_countrymodel.csv')
    write.table(file=fn,
                report,row.names = FALSE,col.names = TRUE,sep=',')
    ##join  country
    tmp3 <- testdf[,.(value=sum(value)),by=.(age,replicate)] #testog
    report3 <- tmp3[,.(median=round(quantile(value,0.5)),
                        LQ=round(quantile(value,0.25)),
                        UQ=round(quantile(value,0.75))),by=.(age)]
    report3 <- report3[order(report3$age),]          #order by age
    cat('3\n')
    fn <- here('disaggregation/R/02children/plots/Report_agemodel.csv')
    write.table(file=fn,
                report3,row.names = FALSE,col.names = TRUE,sep=',')
    ##join yng and old & country
    tmp2 <- testdf[,.(value=sum(value)),by=.(replicate)]
    report2 <- tmp2[,.(median=round(quantile(value,0.5)),
                          LQ=round(quantile(value,0.25)),
                          UQ=round(quantile(value,0.75)))]
    cat('4\n')
    fn <- here('disaggregation/R/02children/plots/Report_model.csv')
    write.table(file=fn,
              report2,row.names = FALSE,col.names = TRUE,sep=',')
    ## prettier reporting
    prettyreport <- data.frame(paste(report2$median,
                                     ' [',report2$LQ,'-',report2$UQ,']',
                                     sep=''))
    names(prettyreport)[1] <- 'median [IQR]'
    cat('5\n')
    fn <- here('disaggregation/R/02children/plots/Report_pretty.csv')
    write.table(file=fn,
                prettyreport,row.names = FALSE,col.names = TRUE,sep=',')
    print(report2)
    return(list(report,report2,report3,report4))
}

LR <- reporting(testof)
save(file=here('disaggregation/R/02children/plots/AllReport.Rdata'),LR)

## make df summing yng and old
testob <- testof[,.(value = sum(value)),
                 by=.(country, replicate, iso3)]
## global by method, young and old
testog <- testof[,.(value = sum(value)),by=.(age,replicate)]
## global by method,
testogb <- testob[,.(value = sum(value)),by=.(replicate)]


scientific_10 <- function(x) {
    parse(text=gsub("e", " %*% 10^", scientific_format()(x)))
}

rme <- function(x,digits=0) round(median(x),digits = digits)
ruq <- function(x,digits=0) round(quantile(x,probs = 0.75),
                                  digits=digits)
rlq <- function(x,digits=0) round(quantile(x,probs = 0.25),
                                  digits=digits)

## --------
pp <- ggplot(data=testogb,aes(x=value)) +
    geom_density(trim=FALSE) +
    xlab('Total paediatric TB incidence (per year)')+
    scale_x_continuous(label=comma,limits=c(0,2.5e6))+
    scale_y_continuous(label=scientific_10)
ggsave(pp,file=here('disaggregation/R/02children/plots/Global_NL.pdf'))

p0 <- ggplot(data=testogb,aes(x=value,y = ..density..)) +
    geom_freqpoly(binwidth=1e5) +
    xlab('Total paediatric TB incidence (per year)') +
    scale_x_continuous(label=comma,limits=c(0,2e6))+
    scale_y_continuous(label=scientific_10)
ggsave(p0,file=here('disaggregation/R/02children/plots/GlobalFP.pdf'))

## country plots
testoc <- testob[,.(median=rme(value),LQ=rlq(value),UQ=ruq(value)),
                 by=.(iso3,country)]

## ------- write out ---------
write.table(testoc,
            file=here('disaggregation/R/02children/plots/paedout.csv'),
            sep=',',col.names=TRUE,row.names=FALSE)

## also write out as K.Rdata
K <- copy(testoc)
K[,median.sd:=(UQ-LQ)/1.36]
K <- merge(K,est[year==estyr,.(iso3,inc,inc.sd)],by='iso3')
K <- merge(K,pop[year==estyr,.(iso3,pop)],by='iso3')
K[,inc.num:=inc * pop / 1e5]
K[,inc.num.sd:=inc.sd * pop / 1e5]
K[,propu15b:=median/inc.num]
K[,propu15b.sd:=propu15b*sqrt(median.sd^2/median^2 +
                              inc.num.sd/inc.num^2)]
K <- K[,.(iso3,propu15b,propu15b.sd)]
save(K,file=here('disaggregation/R/02children/data/K.Rdata'))#child prior


## check against last year
DO <- fread(here('disaggregation/R/02children/data/paedout_2019.csv'))
names(DO)[3:5] <- c('M14','L14','H14')
DO <- merge(DO,testoc,by=c('iso3','country'))
summary(DO$median/DO$M14)
odd <- DO$median/DO$M14
odd <- odd>1.5 | odd<0.5

p1 <- ggplot(data=DO[odd,],aes(x=country,y=median/M14)) +
    geom_point() + coord_flip() + scale_y_sqrt() +
    geom_hline(yintercept=1,col=2) + ylab('Ratio 2020/2019') +
    ggtitle('Countries with > 50% change 2019->2020')
ggsave(p1,file=here('disaggregation/R/02children/plots/Ratio.pdf'),
       width=10)



p <- ggplot(data=DO,aes(x=M14,y=median)) +
    geom_abline(intercept=0,slope=1,linetype=2) +
    geom_errorbar(aes(ymin=LQ,ymax=UQ),width=0,color=2) +
    geom_errorbarh(aes(xmin=L14,xmax=H14),color=2) +
    geom_point() +
    xlab('prior 2020 incidence') +
    ylab('prior 2021 incidence')
p2 <- p + scale_x_continuous(label=comma) +
    scale_y_continuous(label=comma) +
    coord_fixed(ratio=1,xlim=c(0,250e3),ylim=c(0,250e3))
## p2
ggsave(p2,file=here('disaggregation/R/02children/plots/2020CountryCompare.pdf'))

p3 <- p2 + geom_text_repel(aes(x=M14,y=median,label=iso3),
                           hjust=0,vjust=0) +
    scale_x_sqrt(label=comma) + scale_y_sqrt(label=comma) +
    coord_fixed(ratio=1,xlim=c(0,250e3),ylim=c(0,250e3))
ggsave(p3,file=here('disaggregation/R/02children/plots/2020CountryCompareLabeled.pdf'))

p3 <- p + scale_x_log10()+ scale_y_log10()
## p3
ggsave(p3,file=here('disaggregation/R/02children/plots/2020CountryCompareLog.pdf'))

p4 <- p + scale_x_sqrt(label=comma) + scale_y_sqrt(label=comma) +
    coord_fixed(ratio=1,xlim=c(0,250e3),ylim=c(0,250e3))
p4 <- p4  + geom_text(aes(x=M14,y=median,label=iso3),hjust=0,vjust=0)
## p4
ggsave(p4,
       file=here('disaggregation/R/02children/plots/2020CountryCompareSqrt.pdf'),
       height=10,width=10)


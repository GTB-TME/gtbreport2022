## post-processing/plotting utilities


## for parsing outputs
getDRprops <- function(X){
  S <- as.data.table(X)
  knmz <- grep('props',names(S),value=TRUE)
  S <- S[,..knmz]
  S[,id:=1:nrow(S)]
  S <- melt(S,id.vars = 'id')
  S[,patients:=ifelse(grepl('nprops',variable),'new','ret')]
  S[,c('dy','cid','dst'):=tstrsplit(variable,split=',')]
  S[,dy:=gsub('nprops','',dy)]
  S[,dy:=gsub('rprops','',dy)]
  S[,dy:=gsub('\\[','',dy)]
  S[,dy:=as.integer(dy)-1]
  S[,dst:=gsub('\\]','',dst)]
  S[,dst:=ifelse(dst=='1','HRS',
          ifelse(dst=='2','HNR',
          ifelse(dst=='3','RNH',
                 'MDR')))]
  S[,cid:=as.integer(cid)]
  S <- merge(S,ckey,by='cid')
  S[,variable:=NULL]
  S
}

getRRprops <- function(X,J=2){
  if(J==2){
    S <- as.data.table(X)
    knmz <- grep('props',names(S),value=TRUE)
    S <- S[,..knmz]
    S[,id:=1:nrow(S)]
    S <- melt(S,id.vars = 'id')
    S[,patients:=ifelse(grepl('nprops',variable),'new','ret')]
    S[,c('dy','cid','dst'):=tstrsplit(variable,split=',')]
    S[,dy:=gsub('nprops','',dy)]
    S[,dy:=gsub('rprops','',dy)]
    S[,dy:=gsub('\\[','',dy)]
    S[,dy:=as.integer(dy)-1]
    S[,dst:=gsub('\\]','',dst)]
    S[,dst:=ifelse(dst=='1','RS','RR')]
    S[,cid:=as.integer(cid)]
    S <- merge(S,ckey,by='cid')
    S[,variable:=NULL]
  } else if(J==4){
    K <- getDRprops(X)
    K[,dst:=ifelse(dst %in% c('RNH','MDR'),'RR','notRR')]
    S <- K[,.(value=sum(value)),by=.(iso3,dy,id,patients,dst)]
  }
  S
}

makeplotdata <- function(K,maxyr=2021,ui=0.95){
  K[,year:=maxyr-dy]
  KO <- K[dst=='RR',
          .(RR.mid=mean(value)*1e2,
            RR.lo=quantile(value,(1-ui)/2)*1e2,
            RR.hi=quantile(value,1-(1-ui)/2)*1e2
            ),
          by=.(iso3,year,patients,type=dst)]
  merge(KO,regkey[,.(iso3,g_whoregion)],by='iso3',
        all.x=TRUE,all.y=FALSE)
}

makeregionplots <- function(KO,fn){
  for(reg in unique(RPD[,g_whoregion])){
    print(reg)
    lfn <- fn + '_' + reg + '.pdf'
    GP <- ggplot(RPD[g_whoregion==reg],
                 aes(year,RR.mid,
                     ymin=RR.lo,ymax=RR.hi,
                     shape=type,col=patients,lty=type))+
      geom_point()+geom_pointrange()+
      geom_line(data=KO[g_whoregion==reg])+
      geom_ribbon(data=KO[g_whoregion==reg],aes(fill=patients),
                  alpha=0.2,col=NA)+
      scale_y_sqrt()+
      facet_wrap(~iso3,scales='free')+
      rot45 + ylab('Rifampicin resistance (%, square root scale)')
    ggsave(GP,file=lfn,w=12,h=15)
  }
}


##generaters a nearest neighbour structure
## the nearest neighbour structure is included in repo/indata so rest can be run without running this
library(here)
## country list (full)
fn <- here('drtb/data/isoidxL.Rdata') 
if(file.exists(fn)){
  load(fn)
} else {
  load(here('data/tb.rda'))
  isoidxL <- unique(tb$iso3) #217
  save(isoidx,file=here('drtb/data/isoidxL.Rdata'))
}


## === generate spatial structure  ===
library(rworldmap)
library(spdep)
sf::sf_use_s2(FALSE) #switch off spherical geom

## Checking the countries and constructing adjacencies
## Load spatial data and compare with covariates
data(countriesCoarse)                   #load rworldmap data
names(countriesCoarse@data)[30] <- 'iso3' #rename


## Remove Antarctica, Greenland, US Virgin Islands (ATA, GRL, VIR)
WTB <- countriesCoarse
WTB <- WTB[which(!(WTB@data$iso3 %in% c('ATA','GRL','VIR'))),]
WTB <- WTB[which(WTB@data$iso3 %in% isoidx),] #drops BES (Bonnaire) and TKL (Tokelau)
dim(WTB)
tmp <- unique(WTB@data$iso3)
length(tmp) #213
setdiff(isoidx,tmp) #"ANT" "GRL" "SCG" "TKL" (SCG=Serbia & Montenegro)

## use this shorter version
fn <- here('drtb/data/isoidx.Rdata') 
if(file.exists(fn)){
  load(fn)
} else {
  isoidx <- as.character(tmp)
  save(isoidx,file=here('drtb/data/isoidx.Rdata'))
}

## Start making adjacencies
W.nb <- poly2nb(WTB,row.names=as.character(WTB@data$iso3)) #error from version change fixed
W.list <- nb2listw( W.nb,style='B',zero.policy = TRUE)
is.symmetric.nb(W.nb)
## W <- nb2mat(W.nb)
W <- listw2mat(W.list)
str(W)

## look at
plot(WTB)
plot(W.nb,coordinates(WTB),add=TRUE,col=2)


## Generate naive nearest neighbout maps
cs_kn5 <- knn2nb(knearneigh(coordinates(WTB),k=5))


## look at the $k=5$ version
plot(WTB,col='gray');
text(coordinates(WTB),labels=WTB$iso3);
plot(cs_kn5,coordinates(WTB),add=TRUE,col=2);


## Modifying the nearest neighbour linkages
## 
## The naive maps have some features that are undesirable. In particular;
## 
## * they are not connected, which can hamper mixing
## * Russia is not connected to former Soviet republics
##
## To address these issues, we introduce transatlantic links based on flight paths carrying more than 900K per year. We also join up the South Pacific islands WLF and TON to FJI. Finally, we join RUS to EST, LVA, UKR, GEO and AZE, which makes sense both in geographical and political terms.
##
## This function is to carry out these changes and to save the data and plots out.

appendKNN <- function(X,savename){
  ## connect up RUS to FSRs & consider making connected
  toadd <- c('EST','LVA','UKR','GEO','AZE')
  toadd <- which(WTB$iso3 %in% toadd)
  indiso <- which(WTB$iso3=='RUS')
  X[[indiso]] <- c(X[[indiso]],toadd)
  ## transatlantic links see: https://en.wikipedia.org/wiki/Transatlantic_flight#Busiest_transatlantic_routes
  ## those with >900K 
  ## transatlantic links UK, FRA to US 
  toadd <- c('GBR','FRA')
  toadd <- which(WTB$iso3 %in% toadd)
  indiso <- which(WTB$iso3=='USA')
  X[[indiso]] <- c(X[[indiso]],toadd)
  ## UK to CAN
  toadd <- which(WTB$iso3 == 'GBR')
  indiso <- which(WTB$iso3=='CAN')
  X[[indiso]] <- c(X[[indiso]],toadd)
  ## US to ETH
  toadd <- which(WTB$iso3 == 'ETH')
  indiso <- which(WTB$iso3=='USA')
  X[[indiso]] <- c(X[[indiso]],toadd)
  ## BRA to AGO
  toadd <- which(WTB$iso3 == 'BRA')
  indiso <- which(WTB$iso3=='AGO')
  X[[indiso]] <- c(X[[indiso]],toadd)
  ## plot
  pdf(paste0(savename,'.pdf'),w=18,h=12)
  plot(WTB,col='gray')
  plot(X, coordinates(WTB),add=TRUE,col=2)
  text(coordinates(WTB),labels=WTB$iso3)
  text(x=0,y=-60,labels='Not shown: WLF & TON connected to FJI')
  dev.off()
  ## transatlantic links WLF, TON to FJI (don't map)
  toadd <- c('WLF','TON')
  toadd <- which(WTB$iso3 %in% toadd)
  indiso <- which(WTB$iso3=='FJI')
  X[[indiso]] <- c(X[[indiso]],toadd)
  knn <- make.sym.nb(X)
  save(knn,file = paste0(savename,'.Rdata'))
  knn
}

## apply this function:
knnb <- appendKNN(cs_kn5,here('drtb/data/knnb_5'))
W <- nb2mat(knnb,style='B',zero.policy = TRUE)
all(W==t(W))
dim(W)
save(W,file=here('drtb/data/W5.Rdata'))

isoidx <- WTB$iso3
save(isoidx,file=here('drtb/data/isoidx5.Rdata'))


## correction to AFG linkages
load(file=here('drtb/data/W5.Rdata'))
load(file=here('drtb/data/isoidx5.Rdata'))

who <- W[which(isoidx=='AFG'),]
who <- which(who>0)
isoidx[who]

## instead use:
afgn <- c('PAK', 'IND', 'TUR', 'IRN','CHN')
who <- which(isoidx %in% afgn)

W[whoa,] <- W[,whoa] <- 0
W[whoa,who] <- W[who,whoa] <- 1

## check
who <- W[which(isoidx=='AFG'),]
who <- which(who>0)
isoidx[who] #OK

all(W==t(W))
dim(W)
save(W,file=here('drtb/data/W5.Rdata'))

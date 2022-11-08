## calculating for DR report
library(here)
library(glue)
library(data.table)

## utilities
gh <- function(x) glue(here(x))
source(here('drtb/R/utils/ftb.R'))
ssum <- function(x,...) sqrt(sum(x^2,...)) #sqrt sum sqrs
brkt <- function(x,y) paste0(ftb(x),' (',
                             ftb(pmax(0,x-1.96*y)),' - ',
                             ftb(x+1.96*y),')')


## MDR mortality: re-written not to use propagate package
##
M2m <- function(M,
                M.sd,
                p,
                p.sd,
                r = 2.26,
                r.sd = 0.45) {
  ##' m = Mpr
  ##' @param m = MDR-TB mortality, the unknown factor
  ##' @param M = total (HIV- and HIV+) TB mortality
  ##' @param p = overall proportion RR among prevalent TB, combined (among new and retreated) estimate
  ##' @param r = RR or risk of dying from MDR compared to non-MDR TB (proxy for RR deaths),
  ##'     approximated with the OR from meta-analysis 2.26 (1.54-3.33), SE=0.45
  Q <- M * p * r / (1 - p + p * r)
  dQ <- c(p * r / (1 - p + p * r),          #dQ/dM
          M * r / (1 - p + p * r)^2,        #dQ/dp
          M * p * (1-p) / (1 - p + p * r)^2 #dQ/dr
  )
  ## check
  ## v <- c('M','p','r')
  ## calculus::derivative(expression(M * p * r / (1 - p + p * r)),v)
  ## calculus::hessian(expression(M * p * r / (1 - p + p * r)),v)[3,]
  H <- c(
    0,r/(1 - p + p * r) - p * r * (r - 1)/(1 - p + p * r)^2,p/(1 - p + p * r) - p * r * p/(1 - p + p * r)^2,
    r/(1 - p + p * r) - p * r * (r - 1)/(1 - p + p * r)^2,-(M * r * (r - 1)/(1 - p + p * r)^2 + (M * r * (r - 1)/(1 - p + p * r)^2 - M * p * r * (r - 1) * (2 * ((r - 1) * (1 - p + p * r)))/((1 - p + p * r)^2)^2)),M/(1 - p + p * r) - M * r * p/(1 - p + p * r)^2 - ((M * p * (r - 1) + M * p * r)/(1 - p + p * r)^2 - M * p * r * (r - 1) * (2 * (p * (1 - p + p * r)))/((1 - p + p * r)^2)^2),
    p/(1 - p + p * r) - p * r * p/(1 - p + p * r)^2,M/(1 - p + p * r) - M * p * (r - 1)/(1 - p + p * r)^2 - ((M * r * p + M * p * r)/(1 - p + p * r)^2 - M * p * r * p * (2 * ((r - 1) * (1 - p + p * r)))/((1 - p + p * r)^2)^2),-(M * p * p/(1 - p + p * r)^2 + (M * p * p/(1 - p + p * r)^2 - M * p * r * p * (2 * (p * (1 - p + p * r)))/((1 - p + p * r)^2)^2))
  )
  H <- matrix(H,nrow=3,ncol=3)
  H2 <- H^2
  v.sd <- c(M.sd,p.sd,r.sd)
  V <- v.sd^2
  Q.V <- sum(dQ^2*v.sd^2) + 0.5 * t(V) %*% H2 %*% V + 0.25*sum((diag(H) * V)^2)
  Q.sd <- c(sqrt(Q.V))
  out <- c(Q,Q.sd,Q-1.96*Q.sd,Q+1.96*Q.sd)
  names(out) <- c("m", "m.sd", "m.lo", "m.hi")
  return(out)
}



## ===  apply

## load data
load(here('data/global.rda'))
load(here('data/dr.est.rda'))
load(here('drtb/dboutput/db_dr_group.rda'))


## RR proportions
rrglob <- db_dr_group[group_name=='global' & year==2021]

rrglob[,inc.rr:=e_inc_rr_num]
rrglob[,inc.rr.sd:=(e_inc_rr_num_hi-e_inc_rr_num_lo)/3.96]
rrglob[,prop.rr:=inc.rr/inc]
rrglob[,prop.rr.sd:=prop.rr * sqrt((inc.rr.sd/inc.rr)^2+(inc.sd/inc)^2)]
rrglob[,.(prop.rr,prop.rr.sd)]



global.mort.rrF <- M2m(last(global$mort),    #total mortality
                      last(global$mort.sd), #SD in total mortality
                      p = 0.05804381,       #same as 2017 - 2021
                      p.sd = 0.004730819) * last(global$pop) / 1e5

brkt(global.mort.rrF['m'],global.mort.rrF['m.sd'])
## NOTE this is the version to use
## 191 000 (119 000 - 264 000)

## X-chec:
## 2020 report
## 182 000 (range, 113 000–250 000)


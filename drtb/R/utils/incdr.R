## function to calculate DR incidence

## # RR incidence
## $I_\text{rr} = I \left[ (1-f) p_n ((1-r) + r \rho) + fp_r \right]$
## $I$ = TB incidence
## $f$ = proportion retreatment for failure or default
## $p_n$ = prop.rr.new
## $r$ = proportion relapses
## $\rho$ = risk ratio of rr in relapse vs new
## $p_r$ = prop.rr.ret

## version by hand so as not to rely on propagate, which I have trouble compiling
incdr <- function(inc,
                  inc.sd,
                  f,
                  f.sd,
                  pn,
                  pn.sd,
                  r,
                  r.sd,
                  rho,
                  rho.sd,
                  pr,
                  pr.sd,
                  secondorder=TRUE) {

  ## Q <- (inc * ((1 - f) * pn * ((1 - r) + r * rho) + f * pr))
  B <- ((1 - r) + r * rho)
  A <- ((1 - f) * pn * B + f * pr)
  Q <- inc * A
  ## deriv wrt I,f,pn,pr,r,rho
  dQ <- c(A,                          #dQ/dI
          inc * (pr - pn * B),        #dQ/df
          inc * (1-f) * B,            #dQ/dpn
          inc * f,                    #dQ/dpr
          inc * (1-f) * pn * (rho-1), #dQ/dr
          inc * (1-f) * pn * r        #dQ/drho
          )
  ## v <- c('inc','f','pn','pr','r','rho')
  ## NOTE CHECK:
  ## calculus::derivative(expression((inc * ((1 - f) * pn * ((1 - r) + r * rho) + f * pr))),v)
  ## calculus::hessian(expression((inc * ((1 - f) * pn * ((1 - r) + r * rho) + f * pr))),v)[6,]
  H <- c(
 0                               ,pr - pn * ((1 - r) + r * rho),(1 - f) * ((1 - r) + r * rho) ,f   ,(1 - f) * pn * (rho - 1),   (1 - f) * pn * r    ,
 (pr - pn * ((1 - r) + r * rho)), 0                            ,-(inc * ((1 - r) + r * rho))  ,inc ,-(inc * (pn * (rho - 1))), -(inc * (pn * r))   ,
 ((1 - f) * ((1 - r) + r * rho)), -(inc * ((1 - r) + r * rho)) ,0                             ,0   ,inc * ((1 - f) * (rho - 1)),inc * ((1 - f) * r) ,
 f                              , inc                          ,0                             ,0   ,0                          , 0                   ,
 ((1 - f) * pn * (rho - 1))     , -(inc * (pn * (rho - 1)))    ,inc * ((1 - f) * (rho - 1))   ,0   ,0                          ,inc * ((1 - f) * pn),
 ((1 - f) * pn * r)              ,-(inc * (pn * r))            ,inc * ((1 - f) * r)           ,0   ,inc * ((1 - f) * pn)       , 0
         )
  H <- matrix(H,nrow=6,ncol=6)
  H2 <- H^2
  v.sd <- c(inc.sd,f.sd,pn.sd,pr.sd,r.sd,rho.sd)
  V <- v.sd^2
  Q.V <- sum(dQ^2*v.sd^2)
  if(secondorder)
    Q.V <- Q.V + 0.5 * t(V) %*% H2 %*% V + 0.25*sum((diag(H) * V)^2)
  Q.sd <- c(sqrt(Q.V))
  list(Q,Q.sd)
}



## ## TEST with rr.rda
## rr[,c('test.rrinc','test.rrinc.sd'):=
##       incdr(inc,
##             inc.sd,
##             f,
##             f.sd,
##             pn = prop.rr.new,
##             pn.sd = prop.rr.new.sd,
##             r,
##             r.sd,
##             rho = rho.rr,
##             rho.sd = rho.sd.rr,
##             pr = prop.rr.ret,
##             pr.sd = prop.rr.ret.sd,
##             secondorder=TRUE),
##    by=iso3]

## rr[,plot(test.rrinc,inc.rr)];abline(a=0,b=1,col=2)
## rr[,plot(test.rrinc.sd,inc.rr.sd)];abline(a=0,b=1,col=2)

## rr[abs(test.rrinc.sd-inc.rr.sd)>1e-2 & inc.rr>1,.(iso3,inc.rr)]
## rr[abs(test.rrinc.sd-inc.rr.sd)>1e-2 & inc.rr>1,.(iso3,inc.rr)]
## ## iso3    inc.rr
## ## 1:  ASM  2.069613
## ## 2:  CAF  3.745904
## ## 3:  TUV 17.283078


## function for adding incidence to proportions:
##    needs rhofnr.Rdata loaded
##    needs est.rda loaded
addIncidence <- function(D,secondorder=TRUE){
  D[,c('rr','rr.sd'):=.(RR.mid/1e2,(RR.hi-RR.lo)/392)]
  D <- dcast(D,iso3 + year + g_whoregion ~ patients,value.var = c('rr','rr.sd'))
  D <- merge(D,rhofnr,by=c('iso3')) #,'year'
  D <- merge(D,est[,.(iso3,year,inc.num,inc.num.sd=(inc.hi.num-inc.lo.num)/3.92)],
             by=c('iso3','year'),all.x = TRUE,all.y=FALSE)
  ## all RR incidence
  D[,c('inc.rr','inc.rr.sd'):=incdr(inc.num,
                                    inc.num.sd,
                                    f,
                                    f.sd,
                                    pn = rr_new,
                                    pn.sd = rr.sd_new,
                                    r,
                                    r.sd,
                                    rho = rho,
                                    rho.sd = rho.sd,
                                    pr = rr_ret,
                                    pr.sd = rr.sd_ret,
                                    secondorder),
    by=.(iso3,year)]
  ## (inc * ((1 - f) * pn * ((1 - r) + r * rho) + f * pr))
  ## pn = rr_new,
  ## pn.sd = rr.sd_new,
  ## pr = rr_ret,
  ## pr.sd = rr.sd_ret,
  ## RR incidence in new
  D[,inc.new.rr:=inc.num * (1 - f) * rr_new * ((1 - r) + r * rho)]
  xx <- D[,{
    dq <- c((1 - f) * rr_new * ((1 - r) + r * rho), #dI
            - inc.num * rr_new * ((1 - r) + r * rho),   #df
            inc.num * (1 - f) *  ((1 - r) + r * rho),   #drrn
            inc.num * (1 - f) * rr_new * (rho-1),       #dr
            inc.num * (1 - f) * rr_new * (r)            #drho
            );
    v.sd <- c(inc.num.sd,f.sd,rr.sd_new,r.sd,rho.sd);
    q.v <- sum(dq^2*v.sd^2);
    list(inc.new.rr.sd=sqrt(q.v))
  },by=.(iso3,year)]
  D <- merge(D,xx,by=c('iso3','year'))
  ## D[,c('inc.new.rr','inc.new.rr.sd'):=incdr(inc.num,
  ##                                           inc.num.sd,
  ##                                           f,
  ##                                           f.sd,
  ##                                           pn = rr_new,
  ##                                           pn.sd = rr.sd_new,
  ##                                           r,
  ##                                           r.sd,
  ##                                           rho = rho,
  ##                                           rho.sd = rho.sd,
  ##                                           pr = 0,
  ##                                           pr.sd = 0,
  ##                                           secondorder),
  ##   by=.(iso3,year)]
  ## RR incidence in ret
  D[,inc.ret.rr:=inc.num * f * rr_ret]
  D[,inc.ret.rr.sd:=inc.ret.rr * sqrt(
  (inc.num.sd/(inc.num+1e-9))^2+
    (f.sd/f)^2+
    (rr.sd_ret/rr_ret)^2
  )]
  ## D[,c('inc.ret.rr','inc.ret.rr.sd'):=incdr(inc.num,
  ##                                           inc.num.sd,
  ##                                           f,
  ##                                           f.sd,
  ##                                           pn = 0,
  ##                                           pn.sd = 0,
  ##                                           r,
  ##                                           r.sd,
  ##                                           rho = rho,
  ##                                           rho.sd = rho.sd,
  ##                                           pr = rr_ret,
  ##                                           pr.sd = rr.sd_ret,
  ##                                           secondorder),
  ##   by=.(iso3,year)]
  ## prop overall
  D[,rr.prop:=((1 - f) * rr_new * ((1 - r) + r * rho) + f * rr_ret)]
  xx <- D[,{
    B <- ((1 - r) + r * rho);
    dq <- c(
    (rr_ret - rr_new * B),        #dq/df
    (1-f) * B,            #dq/dpn
    f,                    #dq/dpr
    (1-f) * rr_new * (rho-1), #dq/dr
    (1-f) * rr_new * r        #dq/drho
    );
    v.sd <- c(f.sd,rr.sd_new,rr.sd_ret,r.sd,rho.sd);
    q.v <- sum(dq^2*v.sd^2);
    list(rr.prop.sd=sqrt(q.v))
  },by=.(iso3,year)]
  D <- merge(D,xx,by=c('iso3','year'))
  D[,c('rr.prop.lo','rr.prop.hi'):=.(pmax(0,rr.prop-1.96*rr.prop.sd),
                                     pmin(1,rr.prop+1.96*rr.prop.sd)
                                    )]
  ## all incidence in new
  D[,inc.new:=inc.num * (1-f)]
  D[,inc.new.sd:=inc.new * sqrt((f.sd/f)^2+(inc.num.sd/(inc.num+1e-9))^2)]
  ## all incidence in ret
  D[,inc.ret:=inc.num * f ]
  D[,inc.ret.sd:=inc.ret * sqrt((f.sd/f)^2+(inc.num.sd/(inc.num+1e-9))^2)]
  ## return
  D
}

## this calculates PGs props - temporary for use in comparing f/r updates

## add mdr regions
dre <- merge(dre,
             est[year == yr, .(iso3, g.mdr)],
             by = 'iso3', all.x = TRUE)


## === PROPORTIONS RR ===
## --- midpoints ---

## NEW: prop RR, surveillance >= 2017
sel1 <- dre$source.new == 'Surveillance' &
  dre$year.new >= 2017
table(sel1)
dre[sel1, rr.new.Num := rr.new]
dre[sel1, rr.new.Den := r.rlt.new]
dre[sel1, prop.rr.new := rr.new / r.rlt.new, by = iso3]
dre[sel1, test.isbinom(prop.rr.new)] #check

## RET: prop RR, surveillance >= 2017
sel2 <- dre$source.ret == 'Surveillance' &
  dre$year.ret >= 2017 &
  dre$r.rlt.ret > 0 &
  !is.na(dre$rr.ret)
table(sel2)
dre[sel2, rr.ret.Num := rr.ret]
dre[sel2, rr.ret.Den := r.rlt.ret]
dre[sel2, prop.rr.ret := rr.ret / r.rlt.ret, by = iso3]
dre[sel2, test.isbinom(prop.rr.ret)] #check

## NEW: prop RR, surveillance < 2017
sel3 <- dre$source.new == 'Surveillance' &
  dre$year.new < 2017
table(sel3)
dre[sel3, rr.new.Num :=
            as.integer(rowSums(cbind(dr.r.nh.new, mdr.new, xpert.dr.r.new),
                               na.rm = TRUE))]
dre[sel3, rr.new.Den :=
            as.integer(rowSums(cbind(dst.rlt.new, xpert.new),
                               na.rm = TRUE))]
dre[sel3, prop.rr.new := rr.new.Num / rr.new.Den]
dre[sel3 & rr.new.Den == 0, prop.rr.new := NA]
dre[is.nan(prop.rr.new), prop.rr.new := NA]
sel3 <- sel3 & !is.na(dre$prop.rr.new)
dre[sel3, test.isbinom(prop.rr.new)] #check

## RET: prop RR, surveillance < 2017
sel4 <- dre$source.ret == 'Surveillance' &
  dre$year.ret < 2017
table(sel4)
dre[sel4, rr.ret.Num :=
            as.integer(rowSums(cbind(dr.r.nh.ret, mdr.ret, xpert.dr.r.ret),
                               na.rm = TRUE))]
dre[sel4, rr.ret.Den := as.integer(rowSums(cbind(dst.rlt.ret, xpert.ret),
                                           na.rm = TRUE))]
dre[sel4, prop.rr.ret := rr.ret.Num / rr.ret.Den]
dre[sel4 & rr.ret.Den == 0, prop.rr.ret := NA]
dre[is.nan(prop.rr.ret), prop.rr.ret := NA]
sel4 <- sel4 & !is.na(dre$prop.rr.ret)
dre[sel4, test.isbinom(prop.rr.ret)] #check


## --- uncertainty ---
## SDs prop RR, surveillance
## assume binomial errors

## NEW
out <- dre[sel1 | sel3, {
  tmp = cii(rr.new.Den, rr.new.Num)
  list(
    prop.rr.new.sd = tmp$se,
    prop.rr.new.lo = tmp$lower95ci,
    prop.rr.new.hi = tmp$upper95ci
  )
}, by = .(iso3)]
dre <- merge(dre, out, by = 'iso3', all.x = T)
dre[sel1 | sel3, test.isbinom(prop.rr.new.sd)]

## RET
out <- dre[sel2 | sel4, {
  tmp = cii(rr.ret.Den, rr.ret.Num)
  list(
    prop.rr.ret.sd = tmp$se,
    prop.rr.ret.lo = tmp$lower95ci,
    prop.rr.ret.hi = tmp$upper95ci
  )
}, by = .(iso3)]
dre <- merge(dre, out, by = 'iso3', all.x = T)
dre[sel2 | sel4, test.isbinom(prop.rr.ret.sd)]


## NEW zero empiric
(dre[prop.rr.new == 0 | prop.rr.ret == 0,
     .(iso3,
       prop.rr.new,
       prop.rr.new.lo,
       prop.rr.new.hi,
       prop.rr.new.sd)])

## assign midpoint of CI
dre[prop.rr.new.sd == 0,
    prop.rr.new.sd := (prop.rr.new.hi - prop.rr.new.lo) / 3.92]
dre[prop.rr.ret.sd == 0,
    prop.rr.ret.sd := (prop.rr.ret.hi - prop.rr.ret.lo) / 3.92]

## prop RR, survey >=2018
## NEW
sel5 <- dre$source.new != 'Surveillance' &
  dre$year.new >= 2018
table(sel5)
dre[sel5, prop.rr.new := rr.new.pcnt / 100]
dre[sel5, prop.rr.new.lo := rr.new.pcnt.lo / 100]
dre[sel5, prop.rr.new.hi := rr.new.pcnt.hi / 100]
dre[sel5, prop.rr.new.sd := (prop.rr.new.hi - prop.rr.new.lo) / 3.92]

## RET
sel6 <- dre$source.ret != 'Surveillance' &
  dre$year.ret >= 2018
table(sel6)
dre[sel6, prop.rr.ret := rr.ret.pcnt / 100]
dre[sel6, prop.rr.ret.lo := rr.ret.pcnt.lo / 100]
dre[sel6, prop.rr.ret.hi := rr.ret.pcnt.hi / 100]
dre[sel6, prop.rr.ret.sd := (prop.rr.ret.hi - prop.rr.ret.lo) / 3.92]


## prop RR, survey <2018, phenotypic DST

## NEW
sel7 <- dre$source.new != 'Surveillance' &
  dre$year.new < 2018
table(sel7)
dre[sel7 &
    !is.na(rr.new.pcnt),
    prop.rr.new := rr.new.pcnt / 100]
    ## by = 'iso3']
dre[sel7 &
    !is.na(rr.new.pcnt),
    prop.rr.new.sd := (rr.new.pcnt.hi - rr.new.pcnt.lo) / 3.92 / 100]
dre[sel7 &
    !is.na(rr.new.pcnt),
    prop.rr.new.lo := rr.new.pcnt.lo / 100]
dre[sel7 &
    !is.na(rr.new.pcnt),
    prop.rr.new.hi := rr.new.pcnt.hi / 100]

## NEW
dre[sel7 &
    is.na(rr.new.pcnt),
    prop.rr.new :=
      rowSums(cbind(dr.r.nh.new.pcnt, mdr.new.pcnt) / 100,
              na.rm = TRUE)]
dre[sel7 &
    is.na(rr.new.pcnt),
    prop.rr.new.sd :=
      sqrt(
      ((dr.r.nh.new.pcnt.hi - dr.r.nh.new.pcnt.lo) * 1e-2 / 3.92) ^ 2 +
      ((mdr.new.pcnt.hi - mdr.new.pcnt.lo) * 1e-2 / 3.92) ^ 2
      )]
dre[sel7 &
    is.na(prop.rr.new.sd),
    prop.rr.new.sd := (mdr.new.pcnt.hi - mdr.new.pcnt.lo) * 1e-2 / 3.92]

## check:
dre[sel7 &
      !is.na(prop.rr.new) &
      is.na(prop.rr.new.sd), 
    .(iso3,
      source.new,
      year.new,
      dr.r.nh.new.pcnt,
      mdr.new.pcnt,
      prop.rr.new,
      prop.rr.new.sd,
      mdr.new.pcnt.hi,
      mdr.new.pcnt.lo,
      dr.r.nh.new.pcnt.hi,
      dr.r.nh.new.pcnt.lo
      )]

## for 1 case with missing SD, set midpoint NA
dre[sel7 &
    !is.na(prop.rr.new) &
    is.na(prop.rr.new.sd),
    prop.rr.new := NA]



## NEW bug fix (Matteo, July 2018)
## user Xpert if others missing
sel <-
  sel7 &
  (is.na(dre$prop.rr.new.lo) |
     is.na(dre$prop.rr.new.hi)) &
  !is.na(dre$xpert.dr.r.new.pcnt.lo) &
  !is.na(dre$xpert.dr.r.new.pcnt.hi)
table(sel)
dre[sel, prop.rr.new := xpert.dr.r.new.pcnt / 100]
dre[sel, prop.rr.new.lo := xpert.dr.r.new.pcnt.lo / 100]
dre[sel, prop.rr.new.hi := xpert.dr.r.new.pcnt.hi / 100]
dre[sel, prop.rr.new.sd := (prop.rr.new.hi - prop.rr.new.lo) / 3.92]


## RET:
sel8 <- dre$source.ret != 'Surveillance' &
  dre$year.ret < 2018
table(sel8)
dre[sel8 &
    !is.na(rr.ret.pcnt),
    prop.rr.ret := rr.ret.pcnt / 100]
## , by = 'iso3']
dre[sel8 &
    !is.na(rr.ret.pcnt),
    prop.rr.ret.sd := (rr.ret.pcnt.hi - rr.ret.pcnt.lo) * 1e-2 /  3.92]
dre[sel8 &
    !is.na(rr.ret.pcnt),
    prop.rr.ret.lo := rr.ret.pcnt.lo / 100]
dre[sel8 &
    !is.na(rr.ret.pcnt),
    prop.rr.ret.hi := rr.ret.pcnt.hi / 100]


## RET:
dre[sel8 &
    is.na(rr.ret.pcnt),
    prop.rr.ret :=
      rowSums(cbind(dr.r.nh.ret.pcnt, mdr.ret.pcnt) / 100, na.rm=TRUE)]
dre[sel8 &
    is.na(rr.ret.pcnt),
    prop.rr.ret.sd :=
      sqrt(
      ((dr.r.nh.ret.pcnt.hi - dr.r.nh.ret.pcnt.lo) * 1e-2 / 3.92) ^ 2
         +
         ((mdr.ret.pcnt.hi - mdr.ret.pcnt.lo) * 1e-2 / 3.92) ^ 2
      )]
dre[sel8 &
    is.na(prop.rr.ret.sd),
    prop.rr.ret.sd := (mdr.ret.pcnt.hi - mdr.ret.pcnt.lo) * 1e-2 / 3.92]

## check
dre[sel8 &
      !is.na(prop.rr.ret) &
      is.na(prop.rr.ret.sd),
    .(iso3,
      source.ret,
      year.ret,
      dr.r.nh.ret.pcnt,
      mdr.ret.pcnt,
      prop.rr.ret,
      prop.rr.ret.sd,
      mdr.ret.pcnt.hi,
      mdr.ret.pcnt.lo,
      dr.r.nh.ret.pcnt.hi,
      dr.r.nh.ret.pcnt.lo
      )]

## for 1 row with NA except midpoint, set NA
dre[sel8 &
    !is.na(prop.rr.ret) &
    is.na(prop.rr.ret.sd),
    prop.rr.ret := NA]

## RET bug fix (Matteo, July 2018)
## use Xpert if others missing
sel <-
  sel8 &
  (is.na(dre$prop.rr.ret.lo) |
     is.na(dre$prop.rr.ret.hi)) &
  !is.na(dre$xpert.dr.r.ret.pcnt.lo) &
  !is.na(dre$xpert.dr.r.ret.pcnt.hi)
table(sel)
dre[sel, prop.rr.ret := xpert.dr.r.ret.pcnt / 100]
dre[sel, prop.rr.ret.lo := xpert.dr.r.ret.pcnt.lo / 100]
dre[sel, prop.rr.ret.hi := xpert.dr.r.ret.pcnt.hi / 100]
dre[sel, prop.rr.ret.sd := (prop.rr.ret.hi - prop.rr.ret.lo) / 3.92]

## check
dre[sel7 & !is.na(prop.rr.new), test.ispos(prop.rr.new.sd)]
dre[sel8 & !is.na(prop.rr.ret), test.ispos(prop.rr.ret.sd)]


## prop RR, survey <2018, xpert
## NEW
sel9 <- dre$source.new != 'Surveillance' &
  dre$year.new %in% 2013:2017
table(sel9)
table(sel9 & is.na(dre$prop.rr.new))
sel9 <- sel9 & is.na(dre$prop.rr.new)
dre[sel9, prop.rr.new := xpert.dr.r.new.pcnt / 100]
dre[sel9, prop.rr.new.lo := xpert.dr.r.new.pcnt.lo / 100]
dre[sel9, prop.rr.new.hi := xpert.dr.r.new.pcnt.hi / 100]
dre[sel9, prop.rr.new.sd := (prop.rr.new.hi - prop.rr.new.lo) / 3.92]

## RET
sel10 <- sel8 & is.na(dre$prop.rr.ret)
table(sel10)
dre[sel10, prop.rr.ret := xpert.dr.r.ret.pcnt / 100]
dre[sel10, prop.rr.ret.lo := xpert.dr.r.ret.pcnt.lo / 100]
dre[sel10, prop.rr.ret.hi := xpert.dr.r.ret.pcnt.hi / 100]
dre[sel10, prop.rr.ret.sd := (prop.rr.ret.hi - prop.rr.ret.lo) / 3.92]
dre[sel9 & !is.na(prop.rr.new), test.ispos(prop.rr.new.sd)]
dre[sel10 & !is.na(prop.rr.ret), test.ispos(prop.rr.ret.sd)]

## prop RR, missing bounds
## NEW
sel <-
  !is.na(dre$prop.rr.new) &
  (is.na(dre$prop.rr.new.lo) |
   is.na(dre$prop.rr.new.hi))
table(sel)
out <- vlohi(dre$prop.rr.new[sel], dre$prop.rr.new.sd[sel])
dre[sel, prop.rr.new.lo := out[1,]]
dre[sel, prop.rr.new.hi := out[2,]]
dre[!is.na(prop.rr.new), test.bounds(prop.rr.new, prop.rr.new.lo, prop.rr.new.hi)]


sel <-
  !is.na(dre$prop.rr.ret) &
  (is.na(dre$prop.rr.ret.lo) |
   is.na(dre$prop.rr.ret.hi))
table(sel)
out <- vlohi(dre$prop.rr.ret[sel], dre$prop.rr.ret.sd[sel])
dre[sel, prop.rr.ret.lo := out[1,]]
dre[sel, prop.rr.ret.hi := out[2,]]

## check
dre[!is.na(prop.rr.ret), test.bounds(prop.rr.ret, prop.rr.ret.lo, prop.rr.ret.hi)]



## IMPUTATIONS
## pool props by region
suppressWarnings(fit.new <-
                   rma(
                     yi = prop.rr.new,
                     sei = prop.rr.new.sd,
                     data = dre,
                     mods =  ~ g.mdr - 1
                   ))
suppressWarnings(fit.ret <-
                   rma(
                     yi = prop.rr.ret,
                     sei = prop.rr.ret.sd,
                     data = dre,
                     mods =  ~ g.mdr - 1
                   ))

imp.new <- fit.new$b[, 1]
imp.ret <- fit.ret$b[, 1]

imp.new.sd <- fit.new$se
imp.ret.sd <- fit.ret$se
imp.new.ui <- vlohi(unlist(imp.new), imp.new.sd)
imp.ret.ui <- vlohi(unlist(imp.ret), imp.ret.sd)


# check the order of regions in reg is the same as in imp*
#
all.equal(reg, gsub('g.mdr', '', dimnames(fit.new$beta)[[1]]))
all.equal(reg, gsub('g.mdr', '', dimnames(fit.ret$beta)[[1]]))

table(!is.na(dre$prop.rr.new))
table(!is.na(dre$prop.rr.ret))

tmp <- copy(dre)

## not elegant, should set up a proper merge
for (i in 1:length(reg)) {
  sel <- dre$g.mdr == reg[i] & is.na(dre$prop.rr.new)
  dre[sel, prop.rr.new := imp.new[i]]
  dre[sel, prop.rr.new.sd := imp.new.sd[i]]
  dre[sel, prop.rr.new.lo := imp.new.ui[1, i]]
  dre[sel, prop.rr.new.hi := imp.new.ui[2, i]]

  sel <- dre$g.mdr == reg[i] & is.na(dre$prop.rr.ret)
  dre[sel, prop.rr.ret := imp.ret[i]]
  dre[sel, prop.rr.ret.sd := imp.ret.sd[i]]
  dre[sel, prop.rr.ret.lo := imp.ret.ui[1, i]]
  dre[sel, prop.rr.ret.hi := imp.ret.ui[2, i]]
}

## check
dre[, lapply(.SD, test.ispos), .SDcols = c('prop.rr.new',
                                           'prop.rr.new.sd',
                                           'prop.rr.ret',
                                           'prop.rr.ret.sd')]



## prep data, f, r
rr <-
  dre[, .(
    iso3,
    g.whoregion,
    g.mdr,
    year.new,
    source.new,
    all.areas.covered.new,
    surv.quality.new,
    year.ret,
    source.ret,
    all.areas.covered.ret,
    surv.quality.ret,
    prop.rr.new,
    prop.rr.new.sd,
    prop.rr.new.lo,
    prop.rr.new.hi,
    prop.rr.ret,
    prop.rr.ret.sd,
    prop.rr.ret.lo,
    prop.rr.ret.hi
  )]


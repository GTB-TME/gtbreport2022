## function that will identify country iso3 code cn as an index
## and then strip out all data from a stan data input DIN
## if recentOnly=TRUE, it will only remove dps with dy=0
countryomitter <- function(cn,               #country to work on
                           DIN,              #full data object
                           recentOnly=FALSE, #omit specific times only or all data
                           nt=0,             #nt=0 is most recent time
                           verbose=TRUE
                           ){
  doubt <- DIN
  cid <- ckey[iso3==cn,cid] #numeric country id
  cat('Country ',cn,', id=',cid,', nt=',nt,'\n')
  casends <- c("rr","hr","rm","x","dh0","dh1","dh2","dh1r","dh1h")
  for(dtype in c('L','Y')){
    for(pg in c('new','ret')){
      for(dcase in casends){
        ## dtype <- 'L'; pg <- 'new'; dcase <- casends[1] #NOTE for debugging
        caseroot <- glue('{dtype}_{pg}_{dcase}')
        ## print(caseroot)
        idz <- caseroot + '_id'
        tz <- caseroot + '_t'
        nz <- 'N_' + caseroot #N name
        if(cid %in% doubt[[idz]]){
          if(verbose) cat('...found in ',caseroot,'...\n')
          togo <- which(doubt[[idz]]==cid) #which indices to excise
          if(recentOnly){                   #only ditch dy=0
            togo <- togo[which(doubt[[tz]][togo]==nt)]
            if(length(togo)>0)
              if(verbose) cat('......of which 1 has dy=',nt,'\n')
          }
          if(length(togo)>0){                   #safety in case non dy=0
            doubt[[nz]] <- doubt[[nz]] - length(togo) #down length
            doubt[[idz]] <- doubt[[idz]][-togo] #remove from indices
            doubt[[tz]] <- doubt[[tz]][-togo]   #remove from times
            if(dtype=='L'){
              Nz <- caseroot + '_N'
              nc <- ncol(doubt[[Nz]])
              if(is.null(nc)){ #non-matrix form
                doubt[[Nz]] <- doubt[[Nz]][-togo] #remove from counts
              } else {
                doubt[[Nz]] <- doubt[[Nz]][-togo,] #remove from counts
                if(is.null(ncol(doubt[[Nz]]))){
                  doubt[[Nz]] <- matrix(doubt[[Nz]],ncol=nc) #safety
                  attr(doubt[[idz]],'dim') <- 1
                  attr(doubt[[tz]],'dim') <- 1
                }
                ## print(dim(doubt[[Nz]]))
              }
            } else {                             #surveillance
              Mz <- caseroot + '_M'
              Sz <- caseroot + '_S'
              nc <- ncol(doubt[[Mz]])
              if(is.null(nc)){ #non-matrix form
                doubt[[Mz]] <- doubt[[Mz]][-togo]
                doubt[[Sz]] <- doubt[[Sz]][-togo]
              } else { #matrix form
                doubt[[Mz]] <- doubt[[Mz]][-togo,]
                doubt[[Sz]] <- doubt[[Sz]][-togo,]
                if(is.null(ncol(doubt[[Mz]]))){ #safety
                  doubt[[Mz]] <- matrix(doubt[[Mz]],ncol=nc)
                  doubt[[Sz]] <- matrix(doubt[[Sz]],ncol=nc)
                  attr(doubt[[idz]],'dim') <- 1
                  attr(doubt[[tz]],'dim') <- 1
                }
              }
            }
          }
        }
      }
    }
  }
  doubt
}

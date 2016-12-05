## script to re-run nimble model on experiments that previosly failed

## load needed libraries and functions
setwd('~/Dropbox/hawaiiDimensions/metabarcoding/poc')
source('nimble_funs.R')

## load data from different experiments
diffMarkers <- read.csv('clean_diffMarkers.csv', as.is = TRUE)
diffCycles <- read.csv('clean_pcrCycle.csv', as.is = TRUE)
diffTissue <- read.csv('clean_tissue.csv', as.is = TRUE)

## load mcmc output
load('nimble_out.RData')

## helper function to automate re-running
reRunNimble <- function(dat, output, x) {
    d <- .formatData(dat, x)
    fail <- which(output[[2]]$nGewekeFail > 
                      0.1 * sapply(1:dim(d[[3]])[1], 
                                   function(i) ncol(.prepData(d, i)[[3]])) | 
                      is.na(output[[2]]$nGewekeFail))
    
    ## re-run for failures
    for(i in fail) {
        # redo <- .internalRunNimble(d, i, N = 10000, thin = 50, burn = 50)
        redo <- .internalRunNimble(d, i, N = 100, thin = 1, burn = 1)
        output[[2]][i, ] <- redo$summ
        output[[1]][[i]] <- redo$par
    }
    
    return(output)
}


## re-run nimble where needed
diffMarkerOut <- reRunNimble(diffMarkers, diffMarkerOut, 'Primer')
diffCyclesOut <- reRunNimble(diffCycles, diffCyclesOut, 'PCR_cycles')
diffTissueOut <- reRunNimble(diffTissue, diffTissueOut, 'Experiment')

## save it
save(diffMarkerOut, diffCyclesOut, diffTissueOut, file = 'nimble_out.RData')

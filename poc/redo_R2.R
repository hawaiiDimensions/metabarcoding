## script to re-calculate R2 for models already run

## load needed libraries and functions
setwd('~/Dropbox/hawaiiDimensions/metabarcoding/poc')
source('nimble_funs.R')

## load data from different experiments
diffMarkers <- read.csv('clean_diffMarkers.csv', as.is = TRUE)
diffCycles <- read.csv('clean_pcrCycle.csv', as.is = TRUE)
diffTissue <- read.csv('clean_tissue.csv', as.is = TRUE)

## load mcmc output
load('nimble_out.RData')

## helper function to calculate R2
R2redo <- function(dat, output, x) {
    fdat <- .formatData(dat, x)
    out <- sapply(1:nrow(fdat[[1]]), function(i) {
        pdat <- .prepData(fdat, i)
        bayesR2(y = pdat[[3]], n = pdat[[1]], x = pdat[[2]], output[[1]][[i]])
    })
    
    out <- data.frame(t(out))
    return(out)
}

## recalculate R2
diffMarkerOut$summ[, 2:4] <- R2redo(diffMarkers, diffMarkerOut, 'Primer')
diffCyclesOut$summ[, 2:4] <- R2redo(diffCycles, diffCyclesOut, 'PCR_cycles')
diffTissueOut$summ[, 2:4] <- R2redo(diffTissue, diffTissueOut, 'Experiment')

## re-write output
save(diffMarkerOut, diffCyclesOut, diffTissueOut, file = 'nimble_out.RData')

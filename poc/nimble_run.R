## script to run nimble model on all experiments

## load needed libraries and functions
setwd('~/Dropbox/hawaiiDimensions/metabarcoding/poc')
source('nimble_funs.R')

## load data from different experiments
diffMarkers <- read.csv('clean_diffMarkers.csv', as.is=TRUE)
diffCycles <- read.csv('clean_pcrCycle.csv', as.is=TRUE)

## run nimble models
diffMarkerOut <- buildNRunMod(diffMarkers, 'Primer')
diffCyclesOut <- buildNRunMod(diffCycles, 'PCR_cycles')
diffTissueOut <- NULL

save(diffMarkerOut, diffCyclesOut, diffTissueOut, file = 'nimble_out.RData')

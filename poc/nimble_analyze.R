## script to run nimble model on all experiments

## load needed libraries and functions
setwd('~/Dropbox/hawaiiDimensions/metabarcoding/poc')
source('nimble_funs.R')

## load data from different experiments
diffMarkers <- read.csv('clean_diffMarkers.csv', as.is = TRUE)
diffCycles <- read.csv('clean_pcrCycle.csv', as.is = TRUE)
diffTissue <- read.csv('clean_tissue.csv', as.is = TRUE)

## load mcmc output
load('nimble_out.RData')

diffTissueOut[[2]]
effectiveSize(diffTissueOut[[1]][[1]])

plot(diffTissueOut[[1]][[1]][, 1], type = 'l')


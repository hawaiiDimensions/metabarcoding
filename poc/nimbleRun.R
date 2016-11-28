setwd('~/Dropbox/hawaiiDimensions/metabarcoding/poc')
source('runNimble.R')

diffMarkers <- read.csv('clean_diffMarkers.csv', as.is=TRUE)
diffCycles <- read.csv('clean_pcrCycle.csv', as.is=TRUE)

## number or reads for each primer and pool combination (rows:cycles, cols:pools)
totReads <- acast(melt(diffCycles, c('PCR_cycles', 'Pool'), 'total_Reads'), PCR_cycles ~ Pool, 
                  value.var =  'value', max, na.rm = TRUE, fill = 0)

## amount of DNA from each species in each pool (rows:pools, cols:species)
amountDNA <- acast(melt(diffCycles, c('Pool', 'Specimen'), 'amount_DNA'), Pool ~ Specimen, 
                   value.var =  'value', max, na.rm = TRUE, fill = 0)

## number of reads per primer, pool, species combo (dim1:cycles, dim2:pool, dim3:species)
numReads <- acast(melt(diffCycles, c('PCR_cycles', 'Pool', 'Specimen'), 'number_Reads'), 
                  PCR_cycles ~ Pool ~ Specimen, 
                  value.var =  'value', max, na.rm = TRUE, fill = 0)

## run nimble models
diffMarkerOut <- buildNRunMod(diffMarkers, 'Primer')
diffCyclesOut <- buildNRunMod(diffCycles, 'PCR_cycles')

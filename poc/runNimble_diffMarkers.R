library(nimble)
library(coda)
library(plyr)
library(reshape2)
library(parallel)

setwd('~/Dropbox/hawaiiDimensions/metabarcoding/poc')
source('runNimble.R')

diffMarkers <- read.csv('clean_diffMarkers.csv', as.is=TRUE)

## number or reads for each primer and pool combination (rows:primers, cols:pools)
totReads <- acast(melt(diffMarkers, c('Primer', 'Pool'), 'total_Reads'), Primer ~ Pool, 
                  value.var =  'value', max, na.rm = TRUE, fill = 0)

## amount of DNA from each species in each pool (rows:pools, cols:species)
amountDNA <- acast(melt(diffMarkers, c('Pool', 'Specimen'), 'amount_DNA'), Pool ~ Specimen, 
                   value.var =  'value', max, na.rm = TRUE, fill = 0)

## number of reads per primer, pool, species combo (dim1:primer, dim2:pool, dim3:species)
numReads <- acast(melt(diffMarkers, c('Primer', 'Pool', 'Specimen'), 'number_Reads'), 
                  Primer ~ Pool ~ Specimen, 
                  value.var =  'value', max, na.rm = TRUE, fill = 0)

## loop over primers, fitting model to each and calculating: 
## R2 
## effective sample size
## Geweke's convergence test

outDiffMarkers <- lapply(1:nrow(totReads), funciton(i) {
    ## model parameters
    modPar <- runNimble(totReads[i, ], amountDNA, numReads[i, , ], 
                        N = 10000, thin = 50, burn = 100)
    
    ## return R2 and (across all a's) min effective size and Geweke's test
    return(bayesR2(numReads[i, , ], totReads[i, ], modPar), 
           minESS = min(effectiveSize(modPar)), 
           nGewekeFail = sum(abs(geweke.diag(bla)$z) > 1.96))
})


library(nimble)
library(coda)
library(plyr)
library(reshape2)
library(parallel)

setwd('~/Dropbox/hawaiiDimensions/metabarcoding/poc')
source('runNimble.R')

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

## loop over primers, fitting model to each and calculating: 
## R2 
## effective sample size
## Geweke's convergence test

outDiffCycles <- lapply(1:nrow(totReads), funciton(i) {
    ## model parameters
    modPar <- runNimble(totReads[i, ], amountDNA, numReads[i, , ], 
                        N = 10000, thin = 50, burn = 100)
    
    ## return R2 and (across all a's) min effective size and Geweke's test
    return(bayesR2(numReads[i, , ], totReads[i, ], modPar), 
           minESS = min(effectiveSize(modPar)), 
           nGewekeFail = sum(abs(geweke.diag(bla)$z) > 1.96))
})

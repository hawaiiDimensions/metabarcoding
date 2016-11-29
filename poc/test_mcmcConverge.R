## script to test rough convergence of model

setwd('~/Dropbox/hawaiiDimensions/metabarcoding/poc')
source('nimble_funs.R')

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

mod <- mclapply(1:2, mc.cores = 2, FUN = function(i) {
    runNimble(totReads[8, ], amountDNA, numReads[8, , ], N = 10000, thin = 100, burn = 5)
})

save(mod, file = 'test_mcmc.RData')

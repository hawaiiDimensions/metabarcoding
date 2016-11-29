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

effectiveSize(mod[[1]])
effectiveSize(mod[[2]])

plot(mod[[1]][-(1:100), 24], type = 'l')
points(mod[[2]][-(1:100), 24], type = 'l', col = hsv(1, 1, 1, alpha = 0.5))

plot(density(mod[[1]][-(1:100), 38]))
lines(density(mod[[2]][-(1:100), 38]), col = hsv(1, 1, 1, alpha = 0.5))
curve(dexp(x, 0.00001), add = TRUE, col = 'blue')

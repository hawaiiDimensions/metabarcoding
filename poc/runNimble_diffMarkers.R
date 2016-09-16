library(nimble)
library(plyr)
library(reshape2)

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
numReads <- acast(melt(diffMarkers, c('Primer', 'Pool', 'Specimen'), 'number_Reads'), Primer ~ Pool ~ Specimen, 
                  value.var =  'value', max, na.rm = TRUE, fill = 0)

bla <- runNimble(totReads[8, ], amountDNA, numReads[8, , ], N = 10000, thin = 50, burn = 100)

plot(bla[-(1:100), 1], type = 'l')
plot(density(bla[-(1:100), 1]))
acf(bla[-(1:100), 1])

plot(bla[-(1:500), 1], ylim = range(bla[-(1:400), ]), type = 'l')
for(i in 2:ncol(bla)) lines(bla[-(1:400), i])

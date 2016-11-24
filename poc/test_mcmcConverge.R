library(nimble)
library(coda)
library(plyr)
library(reshape2)
library(parallel)

setwd('~/Dropbox/hawaiiDimensions/metabarcoding/poc')
source('runNimble.R')

bla <- runNimble(totReads[8, ], amountDNA, numReads[8, , ], N = 10000, thin = 10, burn = 100)

mean(effectiveSize(bla))

bayesR2(numReads[8, , ], totReads[8, ], bla)

plot(bla[, 1], type = 'l')
plot(density(bla[-(1:100), 1]))

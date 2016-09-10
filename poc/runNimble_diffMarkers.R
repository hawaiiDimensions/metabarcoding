library(nimble)
library(plyr)
library(reshape2)

setwd('~/Dropbox/hawaiiDimensions/metabarcoding/poc')
source('runNimble.R')

diffMarkers <- read.csv('clean_diffMarkers.csv', as.is=TRUE)

totReads <- acast(melt(diffMarkers, c('Primer', 'Pool'), 'total_Reads'), Primer ~ Pool, 
                  value.var =  'value', max, na.rm = TRUE, fill = 0)
amountDNA <- acast(melt(diffMarkers, c('Pool', 'Specimen'), 'amount_DNA'), Pool ~ Specimen, 
                   value.var =  'value', max, na.rm = TRUE, fill = 0)
numReads <- acast(melt(diffMarkers, c('Primer', 'Pool', 'Specimen'), 'number_Reads'), Primer ~ Pool ~ Specimen, 
                  value.var =  'value', max, na.rm = TRUE, fill = 0)

bla <- runNimble(totReads[8, ], amountDNA, numReads[8, , ])

apply(numReads, 1, function(x) sum(x == 0) / prod(dim(x)))

plot(bla[-(1:400), 1], ylim = range(bla[-(1:400), ]), type = 'l')
for(i in 2:ncol(bla)) lines(bla[-(1:400), i])

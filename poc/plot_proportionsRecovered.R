library(plyr)
library(reshape2)

setwd('~/Dropbox/hawaiiDimensions/metabarcoding/poc')

## load data
diffMarkers <- read.csv('clean_diffMarkers.csv', as.is=TRUE)
pcrCycle <- read.csv('clean_pcrCycle.csv', as.is=TRUE)
# tissue <- read.csv(...)

## function to compute proportions
calcProp <- function(dat, splitVar) {
    ddply(dat, splitVar, function(d) {
        prop <- ddply(d, 'Pool', function(dd) {
            length(unique(dd$Specimen[dd$number_Reads > 0])) / 
                length(unique(dat$Specimen[dat$Pool == dd$Pool[1]]))
        })
        
        out <- c(mean(prop$V1), quantile(prop$V1, c(0.025, 0.975)))
        out[2:3] <- out[1] + c(-1, 1) * mean(abs(out[2:3] - out[1])) / 1.96
        
        names(out) <- c('mean', 'sdLo', 'sdHi')
        
        return(out)
    })
}

## function to plot proportion
plotProp <- function(x, col = 'white', ...) {
    plot(x$mean, ylim = c(0, max(x[, -1])),
         xaxt = 'n', pch = 21, bg = col, 
         xlim = c(1, nrow(x)) + c(-1, 1) * 0.05 * nrow(x),
         panel.first = {
             lo <- x$sdLo
             hi <- x$sdHi
             bad <- lo == hi
             lo[bad] <- lo[bad] - .Machine$double.eps^0.2
             hi[bad] <- hi[bad] + .Machine$double.eps^0.2
             arrows(x0 = 1:nrow(x), y0 = lo, y1 = hi,
                              code = 3, length = 0.05, angle = 90)
         }, ...)
    axis(1, at = 1:nrow(x), labels = x[, 1])
}


## vectors of nuclear v. mitochondrial marker names
nuc <- c('18sMach', '18sSSU', '28sMach', 'H3')
names(nuc) <- c('18sM', '18sS', '28sM', 'H3')
mit <- c('12s_', 'CytB_', 'ARFHCO', 'MCOHCO')
names(mit) <- c('12s', 'CytB', 'ARF', 'MCO')

## proportions taxa recovered

propRecoverMarker <- calcProp(diffMarkers, 'Primer')
propRecoverNuc <- propRecoverMarker[propRecoverMarker$Primer %in% nuc, ]
propRecoverNuc$Primer <- names(nuc)[match(propRecoverNuc$Primer, nuc)]
propRecoverMit <- propRecoverMarker[propRecoverMarker$Primer %in% mit, ]
propRecoverMit$Primer <- names(mit)[match(propRecoverMit$Primer, mit)]
propRecoverMarker <- rbind(propRecoverNuc, propRecoverMit)

propRecoverCycle <- calcProp(pcrCycle, 'PCR_cycles')
# propRecoverTissue <- calcProp(tissue, 'Experiment')


## plot it

pdf('ms/fig_propRecover.pdf', width = 7, height = 3)
par(mfrow = c(1, 4), mgp = c(2, 0.75, 0), 
    mar = c(3, 0, 0, 0) + 0.5, oma = c(0, 3.5, 0, 0),
    cex.lab = 1.2)

plotProp(propRecoverNuc, xlab = '', cex = 1.5)
text(2.5, 0, labels = 'Nuclear')
mtext('Marker', side = 1, at = 4.5, line = 2, cex = par('cex')*par('cex.lab'))
mtext('Proportion taxa recovered', side = 2, line = 2, outer = TRUE, 
      cex = par('cex')*par('cex.lab'))

plotProp(propRecoverMit, xlab = '', yaxt = 'n', cex = 1.5)
text(2.5, 0, labels = 'Mitochondrial')

plotProp(propRecoverCycle, xlab = 'Number of cycles', yaxt = 'n', cex = 1.5)

plot(1, xlab = 'Experiment', yaxt = 'n')# plotProp(propRecoverTissue)

dev.off()

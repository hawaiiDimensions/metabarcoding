## script to analyze nimble model on all experiments

## load needed libraries and functions
setwd('~/Dropbox/hawaiiDimensions/metabarcoding/poc')
source('nimble_funs.R')

## load data from different experiments
diffMarkers <- read.csv('clean_diffMarkers.csv', as.is = TRUE)
diffCycles <- read.csv('clean_pcrCycle.csv', as.is = TRUE)
diffTissue <- read.csv('clean_tissue.csv', as.is = TRUE)

## load mcmc output
load('nimble_out.RData')

## vectors of nuclear v. mitochondrial marker names for use in ordering plots
nuc <- c('18sSSU', '18sMach', '28sMach', 'H3')
names(nuc) <- c('18sSSU', '18sM', '28sM', 'H3')
mit <- c('12s_', 'CytB_', 'MCOHCO', 'ARFHCO')
names(mit) <- c('12s', 'CytB', 'MCO', 'ARF')

## ===========
## plotting R2
## ===========

## function to plot R2 values
plotR2 <- function(x, ...) {
    r2mean <- x$R2.mean
    r2ci <- x[, c('R2.ciLo', 'R2.ciHi')]
    
    plot(r2mean, pch = 21, bg = 'white', 
         panel.first = arrows(x0 = 1:length(r2mean), y0 = r2ci[, 1], y1 = r2ci[, 2], 
                              code = 3, length = 0.05, angle = 90), 
         ...)
}

## plotting

pdf('ms/fig_r2.pdf', width = 6, height = 4)
par(mfrow = c(1, 4), oma = c(0, 3, 0, 0) + 0.5, mar = c(6, 0, 3, 0) + 0.5, cex.axis = 1.2)
ylim <- c(0, 1)

## nuc
plotR2(diffMarkerOut[[2]][match(nuc, diffMarkerOut[[2]]$Primer), ], cex = 1.5, 
       xaxt = 'n', xlim = c(0.5, 4.5), ylim = ylim, xlab = '', ylab = '')
axis(1, at = 1:4, labels = names(nuc), las = 2)
mtext(expression('Posterior'~R^2), side = 2, line = 2, outer = TRUE)
mtext('Locus', side = 1, line = 5)
mtext('Nuclear', side = 3, line = 1)

## mit
plotR2(diffMarkerOut[[2]][match(mit, diffMarkerOut[[2]]$Primer), ], cex = 1.5, 
       xaxt = 'n', yaxt = 'n', xlim = c(0.5, 4.5), ylim = ylim, xlab = '', ylab = '')
axis(1, at = 1:4, labels = names(mit), las = 2)
mtext('Locus', side = 1, line = 5)
mtext('Mitochondrial', side = 3, line = 1)

## cycles
plotR2(diffCyclesOut[[2]], cex = 1.5, 
       xaxt = 'n', yaxt = 'n', xlim = c(0.5, 4.5), ylim = ylim, xlab = '', ylab = '')
axis(1, at = 1:4, labels = diffCyclesOut[[2]]$PCR_cycles, las = 2)
mtext('Cycles', side = 1, line = 5)
mtext('PCR', side = 3, line = 1)

## tissue
plotR2(diffTissueOut[[2]], cex = 1.5, 
       xaxt = 'n', yaxt = 'n', xlim = c(0.5, 4.5), ylim = ylim, xlab = '', ylab = '')
axis(1, at = 1:4, labels = diffTissueOut[[2]]$Experiment, las = 2)
mtext('Experiment', side = 1, line = 5)
mtext('Tissue', side = 3, line = 1)

dev.off()


## ==============================
## plotting corrrected read count
## ==============================


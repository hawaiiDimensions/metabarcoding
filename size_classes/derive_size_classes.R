setwd('~/Dropbox/hawaiiDimensions/metabarcoding')

## data from edward's sorting 9/28/16
x <- c(3.5, 5.5, 14, 2.5, 9, 12.5, 4.5, 6, 1, 3.5, 7, 2, 4.5, 4.5, 2.5, 
       3.5, 6, 1.5, 1, 1.5, 4, 6.5, 1, 4, 2.5, 7, 10, 4.5, 6, 1, 4.5, 7, 
       2, 4, 4.5, 3.5, 3.5, 5, 1, 1.5, 3, 5, 1, 2, 6.5, 3.5, 5, 1, 2.5, 
       7.5, 2, 4, 4, 1, 1.5, 4.5, 5, 1, 2.5, 5, 3, 8, 1.5, 5, 2, 4.5, 3.5, 
       1, 1.5, 4.5, 8, 1, 3.5, 8, 4.5, 7, 3, 5, 2, 4, 1.5, 3, 8, 1, 2.5, 
       5.5, 2, 5, 3, 7, 1.5, 3.5, 2, 2.5, 8.5, 1, 4, 7.5, 3.5, 5, 3, 5.5, 
       2, 3.5, 2, 4, 8.5, 1, 1.5, 7, 4, 5.5, 3, 5, 2, 2, 4, 5, 1, 1.5, 3, 
       5, 2.5, 1.5, 2, 4.5, 9.5, 1, 4, 3.5, 6, 4.5, 2, 1.5, 3, 7.5, 1, 2.5)

## proposed size classes
cuts <- c(0, cumsum(round(1.2^(0:7))))

den <- density(x, from = 0, to = max(x))

pdf('fig_sizeClass.pdf', width = 6, height = 4)
par(mar = c(3.5, 3.5, 0, 1) + 0.1, mgp = c(2.5, 1, 0))

plot(den, xlim = range(cuts, x), ylim = c(0, max(1.25*den$y)), 
     panel.first = {
         for(i in 1:(length(cuts) - 1)) {
             rect(xleft = cuts[i], xright = cuts[i+1], 
                  ybottom = par('usr')[3], ytop = par('usr')[4],
                  col = c('gray75', 'gray95')[i %% 2 + 1], border = NA)
         }
     }, 
     axes = FALSE, xlab = 'Length (mm)', ylab = '', main = '')

axis(1, at = cuts, labels = as.character(cuts))
axis(2, at = pretty(0.85*den$y))
mtext('Probability density', side = 2, line = 2.5, at = mean(range(0.85*den$y)))
axis(2, at = par('usr')[4] - 0.1*diff(par('usr')[3:4]), labels = 'counts:', tick = FALSE, las = 2, hadj = 0.5)

text((cuts[-1] + cuts[-length(cuts)]) / 2, rep(par('usr')[4] - 0.1*diff(par('usr')[3:4]), length(cuts) - 1),
     labels = hist(x, breaks = cuts, plot = FALSE)$count)

dev.off()

library(socorro)

setwd('~/Dropbox/hawaiiDimensions/metabarcoding/size_classes')

## data from edward's sorting 
specCount <- read.csv('specimen_counts.csv', as.is = TRUE)
specLength <- read.csv('specimen_lengths.csv', as.is = TRUE)

## edit interval extremes
specCount$intervalMax[nrow(specCount)] <- max(specLength$length)
specCount$intervalMin[2:nrow(specCount)] <- specCount$intervalMax[1:(nrow(specCount) - 1)]

## split 1mm bin counts by bootstrapping measured specimens with probabilities given by counts
p <- specCount$count[cut(specLength$length, c(0, specCount$intervalMax), labels = 1:nrow(specCount))]
breaks <- seq(0, max(specLength$length), by = 0.5)
den <- c(rowSums(replicate(100, {
    l <- sample(specLength$length, size = sum(specCount$count), replace = TRUE, prob = p)
    hist(l, breaks = breaks, plot = FALSE)$counts
})) / 100, 0)

## make a probability density out of that, mimicing an object of class `density'
den <- den / (sum(den) * 0.5)
den <- list(x = breaks, y = den)


## proposed size classes
cuts <- c(0, cumsum(round(1.2^(0:4))))
cuts <- c(cuts, 10, max(specLength$length) + 1)

## counts for each size class
counts <- as.numeric(tapply(specCount$count, cut(specCount$intervalMax, cuts), sum))

## plot it
pdf('fig_sizeClass.pdf', width = 6, height = 4)
par(mar = c(3.5, 3.5, 0, 1) + 0.1, mgp = c(2.5, 1, 0))

plot(den, xlim = range(cuts, den$x), ylim = c(0, max(1.25*den$y)), 
     type = 's', axes = FALSE, xlab = 'Length (mm)', ylab = '', main = '',
     panel.first = {
         for(i in 1:(length(cuts) - 1)) {
             rect(xleft = cuts[i], xright = cuts[i+1], 
                  ybottom = 0, ytop = par('usr')[4],
                  col = c('gray75', 'gray95')[i %% 2 + 1], border = NA)
         }
     })
polygon(c(rep(den$x, each = 2)[-1], max(den$x)+0.5), rep(den$y, each = 2), col = 'black')
abline(v = cuts, col = 'white')

axis(1, at = cuts, 
     labels = as.character(ifelse(cuts == max(cuts), '', cuts)))
axis(1, at = cuts[length(cuts)], labels = expression(infinity))
axis(2, at = pretty(0.85*den$y))
mtext('Probability density', side = 2, line = 2.5, at = mean(range(0.85*den$y)))
axis(2, at = par('usr')[4] - 0.1*diff(par('usr')[3:4]), labels = 'counts:', tick = FALSE, las = 2, hadj = 0.5)

text((cuts[-1] + cuts[-length(cuts)]) / 2, rep(par('usr')[4] - 0.1*diff(par('usr')[3:4]), length(cuts) - 1),
     labels = counts)

dev.off()


## DNA versus length
sizeDNA <- read.table('size_mass_dna.txt', sep = '\t', header=TRUE, as.is = TRUE)

pdf('fig_sizeDNA_log.pdf', width = 5, height = 5)
par(mar = c(3, 3, 0, 0) + 0.5, mgp = c(2.5, 1, 0))

plot(sizeDNA$mm_Length, sizeDNA$mgDNA, log = 'xy', axes = FALSE, 
     xlab = 'Length (mm)', ylab = 'Amount DNA (mg)', xlim = c(0.5, 15), 
     panel.first = {
         for(i in 1:(length(cuts) - 1)) {
             rect(xleft = ifelse(cuts[i] == 0, 10^par('usr')[1], cuts[i]), 
                  xright = cuts[i+1], 
                  ybottom = 10^par('usr')[3], ytop = 10^par('usr')[4],
                  col = c('gray75', 'gray95')[i %% 2 + 1], border = NA)
         }
         abline(v = cuts, col = 'white')
     })

logAxis(1)
logAxis(2)

dev.off()


pdf('fig_sizeDNA_linear.pdf', width = 5, height = 5)
par(mar = c(3, 3, 0, 0) + 0.5, mgp = c(2.5, 1, 0))

plot(sizeDNA$mm_Length, sizeDNA$mgDNA, axes = FALSE, 
     xlab = 'Length (mm)', ylab = 'Amount DNA (mg)', xlim = range(cuts), 
     panel.first = {
         for(i in 1:(length(cuts) - 1)) {
             rect(xleft = cuts[i], xright = cuts[i+1], 
                  ybottom = par('usr')[3], ytop = par('usr')[4],
                  col = c('gray75', 'gray95')[i %% 2 + 1], border = NA)
         }
         abline(v = cuts, col = 'white')
     })

axis(1)
axis(2)

dev.off()

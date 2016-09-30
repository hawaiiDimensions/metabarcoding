setwd('~/Dropbox/hawaiiDimensions/metabarcoding/size_classes')

cuts <- c(0, cumsum(round(1.2^(0:7))))

sizeDNA <- read.table('size_mass_dna.txt', sep = '\t', header=TRUE, as.is = TRUE)
head(sizeDNA)

plot(sizeDNA$mm_Length, sizeDNA$mgDNA, col = as.factor(sizeDNA$Family), 
     pch = as.numeric(as.factor(sizeDNA$Group)), log = 'xy')

abline(v = cuts)

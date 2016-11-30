setwd('~/Dropbox/hawaiiDimensions/metabarcoding/poc')

library(xlsx)
library(reshape2)

## load data
tissue <- read.xlsx('Tissue_Pools_092016_Andy.xlsx', sheetIndex = 1, stringsAsFactors = FALSE)
tissue <- tissue[, (1:ncol(tissue)) <= which(names(tissue) == 'Average_Prop_H_L')]

## correct issue of missing info when no reads obtained:
## for low molecular weight

tissue[is.na(tissue$Prop_Low), 'Sample.3'] <- paste('L', tissue$Sample[is.na(tissue$Prop_Low)], sep = '')
tissue[is.na(tissue$Prop_Low), 'Reads.1'] <- 0
tissue[is.na(tissue$Prop_Low), 'Prop_Low'] <- 0

tapply(tissue$TotalReads.1, tissue$Sample.3, max, na.rm = TRUE)
tapply(tissue$Reads.1, tissue$Sample.3, sum, na.rm = TRUE)


tapply(tissue$TotalReads, tissue$Sample.2, max, na.rm = TRUE)
tapply(tissue$Reads, tissue$Sample.2, sum, na.rm = TRUE)

tissue[is.na(tissue$Prop_Low), 'TotalReads.1'] <- 0

## for hi moleular weight

## combine low and hi molecular weight samples into tiddy format
tissue <- rbind(tissue[, c('Sample', 'Freezer_RT', 'Pool', 'Taxon', 'mg_in_Pool', 
                           'Total_mg_in_pool', 'Sample.2', 'Reads', 'TotalReads')], 
                tissue[, c('Sample', 'Freezer_RT', 'Pool', 'Taxon', 'mg_in_Pool', 
                           'Total_mg_in_pool', 'Sample.3', 'Reads.1', 'TotalReads.1')])






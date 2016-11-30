setwd('~/Dropbox/hawaiiDimensions/metabarcoding/poc')

library(xlsx)
library(reshape2)

## load data
tissue <- read.xlsx('Tissue_Pools_092016_Andy.xlsx', sheetIndex = 1, stringsAsFactors = FALSE)
tissue <- tissue[, (1:ncol(tissue)) <= which(names(tissue) == 'Average_Prop_H_L')]

## correct issue of missing info when no reads obtained:
## for low molecular weight
tissue[is.na(tissue$Prop_Low) | tissue$Prop_Low == 0, 'Sample.3'] <- 
    paste('L', tissue$Sample[is.na(tissue$Prop_Low) | tissue$Prop_Low == 0],
          sep = '')
tissue[is.na(tissue$Prop_Low) | tissue$Prop_Low == 0, 'Reads.1'] <- 0

## for hi moleular weight
tissue[is.na(tissue$Prop_High) | tissue$Prop_High == 0, 'Sample.2'] <- 
    paste('H', tissue$Sample[is.na(tissue$Prop_High) | tissue$Prop_High == 0],
          sep = '')
tissue[is.na(tissue$Prop_High) | tissue$Prop_High == 0, 'Reads'] <- 0

## correct total reads to be only those assigned:
tissue$TotalReads.1 <- tapply(tissue$Reads.1, tissue$Sample.3, sum, 
                              na.rm = TRUE)[tissue$Sample.3]
tissue$TotalReads <- tapply(tissue$Reads, tissue$Sample.2, sum, 
                              na.rm = TRUE)[tissue$Sample.2]


## combine low and hi molecular weight samples into tiddy format
tHi <- tissue[, c('Sample', 'Freezer_RT', 'Taxon', 'mg_in_Pool', 
                  'Total_mg_in_pool', 'Sample.2', 'Reads', 'TotalReads')]
tLo <- tissue[, c('Sample', 'Freezer_RT', 'Taxon', 'mg_in_Pool', 
                  'Total_mg_in_pool', 'Sample.3', 'Reads.1', 'TotalReads.1')]

names(tHi) <- names(tLo) <- c('Pool', 'Freezer_RT', 'Specimen', 'amount_DNA', 
                              'Total_mg_in_pool', 'Experiment', 'number_Reads', 'total_Reads')
tissue <- rbind(tHi, tLo)

## make experiment just refer to high v. low molec weight and freezer v. room tem
tissue$Experiment <- gsub('[0-9]', '', tissue$Experiment)
tissue$Experiment <- paste(substring(tissue$Experiment, 1, 1), 
                           substring(tissue$Experiment, 2, 100), 
                           sep = '_')

## NOTE: amount_DNA = mg of tissue, not extracted DNA
write.csv(tissue, 'clean_tissue.csv', row.names = FALSE)

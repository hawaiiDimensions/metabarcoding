setwd('~/Dropbox/hawaiiDimensions/metabarcoding/poc')

library(xlsx)
library(reshape2)

## load data
tissue <- read.xlsx2('Tissue_Pools_092016_Andy.xlsx', sheetIndex = 1, , stringsAsFactors = FALSE)
head(tissue)
tissue <- tissue[, (1:ncol(tissue)) <= which(names(tissue) == 'Average_Prop_H_L')]

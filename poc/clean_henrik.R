setwd('~/Dropbox/hawaiiDimensions/metabarcoding/poc')

library(xlsx)

diffMarkers <- read.xlsx2('Blast_Results_DifferentMarkersPCRcycles200616.xlsx', sheetName='Blast_Results_DifferentMarkers2')

## DNA content of input is %DNA*15*124.5
setwd('~/Dropbox/hawaiiDimensions/metabarcoding/poc')

library(xlsx)

## load data
diffMarkers <- read.xlsx2('Blast_Results_DifferentMarkersPCRcycles200616.xlsx', sheetName='Blast_Results_DifferentMarkers2', stringsAsFactors=FALSE)
pcrCycle <- read.xlsx2('Blast_Results_DifferentMarkersPCRcycles200616.xlsx', sheetName='Blast_Results_PCRcycleNumber200', stringsAsFactors=FALSE)
pools <- read.csv('DNA_Pools_Sequencing_JAE_052016.csv', as.is = TRUE)

## clean-up pools
colnames(pools) <- sub('\\.', '', sub('.*Pool', 'JAEP', colnames(pools)))
pools <- pools[, -which(colnames(pools) == 'Concetration')]
for(i in grep('JAEP', colnames(pools))) {
    pools[, i] <- pools[, i] * pools$Currentconcentration
}

pools <- acast(melt(pools, id.vars = c('Sample_ID'), measure.vars = grep('JAEP', colnames(pools))), variable ~ Sample_ID)


## clean data
cleanNames <- function(x) {
    out <- gsub('X.', 'percent', names(x))
    out[grep('Reads_with_Blast_hits', out, ignore.case = TRUE)] <- 'number_Reads'
    out[grep('Total_number_of_reads', out, ignore.case = TRUE)] <- 'total_Reads'
    out <- gsub('\\.', '_', out)
    
    return(out)
}

cleanup <- function(x) {
    names(x) <- cleanNames(x)
    
    x$PCR_cycles <- as.numeric(gsub('[^0-9]', '', x$PCR_cycles))
    x$DNA_Template <- as.numeric(gsub('[^0-9]', '', x$DNA_Template))
    x$Primer[grep('ARFHCO', x$Primer)] <- gsub('[0-9]', '', x$Primer[grep('ARFHCO', x$Primer)])
    
    for(i in c('percent_DNA', 'percent_Reads', 'number_Reads', 'total_Reads')) x[, i] <- as.numeric(x[, i])
    
    ## DNA content of input using pool sheet
    x$amount_DNA <- diag(pools[x$Pool, x$Specimen])
    
    return(x)
}

diffMarkers <- cleanup(diffMarkers)
pcrCycle <- cleanup(pcrCycle)

write.table(diffMarkers, file='clean_diffMarkers.csv', sep=',', row.names=FALSE)
write.table(pcrCycle, file='clean_pcrCycle.csv', sep=',', row.names=FALSE)

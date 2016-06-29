setwd('~/Dropbox/hawaiiDimensions/metabarcoding/poc')

library(xlsx)

## load data
diffMarkers <- read.xlsx2('Blast_Results_DifferentMarkersPCRcycles200616.xlsx', sheetName='Blast_Results_DifferentMarkers2', stringsAsFactors=FALSE)
pcrCycle <- read.xlsx2('Blast_Results_DifferentMarkersPCRcycles200616.xlsx', sheetName='Blast_Results_PCRcycleNumber200', stringsAsFactors=FALSE)

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
    
    for(i in c('percent_DNA', 'percent_Reads', 'number_Reads', 'total_Reads')) x[, i] <- as.numeric(x[, i])
    
    ## DNA content of input is %DNA*15*124.5
    x$amount_DNA <- x$percent_DNA * 15 * 124.5
    
    return(x)
}

diffMarkers <- cleanup(diffMarkers)
pcrCycle <- cleanup(pcrCycle)

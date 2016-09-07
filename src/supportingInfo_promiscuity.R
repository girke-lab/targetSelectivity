#!/usr/bin/env Rscript

# (C) 2016 Tyler William H Backman
# Purpose: produce promiscuity probability file for Supporting Information

library(R.utils)

# parse input options
selectivityCountsIndividualFile <- commandArgs(trailingOnly=TRUE)[1]
selectivityCountskClustFile <- commandArgs(trailingOnly=TRUE)[2]
selectivityCountsdomainsFile <- commandArgs(trailingOnly=TRUE)[3]
promiscuityProbabilityFile <- commandArgs(trailingOnly=TRUE)[4]
highlyScreenedCidsFile <- commandArgs(trailingOnly=TRUE)[5]
activeCidsFile <- commandArgs(trailingOnly=TRUE)[6]
outputFilename <- commandArgs(trailingOnly=TRUE)[7]

# test code for running without make:
if(is.na(commandArgs(trailingOnly=TRUE)[1])){
    selectivityCountsIndividualFile <- "working/selectivityCountsIndividual.txt"
    selectivityCountskClustFile <- "working/selectivityCountskClust.txt"
    selectivityCountsdomainsFile <- "working/selectivityCountsdomains.txt"
    promiscuityProbabilityFile <- "working/promiscuityProbability.tab"
    highlyScreenedCidsFile <- "working/highlyScreenedCids.txt"
    activeCidsFile <- "working/activeCids.txt"
    outputFilename <- "working/supportingInfo_promiscuity.tab"
}

# parse input files
selectivityCountsIndividual <- read.table(selectivityCountsIndividualFile)
selectivityCountskClust <- read.table(selectivityCountskClustFile)
selectivityCountsdomains <- read.table(selectivityCountsdomainsFile)
promiscuityProbability <- read.table(promiscuityProbabilityFile)
highlyScreenedCids <- read.table(highlyScreenedCidsFile)[[1]]
activeCids <- read.table(activeCidsFile)[[1]]

# merge everything
mergedTable <- merge(selectivityCountsIndividual, selectivityCountskClust, by.x=1, by.y=1)
mergedTable <- merge(mergedTable, selectivityCountsdomains, by.x=1, by.y=1)
mergedTable <- merge(mergedTable, promiscuityProbability, by.x=1, by.y=1)
mergedTable <- merge(selectivityCountsIndividual, selectivityCountskClust, by.x=1, by.y=1)
mergedTable <- merge(selectivityCountsIndividual, selectivityCountskClust, by.x=1, by.y=1)

# keep only highly screened actives
colnames(mergedTable) <- c("pubChemCid", "targetSelectivity", "clusterSelectivity", "domainSelectivity", "promiscuityProbability")
mergedTable <- mergedTable[mergedTable$pubChemCid %in% highlyScreenedCids,]
mergedTable <- mergedTable[mergedTable$pubChemCid %in% activeCids,]

# sort by decreasing promiscuity
mergedTable <- mergedTable[order(mergedTable$promiscuityProbability, decreasing=T),]

# write out table
write.table(mergedTable, outputFilename, quote=F, sep="\t", row.names=F)

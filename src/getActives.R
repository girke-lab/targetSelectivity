#!/usr/bin/env Rscript

# (C) 2016 Tyler William H Backman
# Purpose: create a text file of all active compounds

library(R.utils)
library(bioassayR)

bioassayDatabaseFile <- commandArgs(trailingOnly=TRUE)[1]
outputCidListFile <- commandArgs(trailingOnly=TRUE)[2]

# test code for running without make:
if(is.na(commandArgs(trailingOnly=TRUE)[1])){
    bioassayDatabaseFile <- "working/bioassayDatabase.sqlite"
    outputCidListFile <- "working/activeCids.txt"
}

database <- connectBioassayDB(bioassayDatabaseFile)
outputCidList <- allCids(database, activesOnly = TRUE)
write.table(outputCidList, file=outputCidListFile, row.names=FALSE, col.names=FALSE)

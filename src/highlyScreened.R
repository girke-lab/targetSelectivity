#!/usr/bin/env Rscript

# (C) 2015 Tyler William H Backman
# Purpose: create a text file of compounds screened at least 
#   minTargets times

library(R.utils)
library(bioassayR)

bioassayDatabaseFile <- commandArgs(trailingOnly=TRUE)[1]
outputCidListFile <- commandArgs(trailingOnly=TRUE)[2]
minTargets <- as.numeric(commandArgs(trailingOnly=TRUE)[3])

# test code for running without make:
if(is.na(commandArgs(trailingOnly=TRUE)[1])){
    bioassayDatabaseFile <- "working/bioassayDatabase.sqlite"
    outputCidListFile <- "working/highlyScreenedCids.txt"
    minTargets <- 10
}

database <- connectBioassayDB(bioassayDatabaseFile)
outputCidList <- screenedAtLeast(database, minTargets, inconclusives=FALSE) 
write.table(outputCidList, file=outputCidListFile, row.names=FALSE, col.names=FALSE)

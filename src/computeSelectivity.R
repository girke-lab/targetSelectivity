#!/usr/bin/env Rscript

# (C) 2016 Tyler William H Backman
# Purpose: compute target selectivity with desired options

library(R.utils)
library(bioassayR)
library(foreach)
library(doMC)

# parse input options
databaseFile <- commandArgs(trailingOnly=TRUE)[1]
highlyScreenedCidsFile <- commandArgs(trailingOnly=TRUE)[2]
outputFile <- commandArgs(trailingOnly=TRUE)[3]
cores <- commandArgs(trailingOnly=TRUE)[4]
category <- commandArgs(trailingOnly=TRUE)[5]
multiTarget <- commandArgs(trailingOnly=TRUE)[6]

# test code for running without make:
if(is.na(commandArgs(trailingOnly=TRUE)[1])){
    databaseFile <- "working/bioassayDatabase.sqlite"
    highlyScreenedCidsFile <- "working/highlyScreenedCids.txt"
    outputFile <- "working/selectivityCounts.txt"
    cores <- 1
    category <- FALSE
    multiTarget <- "keepOne"
}

# parse input files
cids <- read.table(highlyScreenedCidsFile)[[1]]
if(category == "FALSE")
   category <- FALSE

splitCids <- split(cids, as.integer((1:length(cids))/1000))
gc()
registerDoMC(cores=cores)
results <- foreach(i = splitCids, .combine="c", .packages=c("bioassayR")) %dopar% {
    tempDBC <- connectBioassayDB(databaseFile)
    selectivity <- targetSelectivity(tempDBC, i, scoring="total", category=category, multiTarget=multiTarget)
    disconnectBioassayDB(tempDBC)
    selectivity
}

write.table(results, outputFile, col.names=FALSE, row.names=TRUE)

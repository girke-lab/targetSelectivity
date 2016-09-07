#!/usr/bin/env Rscript

# (C) 2016 Tyler William H Backman
# Purpose: compute heavy atom (non hydrogen) count 

library(R.utils)
library(ChemmineR)
library(foreach)
library(doMC)

bioassayCompoundFolder <- commandArgs(trailingOnly=TRUE)[1]
outputFile <- commandArgs(trailingOnly=TRUE)[2]
cores <- commandArgs(trailingOnly=TRUE)[3]

# test code for running without make:
if(is.na(commandArgs(trailingOnly=TRUE)[1])){
    bioassayCompoundFolder <- "working/bioassayCompounds"
    outputFile <- "working/heavycount.txt"
    cores <- 1
}

# parse input files
compoundFileList <- list.files(bioassayCompoundFolder, pattern=".sdf$", full.names = TRUE)

gc()
registerDoMC(cores=cores)
results <- foreach(i = compoundFileList, .combine='c', .packages=c("ChemmineR")) %dopar% {
    compounds <- read.SDFset(i)
    atomCounts <- atomcount(compounds, addH=FALSE)
    moleculeSize <- sapply(atomCounts, function(x) sum(x[names(x) != "H"])) 
    names(moleculeSize) <- sdfid(compounds)
    moleculeSize
}
results <- results[! duplicated(names(results))]

# save output
write.table(results, outputFile, row.names=T, col.names=F)

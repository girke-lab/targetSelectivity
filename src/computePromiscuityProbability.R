#!/usr/bin/env Rscript

# (C) 2016 Tyler William H Backman
# Purpose: compute promiscuity probability for all compounds in a compound vs target matrix

library(R.utils)
library(Matrix)
library(bioassayR)
library(foreach)
library(doMC)

# parse input options
cidsVStargetFile <- commandArgs(trailingOnly=TRUE)[1]
outputFilename <- commandArgs(trailingOnly=TRUE)[2]
cores <- as.numeric(commandArgs(trailingOnly=TRUE)[3])

# test code for running without make:
if(is.na(commandArgs(trailingOnly=TRUE)[1])){
    cidsVStargetFile <- "working/cidsVStargets.RData"
    outputFilename <- "working/promiscuityProbability.tab"
    cores <- 31 
}

# parse input files
load(cidsVStargetFile)
compoundVsTargetMatrix <- results # rows = targets, cols = compounds
rm(results)

splitCols <- split(1:ncol(compoundVsTargetMatrix), 
                   as.integer((1:ncol(compoundVsTargetMatrix))/(ncol(compoundVsTargetMatrix)/cores + 1)))
gc()

registerDoMC(cores=cores)
promiscuityprob <- foreach(cols=splitCols, .combine='c') %dopar% {
   crossReactivityProbability(compoundVsTargetMatrix[,cols], prior=list(hit_ratio_mean=0.01861531, hit_ratio_sd=0.03494048))
}
names(promiscuityprob) <- colnames(compoundVsTargetMatrix)

write.table(promiscuityprob, file=outputFilename, row.names=T, col.names=F)

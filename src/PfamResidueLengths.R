#!/usr/bin/env Rscript

# (C) 2016 Tyler William H Backman
# Purpose: get table of HMM residue lengths

library(R.utils)

# parse input options
pfamstatsFile <- commandArgs(trailingOnly=TRUE)[1]
outputFile <- commandArgs(trailingOnly=TRUE)[2]

# test code for running without make:
if(is.na(commandArgs(trailingOnly=TRUE)[1])){
    pfamstatsFile <- "working/Pfam-A-stats.txt"
    outputFile <- "working/PfamResidueLengths.tab"
}

# parse input options
pfamstats <- read.table(pfamstatsFile)

# create table
domains <- gsub("^(PF\\d*).*", "\\1", pfamstats[,3], perl=TRUE)
domainTable <- cbind(domain=domains, residues=pfamstats[,6])

# save output
write.table(domainTable, outputFile, row.names=FALSE, quote=F)

#!/usr/bin/env Rscript

# (C) 2016 Tyler William H Backman
# Purpose: 
# fix unparseable PAINs files

library(R.utils)
library(ChemmineR)

# parse input options
unparsablePainsFile <- commandArgs(trailingOnly=TRUE)[1]
activeSDFFile <- commandArgs(trailingOnly=TRUE)[2]
outputFile <- commandArgs(trailingOnly=TRUE)[3]

# test code for running without make:
if(is.null(commandArgs(trailingOnly=TRUE)[1])){
    unparsablePainsFile <- "working/activeCompoundsPAINSu.txt"
    activeSDFFile <- "/dev/shm/activeCompounds.sdf"
    outputFile <- "working/activeCompoundsPAINSuFixed.txt"
}

# parse input files
unparsablePains <- as.integer(read.table(unparsablePainsFile)[[1]]) # positions in activeSDFFile starting at 1

# get cids of all of activeSDFFile to subset by unparsablePains
myid <- function(sdfset){
    ids <- sdfid(sdfset)
    names(ids) <- cid(sdfset)
    ids
}
desc <- function(sdfset) {
    cbind(SDFID=myid(sdfset))
}
tempLoc <- tempfile()
system.time(
sdfStream(input=activeSDFFile, output=tempLoc, append=FALSE, fct=desc, Nlines=10000, silent=T)
)
cids <- read.delim(tempLoc, row.names=1)$SDFID
unlink(tempLoc)

painsCids <- cids[unparsablePains]

write.table(painsCids, outputFile, row.names=F, col.names=F)

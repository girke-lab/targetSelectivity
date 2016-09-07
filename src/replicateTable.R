#!/usr/bin/env Rscript

# (C) 2015 Tyler William H Backman
# Purpose: build table of replicate count distribution

# Note: uses lots of ram (~10-20gb) to run table on sparse matrix with real
#   matrix intermediate. Could be made more efficient.

library(R.utils)
library(Matrix)

# parse input options
replicateMatrixFile <- commandArgs(trailingOnly=TRUE)[1]
outputFile <- commandArgs(trailingOnly=TRUE)[2]

# test code for running without make:
if(is.na(commandArgs(trailingOnly=TRUE)[1])){
    replicateMatrixFile <- "working/replicateMatrix.RData"
    outputFile <- "working/replicateTable.csv"
}

# parse input files
load(replicateMatrixFile) # loads replicateMatrix input

# build replicate table
replicateTable <- table(as.matrix(replicateMatrix))

# write out replicate table
write.csv(replicateTable, outputFile, row.names = F)

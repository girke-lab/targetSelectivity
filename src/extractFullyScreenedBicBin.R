#!/usr/bin/env Rscript

# (C) 2016 Tyler William H Backman
# Purpose: perform biclustering to identify largest fully specified submatrix

# Input Format:
#  drug (row) vs cluster (col) object is called "drugComparisonMatrix" with the following values:
#  0 = untested
#  1 = inactive
#  2 = active

library(R.utils)
library(biclust)
library(Matrix)
library(BicBin)

# parse input options
inputMatrixFile <- commandArgs(trailingOnly=TRUE)[1]
outputFile <- commandArgs(trailingOnly=TRUE)[2]
cores <- as.numeric(commandArgs(trailingOnly=TRUE)[3]) 

# test code for running without make:
if(is.na(commandArgs(trailingOnly=TRUE)[1])){
    inputMatrixFile <- "working/cidsVStargetClusters.RData"
    outputFile <- "working/fullyScreened.RData"
    cores <- 1
}

# parse input files
load(inputMatrixFile) # loads "results"
drugTargetSparseMatrix <- results

# create non-sparse matrix
drugTargetMatrix <- as.matrix(drugTargetSparseMatrix)

# convert drugComparisonMatrix to a binary matrix
binaryMatrix <- binarize(drugTargetMatrix, threshold=0) 

# get all compounds screened against over 289 targets
mostScreened <- binaryMatrix[,colSums(binaryMatrix) > 289]
tempDTMatrix <- drugTargetMatrix[,colSums(binaryMatrix) > 289]

# get 100k most highly screened compounds
# sortedCidIndex <- sort(colSums(binaryMatrix), decreasing=T, index.return=T)$ix
# mostScreened <- binaryMatrix[,sortedCidIndex[1:100000]]
# tempDTMatrix <- drugTargetMatrix[,sortedCidIndex[1:100000]]

# compute biclusters
gc()
registerDoMC(cores=cores)
set.seed(123) # get same clusters each time
tries <- cores*1
system.time(
bicluster <- bicBinCluster(binaryMatrix=mostScreened, totalClusters=1, alpha=0.5, beta=0.8, tries=tries, minDensityP1=1)
)

results <- tempDTMatrix[bicluster[[1]]$rows, bicluster[[1]]$cols]
# table(results)
# dim(results)

# remove all-inactive rows and columns
results <- results[,colSums(binarize(results, threshold=1)) > 0]
results <- results[rowSums(binarize(results, threshold=1)) > 0,]

# save results
save(list=c("results"), file=outputFile)

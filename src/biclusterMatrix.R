#!/usr/bin/env Rscript

# (C) 2015 Tyler William H Backman
# Purpose: perform biclustering to identify target/drug groups

# Input Format:
#  drug (row) vs cluster (col) object is called "drugComparisonMatrix" with the following values:
#  0 = untested
#  1 = inactive
#  2 = active
#  3 = untested and annotated
#  4 = inactive and annotated
#  5 = active and annotated 

library(R.utils)
library(biclust)
library(foreach)
library(doMC)
library(doRNG)

# parse input options
drugComparisonMatrixFile <- commandArgs(trailingOnly=TRUE)[1]
bicBinSource <- commandArgs(trailingOnly=TRUE)[2]
outputFile <- commandArgs(trailingOnly=TRUE)[3]
cores <- commandArgs(trailingOnly=TRUE)[4]

# test code for running without make:
if(is.na(commandArgs(trailingOnly=TRUE)[1])){
    drugComparisonMatrixFile <- "working/drugComparisonMatrix.RData"
    bicBinSource <- "src/bicbin.R"
    outputFile <- "working/biClusters.RData_test"
    cores <- 8
}

# parse input files
load(drugComparisonMatrixFile)
source(bicBinSource)

# convert drugComparisonMatrix to a binary matrix
binaryMatrix <- binarize(drugComparisonMatrix, threshold=1) 

# test code: keep only most active rows and cols
# binaryMatrix <- binaryMatrix[,colSums(binaryMatrix) > 30]
# binaryMatrix <- binaryMatrix[rowSums(binaryMatrix) > 10,]

# perform bicbin clustering in parallel
bicBinCluster <- function(binaryMatrix, totalClusters, cores, alpha=0.5, beta=0.5, tries=900){
    M <- nrow(binaryMatrix)
    N <- ncol(binaryMatrix)
    p <- sum(binaryMatrix) / (M*N)
    clusters <- list()
    for(i in 1:totalClusters){
        allResults <- foreach(thisTry = 1:tries) %dorng% {
              res1 <- BicBin(binaryMatrix, alpha, beta, p, proc_genes=TRUE) # by columns
              res2 <- BicBin(binaryMatrix, alpha, beta, p, proc_genes=FALSE) # by rows
              if(res1$score > res2$score)
                return(res1)
              else
                return(res2)
        }
        max <- allResults[[which.max(sapply(allResults, function(x) x$score))]]
        binaryMatrix[max$x == 1, max$y == 1] <- 0 # zero out this clusters region for next round
        clusters[[i]] <- list(rows=row.names(binaryMatrix)[max$x == 1], cols=colnames(binaryMatrix)[max$y == 1], score=max$score)
    }
    return(clusters)
}
gc()
registerDoMC(cores=cores)
set.seed(123) # get same clusters each time
clusterResults <- bicBinCluster(binaryMatrix, totalClusters=30, alpha=0.6, beta=0.6, cores=cores, tries=5000)

# save results
save(list=c("clusterResults"), file=outputFile)

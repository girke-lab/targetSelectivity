#!/usr/bin/env Rscript

# (C) 2016 Tyler William H Backman
# Purpose: produce list of novel drug-target pairs

# Input Format:
#  drug (row) vs cluster (col) object is called "drugComparisonMatrix" with the following values:
#  0 = untested
#  1 = inactive
#  2 = active
#  3 = untested and annotated
#  4 = inactive and annotated
#  5 = active and annotated 

library(R.utils)
library(Matrix)

# parse input options
drugComparisonMatrixFile <- commandArgs(trailingOnly=TRUE)[1]
clusterAnnotationFile <- commandArgs(trailingOnly=TRUE)[2]
outputFile <- commandArgs(trailingOnly=TRUE)[3]

# test code for running without make:
if(is.na(commandArgs(trailingOnly=TRUE)[1])){
    drugComparisonMatrixFile <- "working/drugComparisonMatrix.RData"
    clusterAnnotationFile <- "working/clusterAnnotations.csv"
    outputFile <- "working/supportingInfo_novelTargets.xls"
}

# parse input files
load(drugComparisonMatrixFile) # drugComparisonMatrix
clusterAnnotations <- read.csv(clusterAnnotationFile)

tMatrix <- as(drugComparisonMatrix, "TsparseMatrix")
i <- tMatrix@i
j <- tMatrix@j
x <- tMatrix@x

# get cids and targets for active (#2) positions
activeIndex <- x == 2
cids <- row.names(drugComparisonMatrix)[i + 1][activeIndex]
targetClusters <- colnames(drugComparisonMatrix)[j + 1][activeIndex]

# replace cluster IDs with UniProt IDs
labels <- as.character(clusterAnnotations$accession)
labels[is.na(labels)] <- as.character(clusterAnnotations$protein_name[is.na(labels)])
labels <- gsub("^UniProt_(.*)", "\\1", labels)
targetIds <- labels[match(targetClusters, clusterAnnotations$uniqueClusterIds)]

# write out table
outputTable <- cbind(PubChemCids=cids, UniProtIds=targetIds)
write.table(outputTable, outputFile, quote=F, sep="\t", row.names=F, col.names = T)

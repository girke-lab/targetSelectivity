#!/usr/bin/env Rscript

# (C) 2016 Tyler William H Backman
# Purpose: produce fully screened compound vs target cluster binary matrix

library(R.utils)
library(biclust)

# parse input options
inputMatrixFile <- commandArgs(trailingOnly=TRUE)[1]
clusterAnnotationFile <- commandArgs(trailingOnly=TRUE)[2]
outputFilename <- commandArgs(trailingOnly=TRUE)[3]

# test code for running without make:
if(is.na(commandArgs(trailingOnly=TRUE)[1])){
    inputMatrixFile <- "working/fullyScreened.RData"
    clusterAnnotationFile <- "working/clusterAnnotations.csv"
    outputFilename <- "working/supportingInfo_fullyScreened.tab"
}

# parse input files
load(inputMatrixFile) # loads "results"
compoundTargetMatrix <- t(results)
clusterAnnotations <- read.csv(clusterAnnotationFile)

# make binary
compoundTargetMatrix <- binarize(compoundTargetMatrix, threshold=1)

# replace cluster IDs with UniProt IDs
labels <- as.character(clusterAnnotations$accession)
labels[is.na(labels)] <- as.character(clusterAnnotations$protein_name[is.na(labels)])
labels <- gsub("^UniProt_(.*)", "\\1", labels)
colnames(compoundTargetMatrix) <- labels[match(colnames(compoundTargetMatrix), clusterAnnotations$uniqueClusterIds)]
compoundTargetMatrix <- as.data.frame(compoundTargetMatrix, stringsAsFactors = F)

# write out table
write.table(compoundTargetMatrix, outputFilename, quote=F, sep="\t", row.names=T, col.names = T)

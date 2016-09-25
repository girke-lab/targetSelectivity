#!/usr/bin/env Rscript

# (C) 2016 Tyler William H Backman
# Purpose: produce table enumerating novel drug-target pairs

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
library(xtable)

# parse input options
drugComparisonMatrixFile <- commandArgs(trailingOnly=TRUE)[1]
clusterAnnotationFile <- commandArgs(trailingOnly=TRUE)[2]
outputFile <- commandArgs(trailingOnly=TRUE)[3]

# test code for running without make:
if(is.na(commandArgs(trailingOnly=TRUE)[1])){
    drugComparisonMatrixFile <- "working/drugComparisonMatrix.RData"
    clusterAnnotationFile <- "working/clusterAnnotations.csv"
    outputFile <- "working/drugComparisonTable.tex"
}

# parse input files
load(drugComparisonMatrixFile) # drugComparisonMatrix
clusterAnnotations <- read.csv(clusterAnnotationFile)

# tabulate counts
annotationFrequencies <- table(drugComparisonMatrix)
annotationFrequencies <- as.data.frame(annotationFrequencies)
annotationFrequencies <- cbind(annotationFrequencies[1:3,2], annotationFrequencies[4:6,2])
row.names(annotationFrequencies) <- c("Untested", "Inactive", "Active")
colnames(annotationFrequencies) <- c("Unannotated Targets", "Annotated Targets")
comparisonTable <- xtable(annotationFrequencies, caption="Frequency of bioassay activity values compared to DrugBank annotation")

# make human target only matrix
humanClusterIds <- clusterAnnotations$uniqueClusterIds[grepl("Homo sapiens", clusterAnnotations$organism)]

humanOnlyMatrix <- drugComparisonMatrix[, colnames(drugComparisonMatrix) %in% humanClusterIds]
annotationFrequencies <- table(humanOnlyMatrix)
annotationFrequencies <- as.data.frame(annotationFrequencies)
annotationFrequencies <- cbind(annotationFrequencies[1:3,2], annotationFrequencies[4:6,2])
row.names(annotationFrequencies) <- c("Untested", "Inactive", "Active")
colnames(annotationFrequencies) <- c("Unannotated Targets", "Annotated Targets")
humanComparisonTable <- xtable(annotationFrequencies, caption="Frequency of bioassay activity values compared to DrugBank annotation")

print(comparisonTable, type="latex", file=outputFile, include.rownames=F)
print(humanComparisonTable, type="latex", file=outputFile, include.rownames=F, append=T)

# count number of annotated target columns
sum(colSums(1*(drugComparisonMatrix > 2)) > 0)
sum(colSums(1*(humanOnlyMatrix > 2)) > 0)

# count number of columns with actives but no annotated targets
sum((1*colSums(1*(drugComparisonMatrix == 2)) > 0) * (1*colSums(1*(drugComparisonMatrix > 2)) == 0))
sum((1*colSums(1*(humanOnlyMatrix == 2)) > 0) * (1*colSums(1*(humanOnlyMatrix > 2)) == 0))


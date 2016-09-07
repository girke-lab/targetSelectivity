#!/usr/bin/env Rscript

# (C) 2015 Tyler William H Backman
# Purpose: 
# find GO terms for all annotated protein targets

library(R.utils)
library(biomaRt)

# parse input options
clusterAnnotationsFile <- commandArgs(trailingOnly=TRUE)[1]
outputFile <- commandArgs(trailingOnly=TRUE)[2]

# test code for running without make:
if(is.na(commandArgs(trailingOnly=TRUE)[1])){
    setwd("../")
    clusterAnnotationsFile <- "working/clusterAnnotations.csv"
    outputFile <- "working/clusterGOannotations.csv"
}

# parse input options
inputTargets <- read.csv(clusterAnnotationsFile) 
inputTargets <- as.character(unique(inputTargets$accession[! is.na(inputTargets$accession)]))

# get GO identifiers
unimart <- useMart("unimart")
uniprot <- useDataset("uniprot", mart=unimart)
GOannotations <- getBM(attributes=c("accession", "go_id", "go_name"), filters=c("accession"), mart=uniprot,values=inputTargets)

write.csv(GOannotations, outputFile, row.names=FALSE)

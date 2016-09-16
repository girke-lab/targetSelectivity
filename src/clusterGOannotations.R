#!/usr/bin/env Rscript

# (C) 2015 Tyler William H Backman
# Purpose: 
# find GO terms for all annotated protein targets

library(R.utils)
library(RCurl)
# library(biomaRt)
source("src/uniprotInterface.R")

# parse input options
clusterAnnotationsFile <- commandArgs(trailingOnly=TRUE)[1]
outputFile <- commandArgs(trailingOnly=TRUE)[2]

# test code for running without make:
if(is.na(commandArgs(trailingOnly=TRUE)[1])){
    clusterAnnotationsFile <- "working/clusterAnnotations.csv"
    outputFile <- "working/clusterGOannotations.csv"
}

# parse input options
inputTargets <- read.csv(clusterAnnotationsFile) 
inputTargets <- as.character(unique(inputTargets$accession[! is.na(inputTargets$accession)]))

# get GO identifiers
# unimart <- useMart("unimart")
# uniprot <- useDataset("uniprot", mart=unimart)
# GOannotations <- getBM(attributes=c("accession", "go_id", "go_name"), filters=c("accession"), mart=uniprot,values=inputTargets)

GOannotations <- UniProtAnnotator(inputTargets, columns=c("genes", "organism", "protein names", "go(molecular function)"))

extractTerms <- function(x){
    print(paste("getting terms for UniProt ID:", x[["accession"]], collapse=" "))
    GOresult <- x[["Gene ontology (molecular function)"]]
    splitTerms <- strsplit(GOresult, split=";\\W")[[1]]
    go_id <- gsub("^.*\\[(.*)\\]$", "\\1", splitTerms)
    go_name <- gsub("^(.*)\\W\\[.*\\]$", "\\1", splitTerms)
    return(data.frame(cbind(accession=x[["accession"]], go_id, go_name), stringsAsFactors = F))
}

results <- apply(GOannotations, MARGIN=1, FUN=extractTerms)
results <- do.call("rbind", results)
results$go_name <- paste("F:", results$go_name, sep="")

write.csv(results, outputFile, row.names=FALSE)

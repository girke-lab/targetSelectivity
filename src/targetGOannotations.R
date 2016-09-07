#!/usr/bin/env Rscript

# (C) 2016 Tyler William H Backman
# Purpose: build table of GI->GO term mappings 

library(R.utils)
library(bioassayR)

# parse input options
uniprotGOannotationsFile <- commandArgs(trailingOnly=TRUE)[1]
pubchemDatabaseFile <- commandArgs(trailingOnly=TRUE)[2]
outputFile <- commandArgs(trailingOnly=TRUE)[3]

# test code for running without make:
if(is.na(commandArgs(trailingOnly=TRUE)[1])){
    uniprotGOannotationsFile <- "working/uniprotGOannotations.csv"
    pubchemDatabaseFile <- "working/bioassayDatabase.sqlite"
    outputFile <- "working/targetGOannotations.csv"
}

# parse input options
uniprotGOannotations <- read.csv(uniprotGOannotationsFile, header=FALSE, stringsAsFactors=FALSE)
colnames(uniprotGOannotations) <- c("UniProt", "GO")
pubchemDatabase <- connectBioassayDB(pubchemDatabaseFile)

# get go terms for each target GI in the database by mapping to the first UniProt ID
giList <- allTargets(pubchemDatabase)
GOannotations <- do.call(rbind, lapply(giList, function(gi){
    translatedId <- translateTargetId(pubchemDatabase, gi, "UniProt")[1]
    if(is.na(translatedId)) return(data.frame())
    goTerms <- uniprotGOannotations$GO[uniprotGOannotations$UniProt == translatedId]
    if(length(goTerms) == 0) return(data.frame())
    cbind(gi, goTerms)
}))
colnames(GOannotations) <- c("accession", "go_id")

# note: at time of writing there were terms here for 5968 out of 6480 GI targets

write.csv(GOannotations, outputFile, row.names=FALSE)

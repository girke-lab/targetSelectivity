#!/usr/bin/env Rscript

# (C) 2015 Tyler William H Backman
# Purpose: 
# collapse protein target GO annotations to GO Slim

library(R.utils)
library(GSEABase)
library(plyr)

# parse input options
inputFile <- commandArgs(trailingOnly=TRUE)[1]
goslimFile <- commandArgs(trailingOnly=TRUE)[2]
outputFile <- commandArgs(trailingOnly=TRUE)[3]

# test code for running without make:
if(is.na(commandArgs(trailingOnly=TRUE)[1])){
    inputFile <- "working/clusterGOannotations.csv"
    goslimFile <- "working/goslim_generic.obo"
    outputFile <- "working/clusterGOslimAnnotations.csv"
}

# parse input options
goInput <- read.csv(inputFile)
myCollection <- GOCollection(unique(as.vector(goInput$go_id)))
slim <- getOBOCollection(goslimFile)

# break apart go terms by protein cluster
goClusters <- split(goInput$go_id, goInput$accession)

# get goSlim terms for each protein cluster
# note: try is used because single terms missing from the goslim
#   return an error, and should be skipped
slimClusters <- lapply(goClusters, function(x){
    result <- tryCatch({
        goSlimResults <- goSlim(GOCollection(as.vector(x)), slim, "MF")
        goSlimResults <- goSlimResults[goSlimResults$Count > 0,,drop=FALSE]
        return(cbind(row.names(goSlimResults), as.vector(goSlimResults$Term)))
    }, error = function(e){ return(NA) })
})
slimClusters <- slimClusters[! is.na(slimClusters)]

# reassemble matrix
results <- ldply(slimClusters)
colnames(results) <- colnames(goInput)

# save output
write.csv(results, outputFile, row.names=FALSE)

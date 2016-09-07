#!/usr/bin/env Rscript

# (C) 2015 Tyler William H Backman
# parse kclust output and fix identifier names

library(R.utils)
library(bioassayR)

clusteringResultFolder <- commandArgs(trailingOnly=TRUE)[1]
outputFile <- commandArgs(trailingOnly=TRUE)[2]

# test code for running without make:
if(is.na(commandArgs(trailingOnly=TRUE)[1])){
    clusteringResultFolder <- "working/combinedCluster"
    outputFile <- "working/curatedClusters.txt"
}

# import sequence clusters and reformat
clusters <- read.table(file.path(clusteringResultFolder, "clusters.dmp"), skip=1)
clusterHeaders <- readLines(file.path(clusteringResultFolder, "headers.dmp"))
clusterHeaders <- gsub("^\\d+\\s+(.*)$", "\\1", clusterHeaders, perl = TRUE)
clusterHeaders <- gsub("^>gi\\|(\\d+).*$", "gi_\\1", clusterHeaders, perl = TRUE)
clusterHeaders <- gsub("^>drugbank_target\\|(\\w+)\\s+.*$", "UniProt_\\1", clusterHeaders, perl = TRUE)
clusters <- cbind2(clusterHeaders, clusters[,2])

# export all clusters
write.table(clusters, outputFile, quote = FALSE, row.names = FALSE, col.names = FALSE)

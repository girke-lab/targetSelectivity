#!/usr/bin/env Rscript

# (C) 2016 Tyler William H Backman
# Purpose: produce bicluster lists for supplement

library(R.utils)
library(foreach)
library(Matrix)

# parse input options
biClustersFile <- commandArgs(trailingOnly=TRUE)[1]
clusterAnnotationFile <- commandArgs(trailingOnly=TRUE)[2]
outputFilename1 <- commandArgs(trailingOnly=TRUE)[3]
outputFilename2 <- commandArgs(trailingOnly=TRUE)[4]

# test code for running without make:
if(is.na(commandArgs(trailingOnly=TRUE)[1])){
    biClustersFile <- "working/biClusters.RData"
    clusterAnnotationFile <- "working/clusterAnnotations.csv"
    outputFilename1 <- "working/supportingInfo_biclusterCids.xls"
    outputFilename2 <- "working/supportingInfo_biclusterTargets.xls"
}

# parse input files
load(biClustersFile) # loads clusterResults
clusterAnnotations <- read.csv(clusterAnnotationFile)

# reformat clusterResults into two cols: cid, clusterid
clustersPerCid <- foreach(thisCluster=1:length(clusterResults), .combine='rbind2') %do% {
    cbind(clusterResults[[thisCluster]]$rows, thisCluster)
}
colnames(clustersPerCid) <- c("PubChemCid", "Bicluster")

targetsPerCluster <- foreach(thisCluster=1:length(clusterResults), .combine='rbind2') %do% {
    cbind(clusterResults[[thisCluster]]$cols, thisCluster)
}
colnames(targetsPerCluster) <- c("UniProtID", "Bicluster")

# remove small biclusters
clusterCounts <- table(clustersPerCid[,2])
largeClusters <- names(clusterCounts)[clusterCounts > 1]
clustersPerCid <- clustersPerCid[clustersPerCid[,2] %in% largeClusters,]
targetsPerCluster <- targetsPerCluster[targetsPerCluster[,2] %in% largeClusters,]

# renumber biclusters
newNumbers <- 1:length(unique(clustersPerCid[,2]))
names(newNumbers) <- unique(clustersPerCid[,2])
clustersPerCid[,2] <- as.character(newNumbers[match(clustersPerCid[,2], names(newNumbers))])
targetsPerCluster[,2] <- as.character(newNumbers[match(targetsPerCluster[,2], names(newNumbers))])

# replace cluster IDs with UniProt IDs
labels <- as.character(clusterAnnotations$accession)
labels[is.na(labels)] <- as.character(clusterAnnotations$protein_name[is.na(labels)])
labels <- gsub("^UniProt_(.*)", "\\1", labels)
targetsPerCluster[,1] <- labels[match(targetsPerCluster[,1], clusterAnnotations$uniqueClusterIds)]

# write out table
write.table(clustersPerCid, outputFilename1, quote=F, sep="\t", row.names=F, col.names = T)
write.table(targetsPerCluster, outputFilename2, quote=F, sep="\t", row.names=F, col.names = T)

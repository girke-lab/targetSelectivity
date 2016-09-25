#!/usr/bin/env Rscript

# (C) 2016 Tyler William H Backman
# Purpose: make bipartite graph from output of src/biclusterMatrix.R

library(R.utils)
library(biclust)
library(igraph)
library(reshape2)

# parse input options
drugComparisonMatrixFile <- commandArgs(trailingOnly=TRUE)[1]
biClustersFile <- commandArgs(trailingOnly=TRUE)[2]
outputFile <- commandArgs(trailingOnly=TRUE)[3]

# test code for running without make:
if(is.null(commandArgs(trailingOnly=TRUE)[1])){
    drugComparisonMatrixFile <- "working/drugComparisonMatrix.RData"
    biClustersFile <- "working/biClusters.RData"
    clusterGOslimAnnotationsFile <- "working/clusterGOslimAnnotations.csv"
    clusterAnnotationsFile <- "working/clusterAnnotations.csv"
    biclustergoidsFile <- "working/biclustergoids.RData"
    outputFilename <- "working/drugComparisonGraph.gml"
}

# parse input files
load(drugComparisonMatrixFile) # loads drugComparisonMatrix
load(biClustersFile) # loads clusterResults
clusterGOslimAnnotations <- read.csv(clusterGOslimAnnotationsFile)
clusterAnnotations <- read.csv(clusterAnnotationsFile)
load(biclustergoidsFile) # loads goids

# binarize graph
binaryGraph <- binarize(drugComparisonMatrix, threshold=1)
colnames(binaryGraph) <- paste("t", colnames(binaryGraph), sep="")

# make edge list
edgeList <- melt(binaryGraph)
edgeList <- edgeList[edgeList$value == 1,c(1,2)]

bigraph <- graph.data.frame(edgeList, directed=F)

# label protein targets
vertexNames <- get.vertex.attribute(bigraph, "name")
vertexType <- rep(3, length(vertexNames))
vertexType[grepl("^t", vertexNames)] <- 4
bigraph <- set.vertex.attribute(bigraph, name="Poly", value=as.character(vertexType))

# label biclusters
biclusters <- lapply(clusterResults, function(x){
    x$rows
})
biclusters <- melt(biclusters)
vertexCluster <- rep("none", length(vertexNames))
vertexCluster[match(biclusters[,1], vertexNames)] <- biclusters[,2]
vertexCluster[grepl("^t", vertexNames)] <- "target"
bigraph <- set.vertex.attribute(bigraph, name="cluster", value=vertexCluster)

# add GO terms for targets
targetClusters <- vertexNames[grepl("^t", vertexNames)]
targetClusters <- gsub("^t", "", targetClusters)
mergedClusters <- merge(targetClusters, clusterAnnotations[,1:2], by.x=1, by.y=1, all.x=T, all.y=F)
clusterGOslimAnnotations <- clusterGOslimAnnotations[clusterGOslimAnnotations$go_id %in% goids,]
# clusterGOslimAnnotations <- clusterGOslimAnnotations[clusterGOslimAnnotations$go_name != "molecular_function",]
mergedClusters <- merge(mergedClusters, clusterGOslimAnnotations, by.x=2, by.y=1, all.x=T, all.y=F)
keepTerms <- rev(names(sort(table(clusterGOslimAnnotations$go_name), decreasing=T)))
mergedClusters <- mergedClusters[mergedClusters$go_name %in% keepTerms,]
mergedClusters <- mergedClusters[sort(match(mergedClusters$go_name, keepTerms), index.return=T)$ix,]
mergedClusters <- mergedClusters[! duplicated(mergedClusters$x),]
vertexGOs <- rep("none", length(vertexNames))
vertexGOs[grepl("^t", vertexNames)] <- as.character(mergedClusters$go_name[match(targetClusters, mergedClusters$x)])
bigraph <- set.vertex.attribute(bigraph, name="targetgo", value=vertexGOs)

write_graph(bigraph, outputFilename, format="gml")

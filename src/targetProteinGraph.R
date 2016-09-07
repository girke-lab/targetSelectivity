#!/usr/bin/env Rscript

# (C) 2016 Tyler William H Backman
# Purpose: create target-target interaction graph

library(R.utils)
library(igraph)
library(Matrix)
library(bioassayR)
library(foreach)
library(doMC)
library(biclust)

# parse input options
cidsVStargetFile <- commandArgs(trailingOnly=TRUE)[1]
promiscuityProbabilityFile <- commandArgs(trailingOnly=TRUE)[2]
databaseFile <- commandArgs(trailingOnly=TRUE)[3]
outputFilename <- commandArgs(trailingOnly=TRUE)[4]
cores <- commandArgs(trailingOnly=TRUE)[5]

# test code for running without make:
if(is.na(commandArgs(trailingOnly=TRUE)[1])){
    cidsVStargetFile <- "working/cidsVStargets.RData"
    promiscuityProbabilityFile <- "working/promiscuityProbability.tab"
    databaseFile <- "working/bioassayDatabase.sqlite"
    goannotationsFile <- "working/targetGOslimAnnotations.csv"
    outputFilename <- "working/tpnet.gml"
    cores <- 33
}

# parse input files
load(cidsVStargetFile)
compoundVsTargetMatrix <- results # rows = targets, cols = compounds
rm(results)
promiscuityProbability <- read.table(promiscuityProbabilityFile)
goannotations <- read.csv(goannotationsFile)

# subset matrix and transpose
# compoundVsTargetMatrix <- compoundVsTargetMatrix[1:100, ] # test code
compoundVsTargetMatrix <- t(compoundVsTargetMatrix) # cols = targets, rows = compounds

# remove promiscuous binders
promiscuousCids <- promiscuityProbability[promiscuityProbability[,2] > 0.5, 1]
compoundVsTargetMatrix <- compoundVsTargetMatrix[! row.names(compoundVsTargetMatrix) %in% promiscuousCids,]

# use trinarySimilarity to make target vs target distance matrix
gc()
registerDoMC(cores=cores)
distanceMat <- foreach(col=1:ncol(compoundVsTargetMatrix), .combine='rbind') %dopar% {
    trinarySimilarity(compoundVsTargetMatrix[,col,drop=F], compoundVsTargetMatrix)
}
row.names(distanceMat) <- colnames(distanceMat)

# save computation results
save(list = c("distanceMat"), file = gsub("^(.*).gml$", "\\1.RData", outputFilename))
# load(gsub("^(.*).gml$", "\\1.RData", outputFilename))

# inject 0s at NAs
distanceMat[is.na(distanceMat)] <- 0

# convert to UniProt ids
# database <- connectBioassayDB(databaseFile)
# uniProtIds <- sapply(row.names(distanceMat), function(x) translateTargetId(x, database=database, category="UniProt")[[1]])
# uniProtIds[is.na(uniProtIds)] <- row.names(distanceMat)[is.na(uniProtIds)]
# row.names(distanceMat) <- uniProtIds
# colnames(distanceMat) <- uniProtIds

# binarize with cutoff
binaryDist <- binarize(distanceMat, threshold=0.50)

# remove unconnected nodes and set diagonals to 1
diag(binaryDist) <- 0
zeroNodes <- colSums(binaryDist) == 0
binaryDist <- binaryDist[! zeroNodes, ! zeroNodes]
distanceMat <- distanceMat[! zeroNodes, ! zeroNodes]
diag(binaryDist) <- 1
diag(distanceMat) <- 1

adjgraph <- graph_from_adjacency_matrix(binaryDist, mode="undirected", diag=FALSE)
weightedadjgraph <- graph_from_adjacency_matrix(distanceMat, mode="undirected", diag=FALSE, weighted=TRUE)

# sort annotations so that only the 13 most common GO terms are represented
# and of those, duplicate annotations are resolved by keeping the rarer one
vertexNames <- get.vertex.attribute(adjgraph, "name")
goannotations <- goannotations[goannotations[,3] != "molecular_function",]
keepTerms <- rev(names(sort(table(goannotations[,2]), decreasing=T)[1:11]))
goannotations <- goannotations[goannotations[,1] %in% vertexNames,]
goannotations <- goannotations[goannotations[,2] %in% keepTerms,]
goannotations <- goannotations[sort(match(goannotations[,2], keepTerms), index.return=T)$ix,]

# annotate graph with GO terms
# only works without uniprot conversion!
goannotations <- goannotations[! duplicated(goannotations[,1]),]
vertexGO <- as.character(goannotations[match(vertexNames, goannotations[,1]),3])
adjgraph <- set.vertex.attribute(adjgraph, name="GOslim", value=vertexGO)

write_graph(adjgraph, outputFilename, format="gml")
write_graph(weightedadjgraph, gsub("^(.*).gml$", "\\1_weighted.gml", outputFilename), format="gml")

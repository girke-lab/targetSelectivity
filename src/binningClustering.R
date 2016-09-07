#!/usr/bin/env Rscript

# (C) 2015 Tyler William H Backman
# Purpose: cluster existing FDA approved drugs
#   using single linkage binning clustering

library(R.utils)
library(ChemmineR)

# Imports the cindex() function for computing the Jaccard Index or the Rand Index.
# source("src/clusterIndex.R")

drugbankFile <- commandArgs(trailingOnly=TRUE)[1]
drugbankLinks <- commandArgs(trailingOnly=TRUE)[2]
outputFile <- commandArgs(trailingOnly=TRUE)[3]

# test code for running without make:
if(is.na(commandArgs(trailingOnly=TRUE)[1])){
    drugbankFile <- "working/drugbank.sdf"
    drugbankLinks <- "working/drugbank_links.csv"
    outputFile <- "working/structureClusterDrugs.tab"
}

# parse input files
drugbankCompounds <- read.SDFset(drugbankFile)
drugbankCompounds <- drugbankCompounds[validSDF(drugbankCompounds)]
drugBankTable <- read.csv(drugbankLinks)

# fix cids in drugbank compounds
dbIdsWithCids <- drugBankTable[! is.na(drugBankTable[,c('PubChem.Compound.ID')]),c('DrugBank.ID')]
drugbankCompounds <- drugbankCompounds[datablocktag(drugbankCompounds, "DRUGBANK_ID") %in% dbIdsWithCids]
cidPositions <- match(datablocktag(drugbankCompounds, "DRUGBANK_ID"), drugBankTable[,c('DrugBank.ID')])
cid(drugbankCompounds) <- as.character(drugBankTable[cidPositions,c('PubChem.Compound.ID')])
drugbankCompounds <- drugbankCompounds[! duplicated(cid(drugbankCompounds))]

# compute atom pairs
apset <- sdf2ap(drugbankCompounds)

# exclude any which didn't return atom pairs
apset <- apset[! cid(apset) %in% names(which(sapply(as(apset, "list"), length)==1))] 

# perform hierarchical atom pair clustering
myTempFile <- tempfile()
dummy <- cmp.cluster(db=apset, cutoff=0, save.distances=myTempFile)
load(myTempFile)
unlink(myTempFile)
hc <- hclust(as.dist(distmat), method="average")
hc[["labels"]] <- cid(apset)

# cut hierarchical clustering tree into clusters
# cutoff was chosen by 'common sense' which best merges
# steroids but excludes non-steroids, and merges phenylethylamines but excludes non-phenylethylamines
clusters <- cutree(hc, h=(1-0.30))

write.table(clusters, file=outputFile, quote=FALSE, sep="\t", row.names=TRUE, col.names=FALSE)  

## test code:
#
## get size distribution
#clusterSizes <- table(clusters)
#clusterSizes[clusterSizes > 1]
#max(clusterSizes)
#
## print out largest cluster
#largestCluster <- which.max(clusterSizes)
#largestClusterContents <- names(clusters)[clusters == largestCluster]
#for(i in largestClusterContents){
#    cat(i)
#    cat("\n")
#}
#
## these should cluster together
#clusters[names(clusters) %in% c(6446, 68873, 3000226)]
#clusters[names(clusters) %in% c(4436,5504,681,5826,5906)]
#
## without these
#clusters[names(clusters) %in% c(4463, 68617, 3059)]
#clusters[names(clusters) %in% c(5775,82148,896,3348,4768,4891)]
#
#myRange <- seq(from=0, to=0.6, by=0.01)
#scores <- sapply(as.list(myRange), function(i){
#    clusters <- cutree(hc, h=(1-i))
#
#    together <- clusters[names(clusters) %in% c(4436,5504,681,5826,5906)]
#    togetherTable <- table(together)
#    exclude <- clusters[names(clusters) %in% c(5775,82148,896,3348,4768,4891)]
#    score1 <- max(togetherTable) - sum(exclude %in% names(which.max(togetherTable)))
#
#    together <- clusters[names(clusters) %in% c(6446, 68873, 3000226)]
#    togetherTable <- table(together)
#    exclude <- clusters[names(clusters) %in% c(4463, 68617, 3059)]
#    score2 <- max(togetherTable) - sum(exclude %in% names(which.max(togetherTable)))
#
#    clusterSizes <- table(clusters)
#    clusterSizes[clusterSizes > 1]
#
#    return(cbind(i, sum(score1, score2), max(clusterSizes), sum(clusterSizes[clusterSizes > 1])))
#}) 
#bestScore <- myRange[which.max(scores)]
#max(scores)

# perform clustering with binning clustering
# bins <- cmp.cluster(db=apset, cutoff=0.5)
# binClusters <- bins[,3]
# names(binClusters) <- bins[,1]


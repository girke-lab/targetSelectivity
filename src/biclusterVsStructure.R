#!/usr/bin/env Rscript

# (C) 2015 Tyler William H Backman
# Purpose: compare activity biclusters vs structure and compute jaccard index

library(R.utils)
library(foreach)
library(ChemmineR)
library(clusteval)

# parse input options
biClustersFile <- commandArgs(trailingOnly=TRUE)[1] 
drugBankStructuresFile <- commandArgs(trailingOnly=TRUE)[2]
drugbankLinksFile <- commandArgs(trailingOnly=TRUE)[3]
outputFile <- commandArgs(trailingOnly=TRUE)[4]

# test code for running without make:
if(is.null(commandArgs(trailingOnly=TRUE)[1])){
    biClustersFile <- "working/biClusters.RData"
    drugBankStructuresFile <- "working/drugbank.sdf" 
    drugbankLinksFile <- "working/drugbank_links.csv"
    # outputFile <- "working/biclusterMDS.pdf"
}

# parse input files
load(biClustersFile) # loads clusterResults
drugbankCompounds <- read.SDFset(drugBankStructuresFile)
drugBankTable <- read.csv(drugbankLinksFile)

# reformat clusterResults into two cols: cid, clusterid
clustersPerCid <- foreach(thisCluster=1:length(clusterResults), .combine='rbind2') %do% {
    cbind(clusterResults[[thisCluster]]$rows, thisCluster)
}

targetsPerCluster <- foreach(thisCluster=1:length(clusterResults), .combine='rbind2') %do% {
    cbind(clusterResults[[thisCluster]]$cols, thisCluster)
}

# remove small biclusters
clusterCounts <- table(clustersPerCid[,2])
largeClusters <- names(clusterCounts)[clusterCounts > 1]
clustersPerCid <- clustersPerCid[clustersPerCid[,2] %in% largeClusters,]

# renumber biclusters
newNumbers <- 1:length(unique(clustersPerCid[,2]))
names(newNumbers) <- unique(clustersPerCid[,2])
clustersPerCid[,2] <- as.character(newNumbers[match(clustersPerCid[,2], names(newNumbers))])

# assign each compound to a unique cluster
clustersPerCid <- clustersPerCid[! duplicated(clustersPerCid[,1]),]

# keep only top 12 largest clusters
# topClusters <- names(sort(table(clustersPerCid[,2]), decreasing=T)[1:12])
# clustersPerCid <- clustersPerCid[clustersPerCid[,2] %in% topClusters,]

# fix cids in drugbank compounds
dbIdsWithCids <- drugBankTable[! is.na(drugBankTable[,c('PubChem.Compound.ID')]),c('DrugBank.ID')]
drugbankCompounds <- drugbankCompounds[datablocktag(drugbankCompounds, "DRUGBANK_ID") %in% dbIdsWithCids]
cidPositions <- match(datablocktag(drugbankCompounds, "DRUGBANK_ID"), drugBankTable[,c('DrugBank.ID')])
cid(drugbankCompounds) <- as.character(drugBankTable[cidPositions,c('PubChem.Compound.ID')])
drugbankCompounds <- drugbankCompounds[! duplicated(cid(drugbankCompounds))]

# get structures
biclusterCids <- unique(clustersPerCid[,1])
biclusterCompounds <- drugbankCompounds[cid(drugbankCompounds) %in% biclusterCids]

# Create atom pair distance matrix
apset <- sdf2ap(biclusterCompounds) 

# compute clusters 
myTempFile <- tempfile()
dummy <- cmp.cluster(db=apset, cutoff=0, save.distances=myTempFile)
load(myTempFile)
unlink(myTempFile)
hc <- hclust(as.dist(distmat), method="complete")
hc[["labels"]] <- cid(apset)
clusters <- cutree(hc, k=11)

# compute jaccard coefficient
biclusterLabelsInOrder <- clustersPerCid[match(names(clusters), clustersPerCid[,1]),2]
names(biclusterLabelsInOrder) <- names(clusters)

jaccardSim <- jaccard(clusters, biclusterLabelsInOrder)

# compare to random sample
set.seed(5412)
clusterWeights = table(clusters)/length(unique(clusters))
randomJaccards <- sapply(1:10000, function(x){
    randomStrutures <- sample(1:11, length(biclusterLabelsInOrder), replace=T, prob=clusterWeights)
    jaccard(randomStrutures, biclusterLabelsInOrder)
})
mean(randomJaccards)
sd(randomJaccards)

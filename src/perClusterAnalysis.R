#!/usr/bin/env Rscript

# (C) 2015 Tyler William H Backman
# Purpose: plot and inspect individual drug/target biclusters from output of src/biclusterMatrix.R

library(R.utils)
library(ggplot2)
library(reshape)
library(ChemmineR)
# library(biomaRt)
library(foreach)
library(doMC)
library(Matrix)

# parse input options
drugComparisonMatrixFile <- commandArgs(trailingOnly=TRUE)[1]
biClustersFile <- commandArgs(trailingOnly=TRUE)[2]
curatedClustersFile <- commandArgs(trailingOnly=TRUE)[3]
clusterAnnotationFile <- commandArgs(trailingOnly=TRUE)[4]
drugLinksFile <- commandArgs(trailingOnly=TRUE)[5]
hclustResultFile <- commandArgs(trailingOnly=TRUE)[6]
targetHmmFile <- commandArgs(trailingOnly=TRUE)[7]
PfamClansInputFile <- commandArgs(trailingOnly=TRUE)[8]
# fullActivityMatrix <- commandArgs(trailingOnly=TRUE)[9]
# bioactivityNearestNeighborsFile <- commandArgs(trailingOnly=TRUE)[10]
# structureClustersFile <- commandArgs(trailingOnly=TRUE)[11]
outputFolder <- commandArgs(trailingOnly=TRUE)[9]
cores <- commandArgs(trailingOnly=TRUE)[10]

# test code for running without make:
if(is.na(commandArgs(trailingOnly=TRUE)[1])){
    drugComparisonMatrixFile <- "working/drugComparisonMatrix.RData"
    biClustersFile <- "working/biClusters.RData"
    curatedClustersFile <- "working/curatedClusters.txt"
    clusterAnnotationFile <- "working/clusterAnnotations.csv"
    drugLinksFile <- "working/drugbank_links.csv"
    hclustResultFile <- "../bioassayClusterAnalysis/working/structureClusterDrugs.tab"
    targetHmmFile <- "working/combinedTargetDomainsTwoCols"
    PfamClansInputFile <- "working/Pfam-A.clans.tsv"
    # fullActivityMatrix <- "../bioassayClusterAnalysis/working/activityFps.Rda"
    # bioactivityNearestNeighborsFile <- "../bioassayClusterAnalysis/working/bioactivityNearestNeighbors.Rda"
    # structureClustersFile <- "../bioassayClusterAnalysis/working/structureClusters_20nn.tab"
    outputFolder <- "working/perClusterAnalysis"
    cores <- 1
}

# parse input files
load(drugComparisonMatrixFile) # loads drugComparisonMatrix
load(biClustersFile) # loads clusterResults
curatedClusters <- read.table(curatedClustersFile)
clusterAnnotations <- read.csv(clusterAnnotationFile)
drugLinks <- read.csv(drugLinksFile) 
hclustResult <- read.table(hclustResultFile)
domains <- read.table(targetHmmFile, header = FALSE)
PfamClans <- read.delim(PfamClansInputFile, header=FALSE, na.strings="\\N")
# load(fullActivityMatrix) # loads matrix as "results"
# load(bioactivityNearestNeighborsFile) # loads matrix as "nnm"
# structureClusters <- read.table(structureClustersFile)


# clean up target HMM names
domains[,1] <- gsub("^(PF\\d*).*", "\\1", domains[,1], perl=TRUE)
domains[,2] <- gsub("^drugbank_target\\|(\\w+).*$", "UniProt_\\1", domains[,2], perl = TRUE)
domains[,2] <- gsub("^gi\\|(\\d+).*$", "gi_\\1", domains[,2], perl = TRUE)
colnames(domains) <- c("DOMAIN", "TARGET")

# get cluster stats
clusterStats <- t(sapply(clusterResults, function(x)
    c(length(x$rows), length(x$cols), x$score)))
# clusterStats <- cbind(table(clusterResults$rowClusterIds), table(clusterResults$colClusterIds), c(NA, clusterResults$score))
colnames(clusterStats) <- c("compounds", "target clusters", "score")

# choose only clusters with > 1 compound
largeClusters <- (1:nrow(clusterStats))[clusterStats[,1] > 1]

##################################
# loop over clusters one at a time
##################################
gc()
registerDoMC(cores=cores)
foreach(selectedCluster = largeClusters, .export=ls(envir=globalenv()), .packages=c("ggplot2", "reshape", "ChemmineR", "Matrix")) %dopar% {
# for(selectedCluster in largeClusters){
print(paste("processing bicluster", selectedCluster))
drugIds <- clusterResults[[selectedCluster]]$rows
targetIds <- clusterResults[[selectedCluster]]$cols
subMatrix <- drugComparisonMatrix[as.character(drugIds), as.character(targetIds), drop=F]

################################
# QUANTIFY OVERALL RESULTS
################################
# sort rows and columns of matrix by tanimoto similarity
clusterSort <- function(binaryMatrix){
    if(nrow(binaryMatrix) < 2)
        return(row.names(binaryMatrix))
    fpset <- as(as.matrix(binaryMatrix), "FPset")
    simMA <- sapply(cid(fpset), function(x) fpSim(fpset[x], fpset, sorted=FALSE)) 
    hc <- hclust(as.dist(1-simMA), method="single")
    return(unique(hc$labels[hc$order]))
}
binaryMatrix <- 1*(subMatrix > 1)
sortedClusters <- clusterSort(t(binaryMatrix))
sortedCompounds <- clusterSort(binaryMatrix) 

# use representative protein names for targets
uniProtData <- clusterAnnotations[match(targetIds, clusterAnnotations$uniqueClusterIds),]
uniProtData$protein_name <- gsub("^(.{60}).*$", "\\1...", uniProtData$protein_name)
uniProtData$protein_name <- make.unique(uniProtData$protein_name)
colnames(subMatrix) <- unlist(uniProtData$protein_name) 
sortedClusters <- unlist(uniProtData$protein_name[match(sortedClusters, targetIds)])

# use common names for drugs

row.names(subMatrix) <- as.character(drugLinks$Name[match(row.names(subMatrix), drugLinks$PubChem.Compound.ID)])
sortedCompounds <- as.character(drugLinks$Name[match(sortedCompounds, drugLinks$PubChem.Compound.ID)])

# plot matrix
plotFileName <- file.path(outputFolder, paste("cluster", selectedCluster, ".pdf", sep=""))
drugComparisons <- melt(subMatrix)
colnames(drugComparisons) <- c("drugs", "clusters", "activity")
drugComparisons$drugs <- factor(drugComparisons$drugs, levels=sortedCompounds)
drugComparisons$clusters <- factor(drugComparisons$clusters, levels=sortedClusters)
labels <- c("untested", "inactive", "active", "untested annotated", "inactive annotated", "active annotated")
drugComparisons$activity <- factor(labels[drugComparisons$activity + 1], levels=labels)
# png(filename = plotFileName, width=1024, height=768)
# pdf(file=plotFileName, width=7, height=10)
keyColors <- c("black", "grey", "red", "orange", "blue", "green")
keyColors <- keyColors[as.numeric(row.names(table(subMatrix)))+1]
p <- ggplot(drugComparisons, aes(clusters, drugs)) 
p <- p + geom_tile(aes(fill=activity)) + 
    scale_fill_manual(values=keyColors) + 
    theme(
        axis.text.x=element_text(angle = -45, hjust = 0),
        axis.ticks=element_blank()) +
#         axis.ticks=element_blank(), 
#         axis.text.x = element_blank(),
#         axis.text.y = element_blank()) + 
    labs(
        x="Targets", 
        y="Drugs", 
        title="Drug Target Annotations and Bioactivity", 
        fill="Activity Type")
print(p)
imageWidth <- 3.71 + length(targetIds)*0.235
imageHeight <- 2.75 + length(drugIds)*0.13425
if(imageHeight < 10)
    imageHeight <- 10
if(imageWidth < 7)
    imageWidth <- 7 
ggsave(plotFileName, width=imageWidth, height=imageHeight, limitsize=FALSE)
# dev.off()

# what is the distribution of values
annotationTable <- rep(0,6)
names(annotationTable) <- as.character(0:5)
annotationFrequencies <- table(subMatrix)
annotationTable[names(annotationFrequencies)] <- annotationFrequencies
annotationTable <- cbind(annotationTable[1:3], annotationTable[4:6])
row.names(annotationTable) <- c("Untested", "Inactive", "Active")
colnames(annotationTable) <- c("Unannotated Targets", "Annotated Targets")

# what are the new (unannotated) compounds and targets?
drugCount <- nrow(subMatrix)
targetClusters <- ncol(subMatrix)
clusterTables <- apply(subMatrix, 2, table)
annotatedClusters <- names(clusterTables)[sapply(clusterTables, function(x)
     (!is.na(x["3"]) || !is.na(x["4"]) || !is.na(x["5"])))]
newActiveClusters <- names(clusterTables)[sapply(clusterTables, function(x)
     (!is.na(x["2"]) && is.na(x["3"]) && is.na(x["4"]) && is.na(x["5"])))]
drugTables <- apply(subMatrix, 1, table)
annotatedDrugs <- names(drugTables)[sapply(drugTables, function(x)
     (!is.na(x["3"]) || !is.na(x["4"]) || !is.na(x["5"])))]
newActiveDrugs <- names(drugTables)[sapply(drugTables, function(x)
    (!is.na(x["2"]) && is.na(x["3"]) && is.na(x["4"]) && is.na(x["5"])))]

################################
# INVESTIGATE COMPOUNDS
################################

# compare compounds w/ structural clusters

# get cluster IDs for each drug
hclustResult <- hclustResult[hclustResult[,1] %in% row.names(drugComparisonMatrix),]
drugClusterIds <- hclustResult[match(drugIds, hclustResult[,1]),2]
clusterFreq <- table(drugClusterIds)
clusterFreq <- clusterFreq[order(clusterFreq, decreasing=TRUE)]

clusterSizes <- table(hclustResult[,2])
clusterDistribution <- cbind(clusterSizes[names(clusterFreq)], clusterFreq)
colnames(clusterDistribution) <- c("Structural Cluster Size", "Drugs From This Bicluster")

biggestClusterSize <- clusterFreq[[1]]
biggestClusterName <- names(clusterFreq[1])
biggestClusterTotalSize <- table(hclustResult[,2])[biggestClusterName]

# out of length(drugIds) drugs in this bicluster, biggestClusterSize fall into the same
# structural cluster which consists of biggestClusterTotalSize drugs

################################
# INVESTIGATE TARGETS
################################

# # find top GO terms among targets
# unimart <- useMart("unimart")
# uniprot <- useDataset("uniprot", mart=unimart)
# uniProtAccessions <- unique(uniProtData$accession)
# uniprotGO <- getBM(attributes=c("accession", "go_id", "go_name"), filters=c("accession"), mart=uniprot, values=uniProtAccessions)
# topGO <- table(uniprotGO$go_name)
# topGO <- topGO[order(topGO, decreasing=TRUE)]
# uniprotKeywords <- getBM(attributes=c("accession", "keyword"), filters=c("accession"), mart=uniprot, values=uniProtAccessions)
# topKeywords <- table(uniprotKeywords$keyword)
# topKeywords <- topKeywords[order(topKeywords, decreasing=TRUE)]

# find top Pfam domains among targets
domainsPerCluster <- unlist(lapply(targetIds, function(x){
    allTargetIds <- as.character(unique(curatedClusters[curatedClusters[,2] %in% x,1]))
    allDomains <- domains$DOMAIN[domains$TARGET %in% allTargetIds]
    return(unique(allDomains))
}))
topDomains <- table(domainsPerCluster)
topDomains <- topDomains[order(topDomains, decreasing=TRUE)]
topDomainDescriptions <- as.character(PfamClans[match(names(topDomains), PfamClans[,1]),5])
domainTable <- cbind(paste(names(topDomains), topDomainDescriptions), topDomains)
colnames(domainTable) <- c("Pfam Family", "Frequency")

################################
# INVESTIGATE NEAREST NEIGHBORS 
################################

# # pull fraction of entire distance matrix for this cluster
# drugIdsFullMatrix <- as.matrix(results[row.names(results) %in% drugIds,,drop=F])
# 
# # determine similarity threshold within this cluster
# similarityMatrix <- 1 - dist(drugIdsFullMatrix, method="binary") 
# similarityThreshold <- mean(similarityMatrix) - sd(similarityMatrix) # Z = -1
# 
# # get all nearest neighbors for compounds in this cluster
# clusterNN <- nnm[row.names(nnm) %in% drugIds,,drop=F]
# 
# # count how many are already in this cluster
# uniqueNNs <- unique(row.names(nnm)[clusterNN[,2:dim(clusterNN)[2]]])
# alreadyInClusterCounts <- table(uniqueNNs %in% drugIds)
# if(length(alreadyInClusterCounts) > 1){
#     alreadyInClusterCounts <- table(uniqueNNs %in% drugIds)[["TRUE"]]
# } else {
#     alreadyInClusterCounts <- 0
# }
# notInClusterCIDs <- uniqueNNs[! uniqueNNs %in% drugIds]
# 
# # extract activity matrix for those not in cluster
# newCidsMatrix <- as.matrix(results[row.names(results) %in% notInClusterCIDs,,drop=F])
# 
# # get tanimoto for each new CID against all targets in bicluster
# biclusterActivityFpset <- as(drugIdsFullMatrix, "FPset")
# newCidsActivityFpset <- as(newCidsMatrix, "FPset")
# simMA <- sapply(cid(newCidsActivityFpset), function(x) fpSim(newCidsActivityFpset[x], biclusterActivityFpset, sorted=FALSE, addone=FALSE))
# 
# # apply threshold to similarities
# newCidsMaxSim <- apply(simMA, 2, max)
# thresholdNewCids <- names(newCidsMaxSim)[newCidsMaxSim >= similarityThreshold]
# 
# # investigate structural relationships by looking at full jarvis patrick clustering
# # get structural clusters with > 1 compound from bicluster
# clusterSizesThisBicluster <- table(structureClusters[structureClusters[,1] %in% drugIds, 2])
# clusterSizesThisBicluster <- clusterSizesThisBicluster[order(clusterSizesThisBicluster, decreasing=T)]
# clusterSizesThisBicluster <- clusterSizesThisBicluster[clusterSizesThisBicluster > 1]
# # get structural clusters common with above that contain new cids
# clusterSizesNewCids <- table(structureClusters[structureClusters[,1] %in% thresholdNewCids, 2])
# matchPositions <- match(names(clusterSizesThisBicluster), names(clusterSizesNewCids))
# matchNewCids <- clusterSizesNewCids[matchPositions]
# matchNewCids[is.na(matchNewCids)] <- 0
# # get full cluster sizes
# clusterSizes <- table(structureClusters[structureClusters[,2] %in% names(clusterSizesThisBicluster),2])
# clusterSizes <- clusterSizes[match(names(clusterSizesThisBicluster), names(clusterSizes))]
# # make table
# nnStructureClusterTable <- cbind(clusterSizes, clusterSizesThisBicluster, matchNewCids)
# colnames(nnStructureClusterTable) <- c("Structural Cluster Size", "Drugs From This Bicluster", "New Compounds With Similar Activity")
# 
# # delete objects that we don't need to save
# rm(drugIdsFullMatrix)
# rm(newCidsMatrix)
# rm(biclusterActivityFpset)
# rm(newCidsActivityFpset)

################################
# save results
################################

save(list=c(ls(all.names=T), "clusterStats"), file=file.path(outputFolder, paste("cluster", selectedCluster, ".RData", sep="")))
}

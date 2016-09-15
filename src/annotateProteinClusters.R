#!/usr/bin/env Rscript

# (C) 2015 Tyler William H Backman
# Purpose: annotate kclust clusters via uniprot
#   choose a representative protein target for each cluster with the following priority order:
#   (1) already annotated in DrugBank, (2) human, (3) non-human

library(R.utils)
# library(biomaRt)
source("src/uniprotInterface.R") # replacement for biomaRt
library(bioassayR)

# parse input options
curatedClustersFile <- commandArgs(trailingOnly=TRUE)[1]
databaseFile <- commandArgs(trailingOnly=TRUE)[2]
outputFile <- commandArgs(trailingOnly=TRUE)[3]

# test code for running without make:
if(is.na(commandArgs(trailingOnly=TRUE)[1])){
    curatedClustersFile <- "working/curatedClusters.txt"
    databaseFile <- "working/bioassayDatabase.sqlite"
    outputFile <- "working/clusterAnnotations.csv"
}

# parse input files
curatedClusters <- read.table(curatedClustersFile)
database <- connectBioassayDB(databaseFile)

# for each protein target cluster choose a representative UniProt ID and extract details
uniqueClusterIds <- unique(curatedClusters[,2])
# unimart <- useMart("unimart")
# uniprot <- useDataset("uniprot", mart=unimart)
uniProtData <- sapply(uniqueClusterIds, function(x){
    allTargetIds <- as.character(curatedClusters[curatedClusters[,2] %in% as.character(x),1])

    # find genbank and uniprot ids
    uniProtTargets <- allTargetIds[grep("^UniProt_", allTargetIds)]
    uniProtTargets <- unique(gsub("^UniProt_(.*)", "\\1", uniProtTargets))
    giTargets <- allTargetIds[grep("^gi_", allTargetIds)]
    giTargets <- unique(gsub("^gi_(.*)", "\\1", giTargets))

    # if there are genbank IDs, get them all from the database
    if(length(giTargets) > 0){
        uniProtGITargets <- unique(unlist(lapply(giTargets, translateTargetId, database=database,category="UniProt")))
    } else 
        uniProtGITargets <- c()
    completeTargetList <- unique(c(uniProtTargets, uniProtGITargets))

    # look up uniprot IDs 
    # proteinDetails <- getBM(attributes=c("accession", "gene_name", "organism", "protein_name"), filters=c("accession"), mart=uniprot, values=completeTargetList)
    proteinDetails <- UniProtAnnotator(completeTargetList, useBiomartNames=T)
    
    # choose a representative ID: DrugBank annotated first, then human, otherwise whatever is left
    # if nothing, just return the first identifier from the raw clustering file
    annotatedDetails <- proteinDetails[proteinDetails$accession %in% uniProtTargets,]
    humanDetails <- proteinDetails[grepl("Homo sapiens", proteinDetails$organism),]
    if(nrow(annotatedDetails) != 0)
        return(annotatedDetails[1,]) # first annotated name
    else if(nrow(humanDetails) != 0)
        return(humanDetails[1,]) # first human name
    else if(nrow(proteinDetails) != 0)
        return(proteinDetails[1,]) # first non-human name
    else
        return(t(data.frame(c(NA, NA, NA, allTargetIds[1]), row.names=colnames(proteinDetails))))
})
uniProtData <- t(as.data.frame(matrix(unlist(uniProtData), ncol=ncol(uniProtData)), row.names=row.names(uniProtData)))
if(nrow(uniProtData) != length(uniqueClusterIds))
    stop("ERROR: invalid length of biomaRt returned result")
clusterAnnotations <- cbind(uniqueClusterIds, uniProtData)

write.csv(clusterAnnotations, outputFile, row.names=FALSE)

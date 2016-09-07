#!/usr/bin/env Rscript

# (C) 2015 Tyler William H Backman
# Purpose: perform comparison between bioassay and drugbank targets and store as matrix
#  drug (row) vs cluster (col) object is called "drugComparisonMatrix" with the following values:
#  0 = untested
#  1 = inactive
#  2 = active
#  3 = untested and annotated
#  4 = inactive and annotated
#  5 = active and annotated 

library(R.utils)
library(bioassayR)

# parse input options
bioassayDatabaseFile <- commandArgs(trailingOnly=TRUE)[1]
highlyScreenedCidsFile <- commandArgs(trailingOnly=TRUE)[2]
curatedClustersFile <- commandArgs(trailingOnly=TRUE)[3]
drug_target_uniprot_linksFile <- commandArgs(trailingOnly=TRUE)[4]
drugbank_linksFile <- commandArgs(trailingOnly=TRUE)[5]
outputFile <- commandArgs(trailingOnly=TRUE)[6]

# test code for running without make:
if(is.na(commandArgs(trailingOnly=TRUE)[1])){
    bioassayDatabaseFile <- "working/bioassayDatabase.sqlite"
    highlyScreenedCidsFile <- "working/highlyScreenedCids.txt"
    curatedClustersFile <- "working/curatedClusters.txt"
    drug_target_uniprot_linksFile <- "working/drug_target_uniprot_links.csv"
    drugbank_linksFile <- "working/drugbank_links.csv"
    outputFile <- "working/drugComparisonMatrix.RData"
}

# parse input files
bioassayDatabase <- connectBioassayDB(bioassayDatabaseFile)
highlyScreenedCids <- read.table(highlyScreenedCidsFile)[[1]]
curatedClusters <- read.table(curatedClustersFile)
drug_target_uniprot_links <- read.csv(drug_target_uniprot_linksFile)
drugbank_links <- read.csv(drugbank_linksFile)
mergedDrugBank <- merge(drugbank_links, drug_target_uniprot_links, by = "DrugBank.ID", all.x = T)

# determine list of cids that meet the following criteria:
#   -screened against at least 10 distinct targets
#   -are annotated drugs in drugbank with a listed cid
allDrugs <- unique(drugbank_links$PubChem.Compound.ID[! is.na(drugbank_links$PubChem.Compound.ID)])
highlyScreenedDrugs <- allDrugs[allDrugs %in% highlyScreenedCids]

# build empty matrix
drugComparisonMatrix <- matrix(data = 0, nrow=length(highlyScreenedDrugs), ncol=length(unique(curatedClusters[,2])))
row.names(drugComparisonMatrix) <- highlyScreenedDrugs
colnames(drugComparisonMatrix) <- unique(curatedClusters[,2])

# load active and inactive scores into matrix
for(drug in highlyScreenedDrugs){
    inactives <- row.names(inactiveTargets(bioassayDatabase, drug))
    actives <- row.names(activeTargets(bioassayDatabase, drug))
    inactiveClusters <- unique(curatedClusters[curatedClusters[,1] %in% paste("gi_", inactives, sep=""),2])
    activeClusters <- unique(curatedClusters[curatedClusters[,1] %in% paste("gi_", actives, sep=""),2])
    # note: some targets have both active and inactive scores, for these the
    #  actives will overwrite the inactive score
    drugComparisonMatrix[as.character(drug), as.character(inactiveClusters)] <- 1
    drugComparisonMatrix[as.character(drug), as.character(activeClusters)] <- 2 
}

# load DrugBank target annotations into matrix
for(drug in highlyScreenedDrugs){
    drugTargets <- unique(as.character(mergedDrugBank$UniProt.ID.y[mergedDrugBank$PubChem.Compound.ID %in% drug]))
    drugClusters <- unique(curatedClusters[curatedClusters[,1] %in% paste("UniProt_", drugTargets, sep=""),2])
    drugComparisonMatrix[as.character(drug), as.character(drugClusters)] <- drugComparisonMatrix[as.character(drug), as.character(drugClusters)] + 3
}

# drop "0" columns (targets with no annotations OR tested values)
drugComparisonMatrix <- drugComparisonMatrix[,colSums(drugComparisonMatrix) > 0]

# save resulting matrix
save(list=c("drugComparisonMatrix"), file=outputFile)

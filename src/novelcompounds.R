#!/usr/bin/env Rscript

# (C) 2016 Tyler William H Backman
# Purpose: identify both new compounds targeting the FDA approved space,
# and compounds targeting novel space

library(R.utils)
library(Matrix)

# parse input options
cidsVStargetFile <- commandArgs(trailingOnly=TRUE)[1]
promiscuityProbabilityFile <- commandArgs(trailingOnly=TRUE)[2]
databaseFile <- commandArgs(trailingOnly=TRUE)[3]
outputFilename <- commandArgs(trailingOnly=TRUE)[4]
cores <- commandArgs(trailingOnly=TRUE)[5]

# test code for running without make:
if(is.na(commandArgs(trailingOnly=TRUE)[1])){
    cidsVStargetFile <- "working/cidsVStargets.RData"
    drugbank_linksFile <- "working/drugbank_links.csv"
}

# parse input files
load(cidsVStargetFile)
compoundVsTargetMatrix <- results # rows = targets, cols = compounds
rm(results)
drugbank_links <- read.csv(drugbank_linksFile)
drugCids <- unique(drugbank_links$PubChem.Compound.ID)
drugCids <- drugCids[! is.na(drugCids)]

# make matrix binary
binaryMatrix <- 1*(compoundVsTargetMatrix == 2)

# get list of active targets for FDA approved compounds
FDAapprovedSubset <- binaryMatrix[,colnames(binaryMatrix) %in% drugCids]
fdaTargets <- row.names(FDAapprovedSubset)[rowSums(FDAapprovedSubset) > 0]

# get non-FDA compounds active against FDA approved targets
fdaTargetMatrix <- binaryMatrix[row.names(binaryMatrix) %in% fdaTargets,]
nonDrugMatrix <- fdaTargetMatrix[, ! colnames(fdaTargetMatrix) %in% drugCids]
nonDrugCids <- colnames(nonDrugMatrix)[colSums(nonDrugMatrix) > 0]
# note : nonDrugCids = 486348

# get non-FDA compounds active against FDA approved targets
nonfdaTargetMatrix <- binaryMatrix[! row.names(binaryMatrix) %in% fdaTargets,]
novelTargetCids <- colnames(nonfdaTargetMatrix)[colSums(nonfdaTargetMatrix) > 0]
# note: novelTargetCids = 153402
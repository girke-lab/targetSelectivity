#!/usr/bin/env Rscript

# (C) 2015 Tyler William H Backman
# Purpose: create sparse matrix showing how many times each
#   cid vs target combination is repeated in unique target assays

library(R.utils)
library(bioassayR)
library(Matrix)
library(doMC)
library(foreach)

# parse input options
bioassayDatabaseFile <- commandArgs(trailingOnly=TRUE)[1]
outputFile <- commandArgs(trailingOnly=TRUE)[2]
cores <- commandArgs(trailingOnly=TRUE)[3]

# test code for running without make:
if(is.na(commandArgs(trailingOnly=TRUE)[1])){
    bioassayDatabaseFile <- "working/bioassayDatabase.sqlite"
    outputFile <- "working/replicateMatrix.RData"
    cores <- 1
}

# parse input files
bioassayDatabase <- connectBioassayDB(bioassayDatabaseFile)

# get all targets that were screened at least twice
replicateCounts <- queryBioassayDB(bioassayDatabase, "SELECT COUNT(*) AS replicates, target FROM targets WHERE target_type = 'protein' GROUP BY target")
multiAssayTargets <- replicateCounts$target[replicateCounts$replicates > 1]

# build empty matrix with a row for each cid, and a column for each target
allCids <- queryBioassayDB(bioassayDatabase, "SELECT DISTINCT cid FROM activity")
allCids <- unique(allCids[[1]])
replicateMatrix <- sparseMatrix(
    i = {},
    j = {},
    dims = c(length(allCids), length(multiAssayTargets)),
    dimnames = list(allCids, multiAssayTargets)
)
replicateMatrix <- as(replicateMatrix, "dgCMatrix")

# get a list of assays that contain only a single target
targetCounts <- queryBioassayDB(bioassayDatabase, "SELECT COUNT(*) AS targets, aid FROM targets GROUP BY aid")
singleTargetAssays <- targetCounts$aid[targetCounts$target == 1]

# loop over multi target compounds and get activity matrix for each
disconnectBioassayDB(bioassayDatabase)
gc()
registerDoMC(cores=cores)
results <- foreach(thisTarget = multiAssayTargets, .packages=c("bioassayR")) %dopar% {
    # get activity matrix for this target
    db <- connectBioassayDB(bioassayDatabaseFile)
    aids <- queryBioassayDB(db, paste("SELECT DISTINCT aid FROM targets WHERE target = '", thisTarget, "'", sep=""))[[1]]
    aids <- aids[aids %in% singleTargetAssays]
    if(length(aids) == 0) return(NA) 
    thisTargetAssays <- getAssays(db, aids)
    disconnectBioassayDB(db)
    thisTargetActivityMatrix <- slot(thisTargetAssays, "activity")

    # make activity matrix binary showing actives and inactives as 1
    binaryMatrix <- 1*(thisTargetActivityMatrix > 0)

    # keep only replicates that occur more than once
    replicateCountsByCid <- colSums(binaryMatrix)
    replicateCountsByCid <- replicateCountsByCid[! is.na(replicateCountsByCid)]
    replicateCountsByCid <- replicateCountsByCid[replicateCountsByCid > 1]

    return(replicateCountsByCid)
}

# insert counts into replicateMatrix
targetIndex <- 1:length(multiAssayTargets)
targetIndex <- targetIndex[!is.na(results)]
for(x in targetIndex){
    thisTarget <- multiAssayTargets[[x]]
    replicateCountsByCid <- results[[x]]
    replicateMatrix[as.character(names(replicateCountsByCid)),as.character(thisTarget)] <- replicateCountsByCid
}

# remove zero rows and columns from replicate matrix
replicateMatrix <- replicateMatrix[,colSums(replicateMatrix) > 1, drop = F]
replicateMatrix <- replicateMatrix[rowSums(replicateMatrix) > 1,, drop = F]

# save resulting matrix
save(list=c("replicateMatrix"), file=outputFile)

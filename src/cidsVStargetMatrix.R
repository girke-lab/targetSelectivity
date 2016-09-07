#!/usr/bin/env Rscript

# (C) 2016 Tyler William H Backman
# Purpose: create sparse matrix of highly screened cids
#   vs protein targets, including inactive scores

# targets are clustered if curatedClustersFile is a valid filename,
# individual targets are used if it is "none"

# the option "highlyScreenedCidsFile" specifies a list of compounds
# to use. If this value is "none" all compounds in the database are considered
# instead of just a subset.

library(R.utils)
library(bioassayR)
library(foreach)
library(doMC)

databaseFile <- commandArgs(trailingOnly=TRUE)[1]
highlyScreenedCidsFile <- commandArgs(trailingOnly=TRUE)[2]
curatedClustersFile <- commandArgs(trailingOnly=TRUE)[3]
outputFile <- commandArgs(trailingOnly=TRUE)[4]
cores <- commandArgs(trailingOnly=TRUE)[5]

# test code for running without make:
if(is.na(commandArgs(trailingOnly=TRUE)[1])){
    databaseFile <- "working/bioassayDatabase.sqlite"
    highlyScreenedCidsFile <- "working/highlyScreenedCids.txt"
    curatedClustersFile <- "working/curatedClusters.txt"
    outputFile <- "working/cidsVStargetClusters.RData"
    cores <- 1
}
# parse input files
database <- connectBioassayDB(databaseFile) 

# build list of compounds to screen
if(highlyScreenedCidsFile != "none"){
    cids <- read.table(highlyScreenedCidsFile)[[1]]
} else {
    # note: make these unique because of a few improperly imported assays in the database- also needs fixing
    cids <- unique(as.integer(allCids(database)))
    cids <- cids[cids > 0]
}

# build target translation list
if(curatedClustersFile != "none"){
    curatedClusters <- read.table(curatedClustersFile)
    targetTable <- queryBioassayDB(database, "SELECT aid, target FROM targets WHERE target_type = 'protein'")
    targetClusters <- as.character(curatedClusters[match(paste("gi", targetTable[,2], sep="_"), curatedClusters[,1]),2])
    names(targetClusters) <- as.character(targetTable$aid)
    orderedTargetList <- unique(targetClusters)
} else {
    targetClusters <- FALSE
    orderedTargetList <- allTargets(database)
}

splitCids <- split(cids, as.integer((1:length(cids))/1000))
disconnectBioassayDB(database)
gc()
registerDoMC(cores=cores)
results <- foreach(i = splitCids, .combine='cbind2', .packages=c("bioassayR")) %dopar% {
    tempDBC <- connectBioassayDB(databaseFile) 
    assays <- getBioassaySetByCids(tempDBC, i)
    disconnectBioassayDB(tempDBC)
    perTargetMatrix(assays, inactives = T, assayTargets = targetClusters, targetOrder = orderedTargetList)
}

# remove target rows without any data (activity or inactivity) 
results <- results[rowSums(results) != 0,]

# save output
save(list = c("results"), file = outputFile)

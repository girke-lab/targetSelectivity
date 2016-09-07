#!/usr/bin/env Rscript

# (C) 2015 Tyler William H Backman
# Purpose: extract replicates and investigate one at a time- and output counts

library(R.utils)
library(bioassayR)
library(Matrix)
library(doMC)
library(foreach)

# parse input options
replicateMatrixFile <- commandArgs(trailingOnly=TRUE)[1]
bioassayDatabaseFile <- commandArgs(trailingOnly=TRUE)[2]
outputFile <- commandArgs(trailingOnly=TRUE)[3]
cores <- commandArgs(trailingOnly=TRUE)[4]

# test code for running without make:
if(is.na(commandArgs(trailingOnly=TRUE)[1])){
    replicateMatrixFile <- "working/replicateMatrix.RData"
    bioassayDatabaseFile <- "working/bioassayDatabase.sqlite"
    outputFile <- "working/replicateStats.RData"
    cores <- 2
}

# parse input files
load(replicateMatrixFile) # loads replicateMatrix input
bioassayDatabase <- connectBioassayDB(bioassayDatabaseFile)
cores <- as.numeric(cores)

# make replicate matrix binary (mark all positions tested 2 or more times)
binaryMatrix <- 1*(replicateMatrix > 1)

# get a list of assays that contain only a single target
targetCounts <- queryBioassayDB(bioassayDatabase, "SELECT COUNT(*) AS targets, aid FROM targets GROUP BY aid")
singleTargetAssays <- targetCounts$aid[targetCounts$targets == 1]

# loop over non-zero coordinates in binary matrix and count inactives for each
nonZero <- which(binaryMatrix == 1, arr.ind = T)
indicies <- 1:nrow(nonZero)
indicies <- split(indicies, as.integer((1:length(indicies))/(length(indicies)/cores + 1)))
gc()
registerDoMC(cores=cores)
agreement <- foreach(thisRowSet = indicies, .combine=c, .packages=c("bioassayR")) %dopar% {
    db <- connectBioassayDB(bioassayDatabaseFile)
    result <- foreach(thisRow = thisRowSet, .combine=c, .packages=c("bioassayR")) %do% {
            cid <- row.names(replicateMatrix)[nonZero[thisRow,1]]
            target <- colnames(replicateMatrix)[nonZero[thisRow,2]]
            activity <- queryBioassayDB(db, 
                paste("SELECT * FROM activity NATURAL JOIN targets WHERE cid = '", cid, "' AND target = '", target, "'", sep=""))
            activity <- activity[activity$aid %in% singleTargetAssays,,drop=F] # keep only single target assay data
            activity <- activity[! is.na(activity$activity),,drop=F]
            agreementTable <- table(activity$activity)
            inactiveCount <- agreementTable["0"]
            if(is.na(inactiveCount)) return(0)
            return(inactiveCount)
    }
    disconnectBioassayDB(db)
    return(result)
}

# save(list = c("agreement"), file = "working/agreement.RData")
# load("working/agreement.RData")

# get number of times tested for each pair
timesTested <- replicateMatrix[nonZero]

# investigate pairs
pairCoords <- which(timesTested == 2)
pairInactive <- agreement[pairCoords]
pairInactiveFrequency <- table(pairInactive)

# investigate triplicates
tripleCoords <- which(timesTested == 3)
tripleInactive <- agreement[tripleCoords]
tripleInactiveFrequency <- table(tripleInactive)

# investigate quadruplicates
quadCoords <- which(timesTested == 4)
quadInactive <- agreement[quadCoords]
quadInactiveFrequency <- table(quadInactive)

# write output
replicateStats <- list(pairInactiveFrequency, tripleInactiveFrequency, quadInactiveFrequency)
save(list=c("replicateStats"), file=outputFile)

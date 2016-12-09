#!/usr/bin/env Rscript

# (C) 2016 Tyler William H Backman
# Purpose: compute target selectivity against compounds from an individual domain 

library(R.utils)
library(bioassayR)
library(foreach)
library(doMC)

# parse input options
databaseFile <- commandArgs(trailingOnly=TRUE)[1]
highlyScreenedCidsFile <- commandArgs(trailingOnly=TRUE)[2]
cores <- commandArgs(trailingOnly=TRUE)[3]
drugbank_linksFile <- commandArgs(trailingOnly=TRUE)[4]

# test code for running without make:
if(is.na(commandArgs(trailingOnly=TRUE)[1])){
    databaseFile <- "/dev/shm/bioassayDatabase.sqlite"
    highlyScreenedCidsFile <- "working/highlyScreenedCids.txt"
    cores <- 2
    drugbank_linksFile <- "working/drugbank_links.csv"
}

# parse input files
cids <- read.table(highlyScreenedCidsFile)[[1]]
drugbank_links <- read.csv(drugbank_linksFile)
drugCids <- unique(drugbank_links$PubChem.Compound.ID)
drugCids <- drugCids[! is.na(drugCids)]

# get list of most highly screened domains of FDA approved drugs

domainStats <- function(queryDomain, db){
    # get list of assay ids (aids) whose targets have this domain
    allTargetsWithDomain <- translateTargetId(db, queryDomain, category="GI", fromCategory="domains")
    allAssayTargets <- queryBioassayDB(db, "SELECT aid, target FROM targets WHERE target_type = 'protein'")
    allAssaysWithDomain <- unique(allAssayTargets$aid[allAssayTargets$target %in% allTargetsWithDomain])

    # create compound vs. target matrix
    assays <- getAssays(db, allAssaysWithDomain)
    targetMatrix <- perTargetMatrix(assays)

    # keep only highly screened active compounds
    targetMatrixHs <- targetMatrix[,colSums(targetMatrix > 0) > 9]
    if(ncol(targetMatrixHs) == 0)
        return(NA)
    targetMatrixActives <- targetMatrixHs[,colSums(targetMatrixHs == 2) > 0]
    if(ncol(targetMatrixActives) == 1)
        return(NA)

    # get screening frequency
    screeningFrequency <- colSums(1*(targetMatrixActives > 0))
    drugScreeningFrequency <- screeningFrequency[names(screeningFrequency) %in% drugCids]
    nonDrugScreeningFrequency <- screeningFrequency[! names(screeningFrequency) %in% drugCids]

    # get activity frequency
    activeFrequency <- colSums(1*(targetMatrixActives == 2))
    drugActiveFrequency <- activeFrequency[names(activeFrequency) %in% drugCids]
    nonDrugActiveFrequency <- activeFrequency[! names(activeFrequency) %in% drugCids]
    
    mergedValues <- rbind(
          cbind("drugScreeningFrequency", drugScreeningFrequency),
          cbind("nonDrugScreeningFrequency", nonDrugScreeningFrequency),
          cbind("drugActiveFrequency", drugActiveFrequency),
          cbind("nonDrugActiveFrequency", nonDrugActiveFrequency)
    )
    return(cbind(queryDomain, mergedValues))
}

# test code
allDomains <- queryBioassayDB(myDB, "SELECT distinct identifier from targetTranslations where category = 'domains' limit 10")[[1]]

gc()
registerDoMC(cores=cores)
results <- foreach(i = allDomains, .combine="rbind", .packages=c("bioassayR")) %dopar% {
    print(i)
    tempDBC <- connectBioassayDB(databaseFile)
    result <- domainStats(i, tempDBC)
    disconnectBioassayDB(tempDBC)
    result
}


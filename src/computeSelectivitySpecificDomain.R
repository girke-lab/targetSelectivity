#!/usr/bin/env Rscript

# (C) 2016 Tyler William H Backman
# Purpose: compute target selectivity against compounds from an individual domain 

library(R.utils)
library(bioassayR)
library(foreach)
library(doMC)
library(doRNG)

# parse input options
databaseFile <- commandArgs(trailingOnly=TRUE)[1]
highlyScreenedCidsFile <- commandArgs(trailingOnly=TRUE)[2]
drugbank_linksFile <- commandArgs(trailingOnly=TRUE)[3]
outputFile <- commandArgs(trailingOnly=TRUE)[4]
cores <- commandArgs(trailingOnly=TRUE)[5]

# test code for running without make:
if(is.na(commandArgs(trailingOnly=TRUE)[1])){
    databaseFile <- "/dev/shm/bioassayDatabase.sqlite"
    # databaseFile <- "working/bioassayDatabaseSingleTarget.sqlite"
    highlyScreenedCidsFile <- "working/highlyScreenedCids.txt"
    cores <- 64 
    drugbank_linksFile <- "working/drugbank_links.csv"
    compoundsByDomain <- "working/compoundsByDomain.RData"
    cidsVStargetsFile <- "working/cidsVStargets.RData"
    outputFile <- "working/computeSelectivitySpecificDomain.RData"
}

# parse input files
db <- connectBioassayDB(databaseFile)
cids <- read.table(highlyScreenedCidsFile)[[1]]
drugbank_links <- read.csv(drugbank_linksFile)
drugCids <- unique(drugbank_links$PubChem.Compound.ID)

# load previously saved domainCounts without polluting env
tmp.env <- new.env()
load(compoundsByDomain, envir=tmp.env)
domainCounts <- get("domainCounts", pos=tmp.env)
rm(tmp.env)

# load previously saved domainCounts without polluting env
tmp.env <- new.env()
load(cidsVStargetsFile, envir=tmp.env)
cidsVStargetsMatrix <- get("results", pos=tmp.env)
rm(tmp.env)

# get list of all domains with 2+ targets
totalTargets <- domainCounts[order(domainCounts$totalTargets, decreasing=TRUE),"totalTargets", drop=F]
multiTargetDomains <- row.names(totalTargets)[totalTargets > 1]

# optional code: instead get top 35 domains with most active FDA approved drugs
# multiTargetDomains <- row.names(domainCounts)[order(domainCounts$drugCidCountActives, decreasing=TRUE)[1:35]]

domainStats <- function(queryDomain, db, equalizeMedians=TRUE){
    allTargetsWithDomain <- translateTargetId(db, queryDomain, category="GI", fromCategory="domains")
    matrixSubset <- cidsVStargetsMatrix[rownames(cidsVStargetsMatrix) %in% allTargetsWithDomain,,drop=F]

    # keep only subset with tested compounds
    matrixSubset <- matrixSubset[,colSums(matrixSubset) > 0,drop=F]
    if(ncol(matrixSubset) == 0)
        return(NULL)

    # keep only highly screened active compounds
    targetMatrixHs <- matrixSubset[,colSums(matrixSubset > 0) > 9,drop=F]
    if(ncol(targetMatrixHs) == 0)
        return(NULL)
    targetMatrixActives <- targetMatrixHs[,colSums(targetMatrixHs == 2) > 0,drop=F]
    if(ncol(targetMatrixActives) == 0)
        return(NULL)

    if(sum(colnames(targetMatrixActives) %in% drugCids) == 0)
        return(NULL)
    if(sum(! colnames(targetMatrixActives) %in% drugCids) == 0)
        return(NULL)

    # get screening frequency
    screeningFrequency <- colSums(1*(targetMatrixActives > 0))
    drugScreeningFrequency <- screeningFrequency[names(screeningFrequency) %in% drugCids]
    nonDrugScreeningFrequency <- screeningFrequency[! names(screeningFrequency) %in% drugCids]

    # note: comment this section out to avoid excluding the most highly screened drugs
    # iteratively drop most screened FDA approved drugs until the median screening frequency is the same or lower
    if(equalizeMedians){
        while(median(drugScreeningFrequency) > median(nonDrugScreeningFrequency)){
            mostScreenedDrug <- names(which.max(drugScreeningFrequency))
            activityCol <- targetMatrixActives[,mostScreenedDrug]
            randomRow <- sample(which(activityCol > 0), 1)
            targetMatrixActives[randomRow,mostScreenedDrug] <- 0

            screeningFrequency <- colSums(1*(targetMatrixActives > 0))
            drugScreeningFrequency <- screeningFrequency[names(screeningFrequency) %in% drugCids]
        }
    }

    # get activity frequency
    activeFrequency <- colSums(1*(targetMatrixActives == 2))
    drugActiveFrequency <- activeFrequency[names(activeFrequency) %in% drugCids]
    drugActiveFrequency <- drugActiveFrequency[names(drugActiveFrequency) %in% names(drugScreeningFrequency)]
    nonDrugActiveFrequency <- activeFrequency[! names(activeFrequency) %in% drugCids]
    
    mergedValues <- rbind( 
          cbind("drugScreeningFrequency", drugScreeningFrequency),
          cbind("nonDrugScreeningFrequency", nonDrugScreeningFrequency),
          cbind("drugActiveFrequency", drugActiveFrequency),
          cbind("nonDrugActiveFrequency", nonDrugActiveFrequency)
    )
    return(cbind(queryDomain, mergedValues))
}

gc()
registerDoMC(cores=cores)
set.seed(123)
results <- foreach(i = multiTargetDomains, .combine="rbind", .packages=c("bioassayR")) %dorng% {
    print(i)
    tempDBC <- connectBioassayDB(databaseFile)
    result <- domainStats(i, tempDBC)
    disconnectBioassayDB(tempDBC)
    result
}

results <- data.frame(results, stringsAsFactors=FALSE)
results[,3] <- as.integer(results[,3])
colnames(results) <- c("domain", "category", "frequency")

save(list=c("totalTargets", "results"), file=outputFile)

gc()
registerDoMC(cores=cores)
set.seed(123)
results <- foreach(i = multiTargetDomains, .combine="rbind", .packages=c("bioassayR")) %dorng% {
    print(i)
    tempDBC <- connectBioassayDB(databaseFile)
    result <- domainStats(i, tempDBC, equalizeMedians=FALSE)
    disconnectBioassayDB(tempDBC)
    result
}

results <- data.frame(results, stringsAsFactors=FALSE)
results[,3] <- as.integer(results[,3])
colnames(results) <- c("domain", "category", "frequency")

save(list=c("totalTargets", "results"), file="working/computeSelectivitySpecificDomain_noexclusion.RData")

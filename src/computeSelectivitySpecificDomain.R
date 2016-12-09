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
drugbank_linksFile <- commandArgs(trailingOnly=TRUE)[3]
outputFile <- commandArgs(trailingOnly=TRUE)[4]
cores <- commandArgs(trailingOnly=TRUE)[5]

# test code for running without make:
if(is.na(commandArgs(trailingOnly=TRUE)[1])){
    # databaseFile <- "/dev/shm/bioassayDatabase.sqlite"
    databaseFile <- "working/bioassayDatabaseSingleTarget.sqlite"
    highlyScreenedCidsFile <- "working/highlyScreenedCids.txt"
    cores <- 2
    drugbank_linksFile <- "working/drugbank_links.csv"
    compoundsByDomain <- "working/compoundsByDomain.RData"
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

# get list of most highly active domains of FDA approved drugs
allDomains <- row.names(domainCounts)[order(domainCounts$drugCidCountActives, decreasing=TRUE)[1:20]]

domainStats <- function(queryDomain, db){
    # get list of assay ids (aids) whose targets have this domain
    allTargetsWithDomain <- translateTargetId(db, queryDomain, category="GI", fromCategory="domains")
    allAssayTargets <- queryBioassayDB(db, "SELECT aid, target FROM targets WHERE target_type = 'protein'")
    allAssaysWithDomain <- unique(allAssayTargets$aid[allAssayTargets$target %in% allTargetsWithDomain])
    if(length(allAssaysWithDomain) == 0)
        return(NULL)
    
    # create compound vs. target matrix
    assays <- getAssays(db, allAssaysWithDomain)
    targetMatrix <- perTargetMatrix(assays)

    # keep only highly screened active compounds
    targetMatrixHs <- targetMatrix[,colSums(targetMatrix > 0) > 9,drop=F]
    # targetMatrixHs <- targetMatrix[,colSums(targetMatrix > 0) > 2,drop=F]
    if(ncol(targetMatrixHs) == 0)
        return(NULL)
    targetMatrixActives <- targetMatrixHs[,colSums(targetMatrixHs == 2) > 0,drop=F]
    if(ncol(targetMatrixActives) == 1)
        return(NULL)

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

gc()
registerDoMC(cores=cores)
results <- foreach(i = allDomains, .combine="rbind", .packages=c("bioassayR")) %dopar% {
    print(i)
    tempDBC <- connectBioassayDB(databaseFile)
    result <- domainStats(i, tempDBC)
    disconnectBioassayDB(tempDBC)
    result
}

results <- data.frame(results)
results[,3] <- as.integer(results[,3])
colnames(results) <- c("domain", "category", "frequency")

computeSelectivitySpecificDomain <- results
save(list=c("computeSelectivitySpecificDomain"), file=outputFile)

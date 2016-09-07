#!/usr/bin/env Rscript

# (C) 2016 Tyler William H Backman
# Purpose: compute target selectivity on a per-domain basis 

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
multiTarget <- commandArgs(trailingOnly=TRUE)[6]

# test code for running without make:
if(is.na(commandArgs(trailingOnly=TRUE)[1])){
    databaseFile <- "/dev/shm/bioassayDatabase.sqlite"
    highlyScreenedCidsFile <- "working/highlyScreenedCids.txt"
    drugbank_linksFile <- "working/drugbank_links.csv"
    outputFile <- "working/targetSelectivityByDomain.tab"
    cores <- 2 
    multiTarget <- "keepOne"
}

# parse input files
highlyScreenedCids <- read.table(highlyScreenedCidsFile)[[1]]
drugbank_links <- read.csv(drugbank_linksFile)
drugCids <- unique(drugbank_links$PubChem.Compound.ID)

# get list of all domains in database
database <- connectBioassayDB(databaseFile)
allDomains <- queryBioassayDB(database, "SELECT DISTINCT(identifier) FROM targetTranslations WHERE category = 'domains'")[[1]]

# debug code
# results <- read.table(outputFile)
# missingDomains <- allDomains[! allDomains %in% row.names(results)]
# thisDomain <- missingDomains[2]

# for each domain, compute selectivity of all active targets
gc()
registerDoMC(cores=cores)
results <- foreach(thisDomain = allDomains, .combine="rbind", .packages=c("bioassayR")) %dopar% {
    nullResult <- c(  
        totalDrugs = 0,
        totalOther = 0,
        drugGI = 0,
        otherGI = 0,
        drugkclust = 0,
        otherkclust = 0,
        drugdomain = 0,
        otherdomain = 0 
    )
    nullResult <- t(data.frame(nullResult))
    row.names(nullResult) <- thisDomain
    # thisDomain <- allDomains[1] # TEST CODE
    tempDBC <- connectBioassayDB(databaseFile)
    targetList <- translateTargetId(tempDBC, thisDomain, category="GI", fromCategory="domains") 
    if(length(targetList) < 2)
        if(is.na(targetList))
            return(nullResult)
    activeCompounds <- lapply(targetList, activeAgainst, database=tempDBC)
    activeCids <- unique(unlist(sapply(activeCompounds, function(x) rownames(x))))
    if(is.null(activeCids))
        return(nullResult)
    activeCids <- activeCids[! is.na(activeCids)]
    activeCids <- activeCids[activeCids %in% highlyScreenedCids]
    if(! is.character(activeCids))
        return(nullResult)
    if(length(activeCids) == 0)
        return(nullResult)
    GIselectivity <- targetSelectivity(tempDBC, activeCids, scoring="total", category=FALSE, multiTarget=multiTarget)
    kclustselectivity <- targetSelectivity(tempDBC, activeCids, scoring="total", category="kClust", multiTarget=multiTarget)
    domainselectivity <- targetSelectivity(tempDBC, activeCids, scoring="total", category="domains", multiTarget=multiTarget)
    disconnectBioassayDB(tempDBC)
    result <- c(  
        totalDrugs = sum(activeCids %in% drugCids),
        totalOther = sum(! activeCids %in% drugCids),
        drugGI = median(GIselectivity[names(GIselectivity) %in% drugCids]),
        otherGI = median(GIselectivity[! names(GIselectivity) %in% drugCids]),
        drugkclust = median(kclustselectivity[names(kclustselectivity) %in% drugCids]),
        otherkclust = median(kclustselectivity[! names(kclustselectivity) %in% drugCids]),
        drugdomain = median(domainselectivity[names(domainselectivity) %in% drugCids]),
        otherdomain = median(domainselectivity[! names(domainselectivity) %in% drugCids])
    )
    result <- t(data.frame(result))
    row.names(result) <- thisDomain
    return(result)
}

write.table(results, outputFile, col.names=TRUE, row.names=TRUE)

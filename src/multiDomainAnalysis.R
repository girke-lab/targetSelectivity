#!/usr/bin/env Rscript

# (C) 2016 Tyler William H Backman
# Purpose: create a text file which lists all co-occurances for each domain that occurs in multiple proteins

library(R.utils)
library(bioassayR)
# library(reshape2)
# library(igraph)

# parse input options
databaseFile <- commandArgs(trailingOnly=TRUE)[1]
outputFile <- commandArgs(trailingOnly=TRUE)[2]

# test code for running without make:
if(is.na(commandArgs(trailingOnly=TRUE)[1])){
    # databaseFile <- "/dev/shm/bioassayDatabase.sqlite"
    databaseFile <- "working/bioassayDatabaseSingleTarget.sqlite"
    outputFile <- "working/domainComposition.tab"
}

# parse input files
db <- connectBioassayDB(databaseFile)

# Get domain list for each target
targetList <- allTargets(db)
domainsPerTargetList <- lapply(targetList, translateTargetId, database=db, category="domains")

# keep only multi-domain targets
domainQuantity <- sapply(domainsPerTargetList, length)
domainsPerTargetListMultiples <- domainsPerTargetList[domainQuantity > 1]

# sort domains alphabetically and keep only unique domain orders
domainsSorted <- lapply(domainsPerTargetListMultiples, sort)
uniqueDomainSets <- domainsSorted[! duplicated(sapply(domainsSorted, paste, collapse="_"))]

# for each domain find (1) all domain sets it occurs in, (2) all co-occurance domains
domainList <- unique(unlist(uniqueDomainSets))
domainTableResults <- lapply(domainList, function(queryDomain){
    domainOrders <- uniqueDomainSets[! is.na(sapply(uniqueDomainSets, match, x=queryDomain))]
    allCoOccurances <- unique(unlist(domainOrders))
    allCoOccurances <- allCoOccurances[allCoOccurances != queryDomain]
    allCoOccurances <- sort(allCoOccurances)

    domainOrderStrings <- sapply(domainOrders, paste, collapse="_")
    cbind(queryDomain, domainOrderStrings, paste(allCoOccurances, collapse="_"))
})
outputTable <- do.call(rbind, domainTableResults)
colnames(outputTable) <- c("Pfam_domain", "Domain_compositions", "All_cooccurrences")
write.table(outputTable, outputFile, sep="\t", row.names = FALSE)

# # old code
# domainsSorted <- sapply(domainsPerTargetListMultiples, function(domainList){
#     paste(sort(domainList), collapse="_")    
#     })
# length(unique(domainsSorted))
# 
# # make graph and find number of connected components
# domainsPerTargetTable <- melt(domainsPerTargetList)
# g <- graph_from_data_frame(domainsPerTargetTable, directed = FALSE)
# number <- components(g)$no

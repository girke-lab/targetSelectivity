#!/usr/bin/env Rscript

# (C) 2016 Tyler William H Backman
# Purpose: plot target selectivity against compounds from an individual domain 

library(R.utils)
library(bioassayR)
library(reshape2)
library(igraph)

# parse input options
databaseFile <- commandArgs(trailingOnly=TRUE)[1]
# outputFile <- commandArgs(trailingOnly=TRUE)[2]

# test code for running without make:
if(is.na(commandArgs(trailingOnly=TRUE)[1])){
    # databaseFile <- "/dev/shm/bioassayDatabase.sqlite"
    databaseFile <- "working/bioassayDatabaseSingleTarget.sqlite"
    # outputFile <- "working/multidomain.out"
}

# parse input files
db <- connectBioassayDB(databaseFile)

# count number of multi-domain proteins with identical domains
targetList <- allTargets(db)
domainsPerTargetList <- lapply(targetList, translateTargetId, database=db, category="domains")
domainQuantity <- sapply(domainsPerTargetList, length)
table(domainQuantity)
domainsPerTargetListMultiples <- domainsPerTargetList[domainQuantity > 1]

domainsSorted <- sapply(domainsPerTargetListMultiples, function(domainList){
    paste(sort(domainList), collapse="_")    
    })
length(unique(domainsSorted))

# make graph and find number of connected components
domainsPerTargetTable <- melt(domainsPerTargetList)
g <- graph_from_data_frame(domainsPerTargetTable, directed = FALSE)
number <- components(g)$no
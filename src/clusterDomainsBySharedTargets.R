#!/usr/bin/env Rscript

# (C) 2016 Tyler William H Backman
# Purpose: cluster domains based on the number of shared annotated targets

library(R.utils)
library(bioassayR)
library(reshape2)
library(Matrix)
library(ChemmineR)

# parse input options
databaseFile <- commandArgs(trailingOnly=TRUE)[1]
outputFile <- commandArgs(trailingOnly=TRUE)[2]

# test code for running without make:
if(is.na(commandArgs(trailingOnly=TRUE)[1])){
    # databaseFile <- "/dev/shm/bioassayDatabase.sqlite"
    databaseFile <- "working/bioassayDatabaseSingleTarget.sqlite"
    outputFile <- "working/domainClusters.csv"
}

# parse input files
db <- connectBioassayDB(databaseFile)

# get domains for each target
targetList <- allTargets(db)
domainList <- lapply(targetList, translateTargetId, database=db, category="domains")
names(domainList) <- targetList
domainList <- domainList[! is.na(domainList)]
domains <- melt(domainList)
colnames(domains) <- c("domain", "target")

# make domain vs target matrix
uniqueTargets <- unique(domains$target)
uniqueDomains <- unique(domains$domain)
domainVsTargetMatrix <- sparseMatrix(
    i = match(domains$target, uniqueTargets),
    j = match(domains$domain, uniqueDomains),
    x = 1,
    dims = c(length(uniqueTargets), length(uniqueDomains)),
    dimnames = list(uniqueTargets, uniqueDomains),
    symmetric = FALSE,
    index1 = TRUE)
domainVsTargetMatrix <- t(as.matrix(domainVsTargetMatrix))

# cluster with single linkage clustering 
# by joining any two domains that share at least 70%
# tanimoto similarity across there common set of protein targets
domainVsTargetFP <- as(domainVsTargetMatrix, "FPset") 
# simFun <- function(a, b, c, d){
#     c/min(c + a, c + b)   
# }
# clusters <- cmp.cluster(domainVsTargetFP, cutoff=0.7, method=simFun, quiet=TRUE)
clusters <- cmp.cluster(domainVsTargetFP, cutoff=0.7, method="Tanimoto", quiet=TRUE)
print(length(unique(clusters$CLID_0.7)))

# find a representative domain for each cluster
targetsPerDomainCount <- table(domains$domain)
domainsPerCluster <- tapply(clusters$ids, clusters$CLID_0.7, function(domainsInCluster){
   counts <- targetsPerDomainCount[domainsInCluster]
   maxCounts <- counts[counts == max(counts)]
   baseName <- sort(names(counts))[1]
   return(cbind(domainsInCluster, baseName))
})
renamedClusters <- do.call(rbind, domainsPerCluster)
colnames(renamedClusters) <- c("Domain", "Cluster")

# write out clusters
write.csv(renamedClusters, outputFile, row.names=FALSE)

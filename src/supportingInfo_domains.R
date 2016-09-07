#!/usr/bin/env Rscript

# (C) 2016 Tyler William H Backman
# Purpose: produce list of Pfam domains including median target, cluster, and domain selectivities for FDA approved and non-FDA compounds

library(R.utils)

# parse input options
targetSelectivityByDomainFile <- commandArgs(trailingOnly=TRUE)[1]
outputFilename <- commandArgs(trailingOnly=TRUE)[2]

# test code for running without make:
if(is.na(commandArgs(trailingOnly=TRUE)[1])){
    targetSelectivityByDomainFile <- "working/targetSelectivityByDomain.tab"
    outputFilename <- "working/supportingInfo_domains.tab"
}

# parse input files
targetSelectivityByDomains <- read.table(targetSelectivityByDomainFile)

colnames(targetSelectivityByDomains) <- c("totalFDA", "totalNonFDA", "medianFDAtargetSelectivity",
    "medianNonFDAtargetSelectivity",       "medianFDAclusterSelectivity", "medianNonFDAclusterSelectivity",
    "medianFDAdomainSelectivity", "medianNonFDAdomainSelectivity")
targetSelectivityByDomains <- targetSelectivityByDomains[rowSums(targetSelectivityByDomains[,1:2]) > 0,]

# write out table
write.table(targetSelectivityByDomains, outputFilename, quote=F, sep="\t", row.names=T, col.names = T)

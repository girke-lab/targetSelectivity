#!/usr/bin/env Rscript

# (C) 2016 Tyler William H Backman
# Purpose: 
# collapse pfam GO annotations to GO Slim

library(R.utils)
library(GSEABase)
library(plyr)

# parse input options
pfam2goFile <- commandArgs(trailingOnly=TRUE)[1]
goslimFile <- commandArgs(trailingOnly=TRUE)[2]
outputFile <- commandArgs(trailingOnly=TRUE)[3]

# test code for running without make:
if(is.na(commandArgs(trailingOnly=TRUE)[1])){
    pfam2goFile <- "working/pfam2go"
    goslimFile <- "working/goslim_generic.obo"
    outputFile <- "working/domainGOslimAnnotations.csv"
}

# parse input options
pfam2go <- readLines(pfam2goFile)
slim <- getOBOCollection(goslimFile)

# parse pfam2go
pfam2go <- pfam2go[! grepl("^!", pfam2go)]
pfams <- gsub(".*Pfam:(\\S+).*", "\\1", pfam2go)
goTerms <- gsub(".*(GO:\\S+).*", "\\1", pfam2go)

myCollection <- GOCollection(unique(pfams))

# break apart go terms by domain
goDomains <- split(goTerms, pfams)

# get goSlim terms for each domain
# note: try is used because single terms missing from the goslim
#   return an error, and should be skipped
slimTerms <- lapply(goDomains, function(x){
    result <- tryCatch({
        goSlimResults <- goSlim(GOCollection(as.vector(x)), slim, "MF")
        goSlimResults <- goSlimResults[goSlimResults$Count > 0,,drop=FALSE]
        return(cbind(row.names(goSlimResults), as.vector(goSlimResults$Term)))
    }, error = function(e){ return(NA) })
})
slimTerms <- slimTerms[! is.na(slimTerms)]

# reassemble matrix
results <- ldply(slimTerms)
colnames(results) <- c("domain","go_id","go_name")

# save output
write.csv(results, outputFile, row.names=FALSE)

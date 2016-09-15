# (C) 2016 Tyler W. H. Backman
# Purpose: annotate kclust clusters via UniProt
# Use: UniProtAnnotator(UniProtIDList)

library(RCurl)
library(utils)
.uniprotserverURL <- "http://www.uniprot.org/uniprot/"

.UniProtDownload <- function(accession, columns){
    columns <- URLencode(paste(columns, collapse=","))
    response <- getURL(paste(.uniprotserverURL,
                               "?query=accession%3A",
                               accession,
                               "&columns=",
                               columns,
                               "&format=tab",
                               sep=""))[[1]]
    if(response == "")
        return(NA)
    responsetable <- read.delim(text=response, stringsAsFactors=F, check.names=FALSE)
    responsetable <- cbind(accession=accession, responsetable, stringsAsFactors=FALSE)
    return(responsetable)
}

UniProtAnnotator <- function(accessions, useBiomartNames=FALSE, columns=c("genes", "organism", "protein names")){
    results <- lapply(accessions, .UniProtDownload, columns=columns)
    if(sum(is.na(results)) > 0){
        warning(
            paste(
                "no UniProt annotations for:",
                paste(accessions[is.na(results)], collapse=" "),
                collapse=" "
            )
        )
        results <- results[! is.na(results)]
    }
    results <- do.call(rbind, results)
    if(useBiomartNames){
        # "accession", "gene_name", "organism", "protein_name"
        colnames(results)[colnames(results) == "Entry"] <- "accession"
        colnames(results)[colnames(results) == "Gene names"] <- "gene_name"
        colnames(results)[colnames(results) == "Organism"] <- "organism"
        colnames(results)[colnames(results) == "Protein names"] <- "protein_name"
    }
    return(results)
}

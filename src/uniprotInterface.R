# (C) 2016 Tyler W. H. Backman
# Purpose: annotate kclust clusters via UniProt
# Use: UniProtAnnotator(UniProtIDList)

library(RCurl)
.uniprotserverURL <- "http://www.uniprot.org/uniprot/"

.UniProtDownload <- function(accession){
    response <- getURL(paste(.uniprotserverURL,
                               "?query=accession%3A",
                               accession,
                               "&format=tab"
                               ,sep=""))[[1]]
    if(response == "")
        return(NA)
    lines <- strsplit(response, "\n")[[1]]
    header <- strsplit(lines[1], "\t")[[1]]
    values <- strsplit(lines[2], "\t")[[1]]
    return(t(data.frame(cbind(header, values), row.names=1)))
}

UniProtAnnotator <- function(accessions, useBiomartNames=FALSE){
    results <- lapply(accessions, .UniProtDownload)
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
    results <- as.data.frame(results, row.names=F, stringsAsFactors = F)
    if(useBiomartNames){
        # "accession", "gene_name", "organism", "protein_name"
        colnames(results)[colnames(results) == "Entry"] <- "accession"
        colnames(results)[colnames(results) == "Gene names"] <- "gene_name"
        colnames(results)[colnames(results) == "Organism"] <- "organism"
        colnames(results)[colnames(results) == "Protein names"] <- "protein_name"
    }
    return(results)
}

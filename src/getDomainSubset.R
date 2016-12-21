# (C) 2016 Tyler William H Backman 
# Purpose: takes a list of domains and returns only a representative domain
#   for each co-occurance set

library(igraph)

# read in domain co-occurance table
ds_domainTable <- read.table("working/domainComposition.tab", header=TRUE, stringsAsFactors=FALSE)
ds_domainTable <- ds_domainTable[! duplicated(ds_domainTable[,1]),c(1,3)]

# get domain frequencies
ds_databaseFile <- "working/bioassayDatabaseSingleTarget.sqlite"
ds_db <- connectBioassayDB(ds_databaseFile)
ds_targetList <- allTargets(ds_db)
ds_domainsPerTargetList <- lapply(ds_targetList, translateTargetId, database=db, category="domains")
ds_domainFrequency <- table(unlist(ds_domainsPerTargetList))

# # create a random sample for testing
# queryDomainSet <- sample(unique(domainTable[,1]), 100)
# getUniqueDomainSet(queryDomainSet)

getUniqueDomainSet <- function(queryDomainSet){
    # create edge list with just those edges in queryDomainSet
    domainTableSubset <- ds_domainTable[ds_domainTable[,1] %in% queryDomainSet,]
    edgeList <- apply(domainTableSubset, 1, function(x){
        domain <- x[[1]]
        cooccurrences <- strsplit(x[2], "_")[[1]]
        cooccurrences <- cooccurrences[cooccurrences %in% queryDomainSet]
        if(length(cooccurrences) == 0)
            return(NA)
        cbind(domain, cooccurrences)
    })
    edgeList <- edgeList[! is.na(edgeList)]
    edgeList <- do.call(rbind, edgeList)
    
    # create graph and remove duplicate edges
    g <- graph_from_data_frame(edgeList, directed = FALSE)
    g <- igraph::simplify(g)
    
    # iteratively find maximal cliques
    myComponents <- c()
    # while(length(largest_cliques(g)[[1]]$name) > 1){
    while(vcount(g) > 1){
        largest <- largest_cliques(g)[[1]]
        myComponents <- c(myComponents, list(largest$name))
        g <- delete_vertices(g, largest)
    }
    
    # plot.igraph(g,vertex.size=5, vertex.label.cex=0.2)
    
    # choose for each connected component, the domain which occurs on the most targets
    largestDomains <- sapply(myComponents, function(domainsWithinComponent){
        # domainSizes <- pfstat[pfstat[,1] %in% domainsWithinComponent,]
        # domainSizes[which.max(as.integer(domainSizes[,2])),1]
        names(which.max(ds_domainFrequency[names(ds_domainFrequency) %in% domainsWithinComponent]))
    })
    
    domainsNotInGraph <- queryDomainSet[! queryDomainSet %in% unlist(myComponents)]
    return(c(domainsNotInGraph, largestDomains))
}
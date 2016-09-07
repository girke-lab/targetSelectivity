#!/usr/bin/env Rscript

# (C) 2016 Tyler William H Backman
# Purpose: make heatmap of GO terms vs biclusters
library(R.utils)
library(ggplot2)
library(reshape)
library(GOstats)
library(GSEABase)

# parse input options
perClusterAnalysisFolder <- commandArgs(trailingOnly=TRUE)[1] 
clusterGOannotationsFile <- commandArgs(trailingOnly=TRUE)[2]
clusterAnnotationFile <- commandArgs(trailingOnly=TRUE)[3]
drugComparisonMatrixFile <- commandArgs(trailingOnly=TRUE)[4]
outputFile <- commandArgs(trailingOnly=TRUE)[5]

# test code for running without make:
if(is.null(commandArgs(trailingOnly=TRUE)[1])){
    perClusterAnalysisFolder <- "working/perClusterAnalysis" 
    clusterGOannotationsFile <- "working/clusterGOslimAnnotations.csv"
    clusterAnnotationFile <- "working/clusterAnnotations.csv"
    drugComparisonMatrixFile <- "working/drugComparisonMatrix.RData"
    outputFile <- "working/biclusterGoScatter.pdf"
}

# parse input files
clusterGOannotations <- read.csv(clusterGOannotationsFile)
clusterAnnotations <- read.csv(clusterAnnotationFile)
load(drugComparisonMatrixFile) # loads drugComparisonMatrix object

# loop over biclusters and compile target list
allFiles <- list.files(perClusterAnalysisFolder, pattern="*.RData", full.names=T)
clusterNums <- as.numeric(gsub("^.*?cluster(\\d+).RData$", "\\1", allFiles))
clusterNums <- sort(clusterNums)
uniprotTargets <- lapply(clusterNums, function(cluster){
    load(file.path(perClusterAnalysisFolder, paste("cluster", cluster, ".RData", sep=""))) 
    return(unique(as.character(uniProtData$accession[! is.na(uniProtData$accession)])))
})

# get list of background compounds (all targets in drugComparisonMatrix)
proteinClusterIds <- colnames(drugComparisonMatrix)
allTargetUniprot <- clusterAnnotations[match(proteinClusterIds, clusterAnnotations$uniqueClusterIds),]
allTargetIds <- as.character(unique(allTargetUniprot$accession[! is.na(allTargetUniprot$accession)]))

# get complete GO annotations for background
backgroundGO <- clusterGOannotations[clusterGOannotations$accession %in% allTargetIds,]
goframeData <- data.frame(frame.go_id=backgroundGO$go_id, frame.Evidence="ISA", frame.gene_id=backgroundGO$accession)
goFrame <- GOFrame(goframeData, organism="Other")
goAllFrame <- GOAllFrame(goFrame)
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())

# perform GO enrichment analysis for each bicluster
goResults <- lapply(uniprotTargets, function(genes){
    params <- GSEAGOHyperGParams(name="My Custom GSEA based annot Params",
                                 geneSetCollection=gsc,
                                 geneIds = genes,
                                 universeGeneIds = allTargetIds,
                                 ontology = "MF",
                                 pvalueCutoff = 0.05,
                                 conditional = FALSE,
                                 testDirection = "over")
    over <- hyperGTest(params)
    return(summary(over))
})

# assemble matrix
goids <- unique(unlist(sapply(goResults, function(x) x$GOMFID)))
save(list=c("goids"), file="working/biclustergoids.RData")
goNames <- unique(unlist(sapply(goResults, function(x) x$Term)))
goMatrix <- matrix(data=NA, nrow=length(clusterNums), ncol=length(goNames), dimnames=list(as.character(clusterNums), goNames))
for(i in 1:length(clusterNums)){ 
    goMatrix[as.character(clusterNums[i]),goResults[[i]]$Term] <- goResults[[i]]$Pvalue
}

# sort go terms by abundance 
# colValues <- apply(goMatrix, 2, function(x) sum(!is.na(x)))
colValues <- apply(goMatrix, 2, function(x) min(x[!is.na(x)]))
abundantTerms <- sort(colValues, decreasing=T)
goMatrix <- goMatrix[,names(abundantTerms)]

# shorten long names
# colnames(goMatrix) <- gsub("^(.{40}).*$", "\\1...", colnames(goMatrix))

# plot matrix
colnames(goMatrix) <-  gsub("\\W?activity", "", colnames(goMatrix))
colnames(goMatrix)[colnames(goMatrix) == 
                       "hydrolase, acting on acid anhydrides, in phosphorus-containing anhydrides"] <- "hydrolase, phosphorus anhydrides"
colnames(goMatrix)[colnames(goMatrix) == 
                       "hydrolase, acting on carbon-nitrogen (but not peptide) bonds"] <- "hydrolase, carbon-nitrogen bonds"
colnames(goMatrix)[colnames(goMatrix) == 
                       "transferase, transferring phosphorus-containing groups"] <- "transferase, phosphorus-containing"
rownames(goMatrix) <- as.character(1:length(rownames(goMatrix)))
goTable <- melt(goMatrix)
colnames(goTable) <- c("Bicluster", "GOterm", "Fraction")
goTable$Bicluster <- factor(as.character(goTable$Bicluster), levels=rownames(goMatrix))
goTable$GOterm <- factor(as.character(goTable$GOterm), levels=as.character(colnames(goMatrix)))
goTable <- goTable[! is.na(goTable$Fraction),]
p <- ggplot(goTable, aes(Bicluster, GOterm)) + 
    geom_point(aes(color = Fraction), size=4) +
    scale_colour_gradient(trans = "log10") +
    scale_x_discrete(drop=FALSE) +
    theme( 
        text = element_text(size=15),
        # axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) +
    labs(
        y="Molecular Function GO Slim Term",
        x="Compound-Target Bicluster",
        # title="GO Term Distribution Among Biclusters",
        colour="p-value")
plot(p)
ggsave(outputFile, plot=p, device="pdf", width=9.5, height=7)

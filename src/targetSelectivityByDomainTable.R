#!/usr/bin/env Rscript

# (C) 2016 Tyler William H Backman
# Purpose: compute target selectivity on a per-domain basis 

library(R.utils)
library(xtable)
library(bioassayR)
library(GOstats)
library(GSEABase)
library(reshape)
library(ggplot2)
library(gridExtra)

# parse input options
targetSelectivityByDomainFile <- commandArgs(trailingOnly=TRUE)[1]
pfamClansFile <- commandArgs(trailingOnly=TRUE)[2]
databaseFile <- commandArgs(trailingOnly=TRUE)[3]
pfam2goslimFile <- commandArgs(trailingOnly=TRUE)[4]
humanDomainsTwoCols <- commandArgs(trailingOnly=TRUE)[5]
targetGOslimAnnotationsFile <- commandArgs(trailingOnly=TRUE)[6]
PfamResidueLengthsFile <- commandArgs(trailingOnly=TRUE)[7]
outputFile <- commandArgs(trailingOnly=TRUE)[8]

# test code for running without make:
if(is.null(commandArgs(trailingOnly=TRUE)[1])){
    targetSelectivityByDomainFile <- "working/targetSelectivityByDomain.tab"
    pfamClansFile <- "working/Pfam-A.clans.tsv"
    databaseFile <- "working/bioassayDatabaseSingleTarget.sqlite"
    pfam2goslimFile <- "working/domainGOslimAnnotations.csv"
    humanDomainsTwoCols <- "working/humanDomainsTwoCols"
    targetGOslimAnnotationsFile <- "working/targetGOslimAnnotations.csv"
    PfamResidueLengthsFile <- "working/PfamResidueLengths.tab"
    outputFile <- "working/targetSelectivityByDomainTable.tex"
}

# parse input files
statsTable <- read.table(targetSelectivityByDomainFile)
pfamClans <- read.delim(pfamClansFile, header=FALSE, na.strings="\\N")
database <- connectBioassayDB(databaseFile)
pfam2goslim <- read.csv(pfam2goslimFile)
humanDomains <- read.table(humanDomainsTwoCols, header = FALSE)
targetGOslimAnnotations <- read.csv(targetGOslimAnnotationsFile)
PfamResidueLengths <- read.table(PfamResidueLengthsFile, stringsAsFactors=F, header=T)

# get list of domains with over 99 residues
longDomains <- PfamResidueLengths[PfamResidueLengths[,2] > 99,1]

# clean up human
humanDomains[,1] <- gsub("^(PF\\d*).*", "\\1", humanDomains[,1], perl=TRUE)
colnames(humanDomains) <- c("DOMAIN", "TARGET")
humanDomainsAll <- unique(humanDomains$DOMAIN)

# keep only domains with at least 10 active drugs and other compounds
statsTable <- statsTable[statsTable$totalDrugs > 9,]
statsTable <- statsTable[statsTable$totalOther > 9,]

humanStatsTable <- statsTable[row.names(statsTable) %in% humanDomainsAll,]
humanStatsTable <- humanStatsTable[row.names(humanStatsTable) %in% longDomains,]
topTargets <- humanStatsTable[order(humanStatsTable$drugdomain, decreasing=T),]
topTargets <- topTargets[topTargets$drugdomain > 16.5,]
bottomTargets <- humanStatsTable[order(humanStatsTable$drugdomain, decreasing=F),]
bottomTargets <- bottomTargets[bottomTargets$drugdomain < 6.5,]

mergedTable <- rbind(topTargets, bottomTargets)
mergedTable <- round(mergedTable, 2)

slashTable <- cbind(
      `Domain`=rownames(mergedTable),
      `Individual Targets`=paste(mergedTable$drugGI, mergedTable$otherGI, sep="/"),
      `Distinct Sequence Targets`=paste(mergedTable$drugkclust, mergedTable$otherkclust, sep="/"),
      `Distinct Domain Targets`=paste(mergedTable$drugdomain, mergedTable$otherdomain, sep="/")
)

newNames <- as.character(pfamClans[match(slashTable[,1], pfamClans[,1]),5])
newNames[is.na(newNames)] <- ""
slashTable[,1] <- paste(slashTable[,1], newNames)
# slashTable[,1] <- gsub("^(.{40}).*$", "\\1...", slashTable[,1])

# remove DUFs (domain of unknown function)
domainList <- as.character(unique(slashTable[,1]))
dufList <- domainList[grepl("DUF\\d{4}", domainList)]
slashTable <- slashTable[! slashTable[,1] %in% dufList,]

# keep only one co-occurance for each domain
source("src/getDomainSubset.R")
domainList <- as.character(unique(slashTable[,1]))
domainList <- gsub("\\W.*", "", domainList)
uniques <- getUniqueDomainSet(domainList)
slashTable <- slashTable[gsub("\\W.*", "", slashTable[,1]) %in% uniques,]

xtmp <- xtable(slashTable, caption="Selectivity by domain", label="selectivityByDomain")
print(xtmp, type="latex", file=outputFile, include.rownames=F)

# get reasonable breakpoints for both sets
allValues <- c(statsTable$drugdomain, statsTable$otherdomain)
cuts <- factor(cut(allValues, breaks=8))
ranges <- tapply(allValues, cuts, range)
breakPoints <- sapply(ranges, function(x) x[2])
breakPoints <- c(0, breakPoints)

# make domain selectivity frequency table for FDA approved drugs
x <- statsTable$drugdomain
drugfactor <- factor(cut(x, breaks=breakPoints))
ranges <- tapply(statsTable$drugdomain, drugfactor, range)
ranges <- sapply(ranges, paste, collapse="-")
levels(drugfactor) <- ranges
xoutDrugs <- as.data.frame(table(drugfactor))

# make domain selectivity frequency table for Other compounds
x <- statsTable$otherdomain
otherfactor <- factor(cut(x, breaks=breakPoints))
ranges <- tapply(statsTable$otherdomain, otherfactor, range)
ranges <- sapply(ranges, paste, collapse="-")
levels(otherfactor) <- ranges
xoutOther <- as.data.frame(table(otherfactor))

# merge the two
totalrows <- max(c(nrow(xoutOther), nrow(xoutDrugs)))
emptyMatrix <- matrix(data=0, nrow=totalrows, ncol=2)
emptyMatrix[1:nrow(xoutDrugs),1] <- xoutDrugs[,2]
emptyMatrix[1:nrow(xoutOther),2] <- xoutOther[,2]
minRangesDrug <- rep(1000, totalrows)
minRangesDrug[1:nrow(xoutDrugs)] <- sapply(strsplit(as.character(xoutDrugs[,1]), "-"), function(x) min(as.numeric(x)))
maxRangesDrug <- rep(0, totalrows)
maxRangesDrug[1:nrow(xoutDrugs)] <- sapply(strsplit(as.character(xoutDrugs[,1]), "-"), function(x) max(as.numeric(x)))
minRangesNonDrug <- rep(1000, totalrows)
minRangesNonDrug[1:nrow(xoutOther)] <- sapply(strsplit(as.character(xoutOther[,1]), "-"), function(x) min(as.numeric(x)))
maxRangesNonDrug <- rep(0, totalrows)
maxRangesNonDrug[1:nrow(xoutOther)] <- sapply(strsplit(as.character(xoutOther[,1]), "-"), function(x) min(as.numeric(x)))
row.names(emptyMatrix)  <- paste(mapply(min, minRangesDrug, minRangesNonDrug), mapply(max, maxRangesDrug, maxRangesNonDrug), sep="-")

# translate bin names for later
levels(drugfactor) <- row.names(emptyMatrix)[1:length(levels(drugfactor))]
levels(otherfactor) <- row.names(emptyMatrix)[1:length(levels(otherfactor))]

# xout <- transform(xout, cumFreq = cumsum(Freq), relative = prop.table(Freq))
xtmp <- xtable(data.frame(emptyMatrix), caption="selectivity frequency by domain", label="selectivityByDomainFreq")
digits(xtmp) <- c(0,0,0)
print(xtmp, type="latex", file="working/targetSelectivityByDomainFreq.tex", include.rownames=T)

###################################
# perform GO enrichment analysis on highest and lowest selectivity domains
###################################

# get list of background targets (all targets in all bins)
allDomains <- row.names(statsTable)

# get complete GO annotations for background
backgroundGO <- pfam2goslim[pfam2goslim$domain %in% allDomains,]
goframeData <- data.frame(frame.go_id=backgroundGO$go_id, frame.Evidence="ISA", frame.gene_id=backgroundGO$domain)
goFrame <- GOFrame(goframeData, organism="Other")
goAllFrame <- GOAllFrame(goFrame)
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())

backgroundDomainList <- as.character(unique(backgroundGO$domain))

# get complete GO annotations for all targets
termAbundance <- table(goAllFrame@data[,"go_id"])
termNames <- Term(names(termAbundance))

# perform GO enrichment analysis for each bicluster
goResultFn <- function(domains){
    annotatedDomains <- intersect(domains, backgroundDomainList)
    if(length(annotatedDomains) < 2)
        return(NA)
    params <- GSEAGOHyperGParams(name="My Custom GSEA based annot Params",
                                 geneSetCollection=gsc,
                                 geneIds = annotatedDomains,
                                 universeGeneIds = backgroundDomainList,
                                 ontology = "MF",
                                 pvalueCutoff = 0.05,
                                 conditional = FALSE,
                                 testDirection = "over")
    over <- hyperGTest(params)
    return(summary(over))
}

goDomainsFDA <- split(row.names(statsTable), drugfactor)
goResultsFDA <-  lapply(goDomainsFDA, goResultFn)
goResultsFDA <- goResultsFDA[! is.na(goResultsFDA)]
goResultsFDA <- goResultsFDA[sapply(goResultsFDA, nrow) > 0]
goResultsFDATable <- ldply(goResultsFDA)

goDomainsnonFDA <- split(row.names(statsTable), otherfactor)
goResultsNonFDA <- lapply(goDomainsnonFDA, goResultFn)
goResultsNonFDA <- goResultsNonFDA[! is.na(goResultsNonFDA)]
goResultsNonFDA <- goResultsNonFDA[sapply(goResultsNonFDA, nrow) > 0]
goResultsnonFDATable <- ldply(goResultsNonFDA)

# merge and assemble
goTable <- rbind(cbind(goResultsFDATable, class="FDA Approved Compounds"), cbind(goResultsnonFDATable, class="Non-FDA Compounds"))
# colnames(goTable) <- c("DomainSelectivity", "GOterm", "Fraction", "class")
termOrder <- order(goTable$Pvalue, decreasing=T)
goTable <- goTable[termOrder,]

# fix go term names
goTable$Term <-  gsub("\\W?activity", "", goTable$Term)
goTable$Term <-  gsub("transferring\\W?", "", goTable$Term)
goTable$Term[goTable$Term ==
                   "hydrolase, acting on acid anhydrides, in phosphorus-containing anhydrides"] <- "hydrolase, phosphorus anhydrides"
goTable$Term[goTable$Term ==
                   "hydrolase, acting on carbon-nitrogen (but not peptide) bonds"] <- "hydrolase, carbon-nitrogen bonds"
goTable$Term[goTable$Term ==
                   "transferase, transferring phosphorus-containing groups"] <- "transferase, phosphorus-containing"
# goTable$GOterm <- shortNames[match(goTable$GOterm, names(shortNames))]
# termLabels <- rev(as.character(unique(rev(goTable$GOterm))))
# termLabels <- c(shortNames[! shortNames %in% termLabels], rev(as.character(unique(rev(goTable$GOterm)))))
# termLabels <- termLabels[termLabels != "molecular_function"]
# goTable$GOterm <- factor(goTable$GOterm, levels=termLabels)
# replace terms with real names
# goTable$GOterm <- termNames[match(goTable$GOterm, names(termNames))]
# goTable <- goTable[goTable$GOterm != "molecular_function",]
# goTable$GOterm <- factor(as.character(goTable$GOterm), levels=rev(as.character(unique(rev(goTable$GOterm)))))
# goTable$DomainSelectivity <- factor(as.character(goTable$DomainSelectivity), levels=row.names(emptyMatrix))
levellist <- levels(drugfactor)[levels(drugfactor) %in% as.character(goTable$.id)]
goTable$.id <- factor(as.character(goTable$.id), levels=levellist)
goTable$Term <- factor(as.character(goTable$Term), levels=unique(goTable$Term))

p <- ggplot(goTable, aes(.id, Term)) + 
    geom_point(aes(color = Pvalue), size=4, na.rm=T) +
    facet_grid(. ~ class) +
    # scale_colour_gradient(trans = "log10") +
    scale_x_discrete(drop=FALSE) +
    scale_y_discrete(drop=FALSE) +
    theme( 
        axis.text.x=element_text(angle = -45, hjust = 0),
        text = element_text(size=12),
        legend.position="bottom",
        legend.key.width = unit(2,"cm"),
        # axis.text.x = element_blank(), 
        plot.margin=unit(c(0,1,0,0),"cm"),
        axis.ticks.x = element_blank()) +
    labs(
        y="Molecular Function GO Slim Term",
        x="Pfam Domains Binned by Median Domain\nSelectivity of Active Compounds",
        colour="p-value")
plot(p)
ggsave("working/targetSelectivityByDomain.pdf", plot=p, device="pdf", width=8, height=5)

########################################################
# make plot of slim term abundance in data
########################################################

# get complete GO annotations for all targets
goframeData <- data.frame(frame.go_id=targetGOslimAnnotations$go_id, frame.Evidence="ISA", frame.gene_id=targetGOslimAnnotations$accession)
goFrame <- GOFrame(goframeData, organism="Other")
goAllFrame <- GOAllFrame(goFrame)
termAbundance <- table(goAllFrame@data[,"go_id"])
termNames <- Term(names(termAbundance))

termCountsMelted <- melt(termAbundance)
colnames(termCountsMelted) <- c("GOTerm", "value")
termCountsMelted <- termCountsMelted[match(as.character(goTable$GOMFID), as.character(termCountsMelted$GOTerm)),]
termCountsMelted$GOTerm <- factor(as.character(termCountsMelted$GOTerm), levels=unique(goTable$GOMFID))

# termAbundance <- termAbundance[match(names(abundantTerms), names(termAbundance))]
# termAbundance <- termAbundance[Term(names(termAbundance)) != "molecular_function"]
# termCountsMelted <- melt(termAbundance)
# colnames(termCountsMelted) <- c("GOTerm", "value")
# termCountsMelted <- cbind(termCountsMelted, variable="targets")
# termCountsMelted <- termCountsMelted[!is.na(termCountsMelted$GOTerm),]
# termCountsMelted$GOTerm <- shortNames[match(termCountsMelted$GOTerm, names(shortNames))]
# termCountsMelted$GOTerm <- factor(termCountsMelted$GOTerm, levels=termLabels)

p2 <- ggplot(data=termCountsMelted, aes(x=value, y=GOTerm)) +
    geom_point(size=4) +
    # facet_grid(. ~ variable, scales="free") + 
    guides(colour=guide_legend()) +
    scale_x_continuous(labels=function(n){format(n, scientific = FALSE)}) + 
    # ylab("Molecular Function GO Slim Term") + 
    xlab("Targets in\nPubChem BioAssay") +
    # ggtitle("PubChem Compound GO Slim Term Activity") +
    theme(panel.grid.major.y = element_line(colour = "grey"), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          legend.position="bottom") +
    scale_colour_brewer(palette="Dark2", breaks=c("active", "inactive"))
plot(p2)
ggsave("working/targetSelectivityByDomain_right.pdf", plot=p2, device="pdf", width=2, height=9)

library(gtable)
g <- ggplotGrob(p)
g <- gtable_add_cols(g, unit(4,"cm"))
g <- gtable_add_grob(g, ggplotGrob(p2),
                     t = 1, l=ncol(g), b=4, r=ncol(g))
plot(g)
ggsave("working/targetSelectivityByDomain_both.pdf", plot=g, device="pdf", width=8, height=4.5)


# test this result by working backwards
drugLinksFile <- "working/drugbank_links.csv"
drugLinks <- read.csv(drugLinksFile)
drugCids <- unique(drugLinks$PubChem.Compound.ID)
drugCids <- drugCids[! is.na(drugCids)]
uniProtIds <- unique(as.character(clusterGOannotations[grep("kinase", clusterGOannotations$go_name),1]))
giList <- unique(unlist(lapply(uniProtIds, translateTargetId, database=database, category="GI", fromCategory="UniProt")))
giList <- unique(giList[!is.na(giList)])

activeCids <- lapply(giList, function(x){
    row.names(activeAgainst(database, x))
})
activeCidsUnlisted <- unique(unlist(activeCids))
activeCidsUnlisted <- intersect(activeCidsUnlisted, drugCids)
ts <- targetSelectivity(database, activeCidsUnlisted, category="domains")
median(ts)
# ligase: 16
# structural molecule: 9
# kinase: 11 

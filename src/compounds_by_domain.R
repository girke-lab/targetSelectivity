#!/usr/bin/env Rscript

# (C) 2016 Tyler William H Backman
# Purpose: plot quantity of compounds targeting GO Slim categories

library(R.utils)
library(ggplot2)
library(Matrix)
library(reshape)
library(plyr)
library(scales)
library(GO.db)
library(GSEABase)
library(bioassayR)
library(foreach)
library(doMC)
# library(cowplot)

# parse input options
cidsVStargetFile <- commandArgs(trailingOnly=TRUE)[1]
drugbank_linksFile <- commandArgs(trailingOnly=TRUE)[2]
databaseFile <- commandArgs(trailingOnly=TRUE)[3]
pfamClansFile <- commandArgs(trailingOnly=TRUE)[4]
humanDomainsTwoCols <- commandArgs(trailingOnly=TRUE)[5]
PfamResidueLengthsFile <- commandArgs(trailingOnly=TRUE)[6]
cores <- commandArgs(trailingOnly=TRUE)[7]
outputFile <- commandArgs(trailingOnly=TRUE)[8]

# test code for running without make:
if(is.null(commandArgs(trailingOnly=TRUE)[1])){
    cidsVStargetFile <- "working/cidsVStargets.RData"
    drugbank_linksFile <- "working/drugbank_links.csv"
    databaseFile <- "~/Downloads/pubchem_protein_only.sqlite"
    pfamClansFile <- "working/Pfam-A.clans.tsv"
    humanDomainsTwoCols <- "working/humanDomainsTwoCols"
    PfamResidueLengthsFile <- "working/PfamResidueLengths.tab"
    cores <- 2
    outputFile <- "working/compoundsByDomain.pdf"
}

# parse input files
load(cidsVStargetFile)
compoundVsTargetMatrix <- results # rows = targets, cols = drugs
drugbank_links <- read.csv(drugbank_linksFile)
drugCids <- drugbank_links$PubChem.Compound.ID
drugCids <- unique(drugCids[!is.na(drugCids)])
database <- connectBioassayDB(databaseFile)
pfamClans <- read.delim(pfamClansFile, header=FALSE, na.strings="\\N")
humanDomains <- read.table(humanDomainsTwoCols, header = FALSE)
PfamResidueLengths <- read.table(PfamResidueLengthsFile, stringsAsFactors=F, header=T)

# clean up human
humanDomains[,1] <- gsub("^(PF\\d*).*", "\\1", humanDomains[,1], perl=TRUE)
colnames(humanDomains) <- c("DOMAIN", "TARGET")
humanDomainAbundance <- table(humanDomains$DOMAIN)
allHumanDomains <- unique(humanDomains$DOMAIN)

# get all target domains
targetList <- allTargets(database)
domainList <- lapply(targetList, translateTargetId, database=database, category="domains")
domainList <- unique(unlist(domainList))
domainList <- domainList[! is.na(domainList)]

# iterate over domains and sum compounds targeting each
gc()
registerDoMC(cores=cores)
domainCounts <- foreach(domain=domainList, .combine='rbind', .packages=c("Matrix")) %dopar% {
    print(domain)
    emptyResult <- rep(NA, 5)  
    # translate domain into targets
    GIs <- translateTargetId(database, domain, category="GI", fromCategory="domains")
    if(is.na(GIs))
        return(emptyResult)

    matrixSubset <- compoundVsTargetMatrix[row.names(compoundVsTargetMatrix) %in% GIs,,drop=F]
    if(nrow(matrixSubset) == 0)
        return(emptyResult)
    # matrixSubset <- as.matrix(matrixSubset)
    
    targetingCidCountActives <- sum(Matrix::colSums(matrixSubset == 2) > 0) # active compound sum
    targetingCidCountInactives <- sum((Matrix::colSums(matrixSubset == 2) == 0) & (Matrix::colSums(matrixSubset == 1) > 0)) # inactive compound sum
    drugCidCountActives <- sum(Matrix::colSums(matrixSubset[,colnames(matrixSubset) %in% drugCids,drop=F] == 2) > 0) # active drugs
    drugCidCountInactives <- sum((Matrix::colSums(matrixSubset[,colnames(matrixSubset) %in% drugCids,drop=F] == 2) == 0) 
                                 & (Matrix::colSums(matrixSubset[,colnames(matrixSubset) %in% drugCids,drop=F] == 1) > 0)) # inactive drugs
    nondrugCidCountActives <- targetingCidCountActives - drugCidCountActives # active nondrugs
    nondrugCidCountInactives <- targetingCidCountInactives - drugCidCountInactives # inactive nondrugs
    totalTargets <- nrow(matrixSubset)
    return(cbind(drugCidCountActives, targetingCidCountActives, drugCidCountInactives, targetingCidCountInactives, totalTargets))
}
domainCounts <- as.data.frame(domainCounts)
row.names(domainCounts) <- domainList
save(list = ls(all=TRUE), file = gsub("^(.*).pdf$", "\\1.RData", outputFile))
# load("working/compoundsByDomain.RData")
allDomainCounts <- domainCounts

# add human proteome abundance column
humanDomainAbundanceCol <- humanDomainAbundance[match(row.names(domainCounts), names(humanDomainAbundance))]
humanDomainAbundanceCol[is.na(humanDomainAbundanceCol)] <- 0
names(humanDomainAbundanceCol) <- row.names(domainCounts)
domainCounts <- cbind(domainCounts, humanDomainAbundanceCol)

# remove non-human domains
domainCounts <- domainCounts[domainCounts$humanDomainAbundanceCol > 0,]

# remove domains with under 100 residues
longDomains <- PfamResidueLengths[PfamResidueLengths[,2] > 99,1]
domainCounts <- domainCounts[row.names(domainCounts) %in% longDomains,]

# replace domains with domain names
savedDomainNames <- row.names(domainCounts)
newNames <- as.character(pfamClans[match(savedDomainNames, pfamClans[,1]),5])
newNames[is.na(newNames)] <- ""
row.names(domainCounts) <- paste(savedDomainNames, newNames)
row.names(domainCounts) <- gsub("^\\s*(.*?)\\s*$", "\\1", row.names(domainCounts))

# combine into ggPlot style data.frame 
domainCountsNames <- name_rows(domainCounts)
domainCountsNames <- domainCountsNames[,! colnames(domainCountsNames) %in% "Var1"]
names(domainCountsNames) <- c("active FDA Approved Drugs", "active Non-FDA Compounds", "inactive FDA Approved Drugs", "inactive Non-FDA Compounds", "protein Protein Targets", "human Protein Targets", "GOdomain")
domainCountsNames$GOdomain <- gsub("\\W?domain", "", domainCountsNames$GOdomain)
domainCountsMelted <- melt(domainCountsNames)
domainCountsMelted <- domainCountsMelted[nrow(domainCountsMelted):1,]

# compute order, decreasing enrichment of activity in FDA approved drugs
enrichment <- domainCountsNames$`active FDA Approved Drugs`
enrichmentIndex <- sort(enrichment, decreasing=T, index.return=T)$ix
domainOrder <- domainCountsNames$GOdomain[enrichmentIndex]
domainCountsMelted$GOdomain <- factor(domainCountsMelted$GOdomain, levels=rev(domainOrder))

# split active and fda approved into two columns
domainCountsMelted <- cbind(domainCountsMelted, Activity=domainCountsMelted$variable)
domainCountsMelted$variable <- gsub("^\\w+\\W", "", domainCountsMelted$variable)
domainCountsMelted$variable <- factor(domainCountsMelted$variable, levels=c("FDA Approved Drugs", "Non-FDA Compounds", "Protein Targets"))
domainCountsMelted$Activity <- gsub("^(\\w+)\\W.*", "\\1", domainCountsMelted$Activity)
domainCountsMelted$Activity <- factor(domainCountsMelted$Activity)
domainCountsMelted <- domainCountsMelted[! ((domainCountsMelted$value == 0) & (domainCountsMelted$Activity == "active")),]

# keep only the top and bottom scoring domains
# n <- c(1:25, (length(domainOrder)-25):length(domainOrder))
n <- 1:35
domainCountsMelted <- domainCountsMelted[domainCountsMelted$GOdomain %in% domainOrder[n],]

# rename activity column
colnames(domainCountsMelted)[4] <- "Legend"
levels(domainCountsMelted$Legend) <- c("Active", "Homo sapiens Proteome", "Inactive", "PubChem BioAssay")

p1 <- ggplot(data=domainCountsMelted, aes(x=value, y=GOdomain)) +
    #=coord_trans(x = "log10") +
    geom_point(aes(color = Legend), size=4) + 
    # facet_grid(. ~ variable, scales="fixed") + 
    facet_grid(. ~ variable, scales="free_x") + 
    guides(colour=guide_legend()) +
    # scale_x_continuous(breaks=c(10, 100, 200, 400, 1000, 10000, 100000), labels=function(n){format(n, scientific = FALSE)}) + 
    scale_x_continuous(labels=function(n){format(n, scientific = FALSE)}) +
    ylab("Pfam Domain") + 
    xlab("Quantity") +
    # annotation_logticks(scaled=F, sides="b") + 
    geom_vline(xintercept=0, alpha=0) + 
    theme(
        text = element_text(size=15),
        panel.grid.major.y = element_line(colour = "grey"), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.key = element_blank(),
        # legend.position="bottom") +
        legend.position="none") +
        scale_colour_brewer(palette="Dark2", breaks= c("Active", "Inactive", "Homo sapiens Proteome", "PubChem BioAssay"))
plot(p1) # test code

# save output
ggsave("working/compounds_by_domain_p1.pdf", plot=p1, height=7, width=16)

#####################
# plot top non-fda
#####################

# combine into ggPlot style data.frame 
domainCountsNames <- name_rows(domainCounts)
domainCountsNames <- domainCountsNames[,! colnames(domainCountsNames) %in% "Var1"]
names(domainCountsNames) <- c("active FDA Approved Drugs", "active Non-FDA Compounds", "inactive FDA Approved Drugs", "inactive Non-FDA Compounds", "protein Protein Targets", "human Protein Targets", "GOdomain")
domainCountsNames$GOdomain <- gsub("\\W?domain", "", domainCountsNames$GOdomain)
domainCountsMelted <- melt(domainCountsNames)
domainCountsMelted <- domainCountsMelted[nrow(domainCountsMelted):1,]

enrichment <- domainCountsNames$`active Non-FDA Compounds`*(domainCountsNames$`active FDA Approved Drugs` == 0)
enrichmentIndex <- sort(enrichment, decreasing=T, index.return=T)$ix
domainOrder <- domainCountsNames$GOdomain[enrichmentIndex]
domainCountsMelted$GOdomain <- factor(domainCountsMelted$GOdomain, levels=rev(domainOrder))

# split active and fda approved into two columns
domainCountsMelted <- cbind(domainCountsMelted, Activity=domainCountsMelted$variable)
domainCountsMelted$variable <- gsub("^\\w+\\W", "", domainCountsMelted$variable)
domainCountsMelted$variable <- factor(domainCountsMelted$variable, levels=c("FDA Approved Drugs", "Non-FDA Compounds", "Protein Targets"))
domainCountsMelted$Activity <- gsub("^(\\w+)\\W.*", "\\1", domainCountsMelted$Activity)
domainCountsMelted$Activity <- factor(domainCountsMelted$Activity)
domainCountsMelted <- domainCountsMelted[! ((domainCountsMelted$value == 0) & (domainCountsMelted$Activity == "active")),]

# keep only the top and bottom scoring domains
# n <- c(1:25, (length(domainOrder)-25):length(domainOrder))
n <- 1:35
domainCountsMelted <- domainCountsMelted[domainCountsMelted$GOdomain %in% domainOrder[n],]

# rename activity column
colnames(domainCountsMelted)[4] <- "Legend"
levels(domainCountsMelted$Legend) <- c("Active", "Homo sapiens Proteome", "Inactive", "PubChem BioAssay")

# remove FDA approved and split other into two
domainCountsMelted <- domainCountsMelted[domainCountsMelted$variable != "FDA Approved Drugs",]
variableCol <- as.character(domainCountsMelted$variable)
variableCol[domainCountsMelted$Legend == "Inactive"] <- "Inactive Non-FDA Compounds"
variableCol[domainCountsMelted$Legend == "Active"] <- "Active Non-FDA Compounds"
domainCountsMelted$variable <- as.factor(variableCol)

p2 <- ggplot(data=domainCountsMelted, aes(x=value, y=GOdomain)) +
    #=coord_trans(x = "log10") +
    geom_point(aes(color = Legend), size=4) + 
    # facet_grid(. ~ variable, scales="fixed") + 
    facet_grid(. ~ variable, scales="free_x") + 
    guides(colour=guide_legend()) +
    # scale_x_continuous(breaks=c(10, 100, 200, 400, 1000, 10000, 100000), labels=function(n){format(n, scientific = FALSE)}) + 
    scale_x_continuous(labels=function(n){format(n, scientific = FALSE)}) +
    ylab("Pfam Domain") + 
    xlab("Quantity") +
    # annotation_logticks(scaled=F, sides="b") + 
    geom_vline(xintercept=0, alpha=0) + 
    theme(
        text = element_text(size=15),
        panel.grid.major.y = element_line(colour = "grey"), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.key = element_blank(),
        legend.position="bottom") +
    scale_colour_brewer(palette="Dark2", breaks= c("Active", "Inactive", "Homo sapiens Proteome", "PubChem BioAssay"))
plot(p2) # test code
ggsave("working/compounds_by_domain_p2.pdf", plot=p2, height=8, width=16)

# save output
# gridplot <- plot_grid(p1, p2, labels=c("A", "B"), ncol = 1, nrow = 2)
# save_plot(outputFile, gridplot, base_width=16, base_height=16)
# ggsave(outputFile, plot=p, height=10, width=12)

###############################
# plot domain screening density
###############################

domainStats <- function(domainCounts){
    activeCounts <- domainCounts$drugCidCountActives + domainCounts$targetingCidCountActives
    inactiveCounts <- domainCounts$drugCidCountInactives + domainCounts$targetingCidCountInactives
    activeThresholds <- sapply(1:max(activeCounts), function(x) sum(activeCounts > x))
    inactiveThresholds <- sapply(1:max(inactiveCounts), function(x) sum(inactiveCounts > x))
    
    actives <- cbind(index=1:length(activeThresholds), thr=activeThresholds, var="Active Compounds")
    actives <- data.frame(actives, stringsAsFactors=F)
    actives <- actives[! duplicated(actives$thr),]
    inactives <- cbind(index=1:length(inactiveThresholds), thr=inactiveThresholds, var="Inactive Compounds")
    inactives <- data.frame(inactives, stringsAsFactors=F)
    inactives <- inactives[! duplicated(inactives$thr),]
    melted_data <- rbind(actives, inactives)
    melted_data
}
allDomainStats <- domainStats(allDomainCounts)
humandomainCounts <- allDomainCounts[row.names(allDomainCounts) %in% allHumanDomains,]
hDomainStats <- domainStats(humandomainCounts)

melted_data <- rbind(
    cbind(allDomainStats, label="All Domains"),
    cbind(hDomainStats, label="H. sapiens Domains")
)

# melted_data <- data.frame(melted_data)
# melted_data <- allDomainStats
melted_data[,1] <- as.numeric(melted_data[,1])
melted_data[,2] <- as.numeric(melted_data[,2])

# plot w/ ggplot2
p <- ggplot(melted_data, aes(x=thr, y=index, colour=var)) +
    geom_path(size=1) + 
    facet_grid(. ~ label, scales="fixed", switch="y") + 
    # facet_grid(var ~ ., scales="free", switch="y") + 
    # geom_path(data = filter(melted_data, Function %in% c("Stretched Exponential", "Power Law", "Exponential")), size=1.5) +
    # geom_point(data = filter(melted_data, Function == "Distinct Sequence Targets"), size=2) +
    xlab("Quantity of Pfam Domains Above Threshold") +
    ylab("Screened Compounds") + 
    # scale_y_continuous(labels=function(n){format(n, scientific = FALSE)}) +
    # xlim(0, 1) + 
    scale_colour_brewer(palette="Accent") +
    # coord_trans(y = "log10") +
    annotation_logticks(sides="l") +
    scale_y_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    # coord_trans(y = "log10") +
    # scale_y_continuous(labels=function(n){format(n, scientific = FALSE)}) +
    # scale_x_continuous(breaksr=1:20) + 
    # geom_hline(yintercept=0, alpha=0) +
    theme(
        text = element_text(size=15),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        # axis.title.y=element_blank(),
        legend.key = element_blank(),
        legend.position="bottom",
        legend.title=element_blank()
    ) 
plot(p)
ggsave("working/compounds_by_domain_participation.pdf", plot=p, width=12)

########################
# stats for paper
########################

# Pfam domains represented in PubChem BioAssay that are associated with active compounds
activeCounts <- allDomainCounts$drugCidCountActives + allDomainCounts$targetingCidCountActives
sum(activeCounts > 0)

# total humans domains
length(allHumanDomains)

# fraction of domains with human proteins
activeDomainList <- unique(as.character(row.names(allDomainCounts))[activeCounts > 0])
length(intersect(allHumanDomains, activeDomainList))
100*length(intersect(allHumanDomains, activeDomainList))/length(allHumanDomains)

length(allHumanDomains) - length(intersect(allHumanDomains, domainList))
length(domainList) - length(intersect(allHumanDomains, domainList))

# get fraction of human domains tested in human targets
domainSpeciesTable <- queryBioassayDB(database, "SELECT organism, target, identifier FROM assays NATURAL JOIN targets NATURAL JOIN targetTranslations where category='domains'")
humanDomains <- domainSpeciesTable[domainSpeciesTable$organism == "Homo_sapiens",]
domainsScreenedInHumanTargets <- unique(humanDomains$identifier)
activeDomainsScreenedInHumanTargets <- intersect(activeDomainList, domainsScreenedInHumanTargets)
100*length(intersect(allHumanDomains, activeDomainsScreenedInHumanTargets))/length(allHumanDomains)

# non-human proteins
length(activeDomainList) - length(intersect(allHumanDomains, activeDomainList))

inactiveCounts <- allDomainCounts$drugCidCountInactives + allDomainCounts$targetingCidCountInactives

# actives only
activesOnly <- row.names(allDomainCounts)[(activeCounts > 0) & (inactiveCounts == 0)]
length(activesOnly)
length(intersect(allHumanDomains, activesOnly))

# what is the deal with the all zero rows?
# allZeros <- row.names(domainCounts[(activeCounts + inactiveCounts) == 0,])
# domainList <- lapply(targetList, translateTargetId, database=database, category="domains")
# translateTargetId(database, allZeros[1], category="GI", fromCategory = "domains")

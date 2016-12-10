#!/usr/bin/env Rscript

# (C) 2016 Tyler William H Backman
# Purpose: plot target selectivity against compounds from an individual domain 

library(ggplot2)
library(cowplot)
library(xtable)

# parse input options
inputFile <- commandArgs(trailingOnly=TRUE)[1]
outputFile <- commandArgs(trailingOnly=TRUE)[2]

# test code for running without make:
if(is.na(commandArgs(trailingOnly=TRUE)[1])){
    inputFile <- "working/computeSelectivitySpecificDomain.RData"
    outputFile <- "working/selectivitySpecificDomains.pdf"
}

# parse input files
load(inputFile) # loads results

# make table of number of compounds
splitFreqs <- split(results$category, results$domain)
compoundCounts <- sapply(splitFreqs, table)
domainScreeningCounts <- t(compoundCounts[c("drugScreeningFrequency", "nonDrugScreeningFrequency"),])
domainScreeningCounts <- cbind(domainScreeningCounts, totalTargets[match(row.names(domainScreeningCounts), row.names(totalTargets)),])
colnames(domainScreeningCounts) <- c("Total FDA Approved Compounds", "Total Other Compounds", "Total Protein Targets")
domainScreeningCounts <- domainScreeningCounts[order(domainScreeningCounts[,3], decreasing=TRUE),]

# get names with at least 10 compounds from each
keepDomains <- row.names(domainScreeningCounts)[(domainScreeningCounts[,1] > 9) & (domainScreeningCounts[,2] > 9)]
domainScreeningCounts <- domainScreeningCounts[row.names(domainScreeningCounts) %in% keepDomains,]
results <- results[results$domain %in% keepDomains,]

# eliminate those where median is identical
splitFreqs <- split(results, results$domain)
medianActives <- sapply(splitFreqs, function(x){
    drugMean <- median(x[x$category == "drugActiveFrequency",3])
    nonDrugMean <- median(x[x$category == "nonDrugActiveFrequency",3])
    return(c(drugMean, nonDrugMean))
})
unequalMedianDomains <- colnames(medianActives)[medianActives[1,] != medianActives[2,]]
results <- results[results$domain %in% unequalMedianDomains,]
domainScreeningCounts <- domainScreeningCounts[row.names(domainScreeningCounts) %in% unequalMedianDomains,]

# fix labels
results$category[results$category == "drugActiveFrequency"] <- "FDA Approved"

activeResults <- results[results$category %in% c("drugActiveFrequency", "nonDrugActiveFrequency"),]
activeResults$category[activeResults$category == "drugActiveFrequency"] <- "FDA Approved"
activeResults$category[activeResults$category == "nonDrugActiveFrequency"] <- "Non-FDA"
p1 <- ggplot(activeResults, aes(x=factor(category), y=frequency, fill=category)) +
    geom_boxplot(outlier.shape=NA) + 
    coord_cartesian(ylim=c(0,32)) +
    facet_grid(. ~ factor(domain, levels=row.names(domainScreeningCounts)), scales="fixed") +
    ylab("Active Protein Targets With Domain") +
    xlab("Pfam Domains") +
    scale_fill_brewer(palette="Set1") + 
    theme(
        text = element_text(size=13),
        panel.grid.major.y = element_line(colour = "grey"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line.x = element_blank(),
        axis.line.y = element_line(colour = "black"),
        # axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="bottom",
        # legend.title = element_blank(),
        legend.key = element_blank()
    )
# plot(p1)

inactiveResults <- results[! results$category %in% c("drugActiveFrequency", "nonDrugActiveFrequency"),]
p2 <- ggplot(inactiveResults, aes(x=factor(category), y=frequency, fill=category)) +
    geom_boxplot(outlier.shape=NA) + 
    coord_cartesian(ylim=c(0,67)) +
    facet_grid(. ~ factor(domain, levels=row.names(domainScreeningCounts)), scales="fixed") +
    ylab("Tested Protein Targets With Domain") +
    xlab("Pfam Domains") +
    scale_fill_brewer(palette="Set1") + 
    theme(
        text = element_text(size=13),
        panel.grid.major.y = element_line(colour = "grey"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line.x = element_blank(),
        axis.line.y = element_line(colour = "black"),
        # axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="none",
        # legend.title = element_blank(),
        legend.key = element_blank()
    )
# plot(p2)

gridplot <- plot_grid(p1, p2, labels=c("A", "B"), ncol = 1, nrow = 2)
# plot(gridplot)
save_plot(outputFile, gridplot, base_width=12, base_height=9)

xtmp <- xtable(domainScreeningCounts, caption="Selectivity by domain", label="domainScreeningCounts")
print(xtmp, type="latex", file="working/selectivitySpecificDomains.tex", include.rownames=T)


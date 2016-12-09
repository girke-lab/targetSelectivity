#!/usr/bin/env Rscript

# (C) 2016 Tyler William H Backman
# Purpose: plot target selectivity against compounds from an individual domain 

library(ggplot2)

# parse input options
inputFile <- commandArgs(trailingOnly=TRUE)[1]

# test code for running without make:
if(is.na(commandArgs(trailingOnly=TRUE)[1])){
    inputFile <- "working/computeSelectivitySpecificDomain.RData"
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

activeResults <- results[results$category %in% c("drugActiveFrequency", "nonDrugActiveFrequency"),]
p <- ggplot(activeResults, aes(x=factor(category), y=frequency, fill=category)) +
    geom_boxplot(outlier.shape=NA) + 
    coord_cartesian(ylim=c(0,35)) +
    facet_grid(. ~ factor(domain, levels=row.names(domainScreeningCounts)), scales="fixed") +
    ylab("Active Protein Targets") +
    xlab("Domains") +
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
plot(p)

inactiveResults <- results[! results$category %in% c("drugActiveFrequency", "nonDrugActiveFrequency"),]
p <- ggplot(inactiveResults, aes(x=factor(category), y=frequency, fill=category)) +
    geom_boxplot(outlier.shape=NA) + 
    coord_cartesian(ylim=c(0,70)) +
    facet_grid(. ~ domain, scales="fixed") +
    ylab("Tested Protein Targets") +
    # xlab("Compounds") +
    scale_fill_brewer(palette="Set1") + 
    theme(
        text = element_text(size=13),
        panel.grid.major.y = element_line(colour = "grey"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line.x = element_blank(),
        axis.line.y = element_line(colour = "black"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="bottom",
        # legend.title = element_blank(),
        legend.key = element_blank()
    )
plot(p)


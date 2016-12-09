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
load(inputFile)

activeResults <- results[results$category %in% c("drugActiveFrequency", "nonDrugActiveFrequency"),]
p <- ggplot(activeResults, aes(x=factor(category), y=frequency, fill=category)) +
    geom_boxplot(outlier.shape=NA) + 
    # coord_cartesian(ylim=c(0,2)) +
    facet_grid(. ~ domain, scales="fixed") +
    ylab("Active Protein Targets") +
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

inactiveResults <- results[! results$category %in% c("drugActiveFrequency", "nonDrugActiveFrequency"),]
p <- ggplot(inactiveResults, aes(x=factor(category), y=frequency, fill=category)) +
    geom_boxplot(outlier.shape=NA) + 
    # coord_cartesian(ylim=c(0,2)) +
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
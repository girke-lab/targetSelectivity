#!/usr/bin/env Rscript

# (C) 2016 Tyler William H Backman
# Purpose: 
# plot distribution of target promiscuity (hit ratio)
# for PAINS vs non-PAINS

library(R.utils)
library(ggplot2)
# library(rethinking)
library(bioassayR)
library(grid)
library(gridExtra)

# parse input options
databaseFile <- commandArgs(trailingOnly=TRUE)[1]
drugbank_linksFile <- commandArgs(trailingOnly=TRUE)[2]
highlyScreened <- commandArgs(trailingOnly=TRUE)[3]
activeCidsFile <- commandArgs(trailingOnly=TRUE)[4]
painsFile <- commandArgs(trailingOnly=TRUE)[5]
unparsablePainsFile <- commandArgs(trailingOnly=TRUE)[6]
outputFile <- commandArgs(trailingOnly=TRUE)[7]

# test code for running without make:
if(is.null(commandArgs(trailingOnly=TRUE)[1])){
    # databaseFile <- "~/Downloads/pubchem_protein_only.sqlite"
    databaseFile <- "/dev/shm/bioassayDatabase.sqlite"
    highlyScreened <- "working/highlyScreenedCids.txt"
    activeCidsFile <- "working/activeCids.txt"
    painsFile <- "working/activeCompoundsPAINS.txt"
    unparsablePainsFile <- "working/activeCompoundsPAINSuFixed.txt"
    outputFile <- "working/hitratioDistributionPAINS.pdf"
}

# parse input files
database <- connectBioassayDB(databaseFile)
highlyScreenedcids <- read.table(highlyScreened)[,1]
activeCids <- read.table(activeCidsFile)[[1]]
pains <- as.character(read.table(painsFile)[[1]])
unparsablePains <- read.table(unparsablePainsFile)[[1]]

# collect list of highly screened active pains
highlyScreenedPains <- intersect(pains, highlyScreenedcids)
highlyScreenedActivePains <- intersect(highlyScreenedPains, activeCids)

# get background
highlyScreenedActive <- intersect(highlyScreenedcids, activeCids)
highlyScreenedActiveNonPains <- setdiff(highlyScreenedActive, highlyScreenedActivePains)
highlyScreenedActiveNonPains <- setdiff(highlyScreenedActiveNonPains, unparsablePains)

# options
# totalCompounds <- 10000 # max number of compounds to build distribution from
# totalSamples <- 10*totalCompounds # must be an integer multiple of totalCompounds
totalSamples <- 1000000

# get random samples
getSamples <- function(cid, database, sampleCount=100, prior=list(hit_ratio_mean=0.01861531, hit_ratio_sd=0.03494048)){
    actives <- row.names(activeTargets(database, cid))
    inactives <- row.names(inactiveTargets(database, cid))
    n <- length(actives)
    N <- length(union(actives, inactives))
    alpha <- ((1 - prior$hit_ratio_mean) / prior$hit_ratio_sd ^ 2 - 1 / prior$hit_ratio_mean) * prior$hit_ratio_mean ^ 2
    beta <- alpha * (1 / prior$hit_ratio_mean - 1)
    return(rbeta(sampleCount, n + alpha, N - n + beta))
}

# get samples for pains compounds
# painsSubset <- sample(highlyScreenedActivePains, totalCompounds, replace=FALSE)
samplesPerCompound <- ceiling(totalSamples/length(highlyScreenedActivePains))
set.seed(900)
painsSamples <- unlist(lapply(highlyScreenedActivePains, getSamples, database=database, sampleCount=samplesPerCompound))
painsSamples <- sample(painsSamples, totalSamples, replace=FALSE)

# get samples for nonpromiscuous compounds
# nonPainsSubset <- sample(highlyScreenedActiveNonPains, totalCompounds, replace=FALSE)
samplesPerCompound <- ceiling(totalSamples/length(highlyScreenedActiveNonPains))
nonpainsSamples <- unlist(lapply(highlyScreenedActiveNonPains, getSamples, database=database, sampleCount=samplesPerCompound))
nonpainsSamples <- sample(nonpainsSamples, totalSamples, replace=FALSE)
save(list = ls(all.names=T), file = "working/hitratioDistributionPAINS1M.RData")
# load("working/hitratioDistributionPAINS1M.RData")

# compute Kolmogorov-Smirnov separation
ksResult <- ks.test(painsSamples, nonpainsSamples)
ksResult
ks.test(sample(painsSamples, 100000), sample(nonpainsSamples,100000))

# make ggplot2 compatible object
allSamples <- rbind(cbind.data.frame(Compounds="PAINS", samples=painsSamples), cbind.data.frame(Compounds="NonPAINS", samples=nonpainsSamples))

p0 <- ggplot(allSamples, aes(samples, fill = Compounds, colour = Compounds)) +
    geom_density(alpha = 0.5, adjust=3) +
    xlim(0, 1) +
    coord_cartesian(xlim = c(0, 1), ylim=c(0,55)) +
    ylab("Probability Density") +
    xlab("Hit Ratio (fraction of actives)") +
    scale_fill_manual(values = brewer.pal(8,"Dark2")[5:6]) + 
    scale_colour_manual(values = brewer.pal(8,"Dark2")[5:6]) + 
    theme(
        text = element_text(size=15),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        # legend.key = element_blank(),
        legend.position = "none",
        # axis.text.x=element_blank(),
        axis.title.x=element_blank()
        # axis.title.y=element_blank()
    )
plot(p0)
ggsave("working/hitratioDistributionPAINS_top.pdf", plot=p0, width=12, height=2)

p <- ggplot(allSamples, aes(samples, fill = Compounds, colour = Compounds)) +
    geom_density(alpha = 0.5, adjust=3) +
    xlim(0, 1) +
    coord_cartesian(xlim = c(0, 0.05), ylim=c(0,55)) +
    ylab("Probability Density") +
    xlab("Hit Ratio (fraction of actives)") +
    scale_fill_manual(values = brewer.pal(8,"Dark2")[5:6]) + 
    scale_colour_manual(values = brewer.pal(8,"Dark2")[5:6]) + 
    theme(
        text = element_text(size=15),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        legend.key = element_blank(),
        legend.position="bottom"
    )
plot(p)
# ggsave(outputFile, plot=p, width=12)

p2 <- ggplot(allSamples, aes(samples, fill = Compounds, colour = Compounds)) +
    geom_density(alpha = 0.5, adjust=3) +
    xlim(0, 1) +
    coord_cartesian(xlim = c(0.05,0.35), ylim=c(0,7.5)) +
    # ylab("Probability Density") +
    xlab("Hit Ratio (fraction of actives)") +
    scale_fill_manual(values = brewer.pal(8,"Dark2")[5:6]) + 
    scale_colour_manual(values = brewer.pal(8,"Dark2")[5:6]) + 
    theme(
        text = element_text(size=15),
        axis.title.y=element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        legend.key = element_blank(),
        legend.position="bottom"
    )
plot(p2)
# ggsave("working/hitratioDistribution_0.2_1.pdf", plot=p2, width=12)

p3 <- ggplot(allSamples, aes(samples, fill = Compounds, colour = Compounds)) +
    geom_density(alpha = 0.5, adjust=3) +
    xlim(0, 1) +
    coord_cartesian(xlim = c(0.35,1), ylim=c(0,.06)) +
    # ylab("Probability Density") +
    xlab("Hit Ratio (fraction of actives)") +
    scale_fill_manual(values = brewer.pal(8,"Dark2")[5:6]) + 
    scale_colour_manual(values = brewer.pal(8,"Dark2")[5:6]) + 
    theme(
        text = element_text(size=15),
        axis.title.y=element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        legend.key = element_blank(),
        legend.position="bottom"
    )
plot(p3)

g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}

mylegend<-g_legend(p)

pdf(outputFile, width = 12, height=5.2)
grid.arrange(arrangeGrob(p0), arrangeGrob(p + theme(legend.position="none", plot.margin = unit(c(2,.2,-1,.2), units = "lines"))  + xlab(""),
                                          p2 + theme(legend.position="none", plot.margin = unit(c(2,.2,-1,.2), units = "lines"))  + xlab(""),
                                          p3 + theme(legend.position="none", plot.margin = unit(c(2,.2,-1,.2), units = "lines"))  + xlab(""),
                                          nrow=1, widths=c(0.33, 0.33, 0.33)),
             arrangeGrob(textGrob("Hit Ratio (fraction of actives)", gp=gpar(fontsize=15)), mylegend),
             nrow=3,heights=c(2, 2.5, 0.7))
dev.off()


# compute HDPI
# HPDI(painsSamples, 0.85)
# HPDI(nonpainsSamples, 0.85)

quantile(painsSamples, probs=c(0.85))
quantile(nonpainsSamples, probs=c(0.85))

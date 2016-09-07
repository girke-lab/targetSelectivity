#!/usr/bin/env Rscript

# (C) 2016 Tyler William H Backman
# Purpose: 
# plot distribution of target promiscuity (hit ratio)
# for aggregators and nonaggregators

library(R.utils)
library(ggplot2)
# library(rethinking)
library(bioassayR)
library(grid)
library(gridExtra)
library(RColorBrewer)

# parse input options
databaseFile <- commandArgs(trailingOnly=TRUE)[1]
drugbank_linksFile <- commandArgs(trailingOnly=TRUE)[2]
highlyScreened <- commandArgs(trailingOnly=TRUE)[3]
activeCidsFile <- commandArgs(trailingOnly=TRUE)[4]
outputFile <- commandArgs(trailingOnly=TRUE)[5]

# test code for running without make:
if(is.null(commandArgs(trailingOnly=TRUE)[1])){
    # databaseFile <- "~/Downloads/pubchem_protein_only.sqlite"
    databaseFile <- "/dev/shm/bioassayDatabase.sqlite"
    drugbank_linksFile <- "working/drugbank_links.csv"
    highlyScreened <- "working/highlyScreenedCids.txt"
    activeCidsFile <- "working/activeCids.txt"
    outputFile <- "working/hitratioDistributionAggregators.pdf"
}

# parse input files
database <- connectBioassayDB(databaseFile)
drugbank_links <- read.csv(drugbank_linksFile)
drugCids <- unique(drugbank_links$PubChem.Compound.ID)
highlyScreenedcids <- read.table(highlyScreened)[,1]
activeCids <- read.table(activeCidsFile)[[1]]

# find promsicuous and non-promiscuous compounds
# Compounds that inhibit only in the absence of detergent are considered likely promiscuous aggregators. 
# 584 - with detergent, 585 - without detergent
assays <- getAssays(database, c(584, 585))
activitymatrix <- as.matrix(assays@activity)

# remove all NAs
activitymatrix <- activitymatrix[,!is.na(activitymatrix[1,])]
activitymatrix <- activitymatrix[,!is.na(activitymatrix[2,])]

# get cids
promsicuous <- colnames(activitymatrix)[(activitymatrix[2,] == 2) & (activitymatrix[1,] == 1)]
noninhibitors <- colnames(activitymatrix)[activitymatrix[2,] == 1]
inhibitors <- colnames(activitymatrix)[(activitymatrix[2,] == 2) & (activitymatrix[1,] == 2)]
nonpromsicuous <- union(inhibitors, noninhibitors)

# keep only highly screened compounds
promsicuous <- promsicuous[promsicuous %in% highlyScreenedcids]
nonpromsicuous <- nonpromsicuous[nonpromsicuous %in% highlyScreenedcids]

# keep only those with active scores
promsicuous <- intersect(promsicuous, activeCids)
nonpromsicuous <- intersect(nonpromsicuous, activeCids)

# options
# totalCompounds <- 10000 # max number of compounds to build distribution from
totalSamples <- 1000000 # must be an integer multiple of totalCompounds

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

# get samples for promiscuous compounds
samplesPerCompound <- ceiling(totalSamples/length(promsicuous))
set.seed(900)
promsicuousSamples <- unlist(lapply(promsicuous, getSamples, database=database, sampleCount=samplesPerCompound))
promsicuousSamples <- sample(promsicuousSamples, totalSamples, replace=FALSE)

# get samples for nonpromiscuous compounds
# nonpromsicuousSubset <- sample(nonpromsicuous, totalCompounds, replace=FALSE)
samplesPerCompound <- ceiling(totalSamples/length(nonpromsicuous))
nonpromsicuousSamples <- unlist(lapply(nonpromsicuous, getSamples, database=database, sampleCount=samplesPerCompound))
nonpromsicuousSamples <- sample(nonpromsicuousSamples, totalSamples, replace=FALSE)

save(list = ls(all.names=T), file = "working/hitratioDistributionAggregators1M.RData")
# load("working/hitratioDistributionAggregators1M.RData")

# compute Kolmogorov-Smirnov separation
ksResult <- ks.test(promsicuousSamples, nonpromsicuousSamples)
ksResult
ks.test(sample(promsicuousSamples, 100000), sample(nonpromsicuousSamples,100000))

# make ggplot2 compatible object
allSamples <- rbind(cbind.data.frame(Compounds="Aggregators", samples=promsicuousSamples), cbind.data.frame(Compounds="Nonaggregators", samples=nonpromsicuousSamples))

p0 <- ggplot(allSamples, aes(samples, fill = Compounds, colour = Compounds)) +
    geom_density(alpha = 0.5, adjust=2) +
    xlim(0, 1) +
    coord_cartesian(xlim = c(0, 1), ylim=c(0,65)) +
    ylab("Probability Density") +
    xlab("Hit Ratio (fraction of actives)") +
    # scale_fill_brewer(palette="Dark2") +
    # scale_color_brewer(palette="Dark2") +
    scale_fill_manual(values = brewer.pal(8,"Dark2")[3:4]) + 
    scale_colour_manual(values = brewer.pal(8,"Dark2")[3:4]) + 
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
ggsave("working/hitratioDistributionAggregators_top.pdf", plot=p0, width=12, height=2)

p <- ggplot(allSamples, aes(samples, fill = Compounds, colour = Compounds)) +
    geom_density(alpha = 0.5, adjust=3) +
    xlim(0, 1) +
    coord_cartesian(xlim = c(0, 0.05), ylim=c(0,65)) +
    ylab("Probability Density") +
    xlab("Hit Ratio (fraction of actives)") +
    scale_fill_manual(values = brewer.pal(8,"Dark2")[3:4]) + 
    scale_colour_manual(values = brewer.pal(8,"Dark2")[3:4]) + 
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
    coord_cartesian(xlim = c(0.05,0.35), ylim=c(0,4.2)) +
    # ylab("Probability Density") +
    xlab("Hit Ratio (fraction of actives)") +
    scale_fill_manual(values = brewer.pal(8,"Dark2")[3:4]) + 
    scale_colour_manual(values = brewer.pal(8,"Dark2")[3:4]) + 
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
    coord_cartesian(xlim = c(0.35,1), ylim=c(0,.020)) +
    # ylab("Probability Density") +
    xlab("Hit Ratio (fraction of actives)") +
    scale_fill_manual(values = brewer.pal(8,"Dark2")[3:4]) + 
    scale_colour_manual(values = brewer.pal(8,"Dark2")[3:4]) + 
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



# compute HDPI for drugs and non-drugs
# HPDI(promsicuousSamples, 0.85)
# HPDI(nonpromsicuousSamples, 0.85)

quantile(promsicuousSamples, probs=c(0.85))
quantile(nonpromsicuousSamples, probs=c(0.85))

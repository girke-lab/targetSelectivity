#!/usr/bin/env Rscript

# (C) 2016 Tyler William H Backman
# Purpose: 
# plot distribution of target promiscuity (hit ratio)

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
outputFile <- commandArgs(trailingOnly=TRUE)[5]

# test code for running without make:
if(is.null(commandArgs(trailingOnly=TRUE)[1])){
    databaseFile <- "~/Downloads/bioassayDatabase.sqlite"
    # databaseFile <- "/dev/shm/bioassayDatabase.sqlite"
    drugbank_linksFile <- "working/drugbank_links.csv"
    highlyScreened <- "working/highlyScreenedCids.txt"
    activeCidsFile <- "working/activeCids.txt"
    outputFile <- "working/hitratioDistribution.pdf"
}

# parse input files
database <- connectBioassayDB(databaseFile)
drugbank_links <- read.csv(drugbank_linksFile)
drugCids <- unique(drugbank_links$PubChem.Compound.ID)
highlyScreenedcids <- read.table(highlyScreened)[,1]
activeCids <- read.table(activeCidsFile)[[1]]

# keep only highly screened compounds
highlyScreenedDrugs <- highlyScreenedcids[highlyScreenedcids %in% drugCids]
highlyScreenedNonDrugs <- highlyScreenedcids[! highlyScreenedcids %in% drugCids]

# keep only compounds with at least one active score
highlyScreenedDrugs <- intersect(highlyScreenedDrugs, activeCids)
highlyScreenedNonDrugs <- intersect(highlyScreenedNonDrugs, activeCids)

# options
# totalCompounds <- 10000 # number of non-drug compounds to build distribution from
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

# get drug samples
set.seed(900)
samplesPerCompound <- ceiling(totalSamples/length(highlyScreenedDrugs))
drugSamples <- unlist(lapply(highlyScreenedDrugs, getSamples, database=database, sampleCount=samplesPerCompound))
drugSamples <- sample(drugSamples, totalSamples, replace=FALSE)

# nonDrugSubset <- sample(highlyScreenedNonDrugs, totalCompounds, replace=FALSE)
samplesPerCompound <- ceiling(totalSamples/length(highlyScreenedNonDrugs))
nonDrugSamples <- unlist(lapply(highlyScreenedNonDrugs, getSamples, database=database, sampleCount=samplesPerCompound))
nonDrugSamples <- sample(nonDrugSamples, totalSamples, replace=FALSE)

save(list = ls(all.names=T), file = "working/hitratioDistribution1M.RData")
# load("working/hitratioDistribution1M.RData")

# compute Kolmogorov-Smirnov separation
ksResult <- ks.test(drugSamples, nonDrugSamples)
ksResult # D = 0.50592, p-value < 2.2e-16
# ks.test(sample(drugSamples, 600000), sample(nonDrugSamples,600000))

# make ggplot2 compatible object
allSamples <- rbind(cbind.data.frame(Compounds="FDA Approved", samples=drugSamples), cbind.data.frame(Compounds="Non-FDA", samples=nonDrugSamples))

p0 <- ggplot(allSamples, aes(samples, fill = Compounds, colour = Compounds)) +
    geom_density(alpha = 0.5, adjust=3) +
    xlim(0, 1) +
    coord_cartesian(xlim = c(0, 1), ylim=c(0,50)) +
    ylab("Probability Density") +
    xlab("Hit Ratio (fraction of actives)") +
    scale_fill_brewer(palette="Set1") +
    scale_color_brewer(palette="Set1") +
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
ggsave("working/hitratioDistribution_top.pdf", plot=p0, width=12, height=2)

p <- ggplot(allSamples, aes(samples, fill = Compounds, colour = Compounds)) +
    geom_density(alpha = 0.5, adjust=3) +
    xlim(0, 1) +
    coord_cartesian(xlim = c(0, 0.05), ylim=c(0,55)) +
    ylab("Probability Density") +
    xlab("Hit Ratio (fraction of actives)") +
    scale_fill_brewer(palette="Set1") +
    scale_color_brewer(palette="Set1") +
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
    coord_cartesian(xlim = c(0.05,0.35), ylim=c(0,8)) +
    # ylab("Probability Density") +
    xlab("Hit Ratio (fraction of actives)") +
    scale_fill_brewer(palette="Set1") +
    scale_color_brewer(palette="Set1") +
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
    coord_cartesian(xlim = c(0.35,1), ylim=c(0,.4)) +
    # ylab("Probability Density") +
    xlab("Hit Ratio (fraction of actives)") +
    scale_fill_brewer(palette="Set1") +
    scale_color_brewer(palette="Set1") +
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
# HPDI(drugSamples, 0.90)
# HPDI(nonDrugSamples, 0.90)

quantile(drugSamples, probs=c(0.85)) # 0.2321003
quantile(nonDrugSamples, probs=c(0.85)) # 0.03360872

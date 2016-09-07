#!/usr/bin/env Rscript

# (C) 2015 Tyler William H Backman
# Purpose: make MDS plot of biclusters (color) vs structure (spatial)

library(R.utils)
library(foreach)
library(ChemmineR)
library(ggplot2)
library(gridExtra)
library(grid)

# parse input options
biClustersFile <- commandArgs(trailingOnly=TRUE)[1] 
drugBankStructuresFile <- commandArgs(trailingOnly=TRUE)[2]
drugbankLinksFile <- commandArgs(trailingOnly=TRUE)[3]
outputFile <- commandArgs(trailingOnly=TRUE)[4]

# test code for running without make:
if(is.null(commandArgs(trailingOnly=TRUE)[1])){
    biClustersFile <- "working/biClusters.RData"
    drugBankStructuresFile <- "working/drugbank.sdf" 
    drugbankLinksFile <- "working/drugbank_links.csv"
    outputFile <- "working/biclusterMDS.pdf"
}

# parse input files
load(biClustersFile) # loads clusterResults
drugbankCompounds <- read.SDFset(drugBankStructuresFile)
drugBankTable <- read.csv(drugbankLinksFile)

# reformat clusterResults into two cols: cid, clusterid
clustersPerCid <- foreach(thisCluster=1:length(clusterResults), .combine='rbind2') %do% {
    cbind(clusterResults[[thisCluster]]$rows, thisCluster)
}

# remove small biclusters
clusterCounts <- table(clustersPerCid[,2])
largeClusters <- names(clusterCounts)[clusterCounts > 1]
clustersPerCid <- clustersPerCid[clustersPerCid[,2] %in% largeClusters,]

# renumber biclusters
newNumbers <- 1:length(unique(clustersPerCid[,2]))
names(newNumbers) <- unique(clustersPerCid[,2])
clustersPerCid[,2] <- as.character(newNumbers[match(clustersPerCid[,2], names(newNumbers))])

# assign each compound to a unique cluster
clustersPerCid <- clustersPerCid[! duplicated(clustersPerCid[,1]),]

# keep only top 12 largest clusters
topClusters <- names(sort(table(clustersPerCid[,2]), decreasing=T)[1:12])
clustersPerCid <- clustersPerCid[clustersPerCid[,2] %in% topClusters,]

# fix cids in drugbank compounds
dbIdsWithCids <- drugBankTable[! is.na(drugBankTable[,c('PubChem.Compound.ID')]),c('DrugBank.ID')]
drugbankCompounds <- drugbankCompounds[datablocktag(drugbankCompounds, "DRUGBANK_ID") %in% dbIdsWithCids]
cidPositions <- match(datablocktag(drugbankCompounds, "DRUGBANK_ID"), drugBankTable[,c('DrugBank.ID')])
cid(drugbankCompounds) <- as.character(drugBankTable[cidPositions,c('PubChem.Compound.ID')])
drugbankCompounds <- drugbankCompounds[! duplicated(cid(drugbankCompounds))]

# get structures
biclusterCids <- unique(clustersPerCid[,1])
biclusterCompounds <- drugbankCompounds[cid(drugbankCompounds) %in% biclusterCids]

# Create atom pair distance matrix
apset <- sdf2ap(biclusterCompounds) 

# compute clustering coordinates
clusters <- cmp.cluster(apset, cutoff = 0.5) 
myTempFile <- tempfile(fileext=".pdf")
coords <- cluster.visualize(apset, clusters, size.cutoff=1, quiet = TRUE, non.interactive=myTempFile)
coords <- coords[match(clustersPerCid[,1], rownames(coords)),]
plotdata <- coords[,1:2]
unlink(myTempFile)

# set up ggplot2 frame
dat <- data.frame(xvar=plotdata[,1], yvar=plotdata[,2], bicluster=clustersPerCid[,2])
dat$bicluster <- factor(as.character(dat$bicluster), levels=unique(as.character(dat$bicluster)))

mytheme = theme(plot.margin = unit(c(.2,.2,.2,.2), units = "lines"), 
    text = element_text(size=15),
    axis.text = element_blank(), 
    axis.title = element_blank(), 
    axis.ticks = element_blank(), 
    axis.ticks.length = unit(0, "lines"), 
    axis.ticks.margin = unit(0, "lines"))
mylabs = labs(x = NULL, y = NULL) 

# placeholder plot - prints nothing at all
empty <- ggplot() + geom_point(aes(1, 1), colour = "white") + theme(plot.background = element_blank(), 
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    panel.border = element_blank(), panel.background = element_blank(), axis.title.x = element_blank(), 
    axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), 
    axis.ticks = element_blank())

# scatterplot of x and y variables
scatter <- ggplot(dat, aes(xvar, yvar)) + 
    # geom_point(colour="black", size = 3) +
    geom_point(aes(fill = bicluster), colour="black", pch=21, size=4) + # , shape=19) + 
    scale_fill_brewer(palette="Paired") + 
    xlab("Principal Coordinate 1") + ylab("Principal Coordinate 2") 


scatterTheme <-     theme(plot.margin = unit(c(.2,.2,.2,.2), units = "lines"),
                          axis.ticks = element_blank(), 
                          text = element_text(size=15),
                          # axis.line.x = element_line(colour = "black"),
                          # axis.line.y = element_line(colour = "black"),
                          # panel.grid.major.y = element_blank(),
                          # panel.grid.major.x = element_blank(),
                          panel.background = element_rect(fill = "darkgrey"),
                          panel.grid.minor = element_blank(), 
                          # panel.background = element_blank(), 
                          axis.ticks.length = unit(0, "lines"),
                          legend.position="none",
                          legend.key = element_blank(),
                          axis.text = element_blank())
plot(scatter + scatterTheme)

# marginal density of x - plot on top
plot_top <- ggplot(dat, aes(xvar, fill = bicluster)) + geom_density(alpha = 0.5) + 
    scale_fill_brewer(palette="Paired") + theme(legend.position = "none", panel.background = element_blank()) +
    scale_colour_brewer(palette="Paired") + 
    mytheme + mylabs + theme(plot.margin = unit(c(1,.2,.2,.2), units = "lines"))

# marginal density of y - plot on the right
plot_right <- ggplot(dat, aes(yvar, fill = bicluster)) + geom_density(alpha = 0.5) + 
    coord_flip() + scale_fill_brewer(palette="Paired") + theme(legend.position = "none", panel.background = element_blank()) +
    mytheme + mylabs

# plot with legend
g <- ggplotGrob(scatter + theme(legend.position="bottom", legend.key = element_blank()) + guides(fill=guide_legend(nrow=1,byrow=TRUE)))$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
lheight <- sum(legend$height)
scatter_nolegend <- scatter + scatterTheme

p <- arrangeGrob(
    arrangeGrob(plot_top, empty, scatter_nolegend, plot_right, ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4)),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
plot(p)
ggsave(outputFile, plot=p, width=12, height=8)


#!/usr/bin/env Rscript

# (C) 2016 Tyler William H Backman
# Purpose: 
# plot distribution of target selectivity vs non-H atom count

library(R.utils)
library(ggplot2)
library(cowplot)
library(bioassayR)
library(xtable)

# parse input options
selectivityCountsFile <- commandArgs(trailingOnly=TRUE)[1]
drugbank_linksFile <- commandArgs(trailingOnly=TRUE)[2]
compoundSizeFile <- commandArgs(trailingOnly=TRUE)[3]
highlyScreened <- commandArgs(trailingOnly=TRUE)[4]
activeCidsFile <- commandArgs(trailingOnly=TRUE)[5]
databaseFile <- commandArgs(trailingOnly=TRUE)[6]
outputFile <- commandArgs(trailingOnly=TRUE)[7]

# test code for running without make:
if(is.null(commandArgs(trailingOnly=TRUE)[1])){
    selectivityCountsFile <- "working/selectivityCountskClust.txt"
    drugbank_linksFile <- "working/drugbank_links.csv"
    compoundSizeFile <- "working/heavycount.txt"
    highlyScreened <- "working/highlyScreenedCids.txt"
    activeCidsFile <- "working/activeCids.txt"
    databaseFile <- "~/Downloads/pubchem_protein_only.sqlite"
    outputFile <- "working/plotSelectivityBySize.pdf"
}

# parse input files
database <- connectBioassayDB(databaseFile)
selectivityCounts <- read.table(selectivityCountsFile)
drugbank_links <- read.csv(drugbank_linksFile)
drugCids <- unique(drugbank_links$PubChem.Compound.ID)
compoundSize <- read.table(compoundSizeFile)
highlyScreenedcids <- read.table(highlyScreened)[,1]
activeCids <- read.table(activeCidsFile)[[1]]

# drop compounds that aren't highly screened actives
selectivityCounts <- selectivityCounts[selectivityCounts[,1] %in% activeCids,]
selectivityCounts <- selectivityCounts[selectivityCounts[,1] %in% highlyScreenedcids,]
compoundSize <- compoundSize[compoundSize[,1] %in% highlyScreenedcids,]
compoundSize <- compoundSize[compoundSize[,1] %in% activeCids,]

# build stats on Size (as used in paper)
drugSize <- compoundSize[compoundSize[,1] %in% drugCids,2]
otherSize <- compoundSize[! compoundSize[,1] %in% drugCids,2]
# msvariation table
msvariation <- rbind(
    cbind(mean(drugSize), sd(drugSize), paste(range(drugSize), collapse = "-")),
    cbind(mean(otherSize), sd(otherSize), paste(range(otherSize), collapse = "-"))
)
msvariation
# size of compounds with 5 or more targets
range(compoundSize[compoundSize[,1] %in% selectivityCounts[selectivityCounts[,2] > 3,1],2])
# number of large and small compounds
sum(compoundSize[,2] < 3)
sum(compoundSize[,2] > 190)

# sort compoundSize by size and keep only highly screened actives
compoundSize <- compoundSize[order(compoundSize[,2], decreasing=T),]
compoundSize <- compoundSize[compoundSize[,1] %in% selectivityCounts[,1],]

# get biggest molecules from DrugBank (new version)
# drugbank_links_all <- read.csv("~/Downloads/drug links.csv")
# keep <- 40
# biggestDrugs <- compoundSize[compoundSize[,1] %in% drugbank_links_all$PubChem.Compound.ID,][1:keep,] 
# biggestDrugs <- merge(biggestDrugs, drugbank_links_all, by.x=1, by.y="PubChem.Compound.ID", all.x=F, all.y=F)
# biggestDrugs <- biggestDrugs[order(biggestDrugs[,2], decreasing=T),]

# get biggest molecules from DrugBank (old version)
keep <- 60
biggestDrugs <- compoundSize[compoundSize[,1] %in% drugbank_links$PubChem.Compound.ID,][1:keep,] 
biggestDrugs <- merge(biggestDrugs, drugbank_links, by.x=1, by.y="PubChem.Compound.ID", all.x=F, all.y=F)
biggestDrugs <- biggestDrugs[order(biggestDrugs[,2], decreasing=T),]
biggestDrugs <- biggestDrugs[,c(4, 1, 2)]
colnames(biggestDrugs) <- c("Name", "PubChem CID", "Heavy Atoms")

# get selectivity counts for biggest molecules
totalActive <- targetSelectivity(database, biggestDrugs[,2], scoring="total", category=FALSE)
fractionActive <- targetSelectivity(database, biggestDrugs[,2], scoring="fraction", category=FALSE)
totalTested <- as.integer(totalActive/fractionActive)
biggestDrugs <- cbind(biggestDrugs, totalActive, totalTested)
biggestDrugs <- biggestDrugs[biggestDrugs$totalTested > 9,]

# write latex table to file
xtmp <- xtable(biggestDrugs[1:40,], caption="Largest drugs", label="largestDrugs")
print(xtmp, type="latex", file="working/largestDrugs.tex", include.rownames=F)

# make data.frame for plotting
selectivityVsSize <- merge(selectivityCounts, compoundSize, by.x=1, by.y=1, all.x=F, all.y=F)
colnames(selectivityVsSize) <- c("cid", "selectivity", "Size")

# drop compounds active against > 10 targets
selectivityVsSize <- selectivityVsSize[selectivityVsSize$selectivity < 11,]

selectivityVsSizedrugs <- selectivityVsSize[selectivityVsSize$cid %in% drugCids,]
selectivityVsSizenondrugs <- selectivityVsSize[! selectivityVsSize$cid %in% drugCids,]

plotData <- rbind(
    cbind(selectivityVsSizedrugs, Compounds="FDA Approved"),
    cbind(selectivityVsSizenondrugs, Compounds="Other")
    )

# selectivityVsSize <- selectivityVsSize[1:1000,]

p1 <- ggplot(data=plotData, aes(x=factor(selectivity), y=Size, fill=Compounds)) +
    geom_violin(draw_quantiles=c(0.25, 0.5, 0.75)) +
    # xlab("Active Protein Targets") + 
    ylab("Heavy Atoms") +
    scale_fill_brewer(palette="Set1") +
    # coord_cartesian(ylim=c(0,2000)) +
    theme(
        text = element_text(size=15),
        panel.grid.major.y = element_line(colour = "grey"),
        panel.grid.major.x = element_blank(),
        axis.title.x=element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none",
        legend.key = element_blank()
    )
# plot(p1)
# ggsave(outputFile, plot=p, width=12)


p2 <- ggplot(data=plotData, aes(x=factor(selectivity), y=Size, fill=Compounds)) +
    geom_violin(draw_quantiles=c(0.25, 0.5, 0.75)) +
    xlab("Active Protein Targets") + 
    ylab("Heavy Atoms") +
    scale_fill_brewer(palette="Set1") +
    coord_cartesian(ylim=c(0,50)) +
    theme(
        text = element_text(size=15),
        panel.grid.major.y = element_line(colour = "grey"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="bottom",
        legend.key = element_blank()
    )
# plot(p2)

gridplot <- plot_grid(p1, p2, labels=c("A", "B"), ncol = 1, nrow = 2)
save_plot(outputFile, gridplot, base_width=12, base_height=7)

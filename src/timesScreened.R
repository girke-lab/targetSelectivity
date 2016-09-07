#!/usr/bin/env Rscript

# (C) 2015 Tyler William H Backman
# Purpose: 
# Tabulate the distribution of number of times screened for each compound in the
# database? Drugs vs non-drugs.

library(R.utils)
# library(Matrix)
library(bioassayR)
library(ggplot2)

# parse input options
databaseFile <- commandArgs(trailingOnly=TRUE)[1]
drugbank_linksFile <- commandArgs(trailingOnly=TRUE)[2]
outputFile <- commandArgs(trailingOnly=TRUE)[3]

# test code for running without make:
if(is.na(commandArgs(trailingOnly=TRUE)[1])){
    databaseFile <- "/dev/shm/bioassayDatabase.sqlite"
    # databaseFile <- "~/Downloads/pubchem_protein_only.sqlite"
    drugbank_linksFile <- "working/drugbank_links.csv"
    outputFile <- "working/timesScreened.csv"
}

# parse input files
database <- connectBioassayDB(databaseFile)
drugbank_links <- read.csv(drugbank_linksFile)
drugCids <- unique(drugbank_links$PubChem.Compound.ID)

# compute target screening participation
screenCounts <- queryBioassayDB(database, "SELECT cid, COUNT(DISTINCT target) AS targets FROM activity NATURAL JOIN targets GROUP BY cid")
save(list = c("screenCounts"), file = "working/timesScreened.RData")
# load("working/timesScreened.RData")

drugCounts <- screenCounts$targets[screenCounts$cid %in% drugCids]
otherCounts <- screenCounts$targets[! screenCounts$cid %in% drugCids]

# get highly screened counts
mean(drugCounts[drugCounts > 9])
median(drugCounts[drugCounts > 9])
mean(otherCounts[otherCounts > 9])
median(otherCounts[otherCounts > 9])

# tabulate frequency
cutoffs <- c(1,2,5,10,50,100,200,300,400,500,1000)
rangeLabels <- sapply(1:(length(cutoffs) - 1), function(x){
    if(x == 1) return("1")
    if(x == (length(cutoffs) - 1)) return(paste(cutoffs[x], "+", sep=""))
    paste(cutoffs[x], "-", (cutoffs[x+1] - 1),sep="")
}) 
structureBins <- table(cut(otherCounts, breaks=cutoffs, right=FALSE, labels=rangeLabels))

# tabulate frequency for FDA approved drugs
drugStructureBins <- table(cut(drugCounts, breaks=cutoffs, right=FALSE, labels=rangeLabels))

# format table and save
screenTable <- cbind(drugStructureBins, structureBins)
colnames(screenTable) <- c("FDA Approved Drugs", "Other Compounds")
write.csv(screenTable, outputFile)


# generate screening frequency plot
allCounts <- table(screenCounts$targets)
melted_data <- cbind(targets=names(allCounts), frac=allCounts)
melted_data <- as.data.frame(melted_data, stringsAsFactors = F)
# melted_data$label <- as.factor(melted_data$label)
melted_data$frac <- as.numeric(melted_data$frac)
melted_data$targets <- as.numeric(melted_data$targets)

# remove values above max limit
# melted_data <- melted_data[melted_data$targets < 20,]

p <- ggplot(melted_data, aes(x=targets, y=frac)) +
    geom_path(size=0.8, colour = "darkgreen") +
    # coord_trans(y = "log10") +
    # geom_point(data = filter(melted_data, Function == "Distinct Sequence Targets"), size=2) +
    ylab("Compounds") +
    xlab("Screened Protein Targets") +
    scale_colour_brewer(palette="Dark2") +
    annotation_logticks(sides="l") + 
    scale_y_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    geom_vline(aes(xintercept = 10), colour="grey", linetype = "longdash") + 
    # scale_x_continuous(breaks=1:20) + 
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
ggsave("working/timesScreened.pdf", plot=p, width=12)

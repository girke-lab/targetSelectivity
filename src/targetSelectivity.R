#e!/usr/bin/env Rscript

# (C) 2016 Tyler William H Backman
# Purpose: 
# plot distribution of target selectivity

library(R.utils)
library(ggplot2)
library(reshape)

# parse input options
selectivityCountsIndividualFile <- commandArgs(trailingOnly=TRUE)[1]
selectivityCountskClustFile <- commandArgs(trailingOnly=TRUE)[2]
selectivityCountsDomainsFile <- commandArgs(trailingOnly=TRUE)[3]
drugbank_linksFile <- commandArgs(trailingOnly=TRUE)[4]
outputFile <- commandArgs(trailingOnly=TRUE)[5]

# test code for running without make:
if(is.null(commandArgs(trailingOnly=TRUE)[1])){
    selectivityCountsIndividualFile <- "working/selectivityCountsIndividual.txt"
    selectivityCountskClustFile <- "working/selectivityCountskClust.txt"
    selectivityCountsDomainsFile <- "working/selectivityCountsdomains.txt"
    drugbank_linksFile <- "working/drugbank_links.csv"
    outputFile <- "working/targetSelectivitySequenceClusters.pdf"
}

# parse input files
selectivityCountsIndividual <- read.table(selectivityCountsIndividualFile)
selectivityCountskClust <- read.table(selectivityCountskClustFile)
selectivityCountsDomains <- read.table(selectivityCountsDomainsFile)
drugbank_links <- read.csv(drugbank_linksFile)
drugCids <- unique(drugbank_links$PubChem.Compound.ID)
drugCids <- drugCids[! is.na(drugCids)]

# keep only active compounds (they are already highly screened)
selectivityCountsIndividual <- selectivityCountsIndividual[selectivityCountsIndividual[,2] > 0,]
selectivityCountskClust <- selectivityCountskClust[selectivityCountskClust[,2] > 0,]
selectivityCountsDomains <- selectivityCountsDomains[selectivityCountsDomains[,2] > 0,]

# compute median target selectivities for manuscript
median(selectivityCountsIndividual[selectivityCountsIndividual[,1] %in% drugCids,2])
median(selectivityCountsIndividual[! selectivityCountsIndividual[,1] %in% drugCids,2])

median(selectivityCountskClust[selectivityCountskClust[,1] %in% drugCids,2])
median(selectivityCountskClust[! selectivityCountskClust[,1] %in% drugCids,2])

median(selectivityCountsDomains[selectivityCountsDomains[,1] %in% drugCids,2])
median(selectivityCountsDomains[! selectivityCountsDomains[,1] %in% drugCids,2])

# compute mean selectivities for manuscript
mean(selectivityCountsIndividual[selectivityCountsIndividual[,1] %in% drugCids,2])
sd(selectivityCountsIndividual[selectivityCountsIndividual[,1] %in% drugCids,2])
mean(selectivityCountsIndividual[! selectivityCountsIndividual[,1] %in% drugCids,2])
sd(selectivityCountsIndividual[! selectivityCountsIndividual[,1] %in% drugCids,2])

mean(selectivityCountskClust[selectivityCountskClust[,1] %in% drugCids,2])
sd(selectivityCountskClust[selectivityCountskClust[,1] %in% drugCids,2])
mean(selectivityCountskClust[! selectivityCountskClust[,1] %in% drugCids,2])
sd(selectivityCountskClust[! selectivityCountskClust[,1] %in% drugCids,2])

mean(selectivityCountsDomains[selectivityCountsDomains[,1] %in% drugCids,2])
sd(selectivityCountsDomains[selectivityCountsDomains[,1] %in% drugCids,2])
mean(selectivityCountsDomains[! selectivityCountsDomains[,1] %in% drugCids,2])
sd(selectivityCountsDomains[! selectivityCountsDomains[,1] %in% drugCids,2])


# compute trimmed mean selectivities for manuscript
trimmedCountsI <- selectivityCountsIndividual[selectivityCountsIndividual[,2] < 21,]
mean(trimmedCountsI[trimmedCountsI[,1] %in% drugCids,2])
sd(trimmedCountsI[trimmedCountsI[,1] %in% drugCids,2])
mean(trimmedCountsI[! trimmedCountsI[,1] %in% drugCids,2])
sd(trimmedCountsI[! trimmedCountsI[,1] %in% drugCids,2])

trimmedCountsK <- selectivityCountskClust[selectivityCountskClust[,2] < 21,]
mean(trimmedCountsK[trimmedCountsK[,1] %in% drugCids,2])
sd(trimmedCountsK[trimmedCountsK[,1] %in% drugCids,2])
mean(trimmedCountsK[! trimmedCountsK[,1] %in% drugCids,2])
sd(trimmedCountsK[! trimmedCountsK[,1] %in% drugCids,2])

trimmedCountsD <- selectivityCountsDomains[selectivityCountsDomains[,2] < 21,]
mean(trimmedCountsD[trimmedCountsD[,1] %in% drugCids,2])
sd(trimmedCountsD[trimmedCountsD[,1] %in% drugCids,2])
mean(trimmedCountsD[! trimmedCountsD[,1] %in% drugCids,2])
sd(trimmedCountsD[! trimmedCountsD[,1] %in% drugCids,2])

# ks.test
ks.test(selectivityCountsIndividual[selectivityCountsIndividual[,1] %in% drugCids,2], selectivityCountsIndividual[! selectivityCountsIndividual[,1] %in% drugCids,2])
ks.test(selectivityCountsDomains[selectivityCountsDomains[,1] %in% drugCids,2], selectivityCountsDomains[! selectivityCountsDomains[,1] %in% drugCids,2])

wilcox.test(selectivityCountsIndividual[selectivityCountsIndividual[,1] %in% drugCids,2], selectivityCountsIndividual[! selectivityCountsIndividual[,1] %in% drugCids,2], paired=F, alternative="greater")
wilcox.test(selectivityCountskClust[selectivityCountskClust[,1] %in% drugCids,2], selectivityCountskClust[! selectivityCountskClust[,1] %in% drugCids,2], paired=F, alternative="greater")
wilcox.test(selectivityCountsDomains[selectivityCountsDomains[,1] %in% drugCids,2], selectivityCountsDomains[! selectivityCountsDomains[,1] %in% drugCids,2], paired=F, alternative="greater")

normalizedCounts <- function(x, max=20, drug=T){
    x <- x[x[,2] > 0,]
    if(drug)
        x <- table(x[x[,1] %in% drugCids, 2])
    else
        x <- table(x[! x[,1] %in% drugCids, 2])
    x <- x/sum(x)
    x[1:max]
}

# drug counts
selectivityCountsIndividualDrugs <- normalizedCounts(selectivityCountsIndividual)
selectivityCountskClustDrugs <- normalizedCounts(selectivityCountskClust)
selectivityCountsDomainsDrugs <- normalizedCounts(selectivityCountsDomains)

# non-drug counts
selectivityCountsIndividualNondrugs <- normalizedCounts(selectivityCountsIndividual, drug=F)
selectivityCountskClustNondrugs <- normalizedCounts(selectivityCountskClust, drug=F)
selectivityCountsDomainsNondrugs <- normalizedCounts(selectivityCountsDomains, drug=F)

# make labels
totalDrugs <- sum(selectivityCountsIndividual[,1] %in% drugCids)
totalOther <- sum(! selectivityCountsIndividual[,1] %in% drugCids)
# labels <- c(paste("FDA Approved (", totalDrugs, ")", sep=""), paste("Other (", totalOther, ")", sep=""))
labels <- c("FDA Approved", "Non-FDA")
targetCategories <- c("Target Selectivity", "Cluster Selectivity", "Domain Selectivity")

# combine into ggPlot style data.frame 
drugCounts <- cbind(selectivityCountsIndividualDrugs, selectivityCountskClustDrugs, selectivityCountsDomainsDrugs)
drugCounts <- melt(drugCounts)
drugCounts <- cbind(drugCounts, type=labels[1])
nondrugCounts <- cbind(selectivityCountsIndividualNondrugs, selectivityCountskClustNondrugs, selectivityCountsDomainsNondrugs)
nondrugCounts <- melt(nondrugCounts)
nondrugCounts <- cbind(nondrugCounts, type=labels[2])
allCounts <- rbind(nondrugCounts, drugCounts)
names(allCounts) <- c("targets", "counter", "fraction", "Compounds")
allCounts$targets <- factor(allCounts$targets, levels=1:20)
allCounts$counter <- gsub("selectivityCounts(.*)[ND][a-z]*$", "\\1", as.character(allCounts$counter))
allCounts$counter[allCounts$counter == "Individual"] <- targetCategories[1]
allCounts$counter[allCounts$counter == "kClust"] <- targetCategories[2]
allCounts$counter[allCounts$counter == "Domains"] <- targetCategories[3]
allCounts$counter <- factor(allCounts$counter, levels=targetCategories)
allCounts$Compounds <- factor(allCounts$Compounds, levels=labels)

p <- ggplot(data=allCounts, aes(y=fraction, x=targets)) +
    geom_point(aes(color = Compounds), size=4) + 
    coord_trans(y = "log10") +
    facet_grid(. ~ counter, scales="fixed") + 
    guides(colour=guide_legend()) +
    scale_y_continuous(labels=function(n){format(n, scientific = FALSE)}) + 
    # xlab("Fraction of Total Compounds") + 
    ylab("Fraction of Total Compounds (log10)") + 
    xlab("Active Protein Targets") +
    # ggtitle("Target Selectivity Distribution") +
    theme(
        text = element_text(size=13),
        panel.grid.major.y = element_line(colour = "grey"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        legend.position="bottom",
        legend.key = element_blank()
    ) + annotation_logticks(scaled=F, sides="l") + 
    scale_x_discrete(breaks=c(1, 5, 10, 15, 20)) + 
    scale_y_continuous(breaks=c(0.01, 0.02, 0.10, 0.30)) + 
    scale_colour_brewer(palette="Set1")
plot(p)
ggsave("working/targetSelectivityFlipped.pdf", plot=p, width=12, height=4)

# make box plot of all counts
labels <- rep("FDA Approved", nrow(selectivityCountsIndividual))
labels[! selectivityCountsIndividual[,1] %in% drugCids] <- "non-FDA"
myDataI <- data.frame(compounds=labels, targets=selectivityCountsIndividual[,2], counter="Target Selectivity")

labels <- rep("FDA Approved", nrow(selectivityCountsIndividual))
labels[! selectivityCountskClust[,1] %in% drugCids] <- "non-FDA"
myDataK <- data.frame(compounds=labels, targets=selectivityCountskClust[,2], counter="Cluster Selectivity")

labels <- rep("FDA Approved", nrow(selectivityCountsIndividual))
labels[! selectivityCountsDomains[,1] %in% drugCids] <- "non-FDA"
myDataD <- data.frame(compounds=labels, targets=selectivityCountsDomains[,2], counter="Domain Selectivity")
# boxplot(targets ~ compounds, myData)
myData <- rbind(myDataI, myDataK, myDataD)
    
p <- ggplot(myData, aes(x=factor(compounds), y=targets, fill=compounds)) +
    geom_boxplot(outlier.shape=NA) + 
    coord_cartesian(ylim=c(0,16)) +
    facet_grid(. ~ counter, scales="fixed") +
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
ggsave("working/targetSelectivityDist.pdf", plot=p, width=12, height=4)

#!/usr/bin/env Rscript

# (C) 2016 Tyler William H Backman
# Purpose: 
# Fit functions to target selectivity function
# plot distribution of target selectivity

library(R.utils)
library(ggplot2)
library(reshape)
library(dplyr)
library(ggthemes)

# parse input options
selectivityCountskClustFile <- commandArgs(trailingOnly=TRUE)[1]
drugbank_linksFile <- commandArgs(trailingOnly=TRUE)[2]
outputFile <- commandArgs(trailingOnly=TRUE)[3]

# test code for running without make:
if(is.null(commandArgs(trailingOnly=TRUE)[1])){
    selectivityCountskClustFile <- "working/selectivityCountskClust.txt"
    drugbank_linksFile <- "working/drugbank_links.csv"
    outputFile <- "working/selectivityFit.pdf"
}

# parse input files
selectivityCountskClust <- read.table(selectivityCountskClustFile)
drugbank_links <- read.csv(drugbank_linksFile)
drugCids <- unique(drugbank_links$PubChem.Compound.ID)

# compute normalized drug and non-drug counts
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
selectivityCountskClustDrugs <- normalizedCounts(selectivityCountskClust)
# non-drug counts
selectivityCountskClustNondrugs <- normalizedCounts(selectivityCountskClust, drug=F)

# prepare statistical model
compoundPoints <- selectivityCountskClustNondrugs
# compoundPoints <- selectivityCountskClustDrugs # swap uncommenting with above line to show FDA approved drug data
compoundPoints <- cbind(targets=as.numeric(names(compoundPoints)), fraction=compoundPoints)
compoundPoints <- as.data.frame(compoundPoints)
targets.seq <- seq( from=1 , to=20 , length.out=1000 )

# stretched exponential fit
model <- nls(fraction ~ exp(-(targets * n)^ c), data=compoundPoints, start=list(c=0.75, n = 1))
x0 <- 1/summary(model)$coefficients["n","Estimate"]
c <- summary(model)$coefficients["c","Estimate"]
sexpPoints <- sapply(targets.seq, function(x) exp(-(x / x0)^ c))

# exponential fit
model <- nls(fraction ~ f*exp(-g * targets), data=compoundPoints, start=list(f=0.5, g=0.5))
f <- summary(model)$coefficients["f","Estimate"]
g <- summary(model)$coefficients["g","Estimate"]
epPoints <- sapply(targets.seq, function(x) (f*exp(-g * x)))

# power law fit
model <- nls(fraction ~ a * targets ^ -b, data=compoundPoints, start=list(a=0.5, b=1))
a <- summary(model)$coefficients["a","Estimate"]
b <- summary(model)$coefficients["b","Estimate"]
plPoints <- sapply(targets.seq, function(x) (a * x ^ -b))

# make ggplot2 style data frame
my_data <- data.frame(x = targets.seq, sexpPoints, plPoints, epPoints)
colnames(my_data) <- c("x", "Stretched Exponential", "Power Law", "Exponential")
melted_data <- melt(data = my_data, id.vars = "x")
raw_compoundPoints <- data.frame(x = compoundPoints$targets, ypoints=compoundPoints$fraction)
colnames(raw_compoundPoints) <- c("x", "Distinct Sequence Targets")
melted_compoundPoints <- melt(data = raw_compoundPoints, id.vars = "x")
melted_data <- rbind(melted_data, melted_compoundPoints)
colnames(melted_data) <- c("x", "Function", "value")

# plot fit
p <- ggplot(melted_data, aes(x=x, y=value, colour=Function)) +
  geom_path(data = filter(melted_data, Function %in% c("Stretched Exponential", "Power Law", "Exponential")), size=1.5) +
  geom_point(data = filter(melted_data, Function == "Distinct Sequence Targets"), size=2) +
  ylab("Fraction of Total Compounds") +
  xlab("Active Protein Targets") +
  scale_colour_brewer(palette="Dark2") +
  scale_x_continuous(breaks=1:20) + 
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
  ) # +
  # guides(fill = guide_legend(override.aes = list(size=3,linetype=1)))
# plot(p)

ggsave("working/selectivityFit.pdf", plot=p, width=12)

# compute R^2
coefficient_of_determination <- function(points, data, model){
  SStotal = sum((data - mean(data))^2)
  SSres = sum((data - model)^2)
  return(1 - SSres/SStotal)
}
sexpR2 <- coefficient_of_determination(compoundPoints$targets, compoundPoints$fraction, exp(-(compoundPoints$targets / x0)^ c))
plR2 <- coefficient_of_determination(compoundPoints$targets, compoundPoints$fraction, (a * compoundPoints$targets ^ -b))
epR2 <- coefficient_of_determination(compoundPoints$targets, compoundPoints$fraction, f*exp(-g * compoundPoints$targets))
c(sexpR2, plR2, epR2)

#!/usr/bin/env Rscript

# (C) 2016 Tyler William H Backman
# Purpose: delete all multi-target assays from a bioassayDB database

library(R.utils)
library(bioassayR)

bioassayDatabaseFile <- commandArgs(trailingOnly=TRUE)[1]

# test code for running without make:
if(is.na(commandArgs(trailingOnly=TRUE)[1])){
    bioassayDatabaseFile <- "/dev/shm/bioassayDatabase.sqlite"
}

database <- connectBioassayDB(bioassayDatabaseFile, writeable=T)
targetsPerAssay <- queryBioassayDB(database, "SELECT aid, COUNT(DISTINCT target) FROM targets GROUP BY aid")
duplicateAssays <- targetsPerAssay[targetsPerAssay[,2] > 1,1]
sapply(duplicateAssays, dropBioassay, database=database)
disconnectBioassayDB(database)

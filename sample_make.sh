#!/bin/bash -l

#PBS -j oe
#PBS -l walltime=48:00:00
#PBS -l nodes=1:ppn=64
#PBS -l mem=128gb
#PBS -q highmem

cd $PBS_O_WORKDIR

export cores="$PBS_NP"
export databaseFile="/dev/shm/bioassayDatabase.sqlite"
export drugBankUsername="putDrugBankEmailHere"
export drugBankPassword="putDrugBankPasswordHere"
module load R/3.2.2
module load hmmer/3.1b2
make -e /dev/shm/bioassayDatabase.sqlite
make -e all
rm /dev/shm/bioassayDatabase.sqlite

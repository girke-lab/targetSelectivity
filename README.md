[![DOI](https://zenodo.org/badge/65675614.svg)](https://zenodo.org/badge/latestdoi/65675614)

<a rel="license" href="http://creativecommons.org/licenses/by/3.0/us/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/3.0/us/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/3.0/us/">Creative Commons Attribution 3.0 United States License</a>.

# Druggable Proteome Small Molecule Bioactivity Analysis Pipeline
This repository contains an analysis pipeline for investigating patterns of small molecule bioactivity
and target druggability in large public bioactivity data.
Note that this pipeline downloads the latest data from online databases,
and as a result the output may change over time as these public resources are updated.

## Citation

Backman TW, Evans DS, Girke T.  
Large-scale bioactivity analysis of the small-molecule assayed proteome.  
PLOS ONE. 2017 Feb 8;12(2):e0171413.  
http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0171413

## System requirements
This analysis pipeline is designed to run in parallel on a CentOS linux system. We ran this on the UC Riverside Biocluster system, on a node with 64 AMD Opteron 6376 2.3ghz cpu cores, and 512GB of ram. 

## Use
* Export (with the Bash shell) a viable ramdisk location for the database as shown in sample_make.sh
* Export (with the Bash shell) the number of cpu cores (64 or more recommended) as shown in sample_make.sh
* Install any needed R packages or other software
* Run 'make -e working/<filename>' to create any desired output file described in the makefiles

## Overall folder structure
* **/** - Makefiles, describing all analysis steps
* **/src** - Source code for each analysis step
* **/working** - Temporary folder with results of each analysis step

## This data analysis pipeline requires many R packages, including the following custom R packages
* https://github.com/TylerBackman/BicBin
* https://github.com/TylerBackman/clusteval/tree/biggerInts

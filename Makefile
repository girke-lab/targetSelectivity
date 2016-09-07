################################################################################
# (C) 2016 Tyler W. H. Backman
# Master makefile
################################################################################

# User Options -set up environment
cores = 1 
databaseFile = working/bioassayDatabase.sqlite
# databaseFile = /dev/shm/bioassayDatabase.sqlite # uncomment for linux ramdisk
dataUrl = https://biocluster.ucr.edu/~tbackman
drugBankUsername = putDrugBankEmailHere
drugBankPassword = putDrugBankPasswordHere

##########################################
# download external dependencies 
# Note: all external files should go here
##########################################

# download prebuilt bioassayR PubChem Bioassay database
working/bioassayDatabaseDownloaded.sqlite:
	mkdir -p working
	wget $(dataUrl)/bioassayR/pubchem_protein_only.sqlite -O $@ --no-check-certificate

# download bioassayR protein targets 
working/targets.fasta:
	mkdir -p working
	wget $(dataUrl)/bioassayR/targets.fasta -O $@ --no-check-certificate

# download single SDF file of all active bioassayR compounds
working/bioassayCompounds.sdf:
	mkdir -p working
	wget $(dataUrl)/bioassayR/activeCompounds.sdf -O $@ --no-check-certificate

# download SDF files of all active bioassayR compounds
working/bioassayCompounds:
	mkdir -p $@
	wget $(dataUrl)/bioassayR/activeCompoundsSplit.tgz -O $@.tgz --no-check-certificate
	cd working && tar xfz bioassayCompounds.tgz
	mv working/splitFolder $@
	rm working/bioassayCompounds.tgz

# get structures of FDA approved drugs
working/drugbank.sdf:
	curl -L -o $@.zip -u $(drugBankUsername):$(drugBankPassword) http://www.drugbank.ca/releases/5-0-1/downloads/approved-structures
	unzip $@.zip -d working/
	mv working/structures.sdf $@

# download drugbank target sequences
working/drugbank_targets.fasta:
	mkdir -p working
	curl -L -o $@.zip -u $(drugBankUsername):$(drugBankPassword) http://www.drugbank.ca/releases/5-0-1/downloads/target-approved-polypeptide-sequences
	unzip $@.zip -d working/
	mv working/protein.fasta $@

# download DrugBank FDA Approved External Drug Links
working/drugbank_links.csv:
	curl -L -o $@.zip -u $(drugBankUsername):$(drugBankPassword) http://www.drugbank.ca/releases/5-0-1/downloads/approved-drug-links
	unzip $@.zip -d working/
	mv working/drug\ links.csv $@

# download annotated drugbank targets from DrugBank
working/drug_target_uniprot_links.csv:
	curl -L -o $@.zip -u $(drugBankUsername):$(drugBankPassword) http://www.drugbank.ca/releases/5-0-1/downloads/target-approved-uniprot-links
	unzip $@.zip -d working/
	mv working/uniprot\ links.csv $@

# download Pfam-A HMMs in an HMM library searchable with the hmmscan program 
working/Pfam-A.hmm:
	wget -O $@.gz ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam29.0/Pfam-A.hmm.gz
	gunzip $@.gz
	hmmpress $@

# download a tab separated file containing Pfam-A family and clan information for all Pfam-A families
working/Pfam-A.clans.tsv:
	wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam29.0/Pfam-A.clans.tsv.gz -O $@.gz
	gunzip $@.gz
	dos2unix $@

# download UniProtKB/Swiss-Prot human proteome 
working/UP000005640_9606.fasta:
	wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640_9606.fasta.gz -O $@.gz
	gunzip $@.gz

# download GO annotations for Pfam domains
working/pfam2go:
	wget http://geneontology.org/external2go/pfam2go -O $@

# download Generic GO slim Developed by GO Consortium 
working/goslim_generic.obo:
	wget http://www.geneontology.org/ontology/subsets/goslim_generic.obo -O $@

##########################################
# Organize and process downloaded files
##########################################

# remove all multi-target assays
working/bioassayDatabaseSingleTarget.sqlite: src/singleTargetOnly.R working/bioassayDatabaseDownloaded.sqlite
	cp working/bioassayDatabaseDownloaded.sqlite $@
	$< $@

# create optional memory cached sqlite database if databaseFile path is changed
$(databaseFile): working/bioassayDatabaseSingleTarget.sqlite
	cp -p $< $@

########################################################
# Structurally cluster FDA approved drugs as a reference
########################################################

working/structureClusterDrugs.tab: src/binningClustering.R working/drugbank.sdf working/drugbank_links.csv
	$^ $@

#########################################################
# Cluster protein targets by similarity and annotate them
#########################################################

# combine drugbank targets and bioassay targets
working/combinedTargets.fasta: working/drugbank_targets.fasta working/targets.fasta
	cp $< $@
	cat working/targets.fasta >> $@

# cluster combined protein targets
working/combinedCluster: working/combinedTargets.fasta
	mkdir -p $@
	src/kClust -i $< -d $@ -s 2.93 -M 16000MB

# parse kclust output and fix identifier names
working/curatedClusters.txt: src/curateClusters.R working/combinedCluster
	$^ $@

# annotate kclust clusters via biomaRt
# choose a representative protein target for each cluster with the following priority order:
# (1) already annotated in DrugBank, (2) human, (3) non-human
working/clusterAnnotations.csv: src/annotateProteinClusters.R working/curatedClusters.txt $(databaseFile) 
	$^ $@

# find GO terms for all annotated protein targets
working/clusterGOannotations.csv: src/clusterGOannotations.R working/clusterAnnotations.csv
	$^ $@

# find GO Slim annotations
working/clusterGOslimAnnotations.csv: src/clusterGOslimAnnotations.R working/clusterGOannotations.csv working/goslim_generic.obo
	$^ $@

# download ontology data from UniProt
working/gene_association.goa_uniprot:
	wget ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz -O $@.gz
	gunzip $@.gz

##############################################
# find Pfam domains within all protein targets 
##############################################

# find Pfam domains within all combined targets
working/combinedTargetDomains: working/Pfam-A.hmm working/combinedTargets.fasta
	hmmscan -E 0.01 --domE 0.01 --tblout $@ --cpu $(cores) --noali $^

# simplify domain data into two columns
working/combinedTargetDomainsTwoCols: working/combinedTargetDomains
	awk '{ if (!/^#/) print $$2 " " $$3}' $^ > $@

# get stats on pfam domains
working/Pfam-A-stats.txt: working/Pfam-A.hmm
	hmmstat $< > $@

# get table of HMM residue lengths
working/PfamResidueLengths.tab: src/PfamResidueLengths.R working/Pfam-A-stats.txt
	$^ $@

##############################################
# find Pfam domains within human proteome 
##############################################

# find Pfam domains within all human targets
working/humanDomains: working/Pfam-A.hmm working/UP000005640_9606.fasta
	hmmscan -E 0.01 --domE 0.01 --tblout $@ --cpu $(cores) --noali $^

# simplify domain data into two columns
working/humanDomainsTwoCols: working/humanDomains
	awk '{ if (!/^#/) print $$2 " " $$3}' $^ > $@
	
######################################################
# Summarize bioactivity data
######################################################

# make list of cids screened at least 10 times
working/highlyScreenedCids.txt: src/highlyScreened.R $(databaseFile) 
	$^ $@ 10

# make list of all active compounds 
working/activeCids.txt: src/getActives.R $(databaseFile) 
	$^ $@

# make matrix of highly screened compounds VS protein target clusters
working/cidsVStargetClusters.RData: src/cidsVStargetMatrix.R $(databaseFile) working/highlyScreenedCids.txt working/curatedClusters.txt
	$^ $@ $(cores)

# make matrix of all compounds VS protein targets (not clustered)
working/cidsVStargets.RData: src/cidsVStargetMatrix.R $(databaseFile) 
	$^ none none $@ $(cores)

# use binary biclustering to extract the largest fully screened submatrix
working/fullyScreened.RData: src/extractFullyScreenedBicBin.R working/cidsVStargetClusters.RData
	$^ $@ $(cores)
	
##########################################################
# include per-section makefiles
# Note: these should not depend on one
#       another. Put all shared dependencies in this file.
##########################################################

include Makefile_dataQuality
include Makefile_targetSelectivity
include Makefile_promiscuityModel
include Makefile_annotationAndBiclustering
include Makefile_targetNetwork
include Makefile_SupportingInfo

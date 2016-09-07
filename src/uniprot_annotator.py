#!/usr/bin/python

# (C) 2016 Tyler William H Backman
# Purpose: extract UniProt GO ontologies for each protein target
# in a bioassayR database

import argparse
import os
import sys
import sqlite3
import csv
import re 

parser = \
    argparse.ArgumentParser(description='extract GO terms')
parser.add_argument('-d', '--database', help='bioassayR database',
    required=True)
parser.add_argument('-u', '--uniprot', help='uniprot GO terms',
    required=True)
parser.add_argument('-o', '--outfile', help='output file',
    required=True)
args = vars(parser.parse_args())

# args = {'database': 'working/bioassayDatabase.sqlite',
#        'uniprot': 'working/gene_association.goa_uniprot',
#        'outfile': 'working/targetGOannotations.csv'};

# get uniprot IDs to translate
con = sqlite3.connect(args['database'])
cur = con.cursor()
cur.execute('SELECT DISTINCT identifier FROM targetTranslations WHERE category = "UniProt"')
idList = cur.fetchall()
con.close()
idDict = {}
for id in idList:
    idDict[id[0]] = " "

# parse CSV
uniprotFile = open(args['uniprot'], 'rb')
uniprotIter = csv.reader(uniprotFile, delimiter="\t") 
outputFile = open(args['outfile'], 'wb')
outputIter = csv.writer(outputFile)
for line in uniprotIter:
    if len(line) < 2:
        continue
    if idDict.get(line[1]):
        outputList = [line[1], line[4]]
        outputIter.writerow(outputList)
uniprotFile.close()
outputFile.close()

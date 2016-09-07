#!/usr/bin/python

# (C) 2016 Tyler William H Backman
# Purpose: find PAINS in input compounds

import argparse
from rdkit import Chem
import sys
import string
from rdkit.Chem import FilterCatalog

parser = \
    argparse.ArgumentParser(description='find pains')
parser.add_argument('-i', '--infile', help='input sdf',
    required=True)
parser.add_argument('-o', '--outfile', help='output file',
    required=True)
parser.add_argument('-u', '--unparseable', help='output file for unparsable ids',
    required=True)
args = vars(parser.parse_args())
infile = args['infile']
outfile = args['outfile']
ufile = args['unparseable']

suppl = Chem.SDMolSupplier(infile)

params = FilterCatalog.FilterCatalogParams()
params.AddCatalog(params.FilterCatalogs.PAINS_A)
params.AddCatalog(params.FilterCatalogs.PAINS_B)
params.AddCatalog(params.FilterCatalogs.PAINS_C)
catalog = FilterCatalog.FilterCatalog(params)

f = open(outfile, 'w')
u = open(ufile, 'w')
count = 0
for mol in suppl:
    count += 1
    if mol is None:
        u.write(str(count) + "\n")
        continue
    try:
        if catalog.HasMatch(mol):
            f.write(mol.GetProp("_Name") + "\n")
    except:
        print("error in compound " + mol.GetProp("_Name") + "\n")

f.close()
u.close()

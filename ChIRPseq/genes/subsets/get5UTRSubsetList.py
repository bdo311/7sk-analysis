# get5UTRSubsetList.py
# 2/6/15
# get a list of 5UTR lines (from extract_gene_regions file) based on genes listed in the subset file for TSS/TES

import csv, sys
csv.register_dialect("textdialect",delimiter='\t')

subsets = set()
with open(sys.argv[1],'r') as subsetFile:
	for line in subsetFile:
		subsets.add(line.rstrip())

with open(sys.argv[2],'r') as utrFile, open(sys.argv[3],'w') as subsetOut:
	reader = csv.reader(utrFile, 'textdialect')
	for row in reader:
		name = row[3]
		geneName = name.split('__')[-1]
		if geneName in subsets: subsetOut.write(name+'\n')




# normalizeGroBedgraph.py
# 8/20/14
# normalizes all gro bedgraphs to 10 million reads

import csv, sys
csv.register_dialect("textdialect", delimiter='\t')
csv.register_dialect("spacedialect", delimiter=' ')

ifn = sys.argv[1]
ofn = sys.argv[2]
numReads = int(sys.argv[3])
multFactor = 10000000/float(numReads)

# normalize by scaling everything to a total read density of 750 million
ifile = open(ifn, 'r')
reader = csv.reader(ifile, 'textdialect')

ofile = open(ofn, 'w')
writer = csv.writer(ofile, 'textdialect')

for row in reader:
	#if '_' in row[0]: continue
	if float(row[-1]) == 0: continue
	row[-1] = float(row[-1]) * multFactor
	writer.writerow(row)
ifile.close()
ofile.close()


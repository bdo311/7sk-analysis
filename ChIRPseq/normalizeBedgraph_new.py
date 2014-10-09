# normalizeBedgraph_new.py
# 6/16/14
# normalizes one bedgraph at a time

import csv, sys
csv.register_dialect("textdialect", delimiter='\t')
csv.register_dialect("spacedialect", delimiter=' ')

ifn = sys.argv[1]
ofn = sys.argv[2]

# get biggest allowed sequence sizes
ifile = open("/seq/chromosome/hg19/hg19.sizes", 'r')
reader = csv.reader(ifile, 'textdialect')

chrToLimit = {}
for row in reader:
	chrToLimit[row[0]] = int(row[1])

ifile.close()

# normalize by scaling everything to a total read density of 750 million
ifile = open(ifn, 'r')
reader = csv.reader(ifile, 'textdialect')

sum = 0
for row in reader:
	if '_' in row[0]: continue
	sum += (int(row[2]) - int(row[1])) * float(row[3])

multFactor = 750000000/sum

ifile.seek(0)

ofile = open(ofn, 'w')
writer = csv.writer(ofile, 'textdialect')

for row in reader:
	if '_' in row[0]: continue
	if float(row[-1]) == 0: continue
	if int(row[2]) > chrToLimit[row[0]]: continue
	row[-1] = float(row[-1]) * multFactor
	writer.writerow(row)
ifile.close()
ofile.close()


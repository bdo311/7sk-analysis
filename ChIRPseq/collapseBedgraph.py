# collapseBedgraph.py
# 6/20/14
# deletes redundant rows

import csv, sys
csv.register_dialect("textdialect", delimiter='\t')
csv.register_dialect("spacedialect", delimiter=' ')

# get biggest allowed sequence sizes
ifile = open("/seq/chromosome/mm9/mm9.sizes", 'r')
reader = csv.reader(ifile, 'textdialect')

chrToLimit = {}
for row in reader:
	chrToLimit[row[0]] = int(row[1])

ifile.close()

ifn = sys.argv[1]
ofn = sys.argv[2]

ifile = open(ifn, 'r')
reader = csv.reader(ifile, 'textdialect')

ofile = open(ofn, 'w')
writer = csv.writer(ofile, 'textdialect')

currRow = ['', 0, 0, 0]
for row in reader:
	if '_' in row[0]: continue
	if float(row[-1]) == 0: continue
	if int(row[2]) > chrToLimit[row[0]]: continue

	if float(row[-1]) == float(currRow[-1]) and int(row[1]) <= int(currRow[2]): 
		currRow[2] = row[2]
	else:
		if float(currRow[-1]) > 0: writer.writerow(currRow)
		currRow = row

ifile.close()
ofile.close()

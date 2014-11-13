# getProportionRegions.py
# 5/16/14
# gets length of regions in the genome

import csv, sys, collections
csv.register_dialect("textdialect", delimiter='\t')

fn = sys.argv[1]
prop = collections.defaultdict(lambda: [0, 0])	# [num of regions, length]
total = 0

ifile = open(fn, 'r')
reader = csv.reader(ifile, 'textdialect')

for row in reader:
	length = int(row[3]) - int(row[2])
	if int(row[3]) == 1500000000: continue #for hg19
	if length == 1: continue
	if 'random' in row[1]: continue #for mm9
	prop[row[5]][0] = prop[row[5]][0] + 1
	prop[row[5]][1] = prop[row[5]][1] + length
	total = total + length
	
ifile.close()
print total

ofile = open(fn[:-4] + "_prop.txt", 'w')
writer = csv.writer(ofile, 'textdialect')

header = ["region", "number", "total length", "proportion"]
writer.writerow(header)
for typeReg in prop:
	outputRow = [typeReg, prop[typeReg][0], prop[typeReg][1]]
	outputRow.append(prop[typeReg][1]/float(total))
	writer.writerow(outputRow)
	
ofile.close()
	

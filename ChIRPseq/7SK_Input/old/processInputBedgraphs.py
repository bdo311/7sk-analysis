# File: processInputBedgraphs.py
# 4/23/14
# Scales input bedgraphs by the right amount

import csv, glob
csv.register_dialect("textdialect", delimiter='\t')

filenames = glob.glob("chr*.bedGraph")

scaleFactor = 10000000/xxx

for fn in filenames:
	newfn = fn[:-9] + "_new.bedGraph"
	ifile = open(fn, 'r')
	reader = csv.reader(ifile, 'textdialect')
	ofile = open(newfn, 'w')
	writer = csv.writer(ofile, 'textdialect')
	
	for row in reader:
		outputRow = row[:-1]
		outputRow.append(float(row[-1]) * scaleFactor)
		writer.writerow(outputRow)
		
	ifile.close()
	ofile.close()
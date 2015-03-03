# nanoEnhancers.py
# 4/15/14
# Simply cuts _exp5 in averages to _exp3

import csv
csv.register_dialect("textdialect", delimiter='\t')
files = ["raw_super_exp5.txt", "scaled_super_exp5.txt"]

dir = "/home/raflynn/7SK_ChIRPseq/averages/"
for fn in files:
	ifile = open(dir + fn, 'r')
	reader = csv.reader(ifile, 'textdialect')
	
	ofile = open(dir + fn[:-5] + '3.txt', 'w')
	writer = csv.writer(ofile, 'textdialect')

	for row in reader:
		outputRow = [row[0]]
		outputRow.extend(row[100:301])
		writer.writerow(outputRow)
	ifile.close()
	ofile.close()

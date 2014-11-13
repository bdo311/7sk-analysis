# extendReg.py
# 5/26/14
# extends regular enhancers to +/- 10kb from the center

import csv
csv.register_dialect("textdialect", delimiter='\t')

ifile = open("hg19_H1_reg_enhancers.bed", 'r')
reader = csv.reader(ifile, 'textdialect')

ofile = open("hg19_H1_reg_enhancers_centered10.bed", 'w')
writer = csv.writer(ofile, 'textdialect')

for row in reader:
	ctr = (int(row[1]) + int(row[2]))/2
	row[1] = ctr - 10000
	row[2] = ctr + 10000
	writer.writerow(row)
	
ifile.close()
ofile.close()

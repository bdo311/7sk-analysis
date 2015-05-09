# makeTransPausedBed.py
# 4/1/15

import csv, sys
csv.register_dialect("textdialect",delimiter='\t')

ifile = open(sys.argv[2],'r')
reader = csv.reader(ifile, 'textdialect')
genes = set()
for row in reader:
	genes.add(row[0])
ifile.close()

ifile = open(sys.argv[1],'r')
reader = csv.reader(ifile, 'textdialect')
ofile = open(sys.argv[3],'w')
writer = csv.writer(ofile, 'textdialect')
for row in reader:
	if row[3] in genes: writer.writerow(row)
ifile.close()
ofile.close()
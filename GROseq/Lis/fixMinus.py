# fixMinus.py
# python fixMinus.py <ifile> <ofile>

import csv, sys
csv.register_dialect("textdialect",delimiter='\t')

ifile = open(sys.argv[1],'r')
reader = csv.reader(ifile,'textdialect')
ofile = open(sys.argv[2],'w')
writer = csv.writer(ofile,'textdialect')

reader.next()
reader.next()
for row in reader:
	row[-1]=-float(row[-1])
	writer.writerow(row)

ifile.close()
ofile.close()

# expandEnhancers.py
# 5/27/14
# expands enhancers to n-fold its current size
# usage: python expandEnhancers.py <ifn> <ofn> n

import csv, sys
csv.register_dialect('textdialect', delimiter='\t')

ifn, ofn = sys.argv[1], sys.argv[2]
n = int(sys.argv[3])

ifile = open(ifn, 'r')
reader = csv.reader(ifile, 'textdialect')
ofile = open(ofn, 'w')
writer = csv.writer(ofile, 'textdialect')

for row in reader:
	length = int(row[3]) - int(row[2])
	start = int(row[2]) - (n-1)/2*length
	end = int(row[3]) + (n-1)/2*length
	writer.writerow([row[0], row[1], start, end, row[4], row[5]])
	
ifile.close()
ofile.close()
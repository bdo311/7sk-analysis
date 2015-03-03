# processIndivSE.py
# renames individual human SEs; makes 1kb

import csv, sys, collections
csv.register_dialect("textdialect",delimiter='\t')

ifile = open(sys.argv[1],'r')
reader = csv.reader(ifile, 'textdialect')
ofile = open(sys.argv[2],'w')
writer = csv.writer(ofile, 'textdialect')

enhToCount = collections.defaultdict(lambda: 0)
for row in reader:
	enhName = '__'.join([row[3].split('__')[1],row[1],row[2]])
	ct = enhToCount[enhName] + 1
	enhToCount[enhName] = ct
	enhName = enhName + '__' + str(ct)
	writer.writerow([row[6],row[7],row[8],enhName, 0, '+'])
	
ifile.close()
ofile.close()
	
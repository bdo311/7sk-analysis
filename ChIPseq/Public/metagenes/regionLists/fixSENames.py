# fixSENames.py

import csv, sys
csv.register_dialect("textdialect",delimiter='\t')

ifn=sys.argv[1]
ifile = open(ifn,'r')
reader = csv.reader(ifile, 'textdialect')
ofn=sys.argv[2]
ofile = open(ofn,'w')
writer = csv.writer(ofile, 'textdialect')

writer.writerow(reader.next())
for row in reader:
	name = row[0].split('__')
	name = '__'.join(name[3:7])
	row[0] = name
	writer.writerow(row)

ifile.close()
ofile.close()

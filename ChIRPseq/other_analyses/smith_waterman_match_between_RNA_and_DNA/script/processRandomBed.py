import csv, sys

fn = sys.argv[1]
csv.register_dialect("textdialect", delimiter='\t')

ifile = open(fn, 'r')
reader = csv.reader(ifile, 'textdialect')

ofn = fn[:-4] + '_new.bed'
ofile = open(ofn, 'w')
writer = csv.writer(ofile, 'textdialect')

for row in reader:
	if row[4] != '0' and row[4] != '1': writer.writerow(row)
	
ifile.close()
ofile.close()

import csv, sys
csv.register_dialect("textdialect",delimiter='\t')

ifile = open(sys.argv[1],'r')
reader = csv.reader(ifile, 'textdialect')
ofile = open(sys.argv[2],'w')
writer = csv.writer(ofile, 'textdialect')

alreadyHere = set()

for row in reader:
	id = '__'.join([row[0],row[1],row[2],row[5]])
	if id not in alreadyHere:
		alreadyHere.add(id)
		writer.writerow(row)

ifile.close()
ofile.close()

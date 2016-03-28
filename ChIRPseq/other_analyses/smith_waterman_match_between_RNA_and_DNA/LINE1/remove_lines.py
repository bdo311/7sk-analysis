import csv, sys
csv.register_dialect("textdialect", delimiter='\t')

with open(sys.argv[1],'r') as a, open(sys.argv[2],'r') as b, open(sys.argv[3],'w') as c:
	ar = csv.reader(a, 'textdialect')
	br = csv.reader(b, 'textdialect')
	cw = csv.writer(c, 'textdialect')

	chr_to_max = {}
	for row in br:
		chr_to_max[row[0]] = int(row[1])

	for row in ar:
		start = int(row[1])
		end = int(row[2])
		if row[0] not in chr_to_max: continue
		if chr_to_max[row[0]] > start and chr_to_max[row[0]] > end:
			cw.writerow(row)
			

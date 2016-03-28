import csv, sys, math
csv.register_dialect("textdialect", delimiter='\t')

ifile = open(sys.argv[1], 'r')
reader = csv.reader(ifile, 'textdialect')
ofile = open(sys.argv[2], 'w')
writer = csv.writer(ofile, 'textdialect')

for row in reader:
	val1 = float(row[6])
	val2 = float(row[15])
	if math.fabs(val1-val2)/val1 > 0.5 or math.fabs(val1-val2)/val2 > 0.5: continue
	writer.writerow(row)

ifile.close()
ofile.close()

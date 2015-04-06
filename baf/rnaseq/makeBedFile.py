# makeBedFile.py
# make a good bed file from the Fazzio directional DHS BED file Sarah sent to me on 3/29

import csv, sys
csv.register_dialect('textdialect',delimiter='\t')

ifile = open(sys.argv[1],'rU')
reader = csv.reader(ifile, 'textdialect')
reader.next()

ofile = open(sys.argv[2],'w')
writer = csv.writer(ofile, 'textdialect')

for row in reader:
	name = '__'.join(row)
	strand = row[0]
	info = row[1]
	chr = 'chr' + info.split(':')[0]
	start, end = [int(x) for x in info.split(':')[1].split('-')]
	middle = (start + end)/2
	if strand == '+':
		start = middle
		end = middle + 1000
	else:
		start = middle - 1000
		end = middle
	writer.writerow([chr, start, end, name, 0, strand])
ifile.close()
ofile.close()
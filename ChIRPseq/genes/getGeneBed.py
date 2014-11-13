#getGeneBed.py
#5/8/14
#gets bed file for genes

import csv
csv.register_dialect("textdialect", delimiter='\t')

ifile = open("mm9_refseq_coll.txt", 'r')
reader = csv.reader(ifile, 'textdialect')

ofile = open("mm9_refseq_coll_gene.bed", 'w')
writer = csv.writer(ofile, 'textdialect')

reader.next()
for row in reader:
	chr = row[1]
	if 'random' in chr: continue
		
	writer.writerow([chr, row[3], row[4], row[0], 0, row[2]])

ifile.close()
ofile.close()
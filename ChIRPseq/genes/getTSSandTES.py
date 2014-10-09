# getTSSandTES.py
# 2/25/14
# Takes mm9_refseq_coll.txt and mm9_ensembl_coll.txt and gets TSS and TES 2000-bp windows (1000bp on each side). For purposes of later processing, start and end of TSS and TES windows will be listed in increasing numerical order. Strand will be included so that negative strand coords can be flipped.

import csv

csv.register_dialect("textdialect", delimiter='\t')

ifile = open("hg19_refseq_coll.txt", 'r')
ofile = open("hg19_refseq_startend.txt", 'w')

reader = csv.reader(ifile, 'textdialect')
writer = csv.writer(ofile, 'textdialect')

reader.next()
writer.writerow(['Gene name', 'Chr', 'Strand', 'Start start', 'Start end', 'End start', 'End end', 'Names'])

for row in reader:
	outputRow = row[:3]
	tssStart, tssEnd = int(row[3]) - 1000, int(row[3]) + 1000 #TSS window
	tesStart, tesEnd = int(row[4]) - 1000, int(row[4]) + 1000 #TES window
	outputRow.extend([tssStart, tssEnd, tesStart, tesEnd, row[-1]])

	writer.writerow(outputRow)

ifile.close()
ofile.close()
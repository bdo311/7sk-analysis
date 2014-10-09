# getBed.py
# 5/28/14
# gets bed files from -1000 to +100 of a given gene list

import csv, sys
csv.register_dialect("textdialect", delimiter='\t')

# read in gene files
genefile = open("/home/raflynn/ChIRPseq/genes/mm9_refseq_coll.txt", 'r')
reader = csv.reader(genefile, 'textdialect')

geneToInfo = {}
reader.next()
for row in reader:
	geneToInfo[row[0]] = row[:5]
genefile.close()
	
# write bed file
ifn, ofn = sys.argv[1], sys.argv[2]

ifile = open(ifn, 'r')
reader = csv.reader(ifile, 'textdialect')
ofile = open(ofn, 'w')
writer = csv.writer(ofile, 'textdialect')

for row in reader:
	gene = row[0]
	if gene not in geneToInfo: continue
	info = geneToInfo[gene]
	
	if info[2] == '+': 
		tss = int(info[3])
		start, end = tss - 1000, tss + 100
	else: 
		tss = int(info[4])
		start, end = tss - 100, tss + 1000
		
	writer.writerow([info[1], start, end, info[0], 0, info[2]])
ifile.close()
ofile.close()
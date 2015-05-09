# makeTRDist.py
# makes distal TR file for the start-seq centered TSS's

import csv
csv.register_dialect("textdialect",delimiter='\t')


# get gene ends
ifile = open("mm9_refseq_coll.txt",'r')
reader = csv.reader(ifile, 'textdialect')

geneToEnd = {}
reader.next()
for row in reader:
	if row[2]=="+": geneToEnd[row[0]] = int(row[4])
	else: geneToEnd[row[0]] = int(row[3])
	
ifile.close()

# get gene starts and write file
ifile = open("mm9_tss_startseq_centered_TR_prox.bed",'r')
reader = csv.reader(ifile, 'textdialect')
ofile = open("mm9_tss_startseq_centered_TR_dist.bed",'w')
writer = csv.writer(ofile, 'textdialect')

for row in reader:
	start = int(row[1]) + 150
	end = geneToEnd[row[3]]
	if row[5]=="+": 
		if start+250>end: end = start + 251
		writer.writerow([row[0],start+250,end,row[3],row[4],row[5]])
	else: 
		if end > (start-250): end = start - 251
		writer.writerow([row[0],end,start-250,row[3],row[4],row[5]])
	
ifile.close()
ofile.close()
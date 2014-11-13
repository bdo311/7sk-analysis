# randomSeq.py
# 4/7/14
# gets 75K random 150-bp peaks to use to test the mean SW score with 7SK

import csv, random

csv.register_dialect("textdialect", delimiter='\t')

# reading in sizes
ifile = open("/home/raflynn/7SK_ChIRPseq/genes/mm9.sizes", 'r')
reader = csv.reader(ifile, 'textdialect')

chroms = []
chrToSize = {}
for row in reader:
	if row[0] == 'chrM': continue
	chrToSize[row[0]] = int(row[1])
	chroms.append(row[0])
	
ifile.close()

# making the bed file
ofile = open("/home/raflynn/7SK_ChIRPseq/smith_waterman_match_between_RNA_and_DNA/random.bed", 'w')
writer = csv.writer(ofile, 'textdialect')

for i in range(75000):
	chrom = random.choice(chroms)
	numStart = random.randint(1, chrToSize[chrom] - 200)
	numEnd = numStart + 150
	
	writer.writerow([chrom, numStart, numEnd, 'peak_' + str(i), 0])

ofile.close()
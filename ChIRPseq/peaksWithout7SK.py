# peaksWithout7SK.py
# 3/28/14
# Removes peaks from the MACS peak file that have >= 20-nt overlap with 7SK loci, in order to run the SW script on non-7SK loci

import csv, collections

csv.register_dialect("textdialect", delimiter='\t')

inputName = "/home/raflynn/7SK_ChIRPseq/WT2/combined_peaks.bed"
outputName = "/home/raflynn/7SK_ChIRPseq/WT2/combined_peaks_no7SK.bed"

# get 7SK loci
# ------------

# get 7SK ENSTs
ifile = open("/home/raflynn/7SK_ChIRPseq/genes/mm9_ensembl_genenames.txt", 'r')
reader = csv.reader(ifile, 'textdialect')

enst = set()
for row in reader:
	if row[1] == '7SK': enst.add(row[0])

ifile.close()

# get gene body for 7SK genes
ifile = open("/home/raflynn/7SK_ChIRPseq/genes/mm9_ensembl_coll.txt", 'r')
reader = csv.reader(ifile, 'textdialect')

genes = collections.defaultdict(lambda: []) # by chromosome
reader.next()

for row in reader:
	transcript = row[-1].split(';')[0]
	if transcript in enst: genes[row[1]].append(row)

ifile.close()

# getting rid of peaks
# --------------------

def in_7SK(chrom, start, end):
	listOfLoci = genes[chrom]
	for loc in listOfLoci:
		strand = loc[2]
		if strand == '+': locStart, locEnd = int(loc[3]), int(loc[4])
		else: locStart, locEnd = int(loc[4]), int(loc[3])
		
		if start > locEnd or end < locStart: continue #completely out
		if start < locStart and end - locStart >= 20: return True #left overlap
		if end > locEnd and locEnd - start >= 20: return True #right overlap
		if end > locEnd and locEnd - locStart >= 20: return True #complete overlap
		if end < locEnd and end - start >= 20: return True #completely within	
		
	return False
	
ifile = open(inputName, 'r')
reader = csv.reader(ifile, 'textdialect')

ofile = open(outputName, 'w')
writer = csv.writer(ofile, 'textdialect')

for row in reader:
	chrom = row[0]
	start = int(row[1])
	end = int(row[2])
	
	if not in_7SK(chrom, start, end): writer.writerow(row)
	
ifile.close()
ofile.close()
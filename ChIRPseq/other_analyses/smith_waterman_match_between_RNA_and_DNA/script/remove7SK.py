# remove7SK.py
# 3/28/14
# removes 7SK-overlapping regions from RNA-DNA.best

import csv, collections
csv.register_dialect("textdialect", delimiter='\t')

inputName = "/home/raflynn/7SK_ChIRPseq/smith_waterman_match_between_RNA_and_DNA/Mm_7SK_WT2_best_match/RNA-DNA.best.with_7SK"
outputName = "/home/raflynn/7SK_ChIRPseq/smith_waterman_match_between_RNA_and_DNA/Mm_7SK_WT2_best_match/RNA-DNA.best"

# get 7SK loci
# ------------

# get gene body for 7SK genes
ifile = open("/home/raflynn/7SK_ChIRPseq/smith_waterman_match_between_RNA_and_DNA/sequences/7SK_loci.txt", 'r')
reader = csv.reader(ifile, 'textdialect')

genes = collections.defaultdict(lambda: []) # by chromosome

reader.next()
reader.next()
for row in reader:
	genes[row[5]].append(row)

ifile.close()


# getting rid of peaks
# --------------------

ovl = 10

def in_7SK(chrom, start, end):
	listOfLoci = genes[chrom]
	for loc in listOfLoci:
		locStart, locEnd = int(loc[6]), int(loc[7])
		
		if start > locEnd or end < locStart: continue #completely out
		if start < locStart and end - locStart >= ovl: return True #left overlap
		if end >= locEnd and locEnd - start >= ovl: return True #right overlap
		if end >= locEnd and locEnd - locStart >= ovl: return True #complete overlap
		if end <= locEnd and end - start >= ovl: return True #completely within	
		
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
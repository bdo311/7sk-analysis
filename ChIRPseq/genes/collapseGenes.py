# collapseGenes.py
# 2/24/14
# Takes mm9_refseq.txt and mm9_ensembl.txt and collapses them so that there is only one TSS and TES pair per gene. Uses the TSS and the TES that will give the longest isoform. Will simply label as 'start' and 'end', avoiding the complexity involved with having the TES < TSS. 

import csv
import collections
import Queue

inputFile = "hg19_refseq.bed"
outputFile = "hg19_refseq_coll.txt"

csv.register_dialect("textdialect", delimiter='\t')

# Reading in file, choosing longest isoform for each gene, and adding it to a hash table
# -----------

# each chromosome has a dictionary with name of gene mapping to [tss, tes, strand, names] [0,0,'','']
chrom = collections.defaultdict(lambda: {})

ifile = open(inputFile, 'r')
reader = csv.reader(ifile, 'textdialect')

reader.next()
for row in reader:
	hashMap = chrom[row[1]] #row[1] is chrN
	geneName = row[5]
	strand = row[2]

	tss, tes = int(row[3]), int(row[4])
	
	name = row[0]

	if geneName not in hashMap: # first entry for a gene
		hashMap[geneName] = [tss, tes, strand, name + ';']
	else: # update tss and tes as needed; append name
		info = hashMap[geneName]
		if info[0] > tss: info[0] = tss
		if info[1] < tes: info[1] = tes
		info[3] = info[3] + name + ';'
		hashMap[geneName] = info

ifile.close()

# Sorting the genes by TSS within each chromosome. Uses a priority queue within each chromosome
# -------------

chromSorted = collections.defaultdict(lambda: Queue.PriorityQueue()) #lowest numbers come out first

for chrNum in chrom:
	for gene in chrom[chrNum]:
		tss = chrom[chrNum][gene]
		chromSorted[chrNum].put((tss, gene))

# Printing sorted, collapsed list into a file
# ----------

ofile = open(outputFile, 'w')
writer = csv.writer(ofile, 'textdialect')

header = ["Gene name", "Chr", "Strand", "Start", "End", "Names"]
writer.writerow(header)

for chrNum in chromSorted:
	while not chromSorted[chrNum].empty():
		gene = chromSorted[chrNum].get()[1]
		info = chrom[chrNum][gene]  #[tss, tes, strand, name + ';']

		outputRow = [gene, chrNum, info[2], info[0], info[1], info[3]]
		writer.writerow(outputRow)

ofile.close()


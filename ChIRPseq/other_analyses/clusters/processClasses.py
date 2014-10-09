# processClasses.py
# 5/9/14
# assigns genes to 4 classes, then puts in U1 data
# T = termination defect
# D = divergent tx
# B = both of the above
# N = none

import csv
csv.register_dialect("textdialect", delimiter='\t')

# get a list of all annotations for genes
ifile = open("annotated_genes.txt", 'r')
reader = csv.reader(ifile, 'textdialect')

geneToClass = {}

reader.next()
for row in reader:
	isTerm = (float(row[2]) != 0)
	divTx = (float(row[3]) > 0)
	
	if isTerm and divTx: name = 'B'
	elif isTerm: name = 'T'
	elif divTx: name = 'D'
	else: name = 'N'
	
	geneToClass[row[0]] = name
	
ifile.close()

# get upstream stats
ifile = open("/home/raflynn/ChIRPseq/genes/foundU1_upstream_new.txt", 'r')
reader = csv.reader(ifile, 'textdialect')

geneToUpstream = {}
geneToUpU1 = {}
geneToUpPAS = {}
for row in reader:
	name = row[0]
	length = int(row[1])
	numU1 = int(row[2])
	rateU1 = float(row[3])
	numPAS = int(row[4])
	ratePAS = float(row[5])
	if numPAS == 0: ratio = numU1
	else: ratio = float(numU1)/numPAS
	
	geneToUpstream[row[0]] = [numU1, numPAS, ratio]
	# positions of U1 and PAS
	u1_str = row[6].split(',')
	pas_str = row[7].split(',')
	u1 = []
	for x in u1_str:
		if x != '': u1.append(float(x)/length)
	pas = []
	for x in pas_str:
		if x != '': pas.append(float(x)/length)
		
	geneToUpU1[name] = u1
	geneToUpPAS[name] = pas
ifile.close()

# get gene stats
ifile = open("/home/raflynn/ChIRPseq/genes/foundU1_gene_new.txt", 'r')
reader = csv.reader(ifile, 'textdialect')

geneToGene = {}
geneToU1 = {}
geneToPAS = {}
for row in reader:
	name = row[0]
	length = int(row[1])
	numU1 = int(row[2])
	rateU1 = float(row[3])
	numPAS = int(row[4])
	ratePAS = float(row[5])

	if numPAS == 0: ratio = numU1
	else: ratio = float(numU1)/numPAS
	
	# gene stats
	geneToGene[name] = [numU1, rateU1, numPAS, ratePAS, ratio, length]
	
	# positions of U1 and PAS
	u1_str = row[6].split(',')
	pas_str = row[7].split(',')
	u1 = []
	for x in u1_str:
		if x != '': u1.append(float(x)/length)
	pas = []
	for x in pas_str:
		if x != '': pas.append(float(x)/length)
		
	geneToU1[name] = u1
	geneToPAS[name] = pas
	
ifile.close()

# end file
endFile = open("foundU1_gene_ends_50.txt", 'r')
endReader = csv.reader(endFile, 'textdialect')

endHeader = endReader.next()[2:]
geneToEnds = {}
for row in endReader:
	geneToEnds[row[0]] = row[2:]
	
endFile.close()

# write out files
# u1file = open("gene_up_u1.txt", 'w')
# u1writer = csv.writer(u1file, 'textdialect')
# pasfile = open("gene_up_pas.txt", 'w')
# paswriter = csv.writer(pasfile, 'textdialect')

ofile = open("gene_classes_new_50.txt", 'w')
writer = csv.writer(ofile, 'textdialect')

header = ["gene", "class", "up_u1", "up_pas", "up_ratio", "gene_u1", "gene_u1rate", "gene_pas", "gene_pasrate", "gene_ratio", "length"]
header.extend(endHeader)
writer.writerow(header)
for gene in geneToClass:
	if gene not in geneToGene or gene not in geneToUpstream: continue
	outputRow = [gene, geneToClass[gene]]
	outputRow.extend(geneToUpstream[gene])
	outputRow.extend(geneToGene[gene])
	outputRow.extend(geneToEnds[gene])
	writer.writerow(outputRow)
	
	# u1Row = [gene]
	# u1Row.extend(geneToUpU1[gene])
	# u1writer.writerow(u1Row)
	
	# pasRow = [gene]
	# pasRow.extend(geneToUpPAS[gene])
	# paswriter.writerow(pasRow)
	
ofile.close()
# u1file.close()
# pasfile.close()
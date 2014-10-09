# File: processHiC.py
# 4/26/14
# Converts NM to gene names, and compares contact freq with 7SK WT2 TSS, TES, and gene body ChIRP densities

import csv, collections, numpy
csv.register_dialect("textdialect", delimiter='\t')

inputFile = "7SK_wt_peak_summit.NoTargetGenes.hic.txt" #HiC input file: [NM, freq]

# 1. Read in nm to gene conversion
ifile = open("/home/raflynn/7SK_ChIRPseq/genes/mm9_refseq.txt", 'r')
reader = csv.reader(ifile, 'textdialect')

nmToGene = {}

reader.next()
for row in reader:
	nmToGene[row[0]] = row[-1]
	
ifile.close()

# 2. Match NMs to genes and put contact freqs in a map
ifile = open(inputFile, 'r')
reader = csv.reader(ifile, 'textdialect')

geneToValues = collections.defaultdict(lambda: []) #sometimes multiple NMs per gene

for row in reader:
	# lol edge  cases
	if row == []: continue
	if row[1] == '': continue
	
	# legit edge case
	if row[0] not in nmToGene: continue
	geneToValues[nmToGene[row[0]]].append(float(row[1]))
	
ifile.close()

geneToFreq = {} #average the freqs from all the NMs comprising each gene
for gene in geneToValues:
	geneToFreq[gene] = numpy.mean(geneToValues[gene])
	
# 3. Get length and ChIRP values for each gene, for TSS, TES, and gene body
wtDir = "/home/raflynn/7SK_ChIRPseq/WT2/bins/"

# gets map of gene to [total] except for geneBody, which gets [length of row, length of gene, total]
def getInfo(region):
	fn = wtDir + region + '/allchr.txt'
	ifile = open(fn, 'r')
	reader = csv.reader(ifile, 'textdialect')
	
	geneToInfo = collections.defaultdict(lambda: []) #gene to [length, total]
	for row in reader:
		name = row[0]
		if region=="geneBody": 
			geneToInfo[name].append(int(row[4]) - int(row[3])) #length of gene
			geneToInfo[name].append(int(len(row) - 7)) #number of bins
		geneToInfo[name].append(float(row[6]))
		
	ifile.close()
	return geneToInfo
	
tssMap = getInfo('tss')
geneBodyMap = getInfo('geneBody')
tesMap = getInfo('tes')

# 4. Integrate everything into a file
ofile = open("hiC_with_chirp.txt", 'w')
writer = csv.writer(ofile, 'textdialect')

writer.writerow(["Gene name", "Length", "HiC frequency", "TSS average", "Gene body average", "TES average"])
for gene in geneToFreq:
	if gene not in tssMap: continue
	outputRow = [gene, geneBodyMap[gene][0], geneToFreq[gene]]
	outputRow.extend([tssMap[gene][0]/float(400), geneBodyMap[gene][2]/float(geneBodyMap[gene][1]), tesMap[gene][0]/float(400)])
	writer.writerow(outputRow)
	
ofile.close()
	
		
	

# classifyGenicPeaks.py 
# 7/28/14
# takes in a peak file from macs2 and classifies them by how many are genic (within 1kb of the gene TSS or TSS) and how many are intergenic. Uses peak center for the calculation

import csv, collections, sys, Queue
csv.register_dialect("textdialect", delimiter='\t')

def getGeneStartEnds(fn):
	chrToStartEnds = collections.defaultdict(lambda: [])

	ifile = open(fn, 'r')
	reader = csv.reader(ifile, 'textdialect')

	reader.next()
	for row in reader:
		name = row[0]
		start = int(row[3])
		end = int(row[6])
		chrToStartEnds[row[1]].append([name, start, end])

	return chrToStartEnds

def getInGene(line, chrToStartEnds):
	chrom = line[0]
	peakCtr = (int(line[1]) + int(line[2]))/2
	for gene in chrToStartEnds[chrom]:
		if peakCtr >= gene[1] and peakCtr <= gene[2]: return gene
	return ['Intergenic']

def getGenicClassif(ifn, ofn, chrToStartEnds):
	ifile = open(ifn, 'r')
	reader = csv.reader(ifile, 'textdialect')

	ofile = open(ofn, 'w')
	writer = csv.writer(ofile, 'textdialect')

	numIntergenic = 0
	numGenic = 0
	for row in reader:
		genicStatus = getInGene(row, chrToStartEnds)
		row.extend(genicStatus)
		writer.writerow(row)
		if genicStatus == ['Intergenic']: numIntergenic += 1
		else: numGenic += 1

	ifile.close()
	ofile.close()
	return numGenic, numIntergenic



def main():
	genome, geneAnnot = sys.argv[2], sys.argv[3]
	chrToStartEnds = getGeneStartEnds("/home/raflynn/7SK/ChIRPseq/genes/" + genome + "_" + geneAnnot + "_startend.txt")
	numGenic, numIntergenic = getGenicClassif(sys.argv[1], sys.argv[4], chrToStartEnds)

	print "Genic:", numGenic
	print "Intergenic:", numIntergenic

if __name__ == '__main__':
	main()

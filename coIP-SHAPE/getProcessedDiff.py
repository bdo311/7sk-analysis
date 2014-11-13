# getProcessedDiff.py
# 10/24/14
# compares ARID IP/Input to IgG/Input

import csv, collections, numpy
csv.register_dialect("textdialect", delimiter='\t')

def process(fn):
	ifile = open(fn, 'r')
	reader = csv.reader(ifile, 'textdialect')
	
	treat = collections.defaultdict(lambda: 0)
	input = collections.defaultdict(lambda: 0)
	
	reader.next()
	for row in reader:
		treat[row[2]] += float(row[7])
		input[row[2]] += float(row[8])
		
	ifile.close()
	return treat, input


def main():
	arid = "/arrayAhome/raflynn/7SK/coIP-SHAPE/aridVsInput_allreads/gene_exp.diff"
	igg = "/arrayAhome/raflynn/7SK/coIP-SHAPE/iggVsInput_allreads/gene_exp.diff"
	
	aridA, aridInput = process(arid)
	iggI, iggInput = process(igg)
	
	ofile = open("arid_vs_igg.txt", 'w')
	writer = csv.writer(ofile, 'textdialect')
	
	writer.writerow(['gene', 'arid_arid', 'arid_input', 'igg_igg', 'igg_input', 'lfc'])
	for gene in aridA:
		if aridA[gene] < 0.1 or aridInput[gene] < 0.1 or iggI[gene] < 0.1 or iggInput[gene] < 0.1: continue
		writer.writerow([gene, aridA[gene], aridInput[gene], iggI[gene], iggInput[gene], numpy.log2((aridA[gene]/aridInput[gene])/(iggI[gene]/iggInput[gene]))])
		
	ofile.close()

if __name__ == "__main__":
	main()
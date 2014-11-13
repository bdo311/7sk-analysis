# chirpToExpression.py 
# 7/8/14
# looks for a correlation between 7SK TSS read density and expression level

import csv, collections, numpy
csv.register_dialect("textdialect", delimiter='\t')

def getRefToRNASeq(fn):
	refToRNASeq = {}
	
	ifile = open(fn, 'r')
	reader = csv.reader(ifile, 'textdialect')

	for row in reader:
		if 'refseq' not in row[0]: continue
		refToRNASeq[row[1]] = float(row[2])

	ifile.close()
	return refToRNASeq

def getRefTo7SK(fn, start, stop):
	refToChIRP = {}

	ifile = open(fn, 'r')
	reader = csv.reader(ifile, 'textdialect')

	for row in reader:
		names = row[5].split(';')
		densityBins = [float(x) for x in row[(7 + start):(7 + stop)]]
		density = numpy.sum(densityBins)
		for name in names:
			refToChIRP[name] = density

	ifile.close()
	return refToChIRP	

def main():
	# refseq to rnaseq expression
	refToRNASeq = getRefToRNASeq("/home/raflynn/ChIRPseq/other_analyses/expression/GZ_v65_RNAseq_Expression.txt")

	# chirp genes to 7SK
	refToChIRP = getRefTo7SK("/home/raflynn/ChIRPseq/WT2_new/bins/tss/allchr_sorted.txt",0,400)

	# correlate
	ofile = open("rnaseq_to_7sk.txt", 'w')
	writer = csv.writer(ofile, 'textdialect')
	writer.writerow(['gene', 'rnaseq', 'chirp', 'log_rnaseq', 'log_chirp'])
	for gene in refToRNASeq:
		if gene not in refToChIRP: continue
		expr = refToRNASeq[gene]
		chirp = refToChIRP[gene]

		if expr < 0.5: continue
		outputRow = [gene, expr, chirp, numpy.log2(expr + 1), numpy.log2(chirp + 1)]
		writer.writerow(outputRow)

	ofile.close()



if __name__ == '__main__':
	main()

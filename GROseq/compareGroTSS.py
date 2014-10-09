# compareGroTSS.py
# 9/11/14
# compare sense and divergent GRO reads across different ASO treatments

import os, glob, csv, collections, numpy
from multiprocessing import Pool
csv.register_dialect("textdialect", delimiter='\t')

def rowMean(row):
	return numpy.mean([float(x) for x in row])
	
def processFolder(folder):
	treatment = '_'.join(folder.split('_')[:2])
	fn = "/home/raflynn/GROseq_ASO/" + folder + "/bins/tss/allchr.txt"
	ifile = open(fn, 'r')
	reader = csv.reader(ifile, 'textdialect')
	
	sense = {}
	div = {}
	
	for row in reader:
		geneName = row[0]
		strand = 'plus' if row[2] == '+' else 'minus'
		if strand in folder: #if strand of gene matches strand of bedgraph
			sense[geneName] = rowMean(row[206:306])
		else:
			div[geneName] = rowMean(row[106:206])
			
	return treatment, sense, div

def main():
	folders = glob.glob("GRO_*")
	sense = collections.defaultdict(lambda: {})
	div = collections.defaultdict(lambda: {})
	for folder in folders:
		print folder
		treatment, senseMap, divMap = processFolder(folder)
		sense[treatment].update(senseMap)
		div[treatment].update(divMap)
	
	# write output file
	ofile = open("tss_compared.txt", 'w')
	writer = csv.writer(ofile, 'textdialect')
	
	treatments = sense.keys()
	genes = sense['GRO_12C'].keys()
	
	header = ['']
	header.extend(['sense' + x for x in treatments])
	header.extend(['div' + x for x in treatments])
	writer.writerow(header)
	
	for gene in genes:
		outputRow = [gene]
		for treatment in sense:
			outputRow.append(sense[treatment][gene])
		for treatment in div:
			outputRow.append(div[treatment][gene])
			
		writer.writerow(outputRow)
		
	ofile.close()
	
	
	

if __name__ == "__main__":
	main()
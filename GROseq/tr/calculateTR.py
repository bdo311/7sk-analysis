# calculateTR.py
# 10/18/14

import csv, numpy
import matplotlib as mpl
from scipy.stats import ks_2samp
csv.register_dialect("textdialect", delimiter='\t')

def getHashMap(fn, strand):
	ifile = open(fn, 'r')
	reader = csv.reader(ifile, 'textdialect')
	
	geneToNum = {}
	for row in reader:
		if row[2] != strand: continue #wrong strand
		geneName = row[0] #might be duplicates but ok for now
		#length = int(row[4]) - int(row[3])
		length = 1 #average only
		geneToNum[geneName] = length * float(row[6]) #length * reads/nt average
		
	return geneToNum		
		
def main():

	# calculating traveling ratio
	treatments = ['GRO_' + x + 'comb' for x in ['125', '123', '12C', '65', '63', '6C']]
	tmToTr = {}
	for tm in treatments:
		print tm
		prom_plus = getHashMap(tm + '_plus/bins/promoter/allchr.txt', '+')
		prom_minus = getHashMap(tm + '_minus/bins/promoter/allchr.txt', '-')
		ds_plus = getHashMap(tm + '_plus/bins/downstream/allchr.txt', '+')
		ds_minus = getHashMap(tm + '_minus/bins/downstream/allchr.txt', '-')
		
		tr = []
		for gene in prom_plus:
			if gene not in ds_plus or ds_plus[gene] == 0: continue
			if prom_plus[gene] == 0: continue
			tr.append(prom_plus[gene]/ds_plus[gene])
		for gene in prom_minus:
			if gene not in ds_minus or ds_minus[gene] == 0: continue
			if prom_minus[gene] == 0: continue
			tr.append(prom_minus[gene]/ds_minus[gene])			
		
		tmToTr[tm] = tr
			
	# outputting things
	print "Outputting file"
	ofile = open("tr_comb.txt", 'w')
	writer = csv.writer(ofile, 'textdialect')
	
	for tm in tmToTr:
		outputRow = [tm]
		outputRow.extend(tmToTr[tm])
		writer.writerow(outputRow)
		
	ofile.close()
	
	# making boxplot
	mpl.use('agg')
	import matplotlib.pyplot as plt
	data_to_plot = [[numpy.log10(x + 1) for x in tmToTr[tm]] for tm in treatments]
	fig = plt.figure(1, figsize=(9, 6))
	ax = fig.add_subplot(111)
	ax.set_xticklabels(treatments)
	bp = ax.boxplot(data_to_plot)
	fig.savefig('tr_comb_boxplot.png', bbox_inches='tight')
	
	# ks test
	print "125 vs 12C: ", ks_2samp(tmToTr["GRO_125comb"], tmToTr["GRO_12Ccomb"])[1]
	print "123 vs 12C: ", ks_2samp(tmToTr["GRO_123comb"], tmToTr["GRO_12Ccomb"])[1]
	print "65 vs 6C: ", ks_2samp(tmToTr["GRO_65comb"], tmToTr["GRO_6Ccomb"])[1]
	print "63 vs 6C: ", ks_2samp(tmToTr["GRO_63comb"], tmToTr["GRO_6Ccomb"])[1]
	

if __name__ == '__main__':
	main()
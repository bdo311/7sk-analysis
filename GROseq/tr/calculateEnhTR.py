# calculateEnhTR.py
# 10/29/14

import csv, os
#import numpy
# import matplotlib as mpl
# from scipy.stats import ks_2samp
csv.register_dialect("textdialect", delimiter='\t')

def getHashMap(fn):
	os.chdir("/home/raflynn/7SK/GROseq/")
	ifile = open(fn, 'r')
	reader = csv.reader(ifile, 'textdialect')
	
	geneToNum = {}
	counter = 0
	for row in reader:
		counter = counter + 1
		geneName = str(counter) + row[3] #might be duplicates but ok for now
		#length = int(row[4]) - int(row[3])
		length = 1 #average only
		geneToNum[geneName] = length * float(row[4]) #length * reads/nt average
		
	return geneToNum		
		
def main():

	# calculating traveling ratio
	treatments = ['GRO_' + x + 'comb' for x in ['125', '123', '12C', '65', '63', '6C']]
	tmToTr = {}
	for tm in treatments: #strandedness doesn't matter since txn goes in both directions for reg enhs
		print tm
		ctr_plus = getHashMap(tm + '_plus/bins/regenh_ctr1000/allchr.txt')
		ctr_minus = getHashMap(tm + '_minus/bins/regenh_ctr1000/allchr.txt')
		left_plus = getHashMap(tm + '_plus/bins/regenh_left10000/allchr.txt')
		left_minus = getHashMap(tm + '_minus/bins/regenh_left10000/allchr.txt')
		right_plus = getHashMap(tm + '_plus/bins/regenh_right10000/allchr.txt')
		right_minus = getHashMap(tm + '_minus/bins/regenh_right10000/allchr.txt')
		
		tr = []
		for gene in ctr_plus:
			prox = ctr_plus[gene] + ctr_minus[gene]
			dist = 0.5*(left_plus[gene] + left_minus[gene] + right_plus[gene] + right_minus[gene])
			if prox == 0 or dist == 0: continue
			tr.append(float(prox)/dist)			
		
		tmToTr[tm] = tr
			
	# outputting things
	os.chdir("tr")
	print "Outputting file"
	ofile = open("enh_tr_comb_2000.txt", 'w')
	writer = csv.writer(ofile, 'textdialect')
	
	for tm in tmToTr:
		outputRow = [tm]
		outputRow.extend(tmToTr[tm])
		writer.writerow(outputRow)
		
	ofile.close()
	
	# # making boxplot
	# mpl.use('agg')
	# import matplotlib.pyplot as plt
	# data_to_plot = [[numpy.log2(x + 1) for x in tmToTr[tm]] for tm in treatments]
	# fig = plt.figure(1, figsize=(9, 6))
	# ax = fig.add_subplot(111)
	# ax.set_xticklabels(treatments)
	# bp = ax.boxplot(data_to_plot)
	# fig.savefig('enh_tr_comb_boxplot.png', bbox_inches='tight')
	
	# # ks test
	# print "125 vs 12C: ", ks_2samp(tmToTr["GRO_125comb"], tmToTr["GRO_12Ccomb"])[1]
	# print "123 vs 12C: ", ks_2samp(tmToTr["GRO_123comb"], tmToTr["GRO_12Ccomb"])[1]
	# print "65 vs 6C: ", ks_2samp(tmToTr["GRO_65comb"], tmToTr["GRO_6Ccomb"])[1]
	# print "63 vs 6C: ", ks_2samp(tmToTr["GRO_63comb"], tmToTr["GRO_6Ccomb"])[1]
	

if __name__ == '__main__':
	main()
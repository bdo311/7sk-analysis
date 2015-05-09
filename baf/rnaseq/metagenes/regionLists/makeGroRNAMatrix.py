# makeGroRNAMatrix.py
# 5/7/15

import csv, collections
import numpy as np
csv.register_dialect("textdialect",delimiter='\t')

rnaFolders = ["SRR1265785_rmdup","SRR1265786_rmdup","SRR1265787_rmdup","SRR1265788_rmdup"]
groFolders = ["GRO_12Ccomb_sense","GRO_12Ccomb_antisense","GRO_125comb_sense","GRO_125comb_antisense"]

re = collections.defaultdict(lambda: [])
se = collections.defaultdict(lambda: [])

for rna in rnaFolders:
	ifile = open("/home/raflynn/7SK/baf/rnaseq/metagenes/"+rna+"/bins/RY_enh_centered/allchr_sorted.txt",'r')
	reader = csv.reader(ifile, 'textdialect')
	for row in reader:	
		#re[row[3]].append(float(row[6]))
		re[row[3]].append(np.sum([float(x) for x in row[107:307]]))
	ifile.close()
	
	ifile = open("/home/raflynn/7SK/baf/rnaseq/metagenes/"+rna+"/bins/SE_indiv/allchr_sorted.txt",'r')
	reader = csv.reader(ifile, 'textdialect')
	for row in reader:	
		se[row[3]].append(np.sum([float(x) for x in row[107:307]]))
	ifile.close()
	
for gro in groFolders:
	ifile = open("/arrayAhome/raflynn/7SK/GROseq/Flynn/metagenes/"+gro+"/bins/RY_enh_centered/allchr_sorted.txt",'r')
	reader = csv.reader(ifile, 'textdialect')
	for row in reader:	
		re[row[3]].append(np.sum([float(x) for x in row[107:307]]))
	ifile.close()
	
	ifile = open("/arrayAhome/raflynn/7SK/GROseq/Flynn/metagenes/"+gro+"/bins/SE_indiv/allchr_sorted.txt",'r')
	reader = csv.reader(ifile, 'textdialect')
	for row in reader:	
		se[row[3]].append(np.sum([float(x) for x in row[107:307]]))
	ifile.close()

ofile = open("re_matrix_500.txt",'w')
writer = csv.writer(ofile, 'textdialect')
header = ['name']
header.extend(rnaFolders)
header.extend(groFolders)
writer.writerow(header)
for r in re:
	outputRow = [r]
	outputRow.extend(re[r])
	writer.writerow(outputRow)
ofile.close()

ofile = open("se_matrix_500.txt",'w')
writer = csv.writer(ofile, 'textdialect')
header = ['name']
header.extend(rnaFolders)
header.extend(groFolders)
writer.writerow(header)
for s in se:
	outputRow = [s]
	outputRow.extend(se[s])
	writer.writerow(outputRow)
ofile.close()
# makeCDFMatrix.py
# 3/4/15
# makes matrix of average read density at all enhancers and super enhancers from -500 to +500bp (middle 200 bins)

import csv, collections
import numpy as np
csv.register_dialect("textdialect",delimiter='\t')

reToValues = collections.defaultdict(lambda: [])
seToValues = collections.defaultdict(lambda: [])
tssToValues = collections.defaultdict(lambda: [])

datasets = ["BAF155_3p_merge","BAF155_5p_merge","BAF155_Scr_merge","Input_3p_merge","Input_5p_merge","Input_Scr_merge","Pou5f1_3p_mergeAB","Pou5f1_3p_mergeAD","Pou5f1_5p_mergeAB","Pou5f1_5p_mergeAD","Pou5f1_Scr_mergeAB","Pou5f1_Scr_mergeAD"]
for d in datasets:
	ifile = open("/home/raflynn/7SK/ChIPseq/Flynn/Second_run/metagenes/" + d + "/bins/All_Enhancers/allchr_sorted.txt",'r')
	reader = csv.reader(ifile, 'textdialect')
	
	for row in reader:
		name = '_'.join(row[:3])
		isSE = ('__' in row[3])
		
		values = [float(x) for x in row[107:307]]
		avg = np.mean(values)
		if isSE: seToValues[name].append(avg)
		else: reToValues[name].append(avg)
	
	ifile.close()
	
	ifile = open("/home/raflynn/7SK/ChIPseq/Flynn/Second_run/metagenes/" + d + "/bins/TSS_trans_paused/allchr_sorted.txt",'r')
	reader = csv.reader(ifile, 'textdialect')
	
	for row in reader:
		name = '_'.join(row[:3])
		
		values = [float(x) for x in row[107:307]]
		avg = np.mean(values)
		tssToValues[name].append(avg)
	
	ifile.close()
	
ofile = open("re_matrix.txt",'w')
writer = csv.writer(ofile, 'textdialect')
header = ['']
header.extend(datasets)
writer.writerow(header)
for re in reToValues:
	outputrow = [re]
	outputrow.extend(reToValues[re])
	writer.writerow(outputrow)
ofile.close()

ofile = open("se_matrix.txt",'w')
writer = csv.writer(ofile, 'textdialect')
header = ['']
header.extend(datasets)
writer.writerow(header)
for se in seToValues:
	outputrow = [se]
	outputrow.extend(seToValues[se])
	writer.writerow(outputrow)
ofile.close()	
	
ofile = open("tss_matrix.txt",'w')
writer = csv.writer(ofile, 'textdialect')
header = ['']
header.extend(datasets)
writer.writerow(header)
for tss in tssToValues:
	outputrow = [tss]
	outputrow.extend(tssToValues[tss])
	writer.writerow(outputrow)
ofile.close()	


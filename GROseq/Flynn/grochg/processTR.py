# processTR.py
# 3/1/15

import csv, collections
csv.register_dialect("textdialect",delimiter='\t')

samples = ["GRO_123","GRO_125","GRO_12C","GRO_63","GRO_65","GRO_6C"]
# promoter TR
sampleToTR = {}
for sample in samples:
	geneToInfo = collections.defaultdict(lambda: [0,0])
	promFile = open("/home/raflynn/7SK/GROseq/Flynn/metagenes/" + sample + "comb_sense/bins/TR_prom_trans_paused/allchr_sorted.txt",'r')
	reader = csv.reader(promFile, 'textdialect')
	for row in reader:
		name = row[3]
		geneToInfo[name][0]=float(row[6])
	promFile.close()
	
	dsFile = open("/home/raflynn/7SK/GROseq/Flynn/metagenes/" + sample + "comb_sense/bins/TR_ds_trans_paused/allchr_sorted.txt",'r')
	reader = csv.reader(dsFile, 'textdialect')
	for row in reader:
		name = row[3]
		geneToInfo[name][1]=float(row[6])
	dsFile.close()
	
	trArray = []
	for name in geneToInfo:
		try: tr = geneToInfo[name][0]/geneToInfo[name][1]
		except: tr = 1
		trArray.append(tr)
	sampleToTR[sample] = trArray
	
ofile = open("promoter_tr.txt",'w')
writer = csv.writer(ofile, 'textdialect')
for s in sampleToTR:
	outputrow = [s]
	outputrow.extend(sampleToTR[s])
	writer.writerow(outputrow)
ofile.close()

# RE TR
sampleToTR = {}
for sample in samples:
	geneToInfo = collections.defaultdict(lambda: [0,0,0, 0]) #ctr sense, ctr antisense, left antisense, right sense
	ctrFile = open("/home/raflynn/7SK/GROseq/Flynn/metagenes/" + sample + "comb_sense/bins/RY_enh_TR_center250/allchr_sorted.txt",'r')
	reader = csv.reader(ctrFile, 'textdialect')
	for row in reader:
		name = row[3]
		geneToInfo[name][0]=float(row[6])
	ctrFile.close()
	ctrFile = open("/home/raflynn/7SK/GROseq/Flynn/metagenes/" + sample + "comb_antisense/bins/RY_enh_TR_center250/allchr_sorted.txt",'r')
	reader = csv.reader(ctrFile, 'textdialect')
	for row in reader:
		name = row[3]
		geneToInfo[name][1]=float(row[6])
	ctrFile.close()	
	
	leftFile = open("/home/raflynn/7SK/GROseq/Flynn/metagenes/" + sample + "comb_antisense/bins/RY_enh_TR_left250/allchr_sorted.txt",'r')
	reader = csv.reader(leftFile, 'textdialect')
	for row in reader:
		name = row[3]
		geneToInfo[name][2]=float(row[6])
	leftFile.close()	
	rightFile = open("/home/raflynn/7SK/GROseq/Flynn/metagenes/" + sample + "comb_antisense/bins/RY_enh_TR_right250/allchr_sorted.txt",'r')
	reader = csv.reader(rightFile, 'textdialect')
	for row in reader:
		name = row[3]
		geneToInfo[name][2]=float(row[6])
	rightFile.close()	
	
	trArray = []
	for name in geneToInfo:
		try: 
			ctrAvg = (geneToInfo[name][0] + geneToInfo[name][1])/2
			sideAvg = (geneToInfo[name][2] + geneToInfo[name][3])/2
			tr = ctrAvg/sideAvg
		except: tr = 1
		trArray.append(tr)
	sampleToTR[sample] = trArray
	
ofile = open("re_tr250.txt",'w')
writer = csv.writer(ofile, 'textdialect')
for s in sampleToTR:
	outputrow = [s]
	outputrow.extend(sampleToTR[s])
	writer.writerow(outputrow)
ofile.close()

# SE TR
sampleToTR = {}
for sample in samples:
	geneToInfo = collections.defaultdict(lambda: [0,0,0, 0]) #ctr sense, ctr antisense, left antisense, right sense
	ctrFile = open("/home/raflynn/7SK/GROseq/Flynn/metagenes/" + sample + "comb_sense/bins/SE_indiv_TR_center250/allchr_sorted.txt",'r')
	reader = csv.reader(ctrFile, 'textdialect')
	for row in reader:
		name = row[3]
		geneToInfo[name][0]=float(row[6])
	ctrFile.close()
	ctrFile = open("/home/raflynn/7SK/GROseq/Flynn/metagenes/" + sample + "comb_antisense/bins/SE_indiv_TR_center250/allchr_sorted.txt",'r')
	reader = csv.reader(ctrFile, 'textdialect')
	for row in reader:
		name = row[3]
		geneToInfo[name][1]=float(row[6])
	ctrFile.close()	
	
	leftFile = open("/home/raflynn/7SK/GROseq/Flynn/metagenes/" + sample + "comb_antisense/bins/SE_indiv_TR_left250/allchr_sorted.txt",'r')
	reader = csv.reader(leftFile, 'textdialect')
	for row in reader:
		name = row[3]
		geneToInfo[name][2]=float(row[6])
	leftFile.close()	
	rightFile = open("/home/raflynn/7SK/GROseq/Flynn/metagenes/" + sample + "comb_antisense/bins/SE_indiv_TR_right250/allchr_sorted.txt",'r')
	reader = csv.reader(rightFile, 'textdialect')
	for row in reader:
		name = row[3]
		geneToInfo[name][2]=float(row[6])
	rightFile.close()	
	
	trArray = []
	for name in geneToInfo:
		try: 
			ctrAvg = (geneToInfo[name][0] + geneToInfo[name][1])/2
			sideAvg = (geneToInfo[name][2] + geneToInfo[name][3])/2
			tr = ctrAvg/sideAvg
		except: tr = 1
		trArray.append(tr)
	sampleToTR[sample] = trArray
	
ofile = open("se_tr250.txt",'w')
writer = csv.writer(ofile, 'textdialect')
for s in sampleToTR:
	outputrow = [s]
	outputrow.extend(sampleToTR[s])
	writer.writerow(outputrow)
ofile.close()
	
	

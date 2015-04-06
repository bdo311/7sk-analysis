# processTR.py
# 3/1/15

import csv, collections
csv.register_dialect("textdialect",delimiter='\t')

samples = ["GRO_123","GRO_125","GRO_12C","GRO_63","GRO_65","GRO_6C"]
print "promoter TR"
sampleToTR = {}
sampleToProx = {}
sampleToDist = {}
sampleToNames = []
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
	proxArray = []
	dsArray = []
	for name in geneToInfo:
		try: tr = geneToInfo[name][0]/geneToInfo[name][1]
		except: tr = 1
		trArray.append(tr)
		proxArray.append(geneToInfo[name][0])
		dsArray.append(geneToInfo[name][1])
		if name not in sampleToNames: sampleToNames.append(name)
		
	sampleToTR[sample] = trArray
	sampleToProx[sample] = proxArray
	sampleToDist[sample] = dsArray
	
ofile = open("promoter_tr.txt",'w')
writer = csv.writer(ofile, 'textdialect')
for s in sampleToTR:
	outputrow = [s]
	outputrow.extend(sampleToTR[s])
	writer.writerow(outputrow)
ofile.close()

print "names for individual promoters"
tfile = open("promoter_tr_withnames.txt",'w')
twriter = csv.writer(tfile, 'textdialect')
pfile = open("promoter_prox_withnames.txt",'w')
pwriter = csv.writer(pfile, 'textdialect')
dfile = open("promoter_ds_withnames.txt",'w')
dwriter = csv.writer(dfile, 'textdialect')

samps = sampleToTR.keys()
header = ['name']
header.extend(samps)
twriter.writerow(header)
pwriter.writerow(header)
dwriter.writerow(header)

for i in range(len(sampleToNames)):
	tOutputRow = [sampleToNames[i]]
	pOutputRow = [sampleToNames[i]]
	dOutputRow = [sampleToNames[i]]
	for s in samps:
		tOutputRow.append(sampleToTR[s][i])
		pOutputRow.append(sampleToProx[s][i])
		dOutputRow.append(sampleToDist[s][i])
	twriter.writerow(tOutputRow)
	pwriter.writerow(pOutputRow)
	dwriter.writerow(dOutputRow)
tfile.close()
pfile.close()
dfile.close()

print "RE TR"
sampleToTR = {}
sampleToProx = {}
sampleToDist = {}
sampleToNames = []
for sample in samples:
	geneToInfo = collections.defaultdict(lambda: [0,0,0, 0]) #ctr sense, ctr antisense, left antisense, right sense
	ctrFile = open("/home/raflynn/7SK/GROseq/Flynn/metagenes/" + sample + "comb_sense/bins/RY_enh_TR_center/allchr_sorted.txt",'r')
	reader = csv.reader(ctrFile, 'textdialect')
	for row in reader:
		name = row[3]
		geneToInfo[name][0]=float(row[6])
	ctrFile.close()
	ctrFile = open("/home/raflynn/7SK/GROseq/Flynn/metagenes/" + sample + "comb_antisense/bins/RY_enh_TR_center/allchr_sorted.txt",'r')
	reader = csv.reader(ctrFile, 'textdialect')
	for row in reader:
		name = row[3]
		geneToInfo[name][1]=float(row[6])
	ctrFile.close()	
	
	leftFile = open("/home/raflynn/7SK/GROseq/Flynn/metagenes/" + sample + "comb_antisense/bins/RY_enh_TR_left/allchr_sorted.txt",'r')
	reader = csv.reader(leftFile, 'textdialect')
	for row in reader:
		name = row[3]
		geneToInfo[name][2]=float(row[6])
	leftFile.close()	
	rightFile = open("/home/raflynn/7SK/GROseq/Flynn/metagenes/" + sample + "comb_sense/bins/RY_enh_TR_right/allchr_sorted.txt",'r')
	reader = csv.reader(rightFile, 'textdialect')
	for row in reader:
		name = row[3]
		geneToInfo[name][2]=float(row[6])
	rightFile.close()	
	
	trArray = []
	proxArray = []
	dsArray = []
	for name in geneToInfo:
		try: 
			ctrAvg = (geneToInfo[name][0] + geneToInfo[name][1])
			sideAvg = (geneToInfo[name][2] + geneToInfo[name][3])/2
			tr = ctrAvg/sideAvg
		except: tr = 1
		trArray.append(tr)		
		proxArray.append(ctrAvg)
		dsArray.append(sideAvg)
		if name not in sampleToNames: sampleToNames.append(name)
	sampleToTR[sample] = trArray
	sampleToProx[sample] = proxArray
	sampleToDist[sample] = dsArray
	
ofile = open("re_tr.txt",'w')
writer = csv.writer(ofile, 'textdialect')
for s in sampleToTR:
	outputrow = [s]
	outputrow.extend(sampleToTR[s])
	writer.writerow(outputrow)
ofile.close()


print "names for individual RE"
tfile = open("re_tr_withnames.txt",'w')
twriter = csv.writer(tfile, 'textdialect')
pfile = open("re_prox_withnames.txt",'w')
pwriter = csv.writer(pfile, 'textdialect')
dfile = open("re_ds_withnames.txt",'w')
dwriter = csv.writer(dfile, 'textdialect')

samps = sampleToTR.keys()
header = ['name']
header.extend(samps)
twriter.writerow(header)
pwriter.writerow(header)
dwriter.writerow(header)

for i in range(len(sampleToNames)):
	tOutputRow = [sampleToNames[i]]
	pOutputRow = [sampleToNames[i]]
	dOutputRow = [sampleToNames[i]]
	for s in samps:
		tOutputRow.append(sampleToTR[s][i])
		pOutputRow.append(sampleToProx[s][i])
		dOutputRow.append(sampleToDist[s][i])
	twriter.writerow(tOutputRow)
	pwriter.writerow(pOutputRow)
	dwriter.writerow(dOutputRow)
tfile.close()
pfile.close()
dfile.close()

print "SE TR"
sampleToTR = {}
sampleToProx = {}
sampleToDist = {}
sampleToNames = []
for sample in samples:
	geneToInfo = collections.defaultdict(lambda: [0,0,0, 0]) #ctr sense, ctr antisense, left antisense, right sense
	ctrFile = open("/home/raflynn/7SK/GROseq/Flynn/metagenes/" + sample + "comb_sense/bins/SE_indiv_TR_center/allchr_sorted.txt",'r')
	reader = csv.reader(ctrFile, 'textdialect')
	for row in reader:
		name = row[3]
		geneToInfo[name][0]=float(row[6])
	ctrFile.close()
	ctrFile = open("/home/raflynn/7SK/GROseq/Flynn/metagenes/" + sample + "comb_antisense/bins/SE_indiv_TR_center/allchr_sorted.txt",'r')
	reader = csv.reader(ctrFile, 'textdialect')
	for row in reader:
		name = row[3]
		geneToInfo[name][1]=float(row[6])
	ctrFile.close()	
	
	leftFile = open("/home/raflynn/7SK/GROseq/Flynn/metagenes/" + sample + "comb_antisense/bins/SE_indiv_TR_left/allchr_sorted.txt",'r')
	reader = csv.reader(leftFile, 'textdialect')
	for row in reader:
		name = row[3]
		geneToInfo[name][2]=float(row[6])
	leftFile.close()	
	rightFile = open("/home/raflynn/7SK/GROseq/Flynn/metagenes/" + sample + "comb_sense/bins/SE_indiv_TR_right/allchr_sorted.txt",'r')
	reader = csv.reader(rightFile, 'textdialect')
	for row in reader:
		name = row[3]
		geneToInfo[name][2]=float(row[6])
	rightFile.close()	
	
	trArray = []
	proxArray = []
	dsArray = []
	for name in geneToInfo:
		try: 
			ctrAvg = (geneToInfo[name][0] + geneToInfo[name][1])
			sideAvg = (geneToInfo[name][2] + geneToInfo[name][3])/2
			tr = ctrAvg/sideAvg
		except: tr = 1
		trArray.append(tr)		
		proxArray.append(ctrAvg)
		dsArray.append(sideAvg)
		if name not in sampleToNames: sampleToNames.append(name)
	sampleToTR[sample] = trArray
	sampleToProx[sample] = proxArray
	sampleToDist[sample] = dsArray
	
ofile = open("se_tr.txt",'w')
writer = csv.writer(ofile, 'textdialect')
for s in sampleToTR:
	outputrow = [s]
	outputrow.extend(sampleToTR[s])
	writer.writerow(outputrow)
ofile.close()

print "names for individual SE"

tfile = open("se_tr_withnames.txt",'w')
twriter = csv.writer(tfile, 'textdialect')
pfile = open("se_prox_withnames.txt",'w')
pwriter = csv.writer(pfile, 'textdialect')
dfile = open("se_ds_withnames.txt",'w')
dwriter = csv.writer(dfile, 'textdialect')

samps = sampleToTR.keys()
header = ['name']
header.extend(samps)
twriter.writerow(header)
pwriter.writerow(header)
dwriter.writerow(header)

for i in range(len(sampleToNames)):
	tOutputRow = [sampleToNames[i]]
	pOutputRow = [sampleToNames[i]]
	dOutputRow = [sampleToNames[i]]
	for s in samps:
		tOutputRow.append(sampleToTR[s][i])
		pOutputRow.append(sampleToProx[s][i])
		dOutputRow.append(sampleToDist[s][i])
	twriter.writerow(tOutputRow)
	pwriter.writerow(pOutputRow)
	dwriter.writerow(dOutputRow)
tfile.close()
pfile.close()
dfile.close()
		

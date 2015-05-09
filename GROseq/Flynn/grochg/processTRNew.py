# processTRNew.py
# 5/5/15

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
	promFile = open("/home/raflynn/7SK/GROseq/Flynn/metagenes/" + sample + "comb_sense/bins/TR_prox30_TSS_centered/allchr_sorted.txt",'r')
	reader = csv.reader(promFile, 'textdialect')
	for row in reader:
		name = row[3]
		geneToInfo[name][0]=float(row[6])
	promFile.close()
	
	dsFile = open("/home/raflynn/7SK/GROseq/Flynn/metagenes/" + sample + "comb_sense/bins/TR_dist30_TSS_centered/allchr_sorted.txt",'r')
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
	geneToInfo = collections.defaultdict(lambda: [0,0])
	promFile = open("/home/raflynn/7SK/GROseq/Flynn/metagenes/" + sample + "comb_sense/bins/TR_prox30_RY_enh_centered/allchr_sorted.txt",'r')
	reader = csv.reader(promFile, 'textdialect')
	for row in reader:
		name = row[3]
		geneToInfo[name][0]=float(row[6])
	promFile.close()
	
	dsFile = open("/home/raflynn/7SK/GROseq/Flynn/metagenes/" + sample + "comb_sense/bins/TR_dist30_RY_enh_centered/allchr_sorted.txt",'r')
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
	geneToInfo = collections.defaultdict(lambda: [0,0])
	promFile = open("/home/raflynn/7SK/GROseq/Flynn/metagenes/" + sample + "comb_sense/bins/TR_prox30_SE_indiv_centered/allchr_sorted.txt",'r')
	reader = csv.reader(promFile, 'textdialect')
	for row in reader:
		name = row[3]
		geneToInfo[name][0]=float(row[6])
	promFile.close()
	
	dsFile = open("/home/raflynn/7SK/GROseq/Flynn/metagenes/" + sample + "comb_sense/bins/TR_dist30_SE_indiv_centered/allchr_sorted.txt",'r')
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
		

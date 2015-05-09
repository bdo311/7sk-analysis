# getMaxNtSignal_v2.py
# from the START-seq bedgraphs that have been intersected with enhancers, take the reverse and forward peaks such that R < F and val(R)+val(F) are maximized
# add in the following things:
# peaks must be between -500 to +500bp
# sense peak must have at least signal = 5
# look for best sense peak first in the -500 to +1000 region ???
# then look for antisense peak upstream of it
# 	if antisense peak is the sense peak for another gene, bidirectional
#	otherwise, divergent
# if one does not exist, take best sense peak --> unidirectional gene

import csv, collections, sys, math
csv.register_dialect("textdialect",delimiter='\t')

def getPos(fn):
	ifile = open(fn,'r')
	reader = csv.reader(ifile, 'textdialect')
	tss = collections.defaultdict(lambda: [])
	for row in reader:
		tssName = '__'.join(row[4:])
		row[1:4] = [int(row[1]),int(row[2]),float(row[3])]
		tss[tssName].append(row[:4])
	ifile.close()
	return tss
	
# hash tables of TSS's to all positions with signal
print "Reading in signal"
tss_plus = getPos("startRNA_NELFwt_forward_pro_tss.bedgraph")
tss_minus = getPos("startRNA_NELFwt_reverse_pro_tss.bedgraph")


print "Getting plus TSS's"
tssToBestSense = {}
tssToBestAntisense = {}
tssToClass = {}
tssToInfo = {}
allPlusSense = collections.defaultdict(lambda: []) #all sense peaks found for the plus strand tss's

ofile = open("tss_height_comparison.txt",'w') 
writer = csv.writer(ofile, 'textdialect')

for tss in tss_plus:
	if tss.split('__')[-1] != '+': continue # we only want to look at TSS's for the plus strand
	tssName = tss.split('__')[3]
	
	if tssName[:3]=="Mir" or tssName[:4]=="Snhg" or tssName[:4]=="Snor": continue
	tssToInfo[tssName] = tss
	maxSignal = 0 
	secondMax = 0
	bestSense = []
	for signal in tss_plus[tss]:
		if signal[3] > maxSignal:
			secondMax = maxSignal #update second max
			maxSignal = signal[3]
			bestSense = signal
		elif signal[3] > secondMax:
			secondMax = signal[3]
	
	if maxSignal < 5 or bestSense==[]: continue # sense peak not found/not high enough, or lower limit for sense peaks
	if bestSense in allPlusSense[bestSense[0]]: continue # duplicate sense peak
	tssToBestSense[tssName] = bestSense
	allPlusSense[bestSense[0]].append(bestSense) #add this sense peak to the list of sense peaks from the first strand
	writer.writerow([maxSignal, secondMax, "sense"])

	if tss not in tss_minus: tssToClass[tssName] = "unidirectional"
	else:
		bestAS = []
		maxSignal = 0 #antisense peak can be any size
		secondMax = 0
		for signal in tss_minus[tss]:
			if signal[1] > bestSense[1]: continue #antisense peak has to be upstream of sense peak
			if signal[3] > maxSignal:
				secondMax = maxSignal #update second max
				maxSignal = signal[3]
				bestAS = signal
			elif signal[3] > secondMax:
				secondMax = signal[3]
				
		if bestAS == []: #no best antisense peak found
			tssToClass[tssName] = "unidirectional"
		else:
			tssToClass[tssName] = "divergent"
			tssToBestAntisense[tssName] = bestAS
			writer.writerow([maxSignal, secondMax, "antisense"])
					
print "Getting minus TSS's"
allMinusSense = collections.defaultdict(lambda: []) #all sense peaks found for the minus strand tss's

for tss in tss_minus:
	if tss.split('__')[-1] != '-': continue # we only want to look at TSS's for the minus strand
	tssName = tss.split('__')[3]
	if tssName[:3]=="Mir" or tssName[:4]=="Snhg" or tssName[:4]=="Snor": continue
	tssToInfo[tssName] = tss
	maxSignal = 0 
	secondMax = 0
	bestSense = []
	for signal in tss_minus[tss]:
		if signal[3] > maxSignal:
			secondMax = maxSignal #update second max
			maxSignal = signal[3]
			bestSense = signal
		elif signal[3] > secondMax:
			secondMax = signal[3]
			
	if maxSignal < 5 or bestSense==[]: continue # sense peak not found/not high enough, or lower limit for sense peaks
	if bestSense in allMinusSense[bestSense[0]]: continue # duplicate sense peak
	tssToBestSense[tssName] = bestSense 
	writer.writerow([maxSignal, secondMax, "sense"])

	if tss not in tss_plus: tssToClass[tssName] = "unidirectional"
	else:
		bestAS = []
		maxSignal = 0
		secondMax = 0

		for signal in tss_plus[tss]:
			if signal[1] < bestSense[1]: continue #antisense peak has to be upstream of sense peak
			if signal[3] > maxSignal:
				secondMax = maxSignal #update second max
				maxSignal = signal[3]
				bestAS = signal
			elif signal[3] > secondMax:
				secondMax = signal[3]
		tssToBestAntisense[tss] = bestAS 
		
		# is there an antisense peak?
		if bestAS == []: 
			tssToClass[tssName] = "unidirectional"
		else:
			tssToClass[tssName] = "divergent"
			tssToBestAntisense[tssName] = bestAS
			writer.writerow([maxSignal, secondMax, "antisense"])
			
			#is it bidirectional?
			if bestAS in allPlusSense[bestSense[0]]: 
				for otherTSS in tssToBestSense:
					if tssToBestSense[otherTSS] == bestAS:
						tssToClass[otherTSS] = "bidirectional" # other tss
						tssToClass[tssName] = "bidirectional" # tss
		
ofile.close()

# tssToBestSense = {}
# tssToBestAntisense = {}
# tssToClass = {}
# tssToInfo = {}
ofile = open("startRNA_NELFwt_tss_peaks.txt",'w')
writer = csv.writer(ofile, 'textdialect')

for tss in tssToBestSense:
	outputRow = tssToInfo[tss].split('__')
	if outputRow[-1] == '+': #plus strand
		if tss in tssToBestAntisense: outputRow.extend(tssToBestAntisense[tss])
		else: outputRow.extend(['.']*4)
		outputRow.extend(tssToBestSense[tss])
	else:
		outputRow.extend(tssToBestSense[tss])
		if tss in tssToBestAntisense: outputRow.extend(tssToBestAntisense[tss])
		else: outputRow.extend(['.']*4)
	outputRow.append(tssToClass[tss])
	writer.writerow(outputRow)
		
print "Generating peaks"

ifile = open("startRNA_NELFwt_tss_peaks.txt",'r')
reader = csv.reader(ifile, 'textdialect')
ofile = open("mm9_tss_startseq_centered_1kb.bed",'w')
writer = csv.writer(ofile, 'textdialect')
unfile = open("mm9_tss_startseq_uncentered_1kb.bed",'w')
unwriter = csv.writer(unfile, 'textdialect')

for row in reader:
	strand = row[5]
	try: dist = int(math.fabs(int(row[11])-int(row[7])))
	except: dist = 0
	if strand == '+': outputRow=[row[0],int(row[11])-1000,int(row[11])+1000,row[3],dist,'+']
	else: outputRow=[row[0],int(row[7])-1000,int(row[7])+1000,row[3],dist,'-']
	writer.writerow(outputRow)
	
	outputRow = row[:6]
	unwriter.writerow(outputRow)
	
ifile.close()
ofile.close()


	
	

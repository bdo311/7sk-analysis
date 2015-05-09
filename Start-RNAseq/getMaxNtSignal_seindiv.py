# getMaxNtSignal.py
# from the START-seq bedgraphs that have been intersected with enhancers, take the reverse and forward peaks such that R < F and val(R)+val(F) are maximized

import csv, collections, sys
csv.register_dialect("textdialect",delimiter='\t')

def getPos(fn):
	ifile = open(fn,'r')
	reader = csv.reader(ifile, 'textdialect')
	enh = collections.defaultdict(lambda: [])
	for row in reader:
		enhName = '__'.join(row[4:])
		row[1:4] = [int(row[1]),int(row[2]),float(row[3])]
		enh[enhName].append(row[:4])
	ifile.close()
	return enh
	
# hash tables of enhancer to all positions with signal
fenh = getPos("startRNA_NELFwt_forward_pro_SEindiv.bedgraph")
renh = getPos("startRNA_NELFwt_reverse_pro_SEindiv.bedgraph")

print "Getting reverse and forward starts"

# optimal enhancer reverse and forward starts
ofile = open("startRNA_NELFwt_SEindiv_peaks.txt",'w')
writer = csv.writer(ofile, 'textdialect')

for enh in fenh:
	#print enh
	if enh not in renh: continue
	bestR = []
	bestF = []
	maxSum = 0
	
	for rSignal in renh[enh]:
		for fSignal in fenh[enh]:			
			if rSignal[1] >= fSignal[1]: continue # r must be < f
			sigSum = rSignal[3] + fSignal[3]
			if sigSum > maxSum:
				#print enh, 'r', rSignal
				#print enh, 'f', fSignal
				maxSum = sigSum
				bestR = rSignal
				bestF = fSignal
				
	if maxSum == 0: continue #no reverse signal before forward signal
	outputRow = enh.split('__')
	outputRow.extend(bestR)
	outputRow.extend(bestF)
	#print enh, bestF, bestR
	outputRow.append(bestF[1]-bestR[1])
	writer.writerow(outputRow)
	
ofile.close()
				
			
print "Generating peaks"

ifile = open("startRNA_NELFwt_SEindiv_peaks.txt",'r')
reader = csv.reader(ifile, 'textdialect')
ofile = open("mES_SE_individual_startseq_centered_1kb.bed",'w')
writer = csv.writer(ofile, 'textdialect')
unfile = open("mES_SE_individual_startseq_uncentered_1kb.bed",'w')
unwriter = csv.writer(unfile, 'textdialect')

centers = set()

for row in reader:
	center = (int(row[10]) + int(row[14]))/2 # center between two START-seq peaks
	centerInfo = '_'.join([row[0],row[10],row[14]])
	if centerInfo in centers: continue
	else: centers.add(centerInfo)
	
	outputRow = [row[0],center-1000,center+1000] # -1kb to +1kb from center
	outputRow.append('__'.join(row[:9])) # original 1kb enhancer
	outputRow.append(row[-1]) # distance between peaks
	outputRow.append(row[8]) # strand
	writer.writerow(outputRow)
	
	outputRow = row[:3]
	outputRow.append('__'.join(row[3:7]))
	outputRow.extend(row[7:9])
	unwriter.writerow(outputRow)
	
ifile.close()
ofile.close()


	
	

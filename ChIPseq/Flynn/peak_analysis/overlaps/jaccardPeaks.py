# jaccardPeaks.py
# 1/20/15
# get pairwise jaccard scores

import csv, collections, sys
csv.register_dialect("textdialect", delimiter='\t')


#allSamples=['7sk','BAF155_Scr','BAF155_DMSO1','BAF155_DMSO2','HEXIM1_Scr1','HEXIM1_Scr2','HEXIM1','DDX21']
allSamples = ['BAF155_rep1', 'BAF155_rep2', 'HEXIM1_rep1', 'HEXIM1_rep2', '7SK_WT']

sampleToPeaks = collections.defaultdict(lambda: set())

ifile = open(sys.argv[1],'r')
reader = csv.reader(ifile, 'textdialect')
for row in reader:
	samples = row[6].split(',')
	peak = row[3]
	for s in samples:
		sampleToPeaks[s].add(peak)
		
ifile.close()

def jaccard(s,t):
	set_1 = sampleToPeaks[s]
	set_2 = sampleToPeaks[t]
	n = len(set_1.intersection(set_2))
	return n / float(len(set_1) + len(set_2) - n)
	
ofile = open(sys.argv[2],'w')
writer = csv.writer(ofile, 'textdialect')
header = ['']
header.extend(allSamples)
writer.writerow(header)

for s in allSamples:
	outputRow = [s]
	for t in allSamples:
		outputRow.append(jaccard(s,t))
	writer.writerow(outputRow)
	
ofile.close()


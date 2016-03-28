# collapseOverlapPeaks.py
# following the Vierstra 2014 paper to get definitive peak lists

import csv, os, sys
csv.register_dialect("textdialect", delimiter='\t')

try:
	cutoff = int(sys.argv[1])
except:
	print "Cutoff must be int; exiting."
	exit()
	
ifile = open("MasterList_overlap.bed",'r')
reader = csv.reader(ifile, 'textdialect')

# Pass 1: chr, start, stop for all peaks
peakToLoc = {}
for row in reader:
	try:
		peakToLoc[row[3]] = row[:3]
	except:
		print row

# Pass 2: get only the peaks with the highest q-score
ifile.seek(0)
badPeaks = set()
goodPeaks = {} # peak --> (information, samples with that peak)

for row in reader:
	peakName = row[3]
	if peakName in badPeaks: continue
	
	peaksInInterval = row[5].split('|')[-1].split(';')
	samples = []
	highestQ = 0
	
	# get highest peak
	for peakInfo in peaksInInterval:
		try:
			[sample, peakNum, q] = peakInfo.split('__')
		except:
			print row
		samples.append(sample)
		if float(q) > highestQ: highestQ = float(q)
	for peakInfo in peaksInInterval:
		[sample, peakNum, q] = peakInfo.split('__')
		q = float(q)
		if q == highestQ and q > cutoff:
			goodPeaks[peakInfo] = samples #good peak
		else: badPeaks.add(peakInfo)

		
# Write info for good peaks

ofn = "distinct_peaks_" + str(cutoff) + ".bed"
ofile = open(ofn, 'w')
writer = csv.writer(ofile, 'textdialect')

for peak in goodPeaks:
	loc = peakToLoc[peak]
	samples = goodPeaks[peak]
	outputRow = loc
	outputRow.append(peak)
	outputRow.extend([0, '+'])
	outputRow.append(','.join(set(samples)))
	writer.writerow(outputRow)
	
ofile.close()

os.system("sort -k1,1 -k2,2n " + ofn + " -o " + ofn)
os.system("python jaccardPeaks.py " + ofn + " jaccard_" + str(cutoff) + ".txt")

		
		
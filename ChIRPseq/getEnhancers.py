# getEnhancers.py
# 4/7/14
# Reads in a bedGraph and takes in enhancer annotation file, and gets enhancer and super enhancer read density

import csv
import collections
import math
import multiprocessing
import os

csv.register_dialect("textdialect", delimiter='\t')
currDir = "/home/raflynn/ChIRPseq/WT2_new/"
os.system("mkdir bins/enhancers/")

# Reading in massive bedgraph file: ~20GB memory, 8 minutes CPU time
# ------------------------------------------------------------------

print "Loading bedgraph"

ifile = open(currDir + "combined.bedGraph")
# ifile = open("/home/raflynn/7SK_ChIRPseq/test/mergeMed.bedGraph")

reader = csv.reader(ifile, 'textdialect')

# by chromosome, and by bin
reads = collections.defaultdict(lambda: collections.defaultdict(lambda: [])) 

for row in reader:
	start = int(row[1])
	end = int(row[2])
	score = float(row[3])

	# which bins does each interval go into?
	multipleOf500k = start/500000
	bin1 = multipleOf500k
	bin2 = bin1 - 1

	# result: bins are overlapping
	reads[row[0]][bin1].append([start, end, score])
	reads[row[0]][bin2].append([start, end, score])

ifile.close()

# Reading in enhancer file
# --------------------

print "Loading gene files"

# get enhancers and super enhancers
ifile = open("/home/raflynn/ChIRPseq/genes/mES_enhancers.txt", 'r')
reader = csv.reader(ifile, 'textdialect')

enh = collections.defaultdict(lambda: []) # by chromosome
super = collections.defaultdict(lambda: [])
reader.next()

for row in reader:
	if row[5] == 'YES': super[row[1]].append(row[:6])
	else: enh[row[1]].append(row[:6])

ifile.close()

# Processing bins
# ----------------

# gets the average score within this bin
def getAverageScore(currStart, currEnd, chrom, binNum, rn):
	currBin = reads[chrom][binNum]
	totalLength = currEnd - currStart + 1
	totalScore = 0

	readNumber = rn 			# which read to start at
	for read in currBin[rn:]:
		readStart, readEnd, score = read[0], read[1], read[2]
		#print 'read', readStart, readEnd, score

		if readEnd < currStart: 
			readNumber += 1 	# never read that read again
			continue 			# moves on to next read
		if currEnd < readStart: 
			#print 'totalScore', totalScore
			break 				# end of gene < start of read. No coverage, so score is 0

		if currStart < readStart: currStart = readStart # ignore zeros that are in the read region
		if readEnd >= currEnd: 	# end of gene is within the read
			totalScore += score * (currEnd - currStart + 1)
			#print 'totalScore', totalScore, currEnd - currStart + 1
			break				# don't advance read number because the read extends past end of gene
		else: # end of gene is past the read
			totalScore += score * (readEnd - currStart + 1)
			#print 'totalScore', totalScore, readEnd - currStart + 1
			readNumber += 1 	# advance read number because this read won't be needed anymore

		currStart = readEnd + 1

	return (float(totalScore), readNumber)

# gets the array of bins for each gene
def getBins(chrom, start, end, numBins):
	spacingPerBin = end - start + 1 # one bin for all

	binNum = start/500000

	scores = []
	currStart = start
	readNumber = 0 #index of read, updated each time
	while currStart < end:
		currEnd = currStart + spacingPerBin - 1 # end of my window
		if currEnd > end: currEnd = end # last bin can't go past TES
		
		#print 'curr', currStart, currEnd
		score, readNumber = getAverageScore(currStart, currEnd, chrom, binNum, readNumber) #updates read number also
		scores.append(score)
		currStart = currEnd + 1 # new start of my window

	return scores

# Looping through all genes
# -------------------------

# loop through genes in each chromosome. This does gene bodies
def geneBodyWorker(chrom):	
	ofile = open(currDir + "bins/enhancers/" + chrom + ".txt", 'w')
	writer = csv.writer(ofile, 'textdialect')

	for gene in enh[chrom]:
		if chrom not in reads: continue
		start = int(gene[2])
		end = int(gene[3])

		geneBins = getBins(chrom, start, end, 1)
		outputRow = gene
		outputRow.append(sum(geneBins))
		writer.writerow(outputRow)
	
	for gene in super[chrom]:
		if chrom not in reads: continue
		start = int(gene[2])
		end = int(gene[3])

		geneBins = getBins(chrom, start, end, 1)
		outputRow = gene
		outputRow.append(sum(geneBins))
		writer.writerow(outputRow)

	ofile.close()

# Multiprocessing for each chromosome for gene body
print "Working on gene bodies"

procs = []
for chrom in enh:
	p = multiprocessing.Process(target=geneBodyWorker, args=(chrom,))
	procs.append(p)
	p.start()
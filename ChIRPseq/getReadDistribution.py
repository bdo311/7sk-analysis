# getReadDistribution.py
# 6/19/14; updated to read bedgraph in parallel
# Reads in a bedGraph and takes in whole-genome annotation file, and adds up scores. Very similar to get7SK.py

import re, csv, collections, math, multiprocessing, sys, glob
csv.register_dialect("textdialect", delimiter='\t')
ifolder = sys.argv[1]
	
# 1. Reading in annotations file
# ------------------------------

print "Loading gene files"

ifile = open("/home/raflynn/ChIRPseq/genes/mm9_newannotation.txt", 'r')
reader = csv.reader(ifile, 'textdialect')

genes = collections.defaultdict(lambda: []) # by chromosome
reader.next()

for row in reader:
	genes[row[1]].append(row)

ifile.close()

# 2. Reading in bedgraph file in sets
# -----------------------------------

print "Loading bedgraph"

# Loading reads for each bedgraph, for each chromosome
def getReads(readQueue, chrom, graph):
	ifile = open(graph, 'r')
	reader = csv.reader(ifile, 'textdialect')
	readsForChrom = {}

	for row in reader:
		start = int(row[1])
		end = int(row[2])
		score = float(row[3])

		# which bins does each interval go into?
		multipleOf500k = start/500000
		bin1 = multipleOf500k
		bin2 = bin1 - 1
		
		# result: bins are overlapping
		if bin1 not in readsForChrom.keys(): readsForChrom[bin1] = []
		readsForChrom[bin1].append([start, end, score])
		if bin2 not in readsForChrom.keys(): readsForChrom[bin2] = []
		readsForChrom[bin2].append([start, end, score])

	ifile.close()
	readQueue.put((chrom, readsForChrom))
	print chrom, 'done'

# get the read file for each chromosome
def readBedGraph(chroms):
	manager = multiprocessing.Manager()
	readQueue = manager.Queue()

	# get all chromosomes in a queue
	procs = []
	for chrom in chroms:
		print chrom
		p = multiprocessing.Process(target=getReads, args=(readQueue, chrom, ifolder + '/bedGraphByChr/' + chrom + '.bedGraph'))
		p.start()
		procs.append(p)
	for proc in procs: proc.join()

	# make a large dictionary
	print "Making dictionary"
	reads = {}
	while not readQueue.empty():
		readTuple = readQueue.get()
		reads[readTuple[0]] = readTuple[1]
	return reads

# 3. Processing bins
# ------------------

# gets the average score within this bin
def getAverageScore(chrom, currStart, currEnd, binNum, rn):
	totalLength = currEnd - currStart + 1
	totalScore = 0
	readNumber = rn 			# which read to start at

	if binNum not in reads[chrom]: return 0, readNumber
	readList = reads[chrom][binNum][rn:]
	for read in readList:
		readStart, readEnd, score = read[0], read[1], read[2]

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

	return float(totalScore), readNumber - 10 # as a buffer, reread the prior read

# gets the array of bins for each gene
def getBins(chrom, start, end, numBins, rn):
	spacingPerBin = end - start + 1 # one bin for all

	binNum = start/500000

	scores = []
	currStart = start
	readNumber = rn #index of read, updated each time
	while currStart < end:
		currEnd = currStart + spacingPerBin - 1 # end of my window
		if currEnd > end: currEnd = end # last bin can't go past TES
		score, readNumber = getAverageScore(chrom, currStart, currEnd, binNum, readNumber) #updates read number also
		scores.append(score)
		currStart = currEnd + 1 # new start of my window

	return scores, readNumber

# Looping through all genes
# -------------------------

# loop through genes in each chromosome. This does gene bodies
def geneBodyWorker(chrom, genes):
	ofile = open(ifolder + "/bins/newannot/" + chrom + ".txt", 'w')
	writer = csv.writer(ofile, 'textdialect')

	readNumber = 0
	counter = 0
	for gene in genes[chrom]:
		start = int(gene[2])
		end = int(gene[3])
		if start < 0: continue
		if end > 1e9: continue

		counter += 1
		print counter, readNumber, gene
		geneBins, readNumber = getBins(chrom, start, end, numBins=1, rn=0)
		outputRow = gene
		outputRow.append(sum(geneBins))
		writer.writerow(outputRow)

	ofile.close()

# Multiprocessing for each chromosome for gene body
def processGeneBodies(chroms, genes):
	print "Working on gene bodies"

	procs = []
	for chrom in chroms:
		p = multiprocessing.Process(target=geneBodyWorker, args=(chrom, genes))
		procs.append(p)
		p.start()

	for p in procs:
		p.join()

# 4. run all this on sets of chromosomes
# --------------------------------------

allChroms = ['chr' + str(x) for x in range(1,20)]
allChroms.extend(['chrX', 'chrY', 'chrM'])

for i in range(4):
	chroms = allChroms[(6*i):(6*i+6)]
	# chroms = ['chr12']
	reads = readBedGraph(chroms)
	processGeneBodies(chroms, genes)	
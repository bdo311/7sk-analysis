# makeEmptyFolders_human.py
# 5/23/14; based on makeEmptyFolders.py
# Makes bins and geneBody, tes, and tss folders for each folder

import os, glob, csv, re, collections, math, multiprocessing

csv.register_dialect("textdialect", delimiter='\t')

# Get all folders I'll be operating on. seven_Folders is all the 7SK ChIRPseq experiments. chip_Folders is all the ChIP experiments
# -------------------------------------------------------------------

seven_Folders = ['H1', 'HeLa']
chirpDir = "/home/raflynn/ChIRPseq/"

for folder in seven_Folders:
	os.chdir(chirpDir + folder + '/')
	os.system("mkdir bins inputbins")
	os.chdir(chirpDir + folder + '/bins/')
	os.system("mkdir geneBody tss tes")
	os.chdir(chirpDir + folder + '/inputbins/')
	os.system("mkdir geneBody tss tes")

# Get a map of bin folder to associated bedgraph
# ------------------------------------------

folderToGraph = {}
for folder in seven_Folders:
	binFolder = chirpDir + folder + '/bins/'
	graph = chirpDir + folder + '/combined_norm.bedGraph'
	folderToGraph[binFolder] = graph

	binFolder = chirpDir + folder + '/inputbins/'
	graph = chirpDir + folder + '/input_norm.bedGraph'
	folderToGraph[binFolder] = graph
	

# Reading in gene files
# ---------------------

# Reading in gene body file
print "Loading gene files"

ifile = open(chirpDir + "genes/hg19_refseq_coll.txt", 'r')
reader = csv.reader(ifile, 'textdialect')

genes = collections.defaultdict(lambda: []) # by chromosome
reader.next()

for row in reader:
	genes[row[1]].append(row)

ifile.close()

# Reading in gene TSS/TES file
ifile = open(chirpDir + "genes/hg19_refseq_startend.txt", 'r')
reader = csv.reader(ifile, 'textdialect')

allEnds = collections.defaultdict(lambda: []) # by chromosome

reader.next()
for row in reader:
	allEnds[row[1]].append(row)

ifile.close()

# File processing methods
# -----------------------

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

	return (float(totalScore)/totalLength, readNumber)

# gets the array of bins for each gene
def getBins(chrom, start, end, numBins):
	#spacingPerBin = int(math.ceil((end - start + 1)/6.0)) 
	spacingPerBin = int(math.ceil((end - start)/float(numBins)))

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


# Now time to process files! Each folder will take approximately 20 minutes (with 3 matchings). There are about 36 folders to go through. All this stuff is taken from processMerge.py
# --------------------------

for folder in folderToGraph:
	graph = folderToGraph[folder]

	# Reading in massive bedgraph file: ~20GB memory, 8 minutes CPU time

	print "Loading bedgraph"

	ifile = open(graph, 'r')
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


	# Looping through all genes
	# -------------------------

	# loop through genes in each chromosome. This does gene bodies
	def geneBodyWorker(chrom):	
		ofile = open(folder + "/geneBody/" + chrom + ".txt", 'w')
		writer = csv.writer(ofile, 'textdialect')

		# # testing
		# writer.writerow(['a'])
		# ofile.close()
		# return

		for gene in genes[chrom]:
			if chrom not in reads: continue
			start = int(gene[3])
			end = int(gene[4])

			if end - start < 200 or end - start > 200000: continue #binning doesn't really make sense here so just ignore

			geneBins = getBins(chrom, start, end, 1000)
			strand = gene[2]
			if strand == '-': geneBins = geneBins[::-1] # reading from end to start if antisens

			outputRow = gene
			outputRow.append(sum(geneBins)) # can be sorted later to get total reads (this definition will need to be changed since it needs to be unscaled)
			outputRow.extend(geneBins)
			writer.writerow(outputRow)

		ofile.close()

	# Multiprocessing for each chromosome for gene body
	print "Working on gene bodies"

	procs = []
	for chrom in genes:
		p = multiprocessing.Process(target=geneBodyWorker, args=(chrom,))
		procs.append(p)
		p.start()

	# wait till all are done before moving on
	for p in procs:
		p.join()

	# loop through genes in each chromosome. This does TSS/TES
	def geneEndsWorker(chrom):	
		tssFile = open(folder + "/tss/" + chrom + ".txt", 'w')
		tssWriter = csv.writer(tssFile, 'textdialect')

		tesFile = open(folder + "/tes/" + chrom + ".txt", 'w')
		tesWriter = csv.writer(tesFile, 'textdialect')

		# # testing
		# tssWriter.writerow(['a'])
		# tesWriter.writerow(['a'])
		# tssFile.close()
		# tesFile.close()
		# return

		for gene in allEnds[chrom]:
			if chrom not in reads: continue
			tssStart, tssEnd = int(gene[3]), int(gene[4])
			tesStart, tesEnd = int(gene[5]), int(gene[6])

			# ignore too short or too long RNAs
			start = (tssEnd + tssStart)/2
			end = (tesEnd + tesStart)/2
			if end - start < 200 or end - start > 200000: continue 

			# get bins. 400 bins since 2000/5
			tssBins = getBins(chrom, tssStart, tssEnd, 400)
			tesBins = getBins(chrom, tesStart, tesEnd, 400)

			# if negative, tss becomes tes and vice versa, and gene bins all flip
			strand = gene[2]
			tssRow = [gene[0], gene[1], gene[2]]
			tesRow = [gene[0], gene[1], gene[2]]
			if strand == '-': 
				temp = tesBins
				tesBins = tssBins[::-1]
				tssBins = temp[::-1]
				tssRow.extend([tesStart, tesEnd])
				tesRow.extend([tssStart, tssEnd])
			else:
				tssRow.extend([tssStart, tssEnd])
				tesRow.extend([tesStart, tesEnd])

			tssRow.extend([gene[7], sum(tssBins)])
			tssRow.extend(tssBins)
			tesRow.extend([gene[7], sum(tesBins)])
			tesRow.extend(tesBins)

			tssWriter.writerow(tssRow)
			tesWriter.writerow(tesRow)

		tssFile.close()
		tesFile.close()

	# Multiprocessing for each chromosome for gene TSS/TES
	print "Working on gene TSS and TES"

	for chrom in allEnds:
		p = multiprocessing.Process(target=geneEndsWorker, args=(chrom,))
		p.start()





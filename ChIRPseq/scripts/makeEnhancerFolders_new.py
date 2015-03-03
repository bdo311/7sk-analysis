# makeEnhancerFolders_new.py
# 4/9/14; does regular enhancers and super enhancers
# Makes enhancer bins for each folder

import os, glob, csv, re, collections, math, multiprocessing

csv.register_dialect("textdialect", delimiter='\t')

# Get all folders I'll be operating on. seven_Folders is all the 7SK ChIRPseq experiments. chip_Folders is all the ChIP experiments
# -------------------------------------------------------------------

seven_Folders = []

chirpDir = "/home/raflynn/ChIRPseq/"

for folder in seven_Folders:
	os.chdir(chirpDir + folder + '/bins/')
	os.system("mkdir extreg")

# getting ChIP folders

chipDir = "/home/raflynn/mmchipseq/mes_all/"

chipFiles = glob.glob(chipDir + "ChIP_*.bedGraph")
chipFolders = []
for fileName in chipFiles:
	exp = re.search("ChIP_(.*)_new.bedGraph", fileName).group(1)
	chipFolders.append(exp)
	os.chdir(chipDir + exp + '/')
	os.system("mkdir extreg")

# Get a map of bin folder to associated bedgraph
# ------------------------------------------

folderToGraph = {}

for folder in seven_Folders:
	binFolder = chirpDir + folder + '/bins/'
	graph = chirpDir + folder + '/combined_norm.bedGraph'
	folderToGraph[binFolder] = graph

for folder in chipFolders:
	binFolder = chipDir + folder
	graph = chipDir + 'ChIP_' + folder + '_new.bedGraph'
	folderToGraph[binFolder] = graph

# Reading in enhancer files
# ---------------------

print "Loading enhancer files"

# regular enhancers
ifile = open(chirpDir + "genes/mES_reg_enhancers_centered10.txt", 'r')
reader = csv.reader(ifile, 'textdialect')

reg = collections.defaultdict(lambda: []) # by chromosome
reader.next()

for row in reader:
	reg[row[1]].append(row)

ifile.close()

# super enhancers
# ifile = open(chirpDir + "genes/mES_super_enhancers_expanded.txt", 'r')
# reader = csv.reader(ifile, 'textdialect')

# super = collections.defaultdict(lambda: []) # by chromosome
# reader.next()

# for row in reader:
	# super[row[1]].append(row)

# ifile.close()


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
	print folder
	
	# if len(glob.glob(folder + "/reg/*")) != 0: continue
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


	# Looping through all regular enhancers
	# -------------------------

	# loop through regular enhancers in each chromosome.
	def regWorker(chrom):	
		ofile = open(folder + "/extreg/" + chrom + ".txt", 'w')
		writer = csv.writer(ofile, 'textdialect')

		for enhancer in reg[chrom]:
			if chrom not in reads: continue
			start = int(enhancer[2])
			end = int(enhancer[3])

			enhBins = getBins(chrom, start, end, 400)
			
			outputRow = enhancer
			outputRow.append(sum(enhBins))
			outputRow.extend(enhBins)
			writer.writerow(outputRow)

		ofile.close()

	# Multiprocessing for each chromosome for regular enhancers
	print "Working on regular enhancers"

	procs = []
	for chrom in reg:
		p = multiprocessing.Process(target=regWorker, args=(chrom,))
		procs.append(p)
		p.start()

	# wait till all are done before moving on
	for p in procs:
		p.join()
		
	# # loop through genes in each chromosome. This does super enhancers
	# # def superWorker(chrom):	
		# # ofile = open(folder + "/super/" + chrom + ".txt", 'w')
		# # writer = csv.writer(ofile, 'textdialect')
		
		# # for enhancer in super[chrom]:
			# # if chrom not in reads: continue
			# # start = int(enhancer[2])
			# # end = int(enhancer[3])

			# # enhBins = getBins(chrom, start, end, 1000)
			
			# # outputRow = enhancer
			# # outputRow.append(sum(enhBins))
			# # outputRow.extend(enhBins)
			# # writer.writerow(outputRow)
			
		# # ofile.close()

	# Multiprocessing for each chromosome for super enhancers
	# print "Working on super enhancers"

	# for chrom in super:
		# p = multiprocessing.Process(target=superWorker, args=(chrom,))
		# p.start()

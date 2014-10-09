# makeBins.py
# 5/27/14, based on makeFolders in the past
# Gets bins for all treatments at all regions (genes, enhancers, and 7SK peaks) 

import os, glob, csv, re, collections, math, multiprocessing

csv.register_dialect("textdialect", delimiter='\t')

folderToGraph = {}
chipDir = "/home/raflynn/mmchipseq/mes_all/"
chirpDir = "/home/raflynn/ChIRPseq/"

# -------------------------------------------------------------------
# 1. Get ChIRP folders
# -------------------------------------------------------------------

# seven_Folders = ['H1', 'HeLa']
seven_Folders = ['ActD_new', 'DRB_new', 'Flavo_new', 'JQ11_new', 'JQ14_new', 'WT2_new']

for folder in seven_Folders:
	os.chdir(chirpDir + folder + '/')
	os.system("mkdir bins")
	os.chdir(chirpDir + folder + '/bins/')
	os.system("mkdir geneBody tss tes extsuper extreg extpeaks ctrpeaks")

for folder in seven_Folders:
	binFolder = chirpDir + folder + '/bins/'
	# graph = chirpDir + 'test/mergeMed.bedGraph' #for testing
	graphFolder = chirpDir + folder + '/bedGraphByChr/'
	folderToGraph[binFolder] = graphFolder

# -------------------------------------------------------------------
# 2. Get ChIP folders
# -------------------------------------------------------------------

# chipFiles = glob.glob(chipDir + "ChIP_*.bedGraph")
# chipFolders = []
# for fileName in chipFiles:
	# exp = re.search("ChIP_(.*)_new.bedGraph", fileName).group(1)
	# chipFolders.append(exp)
	# os.chdir(chipDir + exp + '/')
	# os.system("mkdir geneBody tss tes extsuper extreg extpeaks ctrpeaks")
	
# for folder in chipFolders:
	# binFolder = chipDir + folder
	# graph = chipDir + 'ChIP_' + folder + '_new.bedGraph'
	# folderToGraph[binFolder] = graph

# ---------------------
# 3. Reading in gene files
# ---------------------

def getChrToRegion(fn, chrCol, header):
	ifile = open(fn, 'r')
	reader = csv.reader(ifile, 'textdialect')
	
	regions = collections.defaultdict(lambda: []) # by chromosome
	if header: reader.next()
	for row in reader:
		regions[row[chrCol]].append(row)
	ifile.close()
	return regions

print "Loading region files"
genes = getChrToRegion(chirpDir + "genes/mm9_refseq_coll.txt", 1, True)
allEnds = getChrToRegion(chirpDir + "genes/mm9_refseq_startend.txt", 1, True)
#regExt = getChrToRegion(chirpDir + "genes/mES_reg_enhancers_centered10.txt", 1, False)
#superExt = getChrToRegion(chirpDir + "genes/mES_super_enhancers_expanded3.txt", 1, False)
# peaksCtr = getChrToRegion(chirpDir + "WT2_new/peaks/finalpeaks_over200_ctr3.bed", 0, False)
# peaksExt = getChrToRegion(chirpDir + "WT2_new/peaks/finalpeaks_over200_ext3.bed", 0, False)

# peaksCtr_H1 = getChrToRegion(chirpDir + "H1/peaks/finalpeaks_over200_ctr3.bed", 0, False)
# peaksExt_H1 = getChrToRegion(chirpDir + "H1/peaks/finalpeaks_over200_ext3.bed", 0, False)
# peaksCtr_HeLa = getChrToRegion(chirpDir + "HeLa/peaks/finalpeaks_over200_ctr3.bed", 0, False)
# peaksExt_HeLa = getChrToRegion(chirpDir + "HeLa/peaks/finalpeaks_over200_ext3.bed", 0, False)


# -----------------------
# 4. File processing methods
# -----------------------

# gets the average score within this bin
def getAverageScore(currStart, currEnd, chrReads, binNum, rn):
	currBin = chrReads[binNum]
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
def getBins(chrReads, start, end, numBins):
	spacingPerBin = int(math.ceil((end - start)/float(numBins)))
	binNum = start/500000

	scores = []
	currStart = start
	readNumber = 0 #index of read, updated each time
	while currStart < end:
		currEnd = currStart + spacingPerBin - 1 # end of my window
		if currEnd > end: currEnd = end # last bin can't go past TES
		
		score, readNumber = getAverageScore(currStart, currEnd, chrReads, binNum, readNumber) #updates read number also
		scores.append(score)
		currStart = currEnd + 1 # new start of my window

	return scores

# --------------------------
# 5. Now time to process files! (using 4)
# --------------------------

# for each chromosome, get bins corresponding to each region in the chromosome
def regionWorker(regionType, chrom, chrToRegion, graphFolder, startCol, endCol, stranded, strandCol, limitSize, numBins):
	# get the read file for each chromosome
	if not glob.glob(graphFolder + chrom): continue
	fn = graphFolder + chrom + '.txt' 
	chrReads = getReads(fn)

	# get the bins for all regions
	ofile = open(folder + "/" + regionType + "/" + chrom + ".txt", 'w')
	writer = csv.writer(ofile, 'textdialect')

	for region in chrToRegion[chrom]:
		start = int(region[startCol])
		end = int(region[endCol])

		if limitSize:
			if end - start < 200 or end - start > 200000: continue #binning doesn't really make sense here so just ignore

		regionBins = getBins(chrReads, start, end, numBins)
		strand = region[strandCol]
		
		if stranded:
			if strand == '-': regionBins = regionBins[::-1] # reading from end to start if antisens

		outputRow = region
		outputRow.append(sum(regionBins)) # can be sorted later to get total reads (this definition will need to be changed since it needs to be unscaled)
		outputRow.extend(regionBins)
		writer.writerow(outputRow)

	ofile.close()

# Loading reads for each bedgraph, for each chromosome
def getReads(graph):
	ifile = open(graph, 'r')
	reader = csv.reader(ifile, 'textdialect')
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
		reads[bin1].append([start, end, score])
		reads[bin2].append([start, end, score])

	ifile.close()
	return reads
	
def regionProcess(regionType, chrToRegion, graphFolder, startCol, endCol, stranded, strandCol, limitSize, numBins):
	print "Working on " + regionType
	procs = []
	for chrom in chrToRegion:
		p = multiprocessing.Process(target=regionWorker, args=(regionType, chrom, chrToRegion, graphFolder, startCol, endCol, stranded, strandCol, limitSize, numBins))
		procs.append(p)
		p.start()
	for p in procs:
		p.join()	# wait till all are done before moving on
		
for folder in folderToGraph:
	graphFolder = folderToGraph[folder]

	# Reading in massive bedgraph file

	print "Loading bedgraph"
	# reads = getReads(graph)

	# for human only
	if 'H1' in folder:
		peaksCtr = peaksCtr_H1
		peaksExt = peaksExt_H1
	elif 'HeLa' in folder:
		peaksCtr = peaksCtr_HeLa
		peaksExt = peaksExt_HeLa
	
	# Looping through regions: (regionType, chrToRegion, graphFolder, startCol, endCol, stranded, strandCol, limitSize, numBins)
	regionProcess("geneBody", genes, graphFolder, 3, 4, True, 2, True, 1000)
	#regionProcess("extreg", regExt, graphFolder, 2, 3, False, 0, False, 400)
	#regionProcess("extsuper", superExt, graphFolder, 2, 3, False, 0, False, 1000)
	# regionProcess("ctrpeaks", peaksCtr, graphFolder, 1, 2, False, 0, False, 1000)
	# regionProcess("extpeaks", peaksExt, graphFolder, 1, 2, False, 0, False, 1000)
	
	# loop through genes in each chromosome. This does TSS/TES. This is different so can't use the above functions
	def geneEndsWorker(chrom):	
		tssFile = open(folder + "/tss/" + chrom + ".txt", 'w')
		tssWriter = csv.writer(tssFile, 'textdialect')

		tesFile = open(folder + "/tes/" + chrom + ".txt", 'w')
		tesWriter = csv.writer(tesFile, 'textdialect')


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
	procs = []
	for chrom in allEnds:
		p = multiprocessing.Process(target=geneEndsWorker, args=(chrom,))
		procs.append(p)
		p.start()
	for p in procs:
		p.join()	

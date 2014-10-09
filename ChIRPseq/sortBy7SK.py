# sortBy7SK.py
# 3/16/14
# 1. Gets all folders and files
# 2. Sorts chr files based on 7SK genebody

import os, glob, csv, re, multiprocessing

csv.register_dialect("textdialect", delimiter='\t')
chirpDir = "/home/raflynn/7SK_ChIRPseq/"
chipDir = "/home/raflynn/Mm_ChIPseq/mES_sra_data/"

# Get order of 7SK genebody genes
# -------------------------------

ifile = open(chirpDir + "WT2/bins/geneBody/allchr_sorted.txt", 'r')
reader = csv.reader(ifile, 'textdialect')

orderedGenes = []
counter = 0
for row in reader: 			
	counter += 1
	if counter == 8000: break			# we only want the top 8000 genes
	orderedGenes.append(row[0])
	
ifile.close()

# Get all folders I'll be operating on. seven_Folders is all the 7SK ChIRPseq experiments. chip_Folders is all the ChIP experiments
# -------------------------------------------------------------------

seven_Folders = ['ActD', 'DRB', 'Flavo', 'JQ11', 'JQ14', 'WT2']

# getting ChIP folders

chipFiles = glob.glob(chipDir + "ChIP_*.bedGraph")
chipFolders = []
for fileName in chipFiles:
	exp = re.search("ChIP_(.*).bedGraph", fileName).group(1)
	chipFolders.append(exp)

# Get a map of bin folder to associated bedgraph
# ------------------------------------------

folderToGraph = {}

for folder in seven_Folders:
	binFolder = chirpDir + folder + '/bins/'
	graph = chirpDir + folder + '/merge.bedGraph'
	folderToGraph[binFolder] = (graph, folder)

for folder in chipFolders:
	binFolder = chipDir + folder + '/'
	graph = chipDir + 'ChIP_' + folder + '.bedGraph'
	folderToGraph[binFolder] = (graph, folder)

# Make a new sorted file
# ----------------------

def processFile():
	ifile = open("allchr_sorted.txt", 'r')
	reader = csv.reader(ifile, 'textdialect')
	
	# place all rows into a hash table
	geneToData = {}
	for row in reader:
		geneToData[row[0]] = row
	
	ifile.close()
	
	# write new file with rows in order of 7SK order
	ofile = open("allchr_sorted_byWT2.txt", 'w')
	writer = csv.writer(ofile, 'textdialect')
	
	for gene in orderedGenes:
		writer.writerow(geneToData[gene])
		
	ofile.close()

def folderWorker(start, end, folders):
	for i in range(start, end):
		if i >= len(folders): continue
		folder = folders[i]

		folderName = folderToGraph[folder][1]
		print folderName
		
		os.chdir(folder + 'geneBody/')
		processFile()
		os.chdir(folder + 'tss/')
		processFile()
		os.chdir(folder + 'tes/')
		processFile()

folders = folderToGraph.keys()
nproc = 8
numPerProc = len(folders)/8 + 1
for i in range(nproc):
	p = multiprocessing.Process(target=folderWorker, args=(i * numPerProc, (i + 1) * numPerProc, folders))
	# print i * numPerProc, (i + 1) * numPerProc
	p.start()

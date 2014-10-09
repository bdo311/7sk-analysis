# mergeSort.py
# 3/5/14, updated 3/26/14 to process new 7SK ChIRP combined files; updated 4/8/14 to use combined_norm.bedGraph
# 1. Gets all folders and files
# 2. Merges and sorts chr files for all folders
# 3. Runs R on all of them
# 4. Merges genebody, tss, and tes files to get one average file for all

import os, glob, csv, re, multiprocessing

csv.register_dialect("textdialect", delimiter='\t')

# Get all folders I'll be operating on. seven_Folders is all the 7SK ChIRPseq experiments. chip_Folders is all the ChIP experiments
# -------------------------------------------------------------------

seven_Folders = ['ActD', 'DRB', 'Flavo', 'JQ11', 'JQ14', 'WT2']
chirpDir = "/home/raflynn/7SK_ChIRPseq/"

# getting ChIP folders

chipDir = "/home/raflynn/Mm_ChIPseq/mES_sra_data/"

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
	graph = chirpDir + folder + '/combined_norm.bedGraph'
	folderToGraph[binFolder] = (graph, folder)

for folder in chipFolders:
	binFolder = chipDir + folder + '/'
	graph = chipDir + 'ChIP_' + folder + '.bedGraph'
	folderToGraph[binFolder] = (graph, folder)

# Merge and sort
# --------------
# First merge, then sort
# def mergeSort(folder):
	# os.chdir(folder)
	# os.system("cat chr* > allchr.txt")
	# os.system("sort -rn -t $'\t' -k7,7 allchr.txt > allchr_sorted.txt")

# for folder in folderToGraph:
	# mergeSort(folder + 'geneBody/')
	# mergeSort(folder + 'tss/')
	# mergeSort(folder + 'tes/')

# Run R file to get plot and average bins
# ---------------------------------------

scriptFile = "/home/raflynn/7SK_ChIRPseq/visualizeBins_new8000.r"
def folderWorker(start, end, folders):
	for i in range(start, end):
		if i >= len(folders): continue
		folder = folders[i]

		folderName = folderToGraph[folder][1]
		os.chdir(folder + 'geneBody/')
		os.system("Rscript " + scriptFile + " gene " + folderName)
		os.chdir(folder + 'tss/')
		os.system("Rscript " + scriptFile + " tss " + folderName)
		os.chdir(folder + 'tes/')
		os.system("Rscript " + scriptFile + " tes " + folderName)
		# if not glob.glob("*.RData"): os.system("Rscript " + scriptFile + " tes " + folderName)

folders = folderToGraph.keys()
nproc = 8
numPerProc = len(folders)/8 + 1
procs = []
for i in range(nproc):
	p = multiprocessing.Process(target=folderWorker, args=(i * numPerProc, (i + 1) * numPerProc, folders))
	# print i * numPerProc, (i + 1) * numPerProc
	procs.append(p)
	p.start()

for p in procs:
	p.join()
	
# merge all files, for both scaled and raw averages
# -------------------------------------------------

# add chip folders to folderToGraph
chipDir = "/home/raflynn/Mm_ChIPseq/mES_sra_data/"

chipFiles = glob.glob(chipDir + "ChIP_*.bedGraph")
chipFolders = []
for fileName in chipFiles:
	exp = re.search("ChIP_(.*).bedGraph", fileName).group(1)
	chipFolders.append(exp)

for folder in chipFolders:
	binFolder = chipDir + folder + '/'
	graph = chipDir + 'ChIP_' + folder + '.bedGraph'
	folderToGraph[binFolder] = (graph, folder)
	
# take in an avg_*_* file and extract the bin values
def processFile(fileName):
	avgFile = glob.glob(fileName + "*")[0]
	ifile = open(avgFile, 'r')
	reader = csv.reader(ifile, 'textdialect')

	values = []
	reader.next()
	for row in reader:
		values.append(row[1])

	ifile.close()
	return values

# write the mapping to a file
def writeFile(name, mapping):
	ofile = open('/home/raflynn/7SK_ChIRPseq/averages/' + name + '.txt', 'w')
	writer = csv.writer(ofile, 'textdialect')
	for treatment in mapping:
		outputRow = [treatment]
		outputRow.extend(mapping[treatment])
		writer.writerow(outputRow)

	ofile.close()

# 6 hash tables	
rawgene = {}
rawtes = {}
rawtss = {}

scaledgene = {}
scaledtes = {}
scaledtss = {}

# go through each folder and get a bin array
for folder in folderToGraph:
	folderName = folderToGraph[folder][1]

	os.chdir(folder + 'geneBody/')
	rawgene[folderName] = processFile("avgraw_8000")
	scaledgene[folderName] = processFile("avgscaled_8000")

	os.chdir(folder + 'tes/')
	rawtes[folderName] = processFile("avgraw_8000")
	scaledtes[folderName] = processFile("avgscaled_8000")

	os.chdir(folder + 'tss/')
	rawtss[folderName] = processFile("avgraw_8000")
	scaledtss[folderName] = processFile("avgscaled_8000")

# go through each array
writeFile("rawgene_8000", rawgene)
writeFile("rawtss_8000", rawtss)
writeFile("rawtes_8000", rawtes)
writeFile("scaledgene_8000", scaledgene)
writeFile("scaledtss_8000", scaledtss)
writeFile("scaledtes_8000", scaledtes)
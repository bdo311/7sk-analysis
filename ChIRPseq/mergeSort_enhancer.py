# mergeSort_enhancer.py
# 4/8/14; based on mergeSort.py
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
# First merge, then sort, then split into 'no' and 'yes' --> yes is super enhancer
def mergeSort(folder):
	os.chdir(folder)
	os.system("cat chr* > allchr.txt")
	os.system("sort -rn -t $'\t' -k7,7 allchr.txt > allchr_sorted.txt")
	os.system("awk -F \"\t\" '{output=$6\".txt\"; print $0 > output}' allchr_sorted.txt")

for folder in folderToGraph:
	mergeSort(folder + 'enhancers/')

# Run R file to get plot and average bins
# ---------------------------------------

scriptFile = "/home/raflynn/7SK_ChIRPseq/visualizeBins_enhancer.r"
def folderWorker(start, end, folders):
	for i in range(start, end):
		if i >= len(folders): continue
		folder = folders[i]

		folderName = folderToGraph[folder][1]
		os.chdir(folder + 'enhancers/')
		os.system("Rscript " + scriptFile + " gene " + folderName)

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

# 4 hash tables	
raw_super = {}
raw_reg = {}

scaled_super = {}
scaled_reg = {}

# go through each folder and get a bin array
for folder in folderToGraph:
	folderName = folderToGraph[folder][1]

	os.chdir(folder + 'enhancers/')
	raw_super[folderName] = processFile("avgraw_YES")
	raw_reg[folderName] = processFile("avgraw_NO")
	scaled_super[folderName] = processFile("avgscaled_YES")
	scaled_reg[folderName] = processFile("avgscaled_NO")

# go through each array
writeFile("raw_super", raw_super)
writeFile("raw_reg", raw_reg)
writeFile("scaled_super", scaled_super)
writeFile("scaled_reg", scaled_reg)
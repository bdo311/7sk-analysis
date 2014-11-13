# mergeSort_human.py
# 5/23/14, based on mergeSort.py
# 1. Gets all folders and files
# 2. Merges and sorts chr files for all folders
# 3. Runs R on all of them
# 4. Merges genebody, tss, tes files to get one average file for all

import os, glob, csv, re, multiprocessing

csv.register_dialect("textdialect", delimiter='\t')

# Get all folders I'll be operating on. 
# -------------------------------------------------------------------

seven_Folders = ['H1', 'HeLa', 'H1_input', 'HeLa_input']
chirpDir = "/home/raflynn/ChIRPseq/"

# Get a map of bin folder to folder name
# ------------------------------------------

folderToGraph = {}

for folder in seven_Folders:
	binFolder = chirpDir + folder + '/bins/'
	folderToGraph[binFolder] = folder
	
# Merge and sort
# --------------
# First merge, then sort

def mergeSort(folder):
	os.chdir(folder)
	if len(glob.glob("chr*")) != 0: os.system("cat chr* > allchr.txt")
	os.system("sort -rn -t $'\t' -k7,7 allchr.txt > allchr_sorted.txt")

regionTypes = ['geneBody', 'tes', 'tss']
# for folder in folderToGraph:
	# for region in regionTypes:
		# mergeSort(folder + region + '/')

# Run R file to get plot and average bins
# ---------------------------------------

geneScript = "/home/raflynn/ChIRPseq/visualizeBins_new.r"
def folderWorker(start, end, folders):
	for i in range(start, end):
		if i >= len(folders): continue
		folder = folders[i]

		folderName = folderToGraph[folder]
		
		os.chdir(folder + 'geneBody/')
		os.system("Rscript " + geneScript + " gene " + folderName)
		os.chdir(folder + 'tss/')
		os.system("Rscript " + geneScript + " tss " + folderName)
		os.chdir(folder + 'tes/')
		os.system("Rscript " + geneScript + " tes " + folderName)
		
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
	ofile = open('/home/raflynn/ChIRPseq/averages/human/' + name + '.txt', 'w')
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

# go through each folder and get a bin array
for folder in folderToGraph:
	folderName = folderToGraph[folder]

	os.chdir(folder + 'geneBody/')
	rawgene[folderName] = processFile("avgraw")
	os.chdir(folder + 'tes/')
	rawtes[folderName] = processFile("avgraw")	
	os.chdir(folder + 'tss/')
	rawtss[folderName] = processFile("avgraw")

# go through each array
writeFile("rawgene", rawgene)
writeFile("rawtss", rawtss)
writeFile("rawtes", rawtes)
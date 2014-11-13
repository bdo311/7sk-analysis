# mergeSort_all.py
# 4/25/14, based on mergeSort and mergeSort_enhancer_new
# 1. Gets all folders and files
# 2. Merges and sorts chr files for all folders
# 3. Runs R on all of them
# 4. Merges genebody, tss, tes, extsuper, super, reg files to get one average file for all

import os, glob, csv, re, multiprocessing
csv.register_dialect("textdialect", delimiter='\t')

# 1. Get all folders I'll be operating on. 
# ---------------------------------------

folderToGraph = {}
# seven_Folders = ['HeLa']
# seven_Folders = ['50minFP_minus', '50minFP_plus', 'untreated_minus', 'untreated_plus']
seven_Folders = ['7SK_Input', 'ActD_new', 'DRB_new', 'Flavo_new', 'JQ11_new', 'JQ14_new', 'WT2_new']
# seven_Folders = ['12.5minFP_minus', '12.5minFP_plus', '25minFP_minus', '25minFP_plus', '50minFP_minus', '50minFP_plus', '5minFP_minus', '5minFP_plus', 'untreated_minus', 'untreated_plus']
chirpDir = "/home/raflynn/ChIRPseq/"
groDir = chirpDir + 'groseq/'

for folder in seven_Folders:
	binFolder = chirpDir + folder + '/bins/'
	folderToGraph[binFolder] = folder

# getting ChIP folders

# chipDir = "/home/raflynn/mmchipseq/mes_all/"
# chipFiles = glob.glob(chipDir + "ChIP_*.bedGraph")
# chipFolders = []
# for fileName in chipFiles:
	# exp = re.search("ChIP_(.*)_new.bedGraph", fileName).group(1)
	# chipFolders.append(exp)
	
# for folder in chipFolders:
	# binFolder = chipDir + folder + '/'
	# folderToGraph[binFolder] = folder
	
# Merge and sort
# --------------
# First merge, then sort

def flip(folder):
	os.chdir(folder)
	ifile = open("allchr.txt", 'r')
	reader = csv.reader(ifile, 'textdialect')
	ofile = open("allchr_flipped.txt", 'w')
	writer = csv.writer(ofile, 'textdialect')
	
	for row in reader:
		if row[5] == '+': writer.writerow(row)
		else:
			bins = row[7:][::-1]
			outputrow = row[:7]
			outputrow.extend(bins)
			writer.writerow(outputrow)
		
	ifile.close()
	ofile.close()

def mergeSort(folder):
	os.chdir(folder)
	if len(glob.glob("chr*")) != 0: os.system("cat chr* > allchr.txt")
	#flip(folder)
	os.system("sort -rn -t $'\t' -k7,7 allchr.txt > allchr_sorted.txt")

# allRegions = ['CTCF', 'Nanog', 'Sox2', 'REST', 'Oct4', 'Med1', 'p300', 'Stat3']
# allRegions = ['jb']
# allRegions = ['geneBody', 'tes', 'tss']
allRegions = ['top2a']
# enhRegions = ['extreg', 'extsuper', 'reg', 'super']
# for folder in folderToGraph:
# 	for region in allRegions: 
# 		mergeSort(folder + region + '/')
	# for region in enhRegions:
	# 	if 'H1' in folder or 'HeLa' in folder: continue #human data is not paired with enhancer data yet
	# 	mergeSort(folder + region + '/')

# Run R file to get plot and average bins
# ---------------------------------------
lincScript = "/home/raflynn/ChIRPseq/visualizeBins_new_mouse_linc.r"
mouseScript = "/home/raflynn/ChIRPseq/visualizeBins_new_mouse.r"
humanScript = "/home/raflynn/ChIRPseq/visualizeBins_new_human.r"
enhScript = "/home/raflynn/ChIRPseq/visualizeBins_enhancer.r"
def folderWorker(start, end, folders):
	for i in range(start, end):
		if i >= len(folders): continue
		folder = folders[i]

		folderName = folderToGraph[folder]
		isHuman = 'H1' in folderName or 'HeLa' in folderName
		geneScript = humanScript if isHuman else mouseScript

		# for region in allRegions:
		# 	os.chdir(folder + region + '/')
		# 	os.system("Rscript " + geneScript + " gene " + folderName + ' 6')	

		os.chdir(folder + 'top2a/')
		os.system("Rscript " + geneScript + " gene " + folderName + '7')
		# os.chdir(folder + 'jb/')
		# os.system("Rscript " + geneScript + " gene " + folderName + ' 6')
		# os.chdir(folder + 'lincRNA/')
		# os.system("Rscript " + lincScript + " tss " + folderName + ' 8')
		# os.chdir(folder + 'geneBody/')
		# os.system("Rscript " + geneScript + " gene " + folderName + ' 8')
		# os.chdir(folder + 'tss/')
		# os.system("Rscript " + geneScript + " tss " + folderName + ' 8')
		# os.chdir(folder + 'tes/')
		# os.system("Rscript " + geneScript + " tes " + folderName + ' 8')
		# os.chdir(folder + 'ctrpeaks/')
		# os.system("Rscript " + geneScript + " gene " + folderName + ' 8')
		# os.chdir(folder + 'extpeaks/')
		# os.system("Rscript " + geneScript + " gene " + folderName + ' 8')		

		# if not isHuman:
		# 	os.chdir(folder + 'extreg/')
		# 	os.system("Rscript " + enhScript + " super " + folderName + ' 8')
		# 	os.chdir(folder + 'extsuper/')
		# 	os.system("Rscript " + enhScript + " super " + folderName + ' 8')
		# 	os.chdir(folder + 'reg/')
		# 	os.system("Rscript " + enhScript + " super " + folderName + ' 6')
		# 	os.chdir(folder + 'super/')
		# 	os.system("Rscript " + enhScript + " super " + folderName + ' 6')

folders = folderToGraph.keys()
nproc = 8
numPerProc = len(folders)/8 + 1
procs = []
# for i in range(nproc):
# 	p = multiprocessing.Process(target=folderWorker, args=(i * numPerProc, (i + 1) * numPerProc, folders))
# 	print i * numPerProc, (i + 1) * numPerProc
# 	procs.append(p)
# 	p.start()

# for p in procs:
# 	p.join()

# merge all files, for both scaled and raw averages
# -------------------------------------------------
	
# take in an avg_*_* file and extract the bin values
def processFile(fileName):
	avgFile = glob.glob(fileName + "*.txt")[0]
	print avgFile
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
	ofile = open('/home/raflynn/ChIRPseq/averages/' + name + '.txt', 'w')
	writer = csv.writer(ofile, 'textdialect')
	for treatment in mapping:
		outputRow = [treatment]
		outputRow.extend(mapping[treatment])
		writer.writerow(outputRow)

	ofile.close()

# hash tables	
# linc_mouse = {}
# gene_mouse = {}
# tes_mouse = {}
# tss_mouse = {}
# extsuper_mouse = {}
# extreg_mouse = {}
# super_mouse = {}
# reg_mouse = {}
# extpeaks_mouse = {}
# ctrpeaks_mouse = {}
top2a_mouse = {}


# hashes = {}
# for region in allRegions:
# 	hashes[region] = {}

# gene_human = {}
# tes_human = {}
# tss_human = {}
# extpeaks_human = {}
# ctrpeaks_human = {}
#jb_human = {}

# go through each folder and get a bin array
for folder in folderToGraph:
	folderName = folderToGraph[folder]
	#if not glob.glob(folder + 'geneBody/'): continue
	print folderName
	isHuman = 'H1' in folderName or 'HeLa' in folderName
	# os.chdir(folder + 'lincRNA/')
	# linc_mouse[folderName] = processFile("avgraw")
	
	if isHuman:
		os.chdir(folder + 'geneBody/')
		gene_human[folderName] = processFile("avgraw_")
		# os.chdir(folder + 'tes/')
		# tes_human[folderName] = processFile("avgraw_tes_" + folderName)	
		# os.chdir(folder + 'tss/')
		# tss_human[folderName] = processFile("avgraw_tss_" + folderName)
		# os.chdir(folder + 'jb/')
		# jb_human[folderName] = processFile("avgraw_jb_" + folderName)
		# os.chdir(folder + 'extpeaks/')
		# extpeaks_human[folderName] = processFile("avgraw")
		# os.chdir(folder + 'ctrpeaks/')
		# ctrpeaks_human[folderName] = processFile("avgraw")
	else:
		# for region in allRegions:
		# 	os.chdir(folder + region + '/')
		# 	hashes[region][folderName] = processFile("avgraw_")
		os.chdir(folder + 'top2a/')
		top2a_mouse[folderName] = processFile("avgraw_gene_" + folderName)
		# os.chdir(folder + 'geneBody/')
		# gene_mouse[folderName] = processFile("avgraw_gene_" + folderName)
		# os.chdir(folder + 'tes/')
		# tes_mouse[folderName] = processFile("avgraw_tes_" + folderName)
		# os.chdir(folder + 'tss/')
		# tss_mouse[folderName] = processFile("avgraw_tss_" + folderName)
		# os.chdir(folder + 'extsuper/')
		# extsuper_mouse[folderName] = processFile("avgraw")	
		# os.chdir(folder + 'extreg/')
		# extreg_mouse[folderName] = processFile("avgraw")
		# os.chdir(folder + 'extpeaks/')
		# extpeaks_mouse[folderName] = processFile("avgraw")
		# os.chdir(folder + 'ctrpeaks/')
		# ctrpeaks_mouse[folderName] = processFile("avgraw")
		
# go through each array
# for region in allRegions:
# 	writeFile('groseq_' + region, hashes[region])

writeFile("top2a_mouse", top2a_mouse)
# writeFile("linc_mouse", linc_mouse)
# writeFile("gene_mouse", gene_mouse)
# writeFile("tes_mouse", tes_mouse)
# writeFile("tss_mouse", tss_mouse)
# writeFile("extsuper_mouse", extsuper_mouse)
# writeFile("extreg_mouse", extreg_mouse)
# writeFile("extpeaks_mouse", extpeaks_mouse)
# writeFile("ctrpeaks_mouse", ctrpeaks_mouse)

# writeFile("jb_human", jb_human)
# writeFile("gene_human", gene_human)
# writeFile("tes_human", tes_human)
# writeFile("tss_human", tss_human)
# writeFile("extpeaks_human", extpeaks_human)
# writeFile("ctrpeaks_human", ctrpeaks_human)
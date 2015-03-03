# geneToTFMatrix.py
# 4/27/14
# To get a 21000x156 matrix with average density values for all genes

import glob, collections, csv, re

csv.register_dialect("textdialect", delimiter='\t')

# 1. Get all folders I'll be operating on. 
# -------------------------------------------------------------------

seven_Folders = ['7SK_Input', 'ActD', 'Flavo', 'DRB', 'WT2', 'JQ11', 'JQ14']
chirpDir = "/home/raflynn/7SK_ChIRPseq/"

# getting ChIP folders

chipDir = "/home/raflynn/Mm_ChIPseq/mES_sra_data_new/"

chipFiles = glob.glob(chipDir + "ChIP_*.bedGraph")
chipFolders = []
for fileName in chipFiles:
	exp = re.search("ChIP_(.*)_new.bedGraph", fileName).group(1)
	chipFolders.append(exp)

# getting ChIP folders

youngDir = "/home/raflynn/Mm_ChIPseq/mES_young/"

youngFiles = glob.glob(youngDir + "ChIP_*.bedGraph")
youngFolders = []
for fileName in youngFiles:
	exp = re.search("ChIP_(.*)_new.bedGraph", fileName).group(1)
	youngFolders.append(exp)
	
# 2. Get a map of bin folder to folder name
# ------------------------------------------

folderToGraph = {}

for folder in seven_Folders:
	binFolder = chirpDir + folder + '/bins/'
	folderToGraph[binFolder] = folder

for folder in chipFolders:
	binFolder = chipDir + folder + '/'
	folderToGraph[binFolder] = folder
	
for folder in youngFolders:
	binFolder = youngDir + folder + '/'
	folderToGraph[binFolder] = folder
	
# 3. Go into each folder, process the allchr.txt files
# ----------------------------------------------------

# region --> gene --> TF --> average density
tss = collections.defaultdict(lambda: {})
tes = collections.defaultdict(lambda: {})
geneBody = collections.defaultdict(lambda: {})

def addToMap(regionMap, fn, tf, regionType):
	ifile = open(fn, 'r')
	reader = csv.reader(ifile, 'textdialect')
	
	numBins = len(reader.next()) - 7 if regionType == "geneBody" else 400
	ifile.seek(0)
	
	for row in reader:
		regionMap[row[0]][tf] = float(row[6])/numBins
		
	ifile.close()

for folder in folderToGraph:
	folderName = folderToGraph[folder]
	print folderName
	addToMap(tss, folder + "tss/allchr.txt", folderName, 'tss')
	addToMap(tes, folder + "tes/allchr.txt", folderName, 'tes')
	addToMap(geneBody, folder + "geneBody/allchr.txt", folderName, 'geneBody')
	
# 4. Print out information
# ------------------------

header = ['']
header.extend(folderToGraph.values())

tssFile = open("geneToTSS.txt", 'w')
tssWriter = csv.writer(tssFile, 'textdialect')
tssWriter.writerow(header)
for gene in tss:
	outputRow = [gene]
	outputRow.extend([tss[gene][tf] for tf in header[1:]])
	tssWriter.writerow(outputRow)
tssFile.close()

tesFile = open("geneToTES.txt", 'w')
tesWriter = csv.writer(tesFile, 'textdialect')
tesWriter.writerow(header)
for gene in tes:
	outputRow = [gene]
	outputRow.extend([tes[gene][tf] for tf in header[1:]])
	tesWriter.writerow(outputRow)
tesFile.close()

gbFile = open("geneToGB.txt", 'w')
gbWriter = csv.writer(gbFile, 'textdialect')
gbWriter.writerow(header)
for gene in geneBody:
	outputRow = [gene]
	outputRow.extend([geneBody[gene][tf] for tf in header[1:]])
	gbWriter.writerow(outputRow)
gbFile.close()



	
	

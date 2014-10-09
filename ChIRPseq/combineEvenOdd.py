# combineEvenOdd.py
# 3/26/14
# 1. Gets all folders and files in 7SK folder
# 2. Merges even and odd bedgraphs into combined.bedGraph, makes a bw, and uploads it

import os, glob, csv, re, collections, math, multiprocessing

csv.register_dialect("textdialect", delimiter='\t')
csv.register_dialect("spacedialect", delimiter=' ')

# Get all folders I'll be operating on. folders is all the 7SK ChIRPseq experiments.
# -------------------------------------------------------------------
folders = ['ActD', 'DRB', 'Flavo', 'JQ11', 'JQ14', 'WT1', 'WT2']
chirpDir = "/home/raflynn/7SK_ChIRPseq/"

def fromSpaceToTab(inName, outName):
	ifile = open(inName, 'r')
	reader = csv.reader(ifile, 'spacedialect')
	
	ofile = open(outName, 'w')
	writer = csv.writer(ofile, 'textdialect')
	
	for row in reader: writer.writerow(row)
	
	ifile.close()
	ofile.close()

def processFolder(folder):
	# make new subfolders
	# os.chdir(chirpDir + folder + '/')
	# os.system("mkdir bins")
	# os.chdir(chirpDir + folder + '/bins/')
	# os.system("mkdir geneBody tss tes")
	os.chdir(chirpDir + folder)
	
	# print folder, "write even and odd bedgraphs"
	# fromSpaceToTab("even.bedGraph", "even_tabbed.bedGraph")
	# fromSpaceToTab("odd.bedGraph", "odd_tabbed.bedGraph")
	
	# print folder, "make bedgraph with two cols from the even and odd"
	# os.system("unionBedGraphs -i even_tabbed.bedGraph odd_tabbed.bedGraph -header -names Even Odd > combined_unnorm.bedGraph")
		
	# print folder, "combine it with ratio 5"
	# script = chirpDir + "mergeCombined.py"
	# infile = chirpDir + folder + "/combined_unnorm.bedGraph"
	# outfile = chirpDir + folder + "/combined.bedGraph"
	
	# os.system("python  " + script + ' ' + infile + ' ' + outfile)
	
	# ---------------------------------------------------------------------------
	# IN THE FUTURE: NEED TO ADD NORMALIZEBEDGRAPH.PY TO NORM TO 10 MILLION READS
	# ---------------------------------------------------------------------------
	
	print folder, "make a bw file"
	sizeFile = chirpDir + "genes/mm9.sizes"
	os.system("bedGraphToBigWig combined.bedGraph " + sizeFile + " combined.bw")
	
	print folder, "upload to amazon aws"
	outputName = folder + '_combined.bw'
	os.system("../aws put \"x-amz-acl: public-read\" changseq/bdo/chirp/" + outputName + " combined.bw")
	
# Multiprocessing
# ---------------
for folder in folders:
	p = multiprocessing.Process(target=processFolder, args=(folder,))
	p.start()




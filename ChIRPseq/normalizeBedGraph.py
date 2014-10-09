# normalizeBedGraph.py
# 4/8/14, updated 5/22/14 to work on human cells
# Normalizes all ChIRP bedgraphs to a total read density of 750 million, then makes a bw file

import os, glob, csv, re, multiprocessing

csv.register_dialect("textdialect", delimiter='\t')

# Get all folders I'll be operating on. seven_Folders is all the 7SK ChIRPseq experiments. 
# -------------------------------------------------------------------

#seven_Folders = ['ActD_new', 'DRB_new', 'Flavo_new', 'JQ11_new', 'JQ14_new', 'H1', 'HeLa', 'WT2_new']
seven_Folders = ['H1', 'HeLa']
chirpDir = "/home/raflynn/ChIRPseq/"

# normalize by scaling everything to a total read density of 750 million
def normalize(folder):
	os.chdir(chirpDir + folder)
	
	ifile = open("combined_no7SK.bedGraph", 'r')
	reader = csv.reader(ifile, 'textdialect')
	
	sum = 0
	for row in reader:
		if '_' in row[0]: continue
		sum += (int(row[2]) - int(row[1])) * float(row[3])
	
	multFactor = 750000000/sum
	
	ifile.seek(0)
	ofile = open("combined_no7SK_norm.bedGraph", 'w')
	writer = csv.writer(ofile, 'textdialect')
	
	for row in reader:
		if '_' in row[0]: continue
		row[-1] = float(row[-1]) * multFactor
		writer.writerow(row)
		
	ifile.close()
	ofile.close()
	
	# ifile = open("input_tabbed.bedGraph", 'r')
	# reader = csv.reader(ifile, 'textdialect')
	
	# sum = 0
	# for row in reader:
		# sum += (int(row[2]) - int(row[1])) * float(row[3])
	
	# multFactor = 750000000/sum
	
	# ifile.seek(0)
	# ofile = open("input_norm.bedGraph", 'w')
	# writer = csv.writer(ofile, 'textdialect')
	
	# for row in reader:
		# row[-1] = float(row[-1]) * multFactor
		# writer.writerow(row)
		
	# ifile.close()
	# ofile.close()

# normalize, make bw, and upload
def doPipeline(folder):
	normalize(folder)
	os.chdir(chirpDir + folder)	
	# make a bw file
	if 'H' in folder: sizeFile = chirpDir + "genes/hg19.sizes"
	else: sizeFile = chirpDir + "genes/mm9.sizes"
	os.system("bedGraphToBigWig combined_no7SK_norm.bedGraph " + sizeFile + " combined_no7SK_norm.bw")
	#os.system("bedGraphToBigWig input_norm.bedGraph " + sizeFile + " input_norm.bw")
	# upload to amazon aws
	outputName = folder + '_combined_no7SK_norm.bw'
	os.system("../aws put \"x-amz-acl: public-read\" changseq/bdo/chirp/" + outputName + " combined_no7SK_norm.bw")

# multiprocessing
nproc = len(seven_Folders)
procs = []
for i in range(nproc):
	p = multiprocessing.Process(target=doPipeline, args=(seven_Folders[i],))
	procs.append(p)
	p.start()

for p in procs:
	p.join()

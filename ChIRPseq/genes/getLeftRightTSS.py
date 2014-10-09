# getLeftRightTSS.py
# 6/8/14
# get bed files for the left, center, and right TSS's to sue for ATAC-seq

import csv
csv.register_dialect("textdialect", delimiter='\t')

# 1. get left, right, and center genes
def getGeneList(fn):
	ifile = open(fn, 'r')
	reader = csv.reader(ifile, 'textdialect')
	
	genes = set()
	for row in reader:
		genes.add(row[0])
	
	ifile.close()
	return genes

dir = "/home/raflynn/ChIRPseq/WT2_new/bins/tss/clustering/"
left = getGeneList(dir + "left.txt")
center = getGeneList(dir + "center.txt")
right = getGeneList(dir + "right.txt")

# 2. go through the entire TSS file and make separate lists
def writeToFile(row, writer):
	chr = row[1]
	if '_' in chr: return
	name = row[0]
	strand = row[2]
	if strand == '+': start, end = row[3], row[4]
	else: start, end = row[5], row[6]
	
	writer.writerow([chr, start, end, strand, name])

ifile = open("mm9_refseq_startend.txt", 'r')
reader = csv.reader(ifile, 'textdialect')

leftFile = open("mm9_left_start.txt", 'w')
leftWriter = csv.writer(leftFile, 'textdialect')
centerFile = open("mm9_center_start.txt", 'w')
centerWriter = csv.writer(centerFile, 'textdialect')
rightFile = open("mm9_right_start.txt", 'w')
rightWriter = csv.writer(rightFile, 'textdialect')

reader.next()
for row in reader:
	if row[0] in left: writeToFile(row, leftWriter)
	elif row[0] in center: writeToFile(row, centerWriter)
	elif row[0] in right: writeToFile(row, rightWriter)

ifile.close()
leftFile.close()
centerFile.close()
rightFile.close()
		
	
# File: getRegEnhancerCoords.py
# 4/9/14
# Getting +- 2kb for regular enhancers, and getting +- 10% for super enhancers

import csv, os

csv.register_dialect("textdialect", delimiter='\t')
os.chdir("/home/raflynn/7SK_ChIRPseq/genes/")

# REGULAR ENHANCERS
# -----------------
ifile = open("mES_reg_enhancers.txt", 'r')
reader = csv.reader(ifile, 'textdialect')
ofile = open("mES_reg_enhancers_centered.txt", 'w')
writer = csv.writer(ofile, 'textdialect')

for row in reader:
	start, end = int(row[2]), int(row[3])
	center = (end + start)/2
	
	newStart = center - 2000
	newEnd = center + 2000
	
	newRow = row[:2]
	newRow.extend([newStart, newEnd])
	newRow.extend(row[4:6])
	writer.writerow(newRow)
	
ifile.close()
ofile.close()

# SUPER ENHANCERS
# ---------------
ifile = open("mES_super_enhancers.txt", 'r')
reader = csv.reader(ifile, 'textdialect')
ofile = open("mES_super_enhancers_expanded5.txt", 'w')
writer = csv.writer(ofile, 'textdialect')

for row in reader:
	start, end = int(row[2]), int(row[3])
	length = end - start 
	
	newStart = start - (length*2)
	newEnd = end + (length*2)
	
	newRow = row[:2]
	newRow.extend([newStart, newEnd])
	newRow.extend(row[4:6])
	writer.writerow(newRow)
	
ifile.close()
ofile.close()
	

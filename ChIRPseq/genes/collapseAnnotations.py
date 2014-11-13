# collapseAnnotations.py
# 3/16/14
# Collapses mm9.basic.annotation so that it can be used to interrogate reads per type of region

import csv

csv.register_dialect("textdialect", delimiter='\t')

ifile = open("hg19.basic.annotation", 'r')
reader = csv.reader(ifile, 'textdialect')

ofile = open("hg19_newannotation.txt", 'w')
writer = csv.writer(ofile, 'textdialect')

carryOver = False
for row in reader:
	len = int(row[3]) - int(row[2])
	if len == 1: carryOver = True #add 1 to the next one
	else:
		newRow = row[:2]
		if carryOver: newRow.append(int(row[2]) - 1) # extend the row by 1 on the left side
		else: newRow.append(int(row[2]))
		newRow.extend(row[3:])
		writer.writerow(newRow)
		
ifile.close()
ofile.close()

		
		


# mergeCombined.py
# 3/22/14, edited 3/26/14
# Looks at the two-column bedgraph (even and odd) and combines according to this:
# 	If the ratio of even to odd is >x or <1/x, replace value with the smaller of the two values.
#	If either of even or odd is 0, don't write the row
#	Else, average the two.

import csv, collections, math, numpy, sys

csv.register_dialect("textdialect", delimiter='\t')

# from command line
inFileName = sys.argv[1] #"/home/raflynn/7SK_ChIRPseq/WT2/combined_unnorm.bedGraph"
outFileName = sys.argv[2] #"/home/raflynn/7SK_ChIRPseq/WT2/combined2.bedGraph"

# ratio threshold
x = 5

# input file
ifile = open(inFileName, 'r')
reader = csv.reader(ifile, 'textdialect')

# output file
ofile = open(outFileName, 'w')
writer = csv.writer(ofile, 'textdialect')

reader.next()
counter = 0

currRow = ['', 0, 0, 0]
isFirstRow = True
for row in reader:
	counter += 1
	if not counter % 1000000: print counter
	
	even, odd = float(row[3]), float(row[4])
	
	# calculating new score
	if even == 0 or odd == 0: continue
	else:
		ratio = math.floor(even/odd) if even > odd else math.floor(odd/even)
		newscore = numpy.mean([even, odd]) if ratio < x else numpy.amin([even, odd])
		
	# writing new line
	outputRow = row[:3]
	outputRow.append(newscore)
	
	if outputRow[1] == currRow[2] and outputRow[3] == currRow[3]: 
		currRow = [currRow[0], currRow[1], outputRow[2], outputRow[3]]
	else:
		if not isFirstRow: writer.writerow(currRow)
		else: isFirstRow = False
		currRow = outputRow	
	
ifile.close()

	



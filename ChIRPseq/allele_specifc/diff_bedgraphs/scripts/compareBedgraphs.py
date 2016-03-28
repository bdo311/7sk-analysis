# compareBedgraphs.py
# 8/5/14
# compares 2 allele-specific bedgraphs; gets difference and ratio tracks

import csv, sys, math
csv.register_dialect("textdialect", delimiter='\t')
csv.register_dialect("spacedialect", delimiter=' ')

ifn = sys.argv[1]
diff_fn = sys.argv[2]
ratio_fn = sys.argv[3]

# normalize by scaling everything to a total read density of 500 million
ifile = open(ifn, 'r')
reader = csv.reader(ifile, 'textdialect')

diff_file = open(diff_fn, 'w')
diff_writer = csv.writer(diff_file, 'textdialect')

ratio_file = open(ratio_fn, 'w')
ratio_writer = csv.writer(ratio_file, 'textdialect')

pastRow = ['0', 0, 0, 0, 0]
for row in reader:
	g1 = float(row[3])
	g2 = float(row[4])
	
	diff = g1 - g2
	ratio = g1/(g1 + g2) - 0.5
	
	if diff == 0 or ratio == 0: continue
	start, stop = int(row[1]), int(row[2])
	
	if start == pastRow[2] and diff == pastRow[3] and ratio == pastRow[4]: 
		pastRow[2] = stop
	else:
		if pastRow[0] != '0': 
			diff_writer.writerow(pastRow[:4])  # write diff
			pastRow[3] = pastRow[4]
			ratio_writer.writerow(pastRow[:4])  # write ratio
		pastRow = [row[0], start, stop, diff, ratio]  # new row
	
	
if pastRow[0] != '0': 
	diff_writer.writerow(pastRow[:4])  # write diff
	pastRow[3] = pastRow[4]
	ratio_writer.writerow(pastRow[:4])  # write ratio
ifile.close()
diff_file.close()
ratio_file.close()


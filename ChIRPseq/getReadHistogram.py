# getReadHistogram.py
# 3/22/14
# Histogram of reads per base to total bases in those bins

import csv, collections
csv.register_dialect("spacedialect", delimiter=' ')
csv.register_dialect("textdialect", delimiter='\t')

bins = collections.defaultdict(lambda: 0)

# old file
ifile = open("/home/raflynn/7SK_ChIRPseq/WT2/odd.bedGraph", 'r')
reader = csv.reader(ifile, 'spacedialect')

# new file
replfile = open("/home/raflynn/7SK_ChIRPseq/WT2/odd_tabbed.bedGraph", 'w')
writer = csv.writer(replfile, 'textdialect')

counter = 0
for row in reader:
	if counter % 100000 == 0: print counter
	bins[round(float(row[3]), 1)] += (int(row[2]) - int(row[1])) #round (x, 1) rounds to nearest tenth
	counter += 1
	writer.writerow(row) #write bedgraph with tabs

ifile.close()
replfile.close()

# writing histogram

ofile = open("/home/raflynn/7SK_ChIRPseq/WT2/odd_hist.txt", 'w')
writer = csv.writer(ofile, 'textdialect')

for bin in bins:
	row = [bin, bins[bin]]
	writer.writerow(row)
	
ofile.close()
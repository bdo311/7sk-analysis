# combinedStats.py
# 3/22/14
# Looks at the two-column bedgraph (even and odd) and gets a histogram of ratios 

import csv, collections, math

csv.register_dialect("textdialect", delimiter='\t')

ifile = open("/home/raflynn/7SK_ChIRPseq/WT2/combined_unnorm.bedGraph", 'r')
reader = csv.reader(ifile, 'textdialect')

histo = collections.defaultdict(lambda: 0)

reader.next()
counter = 0
for row in reader:
	counter += 1
	if not counter % 1000000: print counter
	
	dist = int(row[2]) - int(row[1]) 
	
	even, odd = float(row[3]), float(row[4])
	
	# infinity
	if even == 0 or odd == 0:
		histo[0] += dist
		continue
	
	# binning
	ratio = math.floor(even/odd) if even > odd else math.floor(odd/even)
	histo[ratio] += dist
	
ifile.close()

# print out file
ofile = open("combined_hist.txt", 'w')
writer = csv.writer(ofile, 'textdialect')
for ratio in histo:
	writer.writerow([ratio, histo[ratio]])
	
ofile.close()

# # only for 0-10
# newhisto = collections.defaultdict(lambda: 0)
# for ratio in histo:
	# if ratio >= 10: newhisto[10] += histo[ratio]
	# else: newhisto[ratio] = histo[ratio]
	
# # print out file
# ofile = open("combined_hist.txt", 'w')
# writer = csv.writer(ofile, 'textdialect')
# for ratio in newhisto:
	# writer.writerow([ratio, newhisto[ratio]])
	
# ofile.close()

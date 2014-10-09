import csv
import collections, os, sys

csv.register_dialect("textdialect", delimiter='\t')

dir = "/home/raflynn/ChIRPseq/" + sys.argv[1] + "/bins/newannot/"
os.chdir(dir)
os.system("cat chr* > allchr.txt")
ifile = open(dir + "allchr.txt", 'r')
reader = csv.reader(ifile, 'textdialect')

regions = collections.defaultdict(lambda: [0,0]) #length, score

for row in reader:
	if int(row[3]) > 1e9: continue
	regions[row[5]][0] += (int(row[3]) - int(row[2]) + 1)
	regions[row[5]][1] += float(row[-1])

ifile.close()
	
ofile = open(dir + "allchr_new_collapsed.txt", 'w')
writer = csv.writer(ofile, 'textdialect')

for region in regions:
	row = [region]
	row.extend(regions[region])
	writer.writerow(row)
	
ofile.close()

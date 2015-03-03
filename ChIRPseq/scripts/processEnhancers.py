# processEnhancers.py
# 7/7/14
# gets super enhancer and enhancer count for WT2

import os, glob, csv, re, multiprocessing
csv.register_dialect("textdialect", delimiter='\t')

os.chdir("/home/raflynn/ChIRPseq/WT2_new/bins/enhancers/")
os.system("cat chr* > allchr.txt")
os.system("sort -rn -t $'\t' -k7,7 allchr.txt > allchr_sorted.txt")

regEnh = [0,0]
superEnh = [0,0]

ifile = open("allchr_sorted.txt", 'r')
reader = csv.reader(ifile, 'textdialect')

for row in reader:
	length = int(row[3]) - int(row[2]) + 1
	reads = float(row[-1])

	if row[-2] == 'YES': 
		superEnh[0] += length
		superEnh[1] += reads
	else:
		regEnh[0] += length
		regEnh[1] += reads

print 'reg', regEnh
print 'super', superEnh

ifile.close()
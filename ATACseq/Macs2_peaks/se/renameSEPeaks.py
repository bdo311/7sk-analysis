# renameSEPeaks.py
# 2/4/14

import csv, collections
csv.register_dialect("textdialect", delimiter='\t')

ifile = open("atac_inSE_noTSS.bed",'r')
reader = csv.reader(ifile, 'textdialect')

ofile = open("atac_inSE_noTSS_BED6_1kb.bed",'w')
writer = csv.writer(ofile, 'textdialect')

seList = collections.defaultdict(lambda: 0)
for row in reader:
	se = '__'.join([row[6],row[4],row[5]])
	seList[se] += 1
	
	se = se + '__' + str(seList[se])

	ctr = (int(row[1])+int(row[2]))/2
	row[1]=ctr-1000
	row[2]=ctr+1000
	writer.writerow([row[0],row[1],row[2],se,0,'+'])
ifile.close()
ofile.close()

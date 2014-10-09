#getUpstreamBed.py
#5/8/14
#gets 2000nt upstream of every gene so that we can get U1/PAS

import csv
csv.register_dialect("textdialect", delimiter='\t')

ifile = open("mm9_refseq_startend.txt", 'r')
reader = csv.reader(ifile, 'textdialect')

ofile = open("mm9_refseq_coll_upstream.bed", 'w')
writer = csv.writer(ofile, 'textdialect')

reader.next()
for row in reader:
	chr = row[1]
	if 'random' in chr: continue
	
	if row[2] == '+':
		upStart = int(row[3]) - 1001
		upEnd = upStart + 2000
	else:
		upStart = int(row[5]) + 1001
		upEnd = upStart + 2000
		
	writer.writerow([chr, upStart, upEnd, row[0], 0, row[2]])

ifile.close()
ofile.close()
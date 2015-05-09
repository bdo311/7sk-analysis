# getSingleNtSignal.py
# from the START-seq bedgraphs that have been intersected with enhancers, take only the nucleotide that has highest signal for that particular enhancer
# later on, the reverse signal has to be at a lower nucleotide than the forward signal, but we can take care of that later??

import csv, collections, sys
csv.register_dialect("textdialect",delimiter='\t')

def getSingleNT(ifile, ofile):
	reader = csv.reader(ifile, 'textdialect')
	writer = csv.writer(ofile, 'textdialect')

	row = reader.next()
	enhancer = '__'.join(row[4:])
	currEnh = enhancer
	currBest = row
	for row in reader:
		# print row
		enhancer = '__'.join(row[4:])
		if enhancer != currEnh: # new enhancer 
			writer.writerow(currBest)
			currBest = row
			
		currEnh = enhancer
		if float(row[3]) > float(currBest[3]):
			currBest = row

	writer.writerow(currBest)
	ifile.close()
	ofile.close()
	
# getting single nt files
# ifile = open("startRNA_NELFwt_forward_pro_RYenh.bedgraph",'r')
# ffile = open("startRNA_NELFwt_forward_pro_RYenh_greatest.bedgraph",'w')
# getSingleNT(ifile, ffile)
# ifile = open("startRNA_NELFwt_reverse_pro_RYenh.bedgraph",'r')
# rfile = open("startRNA_NELFwt_reverse_pro_RYenh_greatest.bedgraph",'w')
# getSingleNT(ifile, rfile)

# combining them
def getEnhSignal(fn):
	ifile = open(fn,'r')
	reader = csv.reader(ifile, 'textdialect')
	enh = {}
	for row in reader:
		enhName = '__'.join(row[4:])
		row[1:4] = [float(x) for x in row[1:4]]
		enh[enhName] = row[:4]
	return enh

fenh = getEnhSignal("startRNA_NELFwt_forward_pro_RYenh_greatest.bedgraph")
renh = getEnhSignal("startRNA_NELFwt_reverse_pro_RYenh_greatest.bedgraph")

#print len(set(fenh.keys()) & set(renh.keys()))
counter = 0
for enh in fenh:
	if enh not in renh: continue
	if renh[enh][1] < fenh[enh][1]: counter += 1
	else: print renh[enh], fenh[enh]
print counter



	
	

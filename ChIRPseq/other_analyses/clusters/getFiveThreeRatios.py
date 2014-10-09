# getFiveThreeRatios.py
# 5/15/14
# get the 5' and 3' U1/PAS ratios

import csv

csv.register_dialect("textdialect", delimiter='\t')

ifile = open("foundU1_gene_new.txt", 'r')
reader = csv.reader(ifile, 'textdialect')

endFile = open("foundU1_gene_ends_50.txt", 'w')
endWriter = csv.writer(endFile, 'textdialect')

def getFiveThree(sites, length):
	five = length/2
	three = length - five
	
	fiveSites = []
	threeSites = []
	for site in sites:
		if site < five: fiveSites.append(site)
		elif site > three: threeSites.append(site)
		
	return fiveSites, threeSites
	
header = ["gene", "length", "5_u1", "5_u1rate", "5_pas", "5_pasrate", "5_ratio", "3_u1", "3_u1rate", "3_pas", "3_pasrate", "3_ratio"]
endWriter.writerow(header)
for row in reader:
	geneName = row[0]
	length = int(row[1])
	if row[6] == '': u1_sites = []
	else: u1_sites = [int(x) for x in row[6].split(',')]
	
	if row[7] == '': pas_sites = []
	else: pas_sites = [int(x) for x in row[7].split(',')]
	
	u1_sites_5, u1_sites_3 = getFiveThree(u1_sites, length)
	pas_sites_5, pas_sites_3 = getFiveThree(pas_sites, length)
	
	outputRow = row[:2]
	
	outputRow.extend([len(u1_sites_5), float(len(u1_sites_5))/(float(length)/4000)])
	num5PAS = len(pas_sites_5) if len(pas_sites_5) > 0 else 1
	outputRow.extend([len(pas_sites_5), float(len(pas_sites_5))/(float(length)/4000)])	
	outputRow.append(float(len(u1_sites_5))/num5PAS)
	
	outputRow.extend([len(u1_sites_3), float(len(u1_sites_3))/(float(length)/4000)])
	num3PAS = len(pas_sites_3) if len(pas_sites_3) > 0 else 1
	outputRow.extend([len(pas_sites_3), float(len(pas_sites_3))/(float(length)/4000)])
	outputRow.append(float(len(u1_sites_3))/num3PAS)
	
	endWriter.writerow(outputRow)
	
	
ifile.close()
endFile.close()

	
	
	
	
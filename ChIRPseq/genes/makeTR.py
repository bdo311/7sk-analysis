# makeEnhTR.py
# make enhancer distal and proximal TR files for start-seq centered enhancers
# edited 5/6: SE will be centered at atacseq peak

import csv
csv.register_dialect("textdialect",delimiter='\t')

ifile = open("mES_SE_individual_1kb_BED6.bed",'r')
reader = csv.reader(ifile, 'textdialect')
pfile = open("mES_SE_individual_atacseq_centered_TR_prox30.bed",'w')
pwriter = csv.writer(pfile, 'textdialect')
dfile = open("mES_SE_individual_atacseq_centered_TR_dist30.bed",'w')
dwriter = csv.writer(dfile, 'textdialect')
for row in reader:
	center = int(row[1]) + 1000
	plus_center = center
	minus_center = center #SE center at ATACseq peak
	
	pwriter.writerow([row[0],plus_center-30,plus_center+300,row[3]+"_plus",row[4],'+'])
	dwriter.writerow([row[0],plus_center+300,plus_center+1300,row[3]+"_plus",row[4],'+'])
	
	pwriter.writerow([row[0],minus_center-300,minus_center+30,row[3]+"_minus",row[4],'-'])
	dwriter.writerow([row[0],minus_center-1300,minus_center-300,row[3]+"_minus",row[4],'-'])

ifile.close()
pfile.close()
dfile.close()

ifile = open("mES_reg_enhancers_RY_startseq_centered_1kb.bed",'r')
reader = csv.reader(ifile, 'textdialect')
pfile = open("mES_reg_enhancers_RY_startseq_centered_TR_prox30.bed",'w')
pwriter = csv.writer(pfile, 'textdialect')
dfile = open("mES_reg_enhancers_RY_startseq_centered_TR_dist30.bed",'w')
dwriter = csv.writer(dfile, 'textdialect')
for row in reader:
	center = int(row[1]) + 1000
	dist = int(row[4])/2
	plus_center = center + dist
	minus_center = center - dist
	
	pwriter.writerow([row[0],plus_center-30,plus_center+300,row[3]+"_plus",row[4],'+'])
	dwriter.writerow([row[0],plus_center+300,plus_center+1300,row[3]+"_plus",row[4],'+'])
	
	pwriter.writerow([row[0],minus_center-300,minus_center+30,row[3]+"_minus",row[4],'-'])
	dwriter.writerow([row[0],minus_center-1300,minus_center-300,row[3]+"_minus",row[4],'-'])


ifile.close()
pfile.close()
dfile.close()


# get gene ends
ifile = open("mm9_refseq_coll.txt",'r')
reader = csv.reader(ifile, 'textdialect')

geneToEnd = {}
reader.next()
for row in reader:
	if row[2]=="+": geneToEnd[row[0]] = int(row[4])
	else: geneToEnd[row[0]] = int(row[3])
	
ifile.close()

# get gene starts and write file
ifile = open("mm9_tss_startseq_centered_1kb.bed",'r')
reader = csv.reader(ifile, 'textdialect')
pfile = open("mm9_tss_startseq_centered_TR_prox30.bed",'w')
pwriter = csv.writer(pfile, 'textdialect')
dfile = open("mm9_tss_startseq_centered_TR_dist30.bed",'w')
dwriter = csv.writer(dfile, 'textdialect')

for row in reader:
	start = int(row[1]) + 1000
	end = geneToEnd[row[3]]
	
	if row[5]=="+":
		pwriter.writerow([row[0],start-30,start+300,row[3],row[4],'+'])
		dwriter.writerow([row[0],start+300,end,row[3],row[4],'+'])
	else:
		pwriter.writerow([row[0],start-300,start+30,row[3],row[4],'-'])
		dwriter.writerow([row[0],end,start-300,row[3],row[4],'-'])
	
ifile.close()
pfile.close()
dfile.close()

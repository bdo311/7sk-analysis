# trim.py
# for trimming pro-seq files
# 5/31/15

import sys, os

if len(sys.argv) != 4: 
	print "Usage: python trim.py <name> <R1.fastq> <R2.fastq>"
	exit
	
name = sys.argv[1]
r1 = sys.argv[2]
r2 = sys.argv[3]

# clip 3' adapter and quality filter, then collapse and trim first 4 nt
r1_coll = r1[:-len(".fastq")]+"_coll.fastq"
cmd = "cat %s | fastx_clipper -n -l33 -Q33 -a AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG | fastq_quality_filter -Q33 -q25 -p80 | python collapse.py > %s"%(r1,r1_coll)
print cmd
os.system(cmd)

# split into files by barcode, and trim first 9 nt
os.system("nohup gzip %s &"%(r1,))

def writeRead(read):
	barcode = read[1][:4]
	if barcode not in barcodeToOFile: return
	barcodeToOFile[barcode].write(read[0])
	barcodeToOFile[barcode].write(read[1][9:])
	barcodeToOFile[barcode].write(read[2])
	barcodeToOFile[barcode].write(read[3][9:])
	
allowedBarcodes = ['CGGA','CCGG','AATG','TTAA','GGTT','TTGT','TGGC','AACC']

barcodeToOFile = {}
for bc in allowedBarcodes:
	barcodeToOFile[bc] = open(r1[:-len(".fastq")]+"_%s.fastq"%(bc,),'w')
	
with open(r1_coll,'r') as ifile:
	read = []
	for line in ifile:
		if len(read)==4:
			writeRead(read)
			read = []
		read.append(line)
	
	writeRead(read) #fencepost
		
for barcode in barcodeToOFile:
	barcodeToOFile[barcode].close()
	
#match with R2
#os.system("nohup gzip %s &"%(r1_coll,))

def writeRead(name, r2file, ofile):
	read = []
	for line in r2file: #iterate over r2file, starting where it last left off
		read.append(line)
		if len(read)==4: #only considering each read after all its info has been read
			r2name = read[0].split(' ')[0]
			if r2name == name: #if read names match
				ofile.write(read[0])
				ofile.write(read[1][:20]+'\n') #only take first 20nt
				ofile.write(read[2])
				ofile.write(read[3][:20]+'\n')
				return
			read = [] #reset
				
for bc in allowedBarcodes:
	with open(r1[:-len(".fastq")]+"_%s.fastq"%(bc,),'r') as ifile, \
	open(r2,'r') as r2file, \
	open(r2[:-len(".fastq")]+"_%s.fastq"%(bc,),'w') as ofile:
		counter = 0
		for line in ifile: #each read in R1 to match
			counter += 1
			if counter % 4 == 1:
				name = line.split(' ')[0]
				writeRead(name, r2file, ofile)
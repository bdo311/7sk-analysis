# trim.py
# for trimming FANTOM5 CAGE files
# 6/9/15

import sys, os

if len(sys.argv) != 4: 
	print "Usage: python trim.py <name> <R1.fastq>"
	exit
	
name = sys.argv[1]
r1 = sys.argv[2]

# clip 3' adapter and quality filter, then collapse
r1_coll = r1[:-len(".fastq")]+"_coll.fastq"
cmd = "cat %s | fastx_clipper -n -l29 -Q33 -a TCGTATGCCGTCTTCT | fastq_quality_filter -q25 -p80 | python collapse.py > %s"%(r1,r1_coll) #no Q33 because the quality files don't require it
print cmd
os.system(cmd)

# split into files by 3nt barcode, and trim first 9 nt
#os.system("nohup gzip %s &"%(r1,))

def writeRead(read):
	barcode = read[1][:3]
	if barcode not in barcodeToOFile: return
	barcodeToOFile[barcode].write(read[0])
	barcodeToOFile[barcode].write(read[1][9:])
	barcodeToOFile[barcode].write(read[2])
	barcodeToOFile[barcode].write(read[3][9:])
	
allowedBarcodes = ['CTT','AGC']

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
	
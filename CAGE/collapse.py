# collapse.py
# Usage: cat in.fastq | python collapse.py > out.fastq
# Use ONLY with trim.py for FANTOM5 CAGE
# 6/10/15

import fileinput, sys

allowedBarcodes = ['CTT','AGC']

alreadySeen = set()
read = []
for line in fileinput.input():
	if len(read)==4:
		if read[1][:3] in allowedBarcodes and read[1] not in alreadySeen: 
			# only take reads with desired barcodes
			# take only the first instance of every unique sequence
			sys.stdout.write(read[0])
			sys.stdout.write(read[1])
			sys.stdout.write(read[2])
			sys.stdout.write(read[3])
			alreadySeen.add(read[1])
		read = []
	read.append(line)
	
#fencepost
if read[1][:3] in allowedBarcodes:
	if read[1] not in alreadySeen: 	# take only the first instance of every unique sequence
		sys.stdout.write(read[0])
		sys.stdout.write(read[1])
		sys.stdout.write(read[2])
		sys.stdout.write(read[3])



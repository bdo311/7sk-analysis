# collapse.py
# Usage: cat in.fastq | python collapse.py > out.fastq
# Use ONLY with trim.py for PRO-seq
# 5/31/15

import fileinput, sys

alreadySeen = set()
read = []
for line in fileinput.input():
	if len(read)==4:
		if read[1] not in alreadySeen: 	# take only the first instance of every unique sequence
			sys.stdout.write(read[0])
			sys.stdout.write(read[1][4:])
			sys.stdout.write(read[2])
			sys.stdout.write(read[3][4:])
			alreadySeen.add(read[1])
		read = []
	read.append(line)
	
#fencepost
if read[1] not in alreadySeen: 	# take only the first instance of every unique sequence
	sys.stdout.write(read[0])
	sys.stdout.write(read[1][4:])
	sys.stdout.write(read[2])
	sys.stdout.write(read[3][4:])



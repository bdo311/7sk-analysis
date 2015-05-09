# trimcapseq.py
# full adapter = TCGTATGCCGTCTTCTGCTTG

import fileinput, sys
counter = 0
trim = 100
for str in fileinput.input():
	counter += 1
	if counter % 4 == 2: 
		trim = str.find("TCGTATGCC")
		if trim == -1:
			if str[-8:]=="TCGTATGC": trim = 92
			elif str[-7:]=="TCGTATG": trim = 93
			elif str[-6:]=="TCGTAT": trim = 94
			elif str[-5:]=="TCGTA": trim = 95	
			elif str[-4:]=="TCGT": trim = 96	
			elif str[-3:]=="TCG": trim = 97
			elif str[-2:]=="TC": trim = 98
			else: trim = 100			
		print str[4:trim]
	elif counter % 4 == 0:
		print str[4:trim]
	else: 
		sys.stdout.write(str)

ifile.close()
	
		
			


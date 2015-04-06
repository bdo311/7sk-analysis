# countChanges_nooutput.py
import csv, collections, math, sys
csv.register_dialect("textdialect",delimiter='\t')

ifile = open(sys.argv[1],'r') # 125p 125n 12Cp 12Cn
reader = csv.reader(ifile, 'textdialect')

groGeomZero = collections.defaultdict(lambda: 0)
groGeomConv = collections.defaultdict(lambda: 0)

counter = 0
for row in reader:
	counter += 1
	if counter % 1000000 == 0: print counter
	
	p5 = float(row[3])
	n5 = float(row[4])
	pc = float(row[5])
	nc = float(row[6])
	
	# if pc = nc = 0:
		# # gro00
	# elif pc = 0:
		# # gro05
	# elif nc = 0:
		# # gro50
	# elif pc < nc: 
		# # gro25
	# elif pc > nc:
		# # gro52
	# else:
		# # gro55
		
	# # for things that have transcription on both, does the ratio get closer to 1? groBoth
	# if pc > 0 and nc > 0:
		# cRatio = pc/nc
		# cRatioOver1 = cRatio if cRatio >= 1 else 1/cRatio
		# if p5*n5==0: groBoth['notconv'] += 1
		# else:
			# fiveRatio = p5/n5
			# fiveRatioOver1 = fiveRatio if fiveRatio >= 1 else 1/fiveRatio
			# if fiveRatioOver1 < cRatioOver1: groBoth['conv'] += 1
			# else: groBoth['notconv'] += 1
	# # for things with transcription on one side, do you get transcription on both? groOne
	# elif pc*nc==0:
		# if p5*n5 != 0: groOne['conv'] += 1
		# else: groOne['notconv'] += 1

	# use geometric mean metric
	ctrl = pc*nc
	five = p5*n5

	if ctrl==0:
		if five==0: groGeomZero['zero'] += 1
		else: 
			groGeomZero['conv'] += 1

	else:

		if five==0: 
			groGeomConv['zero'] += 1

		elif min(p5, n5) < min(pc, nc): 
			groGeomConv['less'] += 1

		elif min(p5, n5) == min(pc, nc): groGeomConv['same'] += 1
		else: 
			groGeomConv['more'] += 1

	
ifile.close()

#print groOne, groBoth
print 'groGeomZero', groGeomZero
print 'groGeomConv', groGeomConv
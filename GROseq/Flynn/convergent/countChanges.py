# countChanges.py
import csv, collections, math
csv.register_dialect("textdialect",delimiter='\t')

ifile = open("gro_100.bed",'r') # 125p 125n 12Cp 12Cn
reader = csv.reader(ifile, 'textdialect')

# # C --> 5pASO
# gro00 = collections.defaultdict(lambda: 0)
# gro05 = collections.defaultdict(lambda: 0)
# gro50 = collections.defaultdict(lambda: 0)
# gro25 = collections.defaultdict(lambda: 0)
# gro52 = collections.defaultdict(lambda: 0)
# gro55 = collections.defaultdict(lambda: 0)

# groBoth = collections.defaultdict(lambda: 0)
# groOne = collections.defaultdict(lambda: 0)


groGeomZero = collections.defaultdict(lambda: 0)
groGeomConv = collections.defaultdict(lambda: 0)

mfile = open("convergent_more.bedGraph",'w')
mwriter = csv.writer(mfile, 'textdialect')

lfile = open("convergent_less.bedGraph",'w')
lwriter = csv.writer(lfile, 'textdialect')

cfile = open("convergent_wt.bedGraph",'w')
cwriter = csv.writer(cfile, 'textdialect')

afile = open("convergent_aso.bedGraph",'w')
awriter = csv.writer(afile, 'textdialect')

cgfile = open("convergent_wt_geom.bedGraph",'w')
cgwriter = csv.writer(cgfile, 'textdialect')

agfile = open("convergent_aso_geom.bedGraph",'w')
agwriter = csv.writer(agfile, 'textdialect')

cmfile = open("convergent_wt_min.bedGraph",'w')
cmwriter = csv.writer(cmfile, 'textdialect')

amfile = open("convergent_aso_min.bedGraph",'w')
amwriter = csv.writer(amfile, 'textdialect')


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
		
	# use minimum metric
	ctrl = min(pc, nc)
	if ctrl > 0:
		outputRow = row[:3]				#
		outputRow.append(ctrl/100)			# convergent WT
		cmwriter.writerow(outputRow)	#		
	five = min(p5, n5)
	if five > 0:
		outputRow = row[:3]				#
		outputRow.append(five/100)			# convergent ASO
		amwriter.writerow(outputRow)	#	
		
	# use geometric mean metric
	ctrl = pc*nc
	five = p5*n5
	if five > 0:
		outputRow = row[:3]			#
		outputRow.append(1)			# convergent ASO
		awriter.writerow(outputRow)	#
		outputRow = row[:3]
		outputRow.append(math.sqrt(five))
		agwriter.writerow(outputRow)
	if ctrl==0:
		if five==0: groGeomZero['zero'] += 1
		else: 
			groGeomZero['conv'] += 1
			outputRow = row[:3]
			outputRow.append(1)
			mwriter.writerow(outputRow)
	else:
		outputRow = row[:3]			#
		outputRow.append(1)			# convergent WT
		cwriter.writerow(outputRow)	#
		outputRow = row[:3]
		outputRow.append(math.sqrt(five))
		cgwriter.writerow(outputRow)
		if five==0: 
			groGeomConv['zero'] += 1
			outputRow = row[:3]
			outputRow.append(1)
			lwriter.writerow(outputRow)
		elif five < ctrl: 
			groGeomConv['less'] += 1
			outputRow = row[:3]
			outputRow.append(1)
			lwriter.writerow(outputRow)
		elif five == ctrl: groGeomConv['same'] += 1
		else: 
			groGeomConv['more'] += 1
			outputRow = row[:3]
			outputRow.append(1)
			mwriter.writerow(outputRow)
	
ifile.close()
mfile.close()
lfile.close()
cfile.close()
afile.close()
cgfile.close()
agfile.close()
cmfile.close()
amfile.close()
#print groOne, groBoth
print 'groGeomZero', groGeomZero
print 'groGeomConv', groGeomConv
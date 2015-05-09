# combineScatterFiles.py

import csv, sys
csv.register_dialect("textdialect",delimiter='\t')

reg = sys.argv[1]

info = []
gp_nums = []
gn_nums = []
cw_nums = []
ca_nums = []

with open(reg+'_gro_wt_pos.txt','r') as gro_pos, open(reg+'_gro_wt_neg.txt','r') as gro_neg, open(reg+'_conv_wt.txt','r') as conv_wt, open(reg+'_conv_aso.txt') as conv_aso:
	gp = csv.reader(gro_pos,'textdialect')
	gn = csv.reader(gro_neg,'textdialect')
	cw = csv.reader(conv_wt,'textdialect')
	ca = csv.reader(conv_aso,'textdialect')

	for row in gp:
		info.append(row[:6])
		gp_nums.append(float(row[6]))
	for row in gn:
		gn_nums.append(float(row[6]))
	for row in cw:
		cw_nums.append(float(row[6]))
	for row in ca:
		ca_nums.append(float(row[6]))

ofile = open(reg+'_info.txt','w')
writer = csv.writer(ofile, 'textdialect')
writer.writerow(['chr','start','stop','name','zero','strand','gro','conv_wt','conv_aso'])
for i in range(len(info)):		
	outputRow = info[i]
	outputRow.extend([gp_nums[i]+gn_nums[i],cw_nums[i],ca_nums[i]])
	writer.writerow(outputRow)
ofile.close()
	
	
		
		

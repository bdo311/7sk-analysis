import csv
import sys
csv.register_dialect("textdialect", delimiter='\t')
import collections

ifile = open(sys.argv[1], 'r')
reader = csv.reader(ifile, 'textdialect')
ofile = open(sys.argv[2], 'w')
writer = csv.writer(ofile, 'textdialect')

def same_row(a, b):
	for i in range(len(a)):
		if a[i] != b[i]: return False
	return True

curr_info = ['', '', '', '', '', '']
overlaps = ['', '', '', '', '']  # [HEXIM, BAF, 7SK, Promoter, Enhancer]

for row in reader:
	info = row[:6]
	if not same_row(curr_info, info):  # previous row done
		outputRow = curr_info
		outputRow.extend(overlaps)
		if curr_info[0] != '': writer.writerow(outputRow)
		curr_info = info
		overlaps = ['', '', '', '', '']

	# for this row
	if row[6] == 'Hexim1': overlaps[0] += row[10] + ','
	elif row[6] == 'BAF': overlaps[1] += row[10] + ','
	elif row[6] == '7SK': overlaps[2] += row[10] + ','
	elif row[6] == 'Promoters': overlaps[3] += row[10] + ','
	elif row[6] == 'Enhancers': overlaps[4] += row[10] + ','

outputRow = curr_info
outputRow.extend(overlaps)
writer.writerow(outputRow)

ifile.close()
ofile.close()

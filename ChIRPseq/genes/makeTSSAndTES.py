# makeTSSAndTES.py
# 9/8/14
# splits up *_startend file into a TSS and TES file
# python makeTSSAndTES.py <ifile> <name>

import csv, sys
csv.register_dialect("textdialect", delimiter='\t')

def main():
	ifile = open(sys.argv[1], 'r')
	reader = csv.reader(ifile, 'textdialect')
	
	tssfile = open(sys.argv[2] + '_tss.txt', 'w')
	tsswriter = csv.writer(tssfile, 'textdialect')
	
	tesfile = open(sys.argv[2] + '_tes.txt', 'w')
	teswriter = csv.writer(tesfile, 'textdialect')
	
	reader.next()
	for row in reader:
		tssrow = row[:5]
		tssrow.append(row[7])
		tesrow = row[:3]
		tesrow.extend(row[5:])
		if row[2] == '+':
			tsswriter.writerow(tssrow)
			teswriter.writerow(tesrow)
		else:
			tsswriter.writerow(tesrow)
			teswriter.writerow(tssrow)	
	
	tssfile.close()
	tesfile.close()
	ifile.close()

if __name__ == '__main__':
	main()
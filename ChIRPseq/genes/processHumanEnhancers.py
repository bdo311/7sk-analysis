# processHumanEnhancers.py
# 8/8/14
# gets human enhancers

import csv, sys
csv.register_dialect("textdialect", delimiter='\t')
csv.register_dialect("commadialect", delimiter=',')

def main():
	ifile = open(sys.argv[1], 'r')
	reader = csv.reader(ifile, 'commadialect')

	rfile = open("hg19_" + sys.argv[2] + '_reg_enhancers.bed', 'w')
	rwriter = csv.writer(rfile, 'textdialect')

	sfile = open("hg19_" + sys.argv[2] + '_super_enhancers.bed', 'w')
	swriter = csv.writer(sfile, 'textdialect') 

	reader.next()
	for row in reader:
		outputRow = [row[1], row[2], row[3], '__'.join([row[0], row[4]])]
		if row[6]=='1': swriter.writerow(outputRow)
		else: rwriter.writerow(outputRow)

	ifile.close()
	rfile.close()
	sfile.close()

if __name__ == '__main__':
	main()

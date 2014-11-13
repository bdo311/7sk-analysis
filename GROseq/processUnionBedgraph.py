# processUnionBedgraph.py
# 10/24/14
# combines the unionbedgraph output from merging two GROseq replicates

import csv, sys
csv.register_dialect("textdialect", delimiter='\t')

def main():
	ifn = sys.argv[1]
	ofn = sys.argv[2]
	
	ifile = open(ifn, 'r')
	reader = csv.reader(ifile, 'textdialect')
	ofile = open(ofn, 'w')
	writer = csv.writer(ofile, 'textdialect')
	
	for row in reader:
		outputrow = row[:3]
		outputrow.append((float(row[3]) + float(row[4]))/2)
		writer.writerow(outputrow)
		
	ifile.close()
	ofile.close()
	
	
if __name__ == '__main__':
	main()
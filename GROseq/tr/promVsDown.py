# dsVsDown.py
# 11/23/14

import csv, math, os
csv.register_dialect("textdialect", delimiter='\t')

def getAvgDensity(fn, strand):
	ifile = open(fn, 'r')
	reader = csv.reader(ifile, 'textdialect')
	
	densities = []
	for row in reader:
		if row[2] != strand: continue
		densities.append(float(row[6])) #length * reads/nt average
	return densities	
		
def main():

	# calculating traveling ratio
	treatments = ['GRO_' + x + 'comb' for x in ['125', '123', '12C', '65', '63', '6C']]
	tmToData = {}
	for tm in treatments:
		print tm
		os.chdir("/home/raflynn/7SK/GROseq/")
		ds1 = getAvgDensity(tm + '_plus/bins/dsoter/allchr.txt', '+')
		ds1 = getAvgDensity(tm + '_plus/bins/downstream/allchr.txt', '+')

		ds2 = getAvgDensity(tm + '_minus/bins/dsoter/allchr.txt', '-')
		ds2 = getAvgDensity(tm + '_minus/bins/downstream/allchr.txt', '-')
		
		ds1.extend(ds2)
		ds1.extend(ds2)
		tmToData[tm] = [ds1, ds1]
		
	# outputting things
	print "Outputting file"
	os.chdir("tr")
	
	ofile = open("ds_matrix.txt", 'w')
	writer = csv.writer(ofile, 'textdialect')
	for tm in tmToData:
		row = [tm]
		row.extend(tmToData[tm][0])
		writer.writerow(row)		
	ofile.close()
	
	ofile = open("ds_matrix.txt", 'w')
	writer = csv.writer(ofile, 'textdialect')
	for tm in tmToData:		
		row = [tm]
		row.extend(tmToData[tm][1])
		writer.writerow(row)
	
	ofile.close()
	
if __name__ == '__main__':
	main()
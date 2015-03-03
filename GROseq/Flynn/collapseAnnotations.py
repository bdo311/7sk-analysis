# collapseAnnotations.py
# collapse HOMER annotations for genomic regions and enhancers/SE's

import csv, collections
csv.register_dialect("textdialect", delimiter='\t')

def processAnnots(tm, i, annots):
	regionToLen = collections.defaultdict(lambda: 0)
	regionToTotalDensity = collections.defaultdict(lambda: 0)
	
	ifile = open("/home/raflynn/7SK/GROseq/" + tm + "/bins/newannot/allchr.txt", 'r')
	reader = csv.reader(ifile, 'textdialect')
	
	for row in reader:
		reg = row[5]
		l = int(row[3]) - int(row[2])
		regionToLen[reg] += l
		den = float(row[7])
		regionToTotalDensity[reg] += den * l
		
	for reg in regionToTotalDensity:
		den = regionToTotalDensity[reg]/regionToLen[reg]
		annots[reg][i] = den
		
	return annots
	
def processEnh(tm, folder, i, annots):
	length = 0
	totalDensity = 0
	
	ifile = open("/home/raflynn/7SK/GROseq/" + tm + "/bins/" + folder + "/allchr.txt", 'r')
	reader = csv.reader(ifile, 'textdialect')
	
	for row in reader:
		l = int(row[2]) - int(row[1])
		length += l
		den = float(row[4])/300 # the number is the sum of all the regionbins
		totalDensity += den * l
		
	den = totalDensity/length
	annots[folder][i] = den
		
	return annots

def main():
	# get list of treatments
	treatments=['12C', '125', '123', '6C', '65', '63']
	treatments_plus = ['GRO_' + tm + 'comb_plus' for tm in treatments]
	treatments_minus = ['GRO_' + tm + 'comb_minus' for tm in treatments]
	treatments = treatments_plus + treatments_minus
	
	# process new annotations
	annots = collections.defaultdict(lambda: [0]*12)
	for i in range(len(treatments)):
		tm=treatments[i]
		print tm
		annots=processAnnots(tm, i, annots)
		
	# process enhancers and super enhancers
	for i in range(len(treatments)):
		tm=treatments[i]
		print tm
		annots=processEnh(tm, "reg", i, annots)
		annots=processEnh(tm, "super", i, annots)

	# write matrix
	ofile = open("annotations_collapsed.txt", 'w')
	writer = csv.writer(ofile, 'textdialect')
	
	header = ['region']
	header.extend(treatments)
	writer.writerow(header)
	for reg in annots:
		outputrow = [reg]
		outputrow.extend(annots[reg])
		writer.writerow(outputrow)
		
	ofile.close()

if __name__ == "__main__":
	main()
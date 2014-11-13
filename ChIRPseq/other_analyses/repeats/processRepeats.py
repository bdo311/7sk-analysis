# processRepeats.py 
# 8/11/14
# for mm9, H1, and HeLa, gets difference in read density between 7SK enriched and input for repeats

import csv, collections, glob, numpy, sys
from scipy import stats
csv.register_dialect("textdialect", delimiter='\t')

def getRepeatMap(folder, org):
	repeatToDensities = collections.defaultdict(lambda: [])

	chroms = glob.glob(folder + '*.txt')
	for chrom in chroms:
		print chrom
		ifile = open(chrom, 'r')
		reader = csv.reader(ifile, 'textdialect')

		for row in reader:
			if org == 'mouse': repeatToDensities['__'.join(row[5:7])].append(float(row[7]))
			else: repeatToDensities[row[4]].append(float(row[7]))

		ifile.close()

	return repeatToDensities

def main():
	chirpDir = "/home/raflynn/ChIRPseq/"
	folder1 = chirpDir + sys.argv[1] + "/bins/repeats/" #target
	folder2 = chirpDir + sys.argv[2] + "/bins/repeats/" #control
	org = sys.argv[3]

	print "Reading repeat bed files"
	repeats1 = getRepeatMap(folder1, org)
	repeats2 = getRepeatMap(folder2, org)

	print "Aggregating data"
	ofile = open(sys.argv[4], 'w')
	writer = csv.writer(ofile, 'textdialect')
	if org == 'mouse': writer.writerow(['repeat_class', 'repeat_family', 'num', '7sk_mean', '7sk_sd', 'input_mean', 'input_sd', 'fc', 'ks_pval'])
	else: writer.writerow(['repeat', 'num', '7sk_mean', '7sk_sd', 'input_mean', 'input_sd', 'fc', 'ks_pval'])
	for repeat in repeats1:
		if repeat not in repeats2: continue
		nums1 = repeats1[repeat]
		nums2 = repeats2[repeat]
		if org == 'mouse': outputRow = [repeat.split('__'), len(nums1)]
		else: outputRow = [repeat, len(nums1)]

		print repeat, len(nums1)
		mean1, sd1 = numpy.mean(nums1), numpy.std(nums1)
		mean2, sd2 = numpy.mean(nums2), numpy.std(nums2)
		outputRow.extend([mean1, sd1, mean2, sd2, mean1/mean2])
		outputRow.append(stats.ks_2samp(nums1, nums2)[1])

		writer.writerow(outputRow)


if __name__ == '__main__':
	main()

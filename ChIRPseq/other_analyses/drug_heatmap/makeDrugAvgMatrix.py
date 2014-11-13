# makeDrugAvgMatrix.py
# 8/11/14
# makes matrix to be used for visualizing averages of 7SK binding with different drug treatments

import csv, Queue, numpy
csv.register_dialect("textdialect", delimiter='\t')

def getGeneToAvg(fn):
	geneToAvg = {}

	ifile = open(fn, 'r')
	reader = csv.reader(ifile, 'textdialect')

	for row in reader:
		nums = [float(x) for x in row[107:307]]
		geneToAvg[row[0]] = numpy.mean(nums)

	ifile.close()

	geneAvgQueue = Queue.PriorityQueue()
	for gene in geneToAvg:
		avg = geneToAvg[gene]
		geneAvgQueue.put((-avg, gene))

	return geneToAvg, geneAvgQueue

def main(): 
	wt, orderedQueue = getGeneToAvg("/home/raflynn/ChIRPseq/WT2_new/bins/tss/allchr_sorted.txt")
	actd = getGeneToAvg("/home/raflynn/ChIRPseq/ActD_new/bins/tss/allchr_sorted.txt")[0]
	flavo = getGeneToAvg("/home/raflynn/ChIRPseq/Flavo_new/bins/tss/allchr_sorted.txt")[0]
	jq1 = getGeneToAvg("/home/raflynn/ChIRPseq/JQ11_new/bins/tss/allchr_sorted.txt")[0]

	ofile = open("tss_avg7sk_drugs_matrix.txt", 'w')
	writer = csv.writer(ofile, 'textdialect')

	writer.writerow(["gene", "wt2", "actd", "flavo", "jq1"])
	counter = 0
	while not orderedQueue.empty():
		counter += 1
		if counter > 5000: break

		gene = orderedQueue.get()[1]
		writer.writerow([gene, wt[gene], actd[gene], flavo[gene], jq1[gene]])

	ofile.close()




if __name__ == '__main__':
	main()

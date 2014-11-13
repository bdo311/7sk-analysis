# getBidrectionalGenes.py
# 5/24/14
# gets a list of bidirectional genes whose TSS's lie inside a certain gap

import csv, sys, collections, Queue
csv.register_dialect("textdialect", delimiter='\t')

# params
inputFN = sys.argv[1]
outputFN = sys.argv[2]
gap = int(sys.argv[3])

# 1. get queues for all chrs
chrToQueue = collections.defaultdict(lambda: Queue.PriorityQueue())

ifile = open(inputFN, 'r')
reader = csv.reader(ifile, 'textdialect')

reader.next()
for row in reader:
	chrom = row[1]
	if '_' in chrom: continue
	tss = int(row[3]) if row[2] == '+' else int(row[4])
	chrToQueue[chrom].put((tss, [row[0], row[2], tss])) #(priority num, [name, strand, tss])

ifile.close()

# 2. find all pairs of bidirectional promoters
def isBidir(left, right):
	tss1, tss2 = left[0], right[0]
	str1, str2 = left[1][1], right[1][1]
	
	if str1 == '-' and str2 == '+' and tss2 - tss1 < gap: return True
	return False
	
ofile = open(outputFN, 'w')
writer = csv.writer(ofile, 'textdialect')

counter = 0
for chrom in chrToQueue:
	queue = chrToQueue[chrom]
	initial = queue.get()
	while not queue.empty():
		next = queue.get()
		if isBidir(initial, next): 
			print initial[1], next[1]	
			outputRow = []
			outputRow.extend(initial[1])
			outputRow.extend(next[1])
			writer.writerow(outputRow)
			counter += 1
		initial = next
		
		
print counter




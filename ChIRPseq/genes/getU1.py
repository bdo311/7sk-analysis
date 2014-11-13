import csv
import time
import re
import Queue
csv.register_dialect('textdialect', delimiter='\t')

# ifile = open("foundU1_mrna_new2.txt", 'r')
# ofile = open("foundU1_mrna_new3.txt", 'w')
ifile = open("foundU1_gene.txt", 'r')
ofile = open("foundU1_gene_new.txt", 'w')
reader = csv.reader(ifile, 'textdialect')
writer = csv.writer(ofile, 'textdialect')

counter = 0
for row in reader:
	counter += 1
	if len(row) < 1: continue
	name = row[0]
	newName = name.split(' ')[0][22:]
	outputRow = [newName]
	outputRow.extend(row[1:])
	writer.writerow(outputRow)

exit()

ifile = open("gene.fasta", 'r')
ofile = open("foundU1_gene.txt", 'w')
# ifile = open("E:/Downloads/allmRNA.fasta", 'r')
# ofile = open("foundU1_mrna_new2.txt", 'w')
writer = csv.writer(ofile, 'textdialect')

U1 = ["GGTAAG", "GGTGAG", "GTGAGT"]
PAS = ['AATAAA', 'ATTAAA']

def processString(currName, currString, writer):
	numBases = len(currString)
	if numBases > 1000000: return
	u1_sites = Queue.PriorityQueue()
	pas_sites = Queue.PriorityQueue()

	for site in U1:
		p = re.compile(site)
		for m in p.finditer(currString):
			pos = m.start() + 1
			u1_sites.put((pos, pos))

	for site in PAS:
		p = re.compile(site)
		for m in p.finditer(currString):
			pos = m.start() + 1
			pas_sites.put((pos, pos))

	u1_arr = []
	while not u1_sites.empty():
		u1_arr.append(str(u1_sites.get()[1]))

	pas_arr = []
	while not pas_sites.empty():
		pas_arr.append(str(pas_sites.get()[1]))

	u1_str = ','.join(u1_arr)
	pas_str = ','.join(pas_arr)

	writer.writerow([currName, numBases, len(u1_arr), len(u1_arr) * 1000/ float(numBases), len(pas_arr), len(pas_arr)* 1000/ float(numBases), u1_str, pas_str])

currName = ''
currString = ''
firstTime = True

for line in ifile:
	if firstTime:
		firstTime = False
		currName = line[1:]
		continue

	if line[0] == '>': 		
		print(currName)
		processString(currName, currString, writer)
		currName = line[1:]
		currString = ''
		continue

	currString = currString + line.upper().rstrip()

# fencepost for the last string
processString(currName, currString, writer)

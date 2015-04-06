# makeFASTA.py
import csv, sys
csv.register_dialect("textdialect",delimiter='\t')

alreadyRead = set()

ifile1 = open(sys.argv[1],'r')
reader1 = csv.reader(ifile1,'textdialect')
reader1.next()

ofile = open(sys.argv[2],'w')
writer = csv.writer(ofile, 'textdialect')

counter = 0
readcount = 0
for row in reader1:
	read = row[0]
	readcount += int(row[1])
	if read not in alreadyRead:
		#alreadyRead.add(read)
		for i in range(int(row[1])):
			counter += 1
			writer.writerow([">" + str(counter)])
			writer.writerow([read])

ifile1.close()
ofile.close()

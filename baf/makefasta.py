# makeFASTA.py
import csv, sys
csv.register_dialect("textdialect",delimiter='\t')

alreadyRead = set()

ifile1 = open(sys.argv[1],'r')
reader1 = csv.reader(ifile1,'textdialect')
ifile2 = open(sys.argv[2],'r')
reader2 = csv.reader(ifile2,'textdialect')
reader1.next()
reader2.next()

ofile = open(sys.argv[3],'w')
writer = csv.writer(ofile, 'textdialect')

counter = 0
readcount = 0
for row in reader1:
	read = row[0]
	#readcount += int(row[1])
	if read not in alreadyRead:
		alreadyRead.add(read)
		counter += 1
		writer.writerow([">" + str(counter)])
		writer.writerow([read])
		
#print readcount
readcount = 0
for row in reader2:
	read = row[0]
	#readcount += int(row[1])
	if read not in alreadyRead:
		alreadyRead.add(read)
		counter += 1
		writer.writerow([">" + str(counter)])
		writer.writerow([read])	
		
#print readcount
ifile1.close()
ifile2.close()
ofile.close()

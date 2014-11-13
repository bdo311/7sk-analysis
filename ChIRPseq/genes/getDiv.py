import csv
csv.register_dialect("textdialect", delimiter='\t')

ifile = open("annotated_genes.txt", 'r')
reader = csv.reader(ifile, 'textdialect')
ofile = open("divandfailed.txt", 'w')

reader.next()
for row in reader:
	if float(row[3]) > 0 and row[2] !='0': ofile.write(row[0] + '\n')

ifile.close()
ofile.close()
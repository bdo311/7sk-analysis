import csv, sys, os

unshuf = sys.argv[1]
shuf = sys.argv[2]

names = []

ifile = open(unshuf, 'r')
for line in ifile:
	if line[0]=='>': names.append(line)
ifile.close()

ifile = open(shuf, 'r')
ofile = open(shuf + '.temp', 'w')
counter = 0
for line in ifile:
	if line[0] == '>': 
		ofile.write(names[counter])
		counter += 1
	else: ofile.write(line)

ifile.close()
ofile.close()

os.system("rm " + shuf)
os.system("mv " + shuf + ".temp " + shuf)


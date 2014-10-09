# processShuffledDNA.py
# 3/7/14
# Usage: python processShuffledDNA.py <template> <new>

import csv, sys, os

direc = "/home/raflynn/7SK_ChIRPseq/smith_waterman_match_between_RNA_and_DNA/Mm_7SK_WT1_best_match/"
templateFileName = sys.argv[1]
newFileName = sys.argv[2]

outputFileName = direc + "output.txt" #will be the new shuffled alignment file

tempFile = open(templateFileName, 'r')
newFile = open(newFileName, 'r')

ofile = open(outputFileName, 'w')

tempLines = tempFile.readlines()
newLines = newFile.readlines()

for i in range(len(tempLines)):
	if i % 5 != 0: ofile.write(newLines[i]) #writing every line that is not the name, from shuffled
	else: ofile.write(tempLines[i]) #writing every line that is the name, from non-shuffled

tempFile.close()
newFile.close()
ofile.close()

os.system("mv " + outputFileName + ' ' + newFileName)

if 'rc' in templateFileName: isRC = 'rc'
else: isRC = 'norc'
os.system("python ../script/sw2bed.py " + newFileName + ' ' + newFileName + '.bed ' + isRC)
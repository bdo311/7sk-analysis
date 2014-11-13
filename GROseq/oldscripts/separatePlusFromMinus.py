# separatePlusFromMinus.py
# 9/11/14
# gro-seq genes need to be separated by plus and minus for cleavage500

import os, glob, csv, collections
from multiprocessing import Pool
csv.register_dialect("textdialect", delimiter='\t')

# write the mapping to a file
def writeFile(name, mapping, direc):
	ofile = open(direc + name + '.txt', 'w')
	writer = csv.writer(ofile, 'textdialect')
	for treatment in mapping:
		outputRow = [treatment]
		outputRow.extend(mapping[treatment])
		writer.writerow(outputRow)

	ofile.close()


# take in an avg_*_* file and extract the bin values
def processFile(fileName, isMinus):
	avgFile = glob.glob(fileName)[0]
	print avgFile
	ifile = open(avgFile, 'r')
	reader = csv.reader(ifile, 'textdialect')

	values = []
	reader.next()
	if isMinus:
		for row in reader:
			values.append(-float(row[1]))
	else:
		for row in reader:
			values.append(row[1])
	ifile.close()
	return values

def processFolder(folder):
	print folder + "/bins/*"
	fullpaths = glob.glob("/home/raflynn/7SK/GROseq/" + folder + "/bins/*")
	print fullpaths
	for fullpath in fullpaths:
		print fullpath
		
		os.chdir(fullpath)
		
		os.system("rm -f allchr_*.txt")
		os.system("awk -F '\t' '{print >> \"allchr_\" $3 \".txt\"}' allchr.txt")
		# os.system("Rscript /home/raflynn/Scripts/metagene_maker/makeMetagenePlot.r " + folder + " cleavage500 6 0 100 cleavage500_+.txt")
		# os.system("Rscript /home/raflynn/Scripts/metagene_maker/makeMetagenePlot.r " + folder + " cleavage500 6 0 100 cleavage500_-.txt")
	
def main():
	folders = glob.glob("GRO_*2_*")

	#p = Pool(3)
	#p.map(processFolder, folders)
	
	#exit()
	name = "groseq_aso2_direc"
	parentDir = "/home/raflynn/7SK/GROseq/"
	regionToFolderAvgs = collections.defaultdict(lambda: {})
	for region in ["tss", "tes", "geneBody"]:
		for direc in ['+', '-']:
			for folder in folders:
				d = parentDir + folder + "/bins/%s/"%region
				print d
				os.chdir(d)
				fn = "allchr_" + direc + '.txt'
				if ('minus' in folder and direc=='-') or ('plus' in folder and direc=='+'):
					regionToFolderAvgs[region + '_' + direc][folder] = processFile(fn, isMinus=False)
				else: 
					regionToFolderAvgs[region + '_' + direc][folder] = processFile(fn, isMinus=True)
			writeFile(name + '_' + region + '_' + direc, regionToFolderAvgs[region + '_' + direc], parentDir + '/averages/')
	
	
	

if __name__ == "__main__":
	main()
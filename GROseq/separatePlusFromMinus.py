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
	avgFile = glob.glob(fileName + "*.txt")[0]
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
	fullpath = folder + "/bins/cleavage500/"
	print fullpath
	os.chdir("/home/raflynn/GROseq_ASO/" + fullpath)
	os.system("rm -f cleavage500*.txt")
	os.system("awk -F '\t' '{print >> \"cleavage500_\" $4 \".txt\"}' allchr.txt")
	os.system("Rscript /home/raflynn/Scripts/metagene_maker/makeMetagenePlot.r " + folder + " cleavage500 6 0 100 cleavage500_+.txt")
	os.system("Rscript /home/raflynn/Scripts/metagene_maker/makeMetagenePlot.r " + folder + " cleavage500 6 0 100 cleavage500_-.txt")
	
def main():
	folders = glob.glob("GRO_*")
	p = Pool(6)
	p.map(processFolder, folders)
	
	name = "groseq_aso_direc"
	parentDir = "/home/raflynn/GROseq_ASO/"
	regionToFolderAvgs = collections.defaultdict(lambda: {})
	for region in ["cleavage500_+", "cleavage500_-"]:
		for folder in folders:
			os.chdir(parentDir + folder + "/bins/cleavage500/")
			fn = "avgraw_" + folder + "_cleavage500_" + region
			if ('minus' in folder and '-' in region) or ('plus' in folder and '+' in region):
				regionToFolderAvgs[region][folder] = processFile(fn, isMinus=False)
			else: 
				regionToFolderAvgs[region][folder] = processFile(fn, isMinus=True)
		writeFile(name + '_' + region, regionToFolderAvgs[region], parentDir + '/averages/')
	
	
	

if __name__ == "__main__":
	main()
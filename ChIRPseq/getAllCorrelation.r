library(data.table)
library(RColorBrewer)

chipDir = "/home/raflynn/Mm_ChIPseq/mES_sra_data/"
chirpDir = "/home/raflynn/7SK_ChIRPseq/"

chips = c("NELF_A", "H3K4me3", "H3K9me3", "H3K27me3", "CDK9", "CTCF", "CTR9", "RNAPII_all")
chirps = c("WT2")

tes = list()
tss = list()
gene = list()

print("Reading in chirps")
for (chirp in chirps) {
	folder = paste(chirpDir, chirp, sep='')
	print(chirp)
	
	print("gene")
	load(paste(folder, '/bins/geneBody/data.RData', sep=''))
	mat = data.extrap[,order(colnames(data.extrap))]
	gene = c(gene, chirp = list(c(mat)))
	remove(data, data.extrap, data.scaled)

	print("tss")
	load(paste(folder, '/bins/tss/data.RData', sep=''))
	mat = data.extrap[,order(colnames(data.extrap))]
	tss = c(tss, chirp = list(c(mat)))
	remove(data, data.extrap, data.scaled)
	
	print("tes")
	load(paste(folder, '/bins/tes/data.RData', sep=''))
	mat = data.extrap[,order(colnames(data.extrap))]
	tes = c(tes, chirp = list(c(mat)))
	remove(data, data.extrap, data.scaled)
}

# Extrapolate bins to exactly 1000
extrap = function(x) {
  num = sum(!is.na(x))
  approx(1:num, x[1:num], n=1000)$y
}

# chips don't have data.extrap apparently. what the hell.
print("Reading chips")
for (chip in chips) {
	folder = paste(chipDir, chip, sep='')
	print(chip)
	
	print("gene")
	load(paste(folder, '/geneBody/data.RData', sep=''))
	data.extrap = apply(data, 2, extrap) # coerce to exactly 1000 bins

	mat = data.extrap[,order(colnames(data.extrap))]
	gene = c(gene, chip = list(c(mat)))
	remove(data, data.extrap, data.scaled)

	print("tss")
	load(paste(folder, '/tss/data.RData', sep=''))
	data.extrap = apply(data, c(1,2), as.numeric) # doesn't need extrapolation, but needs conversion to numeric

	mat = data.extrap[,order(colnames(data.extrap))]
	tss = c(tss, chip = list(c(mat)))
	remove(data, data.extrap, data.scaled)
	
	print("tes")
	load(paste(folder, '/tes/data.RData', sep=''))
	data.extrap = apply(data, c(1,2), as.numeric) # doesn't need extrapolation, but needs conversion to numeric

	mat = data.extrap[,order(colnames(data.extrap))]
	tes = c(tes, chip = list(c(mat)))
	remove(data, data.extrap, data.scaled)
}

allnames = c(chirps,chips)
gene = as.data.frame(gene)
colnames(gene) = allnames
tes = as.data.frame(tes)
colnames(tes) = allnames
tss = as.data.frame(tss)
colnames(tss) = allnames
save.image("allcorr.RData")

colors=colorRampPalette(c("red", "darkred", "black", "goldenrod4", "yellow"))(50)

tssCor = cor(tss)
geneCor = cor(gene)
tesCor = cor(tes)

par(mar=c(10,10,10,10))
pdf("TSS_corr.pdf", height=10, width=10)
heatmap(tssCor, symm=TRUE, col=colors, main="TSS")
dev.off()

pdf("gene_corr.pdf", height=10, width=10)
heatmap(geneCor, symm=TRUE, col=colors, main="Gene body")
dev.off()

pdf("TES_corr.pdf", height=10, width=10)
heatmap(tesCor, symm=TRUE, col=colors, main="TES")
dev.off()
	
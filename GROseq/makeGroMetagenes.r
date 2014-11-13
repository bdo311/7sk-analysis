# gro at tss for wt, ft , and div genes
	
getMatrix_new = function(fn, start, stop) {
	gro12c.sense.raw = read.table(fn, header=FALSE, nrows=5)
	classes <- sapply(gro12c.sense.raw, class)
	gro12c.sense.raw = read.table(fn, header=FALSE, colClasses=classes)
	genenames = gro12c.sense.raw[,1]
	gro12c.sense.raw = t(gro12c.sense.raw[,start:stop])
	colnames(gro12c.sense.raw) = genenames
	return(t(gro12c.sense.raw))
}


bigenes = read.delim("/home/raflynn/7SK/ChIRPseq/genes/mm9_bidir_genes.txt", sep='\t')
bigenes.200 = bigenes[bigenes[,6]-bigenes[,3]<200,]
bidir = c(as.character(bigenes.200[,1]), as.character(bigenes.200[,4]))

annot = read.delim("/home/raflynn/7SK/ChIRPseq/genes/annotated_genes.txt", sep='\t', header=TRUE)
div = as.character(annot[annot$udRNA.lfc!=0,]$gene)
ft = as.character(annot[annot$termination.defect!=0,]$gene)

tss = getMatrix_new("/home/raflynn/7SK/GROseq/combined/GRO_12C_tss_sense.txt", 8, 407)
wt_meta = apply(tss, 2, mean)
div_meta = apply(tss[which(rownames(tss) %in% div),], 2, mean)
ft_meta = apply(tss[which(rownames(tss) %in% ft),], 2, mean)

tss = getMatrix_new("/home/raflynn/7SK/GROseq/combined/GRO_12C_tss_antisense.txt", 8, 407)
wt_meta_as = -apply(tss, 2, mean)
div_meta_as = -apply(tss[which(rownames(tss) %in% div),], 2, mean)
ft_meta_as = -apply(tss[which(rownames(tss) %in% ft),], 2, mean)

pdf("141014_gro12c_tss_div_ft.pdf", width=6, height=6)
plot(wt_meta, type='l', col='red', lwd=2, ylim=c(-15,30), xaxt='n', xlab="Position relative to TSS (bp)", main="GROseq signal at TSS", ylab="GROseq reads per nt")
lines(div_meta, col="green", lwd=2)
lines(ft_meta, col="gray", lwd=2)
lines(wt_meta_as, col="red", lwd=2)
lines(div_meta_as, col="green", lwd=2)
lines(ft_meta_as, col="gray", lwd=2)
legend("topleft", legend=c("All TSS", "Div. Txn. TSS", "Fail Term. TSS"), lty=1, lwd=2, col=c("red", "green", "gray"))
axis(1, at=c(0,100,200,300,400), labels=c("-1000", "-500", 0, 500, 1000))
abline(h=0)
dev.off()


# all-gro
pdf("141013_allgenes_gro_overall_metagene.pdf", width=9, height=6)
print("tss")
tss = getMatrix_new("/home/raflynn/7SK/GROseq/combined/GRO_12C_tss_sense.txt", 8, 207)
print("tes")
tes = getMatrix_new("/home/raflynn/7SK/GROseq/combined/GRO_12C_tes_sense.txt", 208, 407)
print("gb")
geneBody = read.delim("/home/raflynn/7SK/GROseq/combined/avgraw_GRO_12C_geneBody_GRO_12C_geneBody_sense.txt", colClasses=c("character", "numeric"), header=TRUE)[,1]
tss = apply(tss, 2, mean)
tes = apply(tes, 2, mean)
# geneBody = approx(1:400, geneBody[1:400], n=200)$y

all = c(tss, geneBody, tes)
plot(all, type='l', xaxt='n', xlab=NA, ylab="GROseq read density per nt", ylim=c(-5, 7), lwd=2)
axis(1, at=c(0, 100, 200, 600, 700, 800), labels=c(-1000, -500, "TSS", "TES", +500, +1000))

counter=1
for(treatment in c("125", "123")) {
	print(treatment)
	counter=counter + 1
	print("tss")
	tss = getMatrix_new(paste("/home/raflynn/7SK/GROseq/combined/GRO_",treatment,"_tss_sense.txt",sep=''), 8, 207)
	print("tes")
	tes = getMatrix_new(paste("/home/raflynn/7SK/GROseq/combined/GRO_",treatment,"_tes_sense.txt",sep=''), 208, 407)
	print("gb")
	geneBody = read.delim(paste("/home/raflynn/7SK/GROseq/combined/avgraw_GRO_", treatment, "_geneBody_GRO_", treatment, "_geneBody_sense.txt", sep=''), colClasses=c("character", "numeric"), header=TRUE)[,1]
	tss = apply(tss, 2, mean)
	tes = apply(tes, 2, mean)
	# geneBody = approx(1:400, geneBody[1:400], n=200)$y

	all = c(tss, geneBody, tes)
	lines(all, type='l', col=counter, lwd=2)
}

counter = 0
for(treatment in c("12C", "125", "123")) {
	print(treatment)
	counter=counter + 1
	print("tss")
	tss = getMatrix_new(paste("/home/raflynn/7SK/GROseq/combined/GRO_",treatment,"_tss_antisense.txt",sep=''), 8, 207)
	print("tes")
	tes = getMatrix_new(paste("/home/raflynn/7SK/GROseq/combined/GRO_",treatment,"_tes_antisense.txt",sep=''), 208, 407)
	print("gb")
	geneBody = read.delim(paste("/home/raflynn/7SK/GROseq/combined/avgraw_GRO_", treatment, "_geneBody_GRO_", treatment, "_geneBody_antisense.txt", sep=''), colClasses=c("character", "numeric"), header=TRUE)[,1]
	tss = apply(tss, 2, mean)
	tes = apply(tes, 2, mean)
	# geneBody = approx(1:400, geneBody[1:400], n=200)$y

	all = c(tss, geneBody, tes)
	all = -all
	lines(all, type='l', col=counter, lwd=2)
}
legend("topright", legend=c("12hr Control", "12hr 5' ASO", "12hr 3' ASO"), lty=1, lwd=2, col=c("black", "red", "green"))
abline(h=0)
abline(v=c(200,600),col="lightblue")
dev.off()
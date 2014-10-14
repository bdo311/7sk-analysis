# makeOverallMetagene.r
# 9/8/14
# uses TSS/TES/gene body files to draw overall metagenes

# Extrapolate bins to exactly the number we need
extrap = function(x, numBins) {
  num = sum(!is.na(x))
  if (num==numBins) return(x)
  approx(1:num, x[1:num], n=numBins)$y
}

# get matrix for files
getMatrix = function(fn, n) {
	no_col = count.fields(fn, sep = "\t") 
	data = t(read.table(fn, sep='\t', fill=TRUE, col.names=1:max(no_col))) #can't use fread

	colnames(data) = data[1,]
	data = data[8:nrow(data),]
	data = apply(data, c(1,2), as.numeric)
	apply(data, 2, extrap, numBins=n) # coerce to exactly some number of bins
}

# pdf("140908_allgenes_overall_metagene.pdf", width=9, height=4.5)
# print("tss")
# tss = getMatrix("/home/raflynn/ChIRPseq/WT2_new/bins/tss/allchr.txt", n=400)
# print("tes")
# tes = getMatrix("/home/raflynn/ChIRPseq/WT2_new/bins/tes/allchr.txt", n=400)
# print("gb")
# geneBody = getMatrix("/home/raflynn/ChIRPseq/WT2_new/bins/geneBody/allchr.txt", n=200)
# tss = tss[1:200,which(colnames(tss) %in% colnames(geneBody))]
# tes = tes[201:400,which(colnames(tss) %in% colnames(geneBody))]
# metagene = rbind(tss, geneBody, tes)

# save(tss, tes, geneBody, metagene, file="data.RData")
# # load("data.RData") #tss, tes, and geneBody matrices

# all = apply(metagene, 1, mean)
# plot(all, type='l', xaxt='n', xlab=NA, ylab="7SK read density per nt", ylim=c(0.3, 1.5), lwd=2)
# axis(1, at=c(0, 100, 200, 400, 500, 600), labels=c(-1000, -500, "TSS", "TES", +500, +1000))

# counter=1
# for(treatment in c("ActD_new", "JQ11_new", "Flavo_new", "DRB_new", "7SK_Input")) {
	# print(treatment)
	# counter=counter + 1
	# print("tss")
	# tss = getMatrix(paste("/home/raflynn/ChIRPseq/", treatment, "/bins/tss/allchr.txt", sep=''), n=400)
	# print("tes")
	# tes = getMatrix(paste("/home/raflynn/ChIRPseq/", treatment, "/bins/tes/allchr.txt", sep=''), n=400)
	# print("gb")
	# geneBody = getMatrix(paste("/home/raflynn/ChIRPseq/", treatment, "/bins/geneBody/allchr.txt", sep=''), n=200)
	# tss = tss[1:200,which(colnames(tss) %in% colnames(geneBody))]
	# tes = tes[201:400,which(colnames(tss) %in% colnames(geneBody))]
	# metagene = rbind(tss, geneBody, tes)

	# # save(tss, tes, geneBody, metagene, file="data.RData")
	# # load("data.RData") #tss, tes, and geneBody matrices

	# all = apply(metagene, 1, mean)
	# lines(all, type='l', col=counter, lwd=2)
# }
# legend("topright", legend=c("WT2_new", "ActD_new", "JQ11_new", "Flavo_new", "DRB_new", "7SK_Input"), lty=1, lwd=2, col=1:6)

# dev.off()
	
getMatrix_new = function(fn, start, stop) {
	gro12c.sense.raw = read.table(fn, header=FALSE, nrows=5)
	classes <- sapply(gro12c.sense.raw, class)
	gro12c.sense.raw = read.table(fn, header=FALSE, colClasses=classes)
	genenames = gro12c.sense.raw[,1]
	gro12c.sense.raw = t(gro12c.sense.raw[,start:stop])
	colnames(gro12c.sense.raw) = genenames
	return(t(gro12c.sense.raw))
}

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

## this script is to look at directionality of CTCF and other sites around enhancers

cnet_s1 = read.delim("/home/raflynn/Mitosis/NETseq/metagenes/cNET_HeLa_Rep1_sense/bins/Churchman_CTCF_all/allchr_sorted.txt", header=FALSE, colClasses=c(rep("character",6),rep("numeric",401)))
rownames(cnet_s1) = apply(cnet_s1, 1, function(x) paste(x[1], x[2], x[3], sep='_'))
cnet_as1 = read.delim("/home/raflynn/Mitosis/NETseq/metagenes/cNET_HeLa_Rep1_antisense/bins/Churchman_CTCF_all/allchr_sorted.txt", header=FALSE, colClasses=c(rep("character",6),rep("numeric",401)))
rownames(cnet_as1) = apply(cnet_as1, 1, function(x) paste(x[1], x[2], x[3], sep='_'))
cnet_as1 = cnet_as1[rownames(cnet_s1),]
cnet_s2 = read.delim("/home/raflynn/Mitosis/NETseq/metagenes/cNET_HeLa_Rep2_sense/bins/Churchman_CTCF_all/allchr_sorted.txt", header=FALSE, colClasses=c(rep("character",6),rep("numeric",401)))
rownames(cnet_s2) = apply(cnet_s2, 1, function(x) paste(x[1], x[2], x[3], sep='_'))
cnet_as2 = cnet_s2[rownames(cnet_s1),]
cnet_as2 = read.delim("/home/raflynn/Mitosis/NETseq/metagenes/cNET_HeLa_Rep2_antisense/bins/Churchman_CTCF_all/allchr_sorted.txt", header=FALSE, colClasses=c(rep("character",6),rep("numeric",401)))
rownames(cnet_as2) = apply(cnet_as2, 1, function(x) paste(x[1], x[2], x[3], sep='_'))
cnet_as2 = cnet_as2[rownames(cnet_s1),]

gro_s = read.delim("/home/raflynn/Mitosis/NETseq/metagenes/GRO_HeLa_sense/bins/Churchman_CTCF_all/allchr_sorted.txt", header=FALSE, colClasses=c(rep("character",6),rep("numeric",401)))
rownames(gro_s) = apply(gro_s, 1, function(x) paste(x[1], x[2], x[3], sep='_'))
gro_s = gro_s[rownames(cnet_s1),]
gro_as = read.delim("/home/raflynn/Mitosis/NETseq/metagenes/GRO_HeLa_antisense/bins/Churchman_CTCF_all/allchr_sorted.txt", header=FALSE, colClasses=c(rep("character",6),rep("numeric",401)))
rownames(gro_as) = apply(gro_as, 1, function(x) paste(x[1], x[2], x[3], sep='_'))
gro_as = gro_as[rownames(cnet_s1),]

cnet_all = cbind(
apply(cnet_s1[,208:307],1,sum),
apply(cnet_as1[,108:207],1,sum),
apply(cnet_s2[,208:307],1,sum),
apply(cnet_as2[,108:207],1,sum))
rownames(cnet_all) = rownames(cnet_s1)
colnames(cnet_all) = c("s1", "as1", "s2", "as2")
cnet_all = as.data.frame(cnet_all)
cnet_all$sc = with(cnet_all, (s1 + s2)/2)
cnet_all$asc = with(cnet_all, (as1 + as2)/2)

cnet_all2 = cnet_all
cnet_all2[cnet_all2 < 10] = 10  # average of 0.1

gro_all = cbind(
apply(gro_s[,208:307],1,sum),
apply(gro_as[,108:207],1,sum))
rownames(gro_all) = rownames(cnet_s1)
colnames(gro_all) = c("s", "as")
gro_all = as.data.frame(gro_all)

gro_all2 = gro_all
gro_all2[gro_all2 < 0.1] = 0.1  # average of 0.001

## Churchman NETseq and GROseq plots
pdf("hg19_ctcf_gro_cnet_plots.pdf", width=8, height=11)
par(mfrow=c(3,2))
with(cnet_all, plot(log2(s1/as1), log2(s2/as2), pch=19, cex=0.2, main="Directionality at CTCF sites,\nChurchman NETseq (HeLa) R1 vs R2"))
text(-6, 10, label="r = 0.785")
abline(h=0, v=0, lty=2, lwd=2, col="lightblue")
with(cnet_all, plot(log2(sc/asc), log2(sc+asc), pch=19, cex=0.2, main="Churchman NETseq at CTCF sites", xlim=c(-15,15)))
with(gro_all, plot(log2(s/as), log2(s+as), pch=19, cex=0.2, main="Directionality at CTCF sites:\nGROseq (HeLa)"))
plot(log2(cnet_all$sc/cnet_all$asc), log2(gro_all$s/gro_all$as), pch=19, cex=0.2,
main="Directionality at CTCF sites:\nGROseq vs Churchman NETseq",
xlab="Churchman NETseq directionality", ylab="GROseq directionality")
text(-9, 9, label="r = 0.694")
abline(h=0, v=0, lty=2, lwd=2, col="lightblue")
greater_than = sapply(0:10, function(x) sum(log2(gro_all2$s/gro_all2$as) > x))
less_than = sapply(0:10, function(x) sum(log2(gro_all2$s/gro_all2$as) < -x))
plot(greater_than, type='b', col='blue', xlab="Absolute value of directionality", ylab="Number of CTCF Sites")
points(less_than, type='b', col="orange")
legend("topright", legend=c("Sense biased", "Antisense biased"), lty=2, lwd=1, col=c("blue", "orange"))
title("Number of CTCF sites biased\n in sense or antisense orientation")
dev.off()
cor(log2(cnet_all2$sc/cnet_all2$asc), log2(gro_all2$s/gro_all2$as))
cor(log2(cnet_all2$s1/cnet_all2$as1), log2(cnet_all2$s2/cnet_all2$as2))


## What happens when we filter to only include ChromHMM sites?
pnet_s1 = read.delim("/home/raflynn/Mitosis/NETseq/metagenes/pNET_allPhos_rep1_sense/bins/Churchman_CTCF_ChromHMMonly/allchr_sorted.txt", header=FALSE, colClasses=c(rep("character",6),rep("numeric",401)))
rownames(pnet_s1) = apply(pnet_s1, 1, function(x) paste(x[1], x[2], x[3], sep='_'))
pnet_as1 = read.delim("/home/raflynn/Mitosis/NETseq/metagenes/pNET_allPhos_rep1_antisense/bins/Churchman_CTCF_ChromHMMonly/allchr_sorted.txt", header=FALSE, colClasses=c(rep("character",6),rep("numeric",401)))
rownames(pnet_as1) = apply(pnet_as1, 1, function(x) paste(x[1], x[2], x[3], sep='_'))
pnet_as1 = pnet_as1[rownames(pnet_s1),]

normalize = function(x) {
z = x - min(x)
z/max(z)
}

ma <- function(arr, n){
  res = arr
  for(i in 1:length(arr)){
	res[i] = ifelse(i < n, mean(arr[1:i]), mean(arr[(i-n):i]))
  }
  res
}

pnet_s1_hi = pnet_s1[apply(pnet_s1[,108:307],1,sum)>30,][108:307]
pnet_as1_hi = pnet_as1[apply(pnet_as1[,108:307],1,sum)>30,][108:307]
pnet_s1_norm = t(apply(pnet_s1_hi, 1, normalize))
pnet_as1_norm = t(apply(pnet_as1_hi, 1, normalize))
pnet_s1_norm[!is.finite(pnet_s1_norm)] = 0
pnet_as1_norm[!is.finite(pnet_as1_norm)] = 0

psn_finite = pnet_s1_norm[apply(pnet_s1_norm,1,sum)>0,]
pan_finite = pnet_as1_norm[apply(pnet_as1_norm,1,sum)>0,]

pdf("psn.pdf", width=7, height=7)
for (i in 3:9) {
	km = kmeans(psn_finite, centers=i, nstart=5, iter.max=30)
	par(mfrow=c(3,3))
	for (j in 1:i) {
		plot(ma(apply(psn_finite[km$cluster==j,],2,mean),5), type='l', ylab="Density", main=sum(km$cluster==j))
		abline(v=100, lty=2, col="green")
	}
}
dev.off()

pdf("pan.pdf", width=7, height=7)
for (i in 3:9) {
	km = kmeans(pan_finite, centers=i, nstart=5, iter.max=30)
	par(mfrow=c(3,3))
	for (j in 1:i) {
		plot(ma(apply(pan_finite[km$cluster==j,],2,mean),5), type='l', ylab="Density", main=sum(km$cluster==j))
		abline(v=100, lty=2, col="green")
	}
}
dev.off()

require(RColorBrewer)
require(pheatmap)
colorRamp = colorRampPalette(c("white", "firebrick"))
pdf("pnet_sense.pdf", width=8, height=8)
pheatmap(psn_finite, cluster_rows=TRUE, cluster_cols=FALSE, scale="none", col=colorRamp(50), show_rownames = F,show_colnames = F)
dev.off()

pheatmap(psn_finite[order(apply(psn_finite, 1, which.max)),],cluster_rows=FALSE, cluster_cols=FALSE, scale="none", col=colorRamp(50), show_rowname=FALSE, show_colnames=FALSE)

pheatmap(pan_finite[order(apply(pan_finite, 1, which.max), decreasing=TRUE),],cluster_rows=FALSE, cluster_cols=FALSE, scale="none", col=colorRamp(50), show_rowname=FALSE, show_colnames=FALSE)

## set up data
cnet_s1 = read.delim("/home/raflynn/Mitosis/NETseq/metagenes/cNET_HeLa_Rep1_sense/bins/HeLa_reg_enhancers/allchr_sorted.txt", header=FALSE, colClasses=c(rep("character",6),rep("numeric",401)))
rownames(cnet_s1) = apply(cnet_s1, 1, function(x) paste(x[1], x[2], x[3], sep='_'))
cnet_as1 = read.delim("/home/raflynn/Mitosis/NETseq/metagenes/cNET_HeLa_Rep1_antisense/bins/HeLa_reg_enhancers/allchr_sorted.txt", header=FALSE, colClasses=c(rep("character",6),rep("numeric",401)))
rownames(cnet_as1) = apply(cnet_as1, 1, function(x) paste(x[1], x[2], x[3], sep='_'))
cnet_as1 = cnet_as1[rownames(cnet_s1),]
cnet_s2 = read.delim("/home/raflynn/Mitosis/NETseq/metagenes/cNET_HeLa_Rep2_sense/bins/HeLa_reg_enhancers/allchr_sorted.txt", header=FALSE, colClasses=c(rep("character",6),rep("numeric",401)))
rownames(cnet_s2) = apply(cnet_s2, 1, function(x) paste(x[1], x[2], x[3], sep='_'))
cnet_as2 = cnet_s2[rownames(cnet_s1),]
cnet_as2 = read.delim("/home/raflynn/Mitosis/NETseq/metagenes/cNET_HeLa_Rep2_antisense/bins/HeLa_reg_enhancers/allchr_sorted.txt", header=FALSE, colClasses=c(rep("character",6),rep("numeric",401)))
rownames(cnet_as2) = apply(cnet_as2, 1, function(x) paste(x[1], x[2], x[3], sep='_'))
cnet_as2 = cnet_as2[rownames(cnet_s1),]

gro_s = read.delim("/home/raflynn/Mitosis/NETseq/metagenes/GRO_HeLa_sense/bins/HeLa_reg_enhancers/allchr_sorted.txt", header=FALSE, colClasses=c(rep("character",6),rep("numeric",401)))
rownames(gro_s) = apply(gro_s, 1, function(x) paste(x[1], x[2], x[3], sep='_'))
gro_s = gro_s[rownames(cnet_s1),]
gro_as = read.delim("/home/raflynn/Mitosis/NETseq/metagenes/GRO_HeLa_antisense/bins/HeLa_reg_enhancers/allchr_sorted.txt", header=FALSE, colClasses=c(rep("character",6),rep("numeric",401)))
rownames(gro_as) = apply(gro_as, 1, function(x) paste(x[1], x[2], x[3], sep='_'))
gro_as = gro_as[rownames(cnet_s1),]

chirp_7sk = read.delim("/home/raflynn/Mitosis/NETseq/metagenes/7SK/bins/HeLa_reg_enhancers/allchr_sorted.txt", header=FALSE, colClasses=c(rep("character",6),rep("numeric",401)))
rownames(chirp_7sk) = apply(chirp_7sk, 1, function(x) paste(x[1], x[2], x[3], sep='_'))
chirp_7sk = chirp_7sk[rownames(cnet_s1),]

## aggregate data
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
gro_all2[gro_all2 < 0.5] = 0.5  # average of 0.005

gro_all$chirp7sk_s = apply(chirp_7sk[,208:307],1,sum)
gro_all$chirp7sk_as = apply(chirp_7sk[,108:207],1,sum)
gro_all2$chirp7sk_s = gro_all$chirp7sk_s
gro_all2$chirp7sk_as = gro_all$chirp7sk_as
gro_all2$chirp7sk_s[gro_all2$chirp7sk_s < 30] = 30  # average of 0.3
gro_all2$chirp7sk_as[gro_all2$chirp7sk_as < 30] = 30

## Churchman NETseq and GROseq comparison plots, + comparison with 7SK
pdf("hg19_hela_enhancers_gro_cnet_plots.pdf", width=7, height=10)
par(mfrow=c(3,2))
with(cnet_all, plot(log2(s1/as1), log2(s2/as2), pch=19, cex=0.2, main="Directionality at HeLa enhancers,\nChurchman NETseq (HeLa) R1 vs R2"))
text(-6, 10, label="r = 0.701")
abline(h=0, v=0, lty=2, lwd=2, col="lightblue")
with(cnet_all, plot(log2(sc/asc), log2(sc+asc), pch=19, cex=0.2, main="Churchman NETseq at HeLa enhancers", xlim=c(-15,15)))
with(gro_all, plot(log2(s/as), log2(s+as), pch=19, cex=0.2, main="Directionality at HeLa enhancers:\nGROseq (HeLa)"))
plot(log2(cnet_all$sc/cnet_all$asc), log2(gro_all$s/gro_all$as), pch=19, cex=0.2,
main="Directionality at HeLa enhancers:\nGROseq vs Churchman NETseq",
xlab="Churchman NETseq directionality", ylab="GROseq directionality")
text(-9, 9, label="r = 0.611")
abline(h=0, v=0, lty=2, lwd=2, col="lightblue")
with(gro_all, plot(log2(chirp7sk_s/chirp7sk_as), log2(chirp7sk_s+chirp7sk_as), pch=19, cex=0.2, main="Directionality at HeLa enhancers:\n7SK (HeLa)"))
with(gro_all, plot(log2(s/as), log2(chirp7sk_s/chirp7sk_as), pch=19, cex=0.1, ylim=c(-5,5),
main="GROseq vs 7SK directionality,\nHeLa enhancers", xlab="GROseq", ylab="7SK ChIRPseq"))
abline(h=0, v=0, lty=2, lwd=2, col="lightblue")
text(-7, 4, "r = 0.064,\np = 8.8e-15")
dev.off()
cor(log2(cnet_all2$sc/cnet_all2$asc), log2(gro_all2$s/gro_all2$as))
cor(log2(cnet_all2$s1/cnet_all2$as1), log2(cnet_all2$s2/cnet_all2$as2))
cor(log2(gro_all2$s/gro_all2$as), log2(gro_all2$chirp7sk_s/gro_all2$chirp7sk_as)) 

## let's work with the GROseq. Metagene directionality. THIS IS DIRECTLY COPIED FROM mm9_enhancer_directionality.r
top1k = with(gro_all, gro_all[order(s + as, decreasing=TRUE),])
top1k$dir = (with(top1k, s > as) - 0.5) * 2
sense1_top1k = gro_s[rownames(top1k),]
as1_top1k = gro_as[rownames(top1k),]

more_top1k = matrix(nrow=nrow(sense1_top1k), ncol=400)
for (i in seq(1, nrow(sense1_top1k))) {
	if (top1k[i,]$dir == 1) {
		more_top1k[i,] = as.numeric(sense1_top1k[i,8:407]);
	} else {
		more_top1k[i,] = rev(as.numeric(as1_top1k[i,8:407]));
	}
}

less_top1k = matrix(nrow=nrow(sense1_top1k), ncol=400)
for (i in seq(1, nrow(sense1_top1k))) {
	if (top1k[i,]$dir == 1) {
		less_top1k[i,] = as.numeric(as1_top1k[i,8:407]);
	} else {
		less_top1k[i,] = rev(as.numeric(sense1_top1k[i,8:407]));
	}
}

chirp7sk_top1k = chirp_7sk[rownames(top1k),]
dir7sk_top1k = matrix(nrow=nrow(sense1_top1k), ncol=400)
for (i in seq(1, nrow(chirp7sk_top1k))) {
	if (top1k[i,]$dir == 1) {
		dir7sk_top1k[i,] = as.numeric(chirp7sk_top1k[i,8:407]);
	} else {
		dir7sk_top1k[i,] = rev(as.numeric(chirp7sk_top1k[i,8:407]));
	}
}

ma <- function(arr, n){
  res = arr
  for(i in 1:length(arr)){
	res[i] = ifelse(i < n, mean(arr[1:i]), mean(arr[(i-n):i]))
  }
  res
}

pdf("hg19_hela_enhancer_GRO_directional_metagenes.pdf", width=9, height=4)
par(mfrow=c(1,2))
plot(ma(apply(more_top1k[1:1000,], 2, mean),5), type='l', ylim=c(-1,2),
ylab="GROseq signal", xlab="Position", xaxt='n')
lines(-ma(apply(less_top1k[1:1000,], 2, mean),5))
abline(h=0)
abline(v=200, lty=2, col="lightblue")
axis(1, at=seq(0, 400, 100), labels=seq(-1000,1000,500))
title("GROseq, top 1000 enhancers, oriented")

plot(apply(dir7sk_top1k[1:1000,], 2, mean), type='l',
ylab="7SK ChIRPseq signal", xlab="Position", xaxt='n')
abline(v=200, lty=2, col="lightblue")
title("7SK, top 1000 enhancers, oriented")
axis(1, at=seq(0, 400, 100), labels=seq(-1000,1000,500))
dev.off()

### Ranking by directionality (GROseq only)
topdir = with(gro_all2, gro_all2[order(abs(log2(s/as)), decreasing=TRUE),])
enhancer_info = gro_s[rownames(topdir),1:6]
enhancer_info[,2] = as.numeric(enhancer_info[,2]) + 1000
enhancer_info[,3] = as.numeric(enhancer_info[,3]) - 999
enhancer_info[,5] = with(topdir, abs(log2(s/as)))
enhancer_info[,6] = sapply(with(topdir, sign(log2(s/as))), function(x) ifelse(x>0, "+", "-"))
write.table(enhancer_info, "hg19_hela_enhancer_rankedbydir.bed", sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

### Closest genes??? (GRO)
system("sort -k1,1 -k2,2n hg19_hela_enhancer_rankedbydir.bed > hg19_hela_enhancer_rankedbydir_sorted.bed")
system("sort -k1,1 -k2,2n ~/7SK/ChIRPseq/genes/hg19_refseq_tss_BED6.txt > hg19_refseq_tss_BED6_sorted.bed")
system("bedtools closest -a hg19_hela_enhancer_rankedbydir_sorted.bed -b hg19_refseq_tss_BED6_sorted.bed -t first > hg19_hela_enhancer_rankedbydir_closesttss.txt")

data = read.delim("hg19_hela_enhancer_rankedbydir_closesttss.txt", header=FALSE)
data = data[order(data$V5, decreasing=TRUE),]

data$tss_to_right = ((data$V9 + data$V8)/2 - data$V2) > 0
data$enhancer_sense = (data$V6 == "+")
data$tss_downstream = (data$tss_to_right == data$enhancer_sense)
data$tss_upstream = (data$tss_to_right != data$enhancer_sense)
data$same_strand = (data$V6 == data$V12)
data$tss_same_down = (data$tss_downstream & data$same_strand)
data$tss_same_up = (data$tss_upstream & data$same_strand)
data$dist = abs((data$V9 + data$V8)/2 - data$V2)
rownames(data) = apply(data, 1, function(x) paste(x[1], as.numeric(x[2]) - 1000, as.numeric(x[3]) + 999, sep='_'))
data$amtgro = gro_all[rownames(data),]$s + gro_all[rownames(data),]$as

## Same strand as nearest gene - dependence on directionality index
data$amtgro2 = data$amtgro
data$amtgro2[data$amtgro2 < 1] = 1
pdf("hg19_hela_enhancer_GRO_directionality_withtss.pdf", width=15, height=4)
par(mfrow=c(1,4))
plot(ma(data$same_strand, 1000), type='l',
xlab="Enhancers, ordered", ylab="Probability (moving average) of same orientation")
lines(ma(data[order(data$amtgro, decreasing=TRUE),]$same_strand, 1000), col='red')
title("Probability that nearest TSS is oriented\n the same as enhancer")
legend("topright", legend=c("Ordered by Directionality Index", "Ordered by absolute GROseq signal"), 
col=c("black","red"), lty=1, lwd=2)
abline(v=10123, lty=2, col="lightblue")
abline(h=0.5, col="lightblue")

## Is the nearest gene upstream? And dependence on DI
plot(ma(data$tss_upstream, 1000), type='l',
xlab="Enhancers, ordered", ylab="Probability (moving average) of TSS upstream", ylim=c(0.47, 0.85))
lines(ma(data[order(data$amtgro, decreasing=TRUE),]$tss_upstream, 1000), col='red')
title("Probability that nearest TSS is upstream \nof oriented enhancer")
legend("topright", legend=c("Ordered by Directionality Index", "Ordered by absolute GROseq signal"), 
col=c("black","red"), lty=1, lwd=2)
abline(v=10123, lty=2, col="lightblue")
abline(h=0.5, col="lightblue")

## Upstream AND same strand? Dependence on DI
plot(ma(data$tss_same_up, 1000), type='l',
xlab="Enhancers, ordered", ylab="Probability (moving average) of same orientation and TSS upstream",
ylim=c(0.2, 0.7))
lines(ma(data[order(data$amtgro, decreasing=TRUE),]$tss_same_up, 1000), col='red')
title("Probability of both happening")
legend("topright", legend=c("Ordered by Directionality Index", "Ordered by absolute GROseq signal"), 
col=c("black","red"), lty=1, lwd=2)
abline(v=10123, lty=2, col="lightblue")
abline(h=0.25, col="lightblue")

## Distance to nearest gene - dependent on DI or just absolute txn?
plot(ma(log10(data$dist + 1), 1000), type='l',
xlab="Enhancers, ordered", ylab="log2 Distance (moving average) to nearest TSS", ylim=c(2.5,6))
lines(ma(log10(data[order(data$amtgro, decreasing=TRUE),]$dist + 1), 1000), col='red')
legend("topleft", legend=c("Ordered by Directionality Index", "Ordered by absolute GROseq signal"), 
col=c("black","red"), lty=1, lwd=2)
abline(v=10123, lty=2, col="lightblue")
title("Distance of nearest TSS to enhancer")
dev.off()
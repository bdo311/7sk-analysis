# this set of scripts is to explore directionality of mm9 enhancer GROseq, PROseq, and 7SK ChIRPseq signal

## Set up data

sense1 = read.delim("/arrayAhome/raflynn/7SK/GROseq/Flynn/metagenes/GRO_Scr1_norm_sense/bins/RY_enh_centered/allchr_sorted.txt", header=FALSE, colClasses=c(rep("character",6),rep("numeric",401)))
rownames(sense1) = sense1[,4]
sense2 = read.delim("/arrayAhome/raflynn/7SK/GROseq/Flynn/metagenes/GRO_Scr2_norm_sense/bins/RY_enh_centered/allchr_sorted.txt", header=FALSE, colClasses=c(rep("character",6),rep("numeric",401)))
rownames(sense2) = sense2[,4]
sense2 = sense2[rownames(sense1),]
sense3 = read.delim("/arrayAhome/raflynn/7SK/GROseq/Flynn/metagenes/GRO_Scr3_norm_sense/bins/RY_enh_centered/allchr_sorted.txt", header=FALSE, colClasses=c(rep("character",6),rep("numeric",401)))
rownames(sense3) = sense3[,4]
sense3 = sense3[rownames(sense1),]
sense4 = read.delim("/arrayAhome/raflynn/7SK/GROseq/Flynn/metagenes/GRO_Scr4_norm_sense/bins/RY_enh_centered/allchr_sorted.txt", header=FALSE, colClasses=c(rep("character",6),rep("numeric",401)))
rownames(sense4) = sense4[,4]
sense4 = sense4[rownames(sense1),]
sensec = read.delim("/arrayAhome/raflynn/7SK/GROseq/Flynn/metagenes/GRO_Scr_comb_norm_sense/bins/RY_enh_centered/allchr_sorted.txt", header=FALSE, colClasses=c(rep("character",6),rep("numeric",401)))
rownames(sensec) = sensec[,4]
sensec = sensec[rownames(sense1),]
as1 = read.delim("/arrayAhome/raflynn/7SK/GROseq/Flynn/metagenes/GRO_Scr1_norm_antisense/bins/RY_enh_centered/allchr_sorted.txt", header=FALSE, colClasses=c(rep("character",6),rep("numeric",401)))
rownames(as1) = as1[,4]
as1 = as1[rownames(sense1),]
as2 = read.delim("/arrayAhome/raflynn/7SK/GROseq/Flynn/metagenes/GRO_Scr2_norm_antisense/bins/RY_enh_centered/allchr_sorted.txt", header=FALSE, colClasses=c(rep("character",6),rep("numeric",401)))
rownames(as2) = as2[,4]
as2 = as2[rownames(sense1),]
as3 = read.delim("/arrayAhome/raflynn/7SK/GROseq/Flynn/metagenes/GRO_Scr3_norm_antisense/bins/RY_enh_centered/allchr_sorted.txt", header=FALSE, colClasses=c(rep("character",6),rep("numeric",401)))
rownames(as3) = as3[,4]
as3 = as3[rownames(sense1),]
as4 = read.delim("/arrayAhome/raflynn/7SK/GROseq/Flynn/metagenes/GRO_Scr4_norm_antisense/bins/RY_enh_centered/allchr_sorted.txt", header=FALSE, colClasses=c(rep("character",6),rep("numeric",401)))
rownames(as4) = as4[,4]
as4 = as4[rownames(sense1),]
asc = read.delim("/arrayAhome/raflynn/7SK/GROseq/Flynn/metagenes/GRO_Scr_comb_norm_antisense/bins/RY_enh_centered/allchr_sorted.txt", header=FALSE, colClasses=c(rep("character",6),rep("numeric",401)))
rownames(asc) = asc[,4]
asc = asc[rownames(sense1),]

s1p = read.delim("/arrayAhome/raflynn/7SK/PROseq/metagenes/PRO_Rep1_sense/bins/RY_enh_centered/allchr_sorted.txt", header=FALSE, colClasses=c(rep("character",6),rep("numeric",401)))
rownames(s1p) = s1p[,4]
s1p = s1p[rownames(sense1),]
s2p = read.delim("/arrayAhome/raflynn/7SK/PROseq/metagenes/PRO_Rep2_sense/bins/RY_enh_centered/allchr_sorted.txt", header=FALSE, colClasses=c(rep("character",6),rep("numeric",401)))
rownames(s2p) = s2p[,4]
s2p = s2p[rownames(sense1),]
as1p = read.delim("/arrayAhome/raflynn/7SK/PROseq/metagenes/PRO_Rep1_antisense/bins/RY_enh_centered/allchr_sorted.txt", header=FALSE, colClasses=c(rep("character",6),rep("numeric",401)))
rownames(as1p) = as1p[,4]
as1p = as1p[rownames(sense1),]
as2p = read.delim("/arrayAhome/raflynn/7SK/PROseq/metagenes/PRO_Rep2_antisense/bins/RY_enh_centered/allchr_sorted.txt", header=FALSE, colClasses=c(rep("character",6),rep("numeric",401)))
rownames(as2p) = as2p[,4]
as2p = as2p[rownames(sense1),]

chirp_7sk = read.delim("/arrayAhome/raflynn/7SK/ChIRPseq/metagenes/7SK_mES_WT/bins/RY_enh_centered/allchr_sorted.txt", header=FALSE, colClasses=c(rep("character",6),rep("numeric",401)))
rownames(chirp_7sk) = chirp_7sk[,4]
chirp_7sk = chirp_7sk[rownames(sense1),]

## Aggregate data

all = cbind(
apply(sense1[,208:307],1,sum),
apply(as1[,108:207],1,sum),
apply(sense2[,208:307],1,sum),
apply(as2[,108:207],1,sum),
apply(sense3[,208:307],1,sum),
apply(as3[,108:207],1,sum),
apply(sense4[,208:307],1,sum),
apply(as4[,108:207],1,sum),
apply(sensec[,208:307],1,sum),
apply(asc[,108:207],1,sum),
apply(chirp_7sk[,208:307],1,sum),
apply(chirp_7sk[,108:207],1,sum))
rownames(all) = rownames(sense1)
colnames(all) = c("s1", "as1", "s2", "as2", "s3", "as3", "s4", "as4", "sc", "asc", "s7sk", "as7sk")
all = as.data.frame(all)
all2 = all
all2[all2 < 40] = 40  # corresponds to read density = 0.4

allpro = cbind(
apply(s1p[,208:307],1,sum),
apply(as1p[,108:207],1,sum),
apply(s2p[,208:307],1,sum),
apply(as2p[,108:207],1,sum))
rownames(allpro) = rownames(sense1)
colnames(allpro) = c("s1p", "as1p", "s2p", "as2p")
allpro = as.data.frame(allpro)
allpro2 = allpro
allpro2[allpro2 < 0.2] = 0.2 # corresponds to read density 0.002

#################
### Enhancers ###
#################

## Paired Plots

pdf("mm9_enhancer_GROPRO_directionality_reps.pdf", width=7, height=7)
with(all, pairs(~log2(s1/as1)+log2(s2/as2)+log2(s3/as3)+log2(s4/as4), pch=19, cex=0.2,
main="GROseq replicates (enhancer directionality)"))
with(allpro, pairs(~log2(s1p/as1p)+log2(s2p/as2p), pch=19, cex=0.2,
main="PROseq replicates (enhancer directionality)"))
plot(log2(all$sc/all$asc), log2((allpro$s1p/allpro$as1p + allpro$s2p/allpro$as2p)/2), pch=19, cex=0.2,
main="GROseq vs PROseq directionality (enhancers)")
abline(v=0, h=0, lty=2, lwd=2, col="lightblue")
dev.off()
with(all2, cor(data.frame(log2(s1/as1),log2(s2/as2),log2(s3/as3),log2(s4/as4))))
with(allpro2, cor(data.frame(log2(s1p/as1p),log2(s2p/as2p)))
cor(log2(all2$sc/all2$asc), log2((allpro2$s1p/allpro2$as1p + allpro2$s2p/allpro2$as2p)/2))

## Heatmaps for all 4 replicates (GROseq only)

chg = with(all, cbind(log2(s1/as1), log2(s2/as2), log2(s3/as3), log2(s4/as4)))
#chg2 = chg[log2(apply(all, 1, sum))>10,]
chg2 = chg
chg2 = chg2[is.finite(apply(chg2, 1, sum)),]
chg2[chg2 > 5] = 5
chg2[chg2 < -5] = -5

require(RColorBrewer)
colorRamp = colorRampPalette(c("blue", "yellow", "red"))
require(pheatmap)
pdf("heatmap_rel.pdf", width=3, height=7)
pheatmap(chg2, cluster_rows=TRUE, cluster_cols=FALSE, scale="none", col=colorRamp(50))
dev.off()

sl = function(x) {
	return(sign(x) * log2(abs(x)))
}

chg3 = with(all, cbind(sl(s1 - as1), sl(s2 - as2), sl(s3 - as3), sl(s4 - as4)))
chg3[is.nan(chg3)] = 0
chg2[chg2 > 10] = 10
chg2[chg2 < -10] = -10
colorRamp = colorRampPalette(c("blue", "yellow", "red"))
require(pheatmap)
pdf("heatmap_abs.pdf", width=3, height=7)
pheatmap(chg3, cluster_rows=TRUE, cluster_cols=FALSE, scale="none", col=colorRamp(50))
dev.off()

## Directionality plot (-100% to +100%), V-plot, 7SK

pdf("mm9_enhancer_GROPRO_directionality_plots.pdf", width=7, height=9)
par(mfrow=c(3,2))

top1k = with(all, all[order(sc + asc, decreasing=TRUE),])
top1k$dir = (with(top1k, sc > asc) - 0.5) * 2
top1k$pct_dir = with(top1k, (sc - asc)/(sc + asc))
p500 = with(top1k, pct_dir[is.finite(pct_dir)][1:500])
p1000 = with(top1k, pct_dir[is.finite(pct_dir)][1:1000])
p2000 = with(top1k, pct_dir[is.finite(pct_dir)][1:2000])
pAll = with(top1k, pct_dir[is.finite(pct_dir)])

plot(density(p500), xlim=c(-1.2, 1.2), ylim=c(0,1), main="Directionality of enhancers",
xlab="Directionality", ylab="Density")
lines(density(p1000), col="red")
lines(density(p2000), col="green")
lines(density(pAll), col="blue")
legend("topright", legend=c("Top 500", "Top 1000", "Top 2000", "All"), 
col=c("black", "red", "green", "blue"), lty=1)

with(all, plot(log2(sc/asc), log2(sc+asc), pch=19, cex=0.3))
title("GROseq directionality vs signal")
with(allpro, plot(log2((s1p/as1p + s2p/as2p)/2), log2((s1p + as1p + s2p + as2p)/2), pch=19, cex=0.3, xlim=c(-10,10)))
title("PROseq directionality vs signal")
with(all, plot(log2(s7sk/as7sk), log2(s7sk + as7sk), pch=19,cex=0.3))
title("7SK directionality vs signal")
with(top1k[1:1000,], plot(log2(sc/asc), log2(s7sk/as7sk), pch=19,cex=0.3))
abline(h=0, lty=2, col="lightblue")
abline(v=0, lty=2, col="lightblue")
title("GROseq vs 7SK directionality, top 1000 enhancers")
text(x=4, y=4, labels="r=0.178")
text(x=4, y=3.5, labels="P=3e-6")
dev.off()

zzz = with(top1k[1:1000,], cbind(log2(s1/as1), log2(s7sk/as7sk)))
zzz = zzz[is.finite(apply(zzz,1,sum)),]
cor.test(zzz[,1], zzz[,2]) # r = 0.178, p = 3.0e-6

## Metagenes with enhancers oriented correctly (GROseq only)

sense1_top1k = sense1[rownames(top1k),]
as1_top1k = as1[rownames(top1k),]

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

pdf("mm9_enhancer_GRO_directional_metagenes.pdf", width=9, height=4)
par(mfrow=c(1,2))
plot(apply(more_top1k[1:1000,], 2, mean), type='l', ylim=c(-4, 12),
ylab="GROseq signal", xlab="Position", xaxt='n')
lines(-apply(less_top1k[1:1000,], 2, mean))
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

plot(log2(all2$sc/all2$asc), log2(all2$sc + all2$asc))
sum(all2$sc == 40 & all2$asc == 40)
topdir = with(all2, all2[order(abs(log2(sc/asc)), decreasing=TRUE),])
enhancer_info = sense1[rownames(topdir),1:6]
enhancer_info[,2] = as.numeric(enhancer_info[,2]) + 1000
enhancer_info[,3] = as.numeric(enhancer_info[,3]) - 999
enhancer_info[,5] = with(topdir, abs(log2(sc/asc)))
enhancer_info[,6] = sapply(with(topdir, sign(log2(sc/asc))), function(x) ifelse(x>0, "+", "-"))
write.table(enhancer_info, "mm9_enhancer_rankedbydir.bed", sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

### Closest genes??? (GRO)
system("sort -k1,1 -k2,2n mm9_enhancer_rankedbydir.bed > mm9_enhancer_rankedbydir_sorted.bed")
system("bedtools closest -a mm9_enhancer_rankedbydir_sorted.bed -b mm9_tss_startseq_centered_100_sorted.bed -t first > mm9_enhancer_rankedbydir_closesttss.txt")

data = read.delim("mm9_enhancer_rankedbydir_closesttss.txt", header=FALSE)
data = data[order(data$V5, decreasing=TRUE),]

data$tss_to_right = ((data$V9 + data$V8)/2 - data$V2) > 0
data$enhancer_sense = (data$V6 == "+")
data$tss_downstream = (data$tss_to_right == data$enhancer_sense)
data$tss_upstream = (data$tss_to_right != data$enhancer_sense)
data$same_strand = (data$V6 == data$V12)
data$tss_same_down = (data$tss_downstream & data$same_strand)
data$tss_same_up = (data$tss_upstream & data$same_strand)
data$dist = abs((data$V9 + data$V8)/2 - data$V2)
rownames(data) = data[,4]
data$amtgro = all[rownames(data),]$sc + all[rownames(data),]$asc

ma <- function(arr, n){
  res = arr
  for(i in 1:length(arr)){
	res[i] = ifelse(i < n, mean(arr[1:i]), mean(arr[(i-n):i]))
  }
  res
}

## Same strand as nearest gene - dependence on directionality index
data$amtgro2 = data$amtgro
data$amtgro2[data$amtgro2 < 80] = 80
pdf("mm9_enhancer_GRO_directionality_withtss.pdf", width=15, height=4)
par(mfrow=c(1,4))
plot(ma(data$same_strand, 500), type='l',
xlab="Enhancers, ordered", ylab="Probability (moving average) of same orientation")
lines(ma(data[order(data$amtgro, decreasing=TRUE),]$same_strand, 500), col='red')
title("Probability that nearest TSS is oriented\n the same as enhancer")
legend("topright", legend=c("Ordered by Directionality Index", "Ordered by absolute GROseq signal"), 
col=c("black","red"), lty=1, lwd=2)
abline(v=3298, lty=2, col="lightblue")
abline(h=0.5, col="lightblue")

## Is the nearest gene upstream? And dependence on DI
plot(ma(data$tss_upstream, 500), type='l',
xlab="Enhancers, ordered", ylab="Probability (moving average) of TSS upstream", ylim=c(0.47, 0.75))
lines(ma(data[order(data$amtgro, decreasing=TRUE),]$tss_upstream, 500), col='red')
title("Probability that nearest TSS is upstream \nof oriented enhancer")
legend("topright", legend=c("Ordered by Directionality Index", "Ordered by absolute GROseq signal"), 
col=c("black","red"), lty=1, lwd=2)
abline(v=3298, lty=2, col="lightblue")
abline(h=0.5, col="lightblue")

## Upstream AND same strand? Dependence on DI
plot(ma(data$tss_same_up, 500), type='l',
xlab="Enhancers, ordered", ylab="Probability (moving average) of same orientation and TSS upstream",
ylim=c(0.2, 0.7))
lines(ma(data[order(data$amtgro, decreasing=TRUE),]$tss_same_up, 500), col='red')
title("Probability of both happening")
legend("topright", legend=c("Ordered by Directionality Index", "Ordered by absolute GROseq signal"), 
col=c("black","red"), lty=1, lwd=2)
abline(v=3298, lty=2, col="lightblue")
abline(h=0.25, col="lightblue")

## Distance to nearest gene - dependent on DI or just absolute txn?
plot(ma(log10(data$dist + 1), 500), type='l',
xlab="Enhancers, ordered", ylab="log2 Distance (moving average) to nearest TSS", ylim=c(3,6))
lines(ma(log10(data[order(data$amtgro, decreasing=TRUE),]$dist + 1), 500), col='red')
legend("topleft", legend=c("Ordered by Directionality Index", "Ordered by absolute GROseq signal"), 
col=c("black","red"), lty=1, lwd=2)
abline(v=3298, lty=2, col="lightblue")
title("Distance of nearest TSS to enhancer")
dev.off()

## linear modeling all of the effects
with(data, plot(V5, log2(dist)))
with(data, plot(log2(amtgro), log2(dist)))
with(data, plot(log2(amtgro), log2(dist), pch=19, cex=0.2))
with(data, plot(V5, log2(dist), pch=19, cex=0.2))
with(data, plot(log2(dist), V5, pch=19, cex=0.2))
with(data, plot(log2(dist), log2(amtgro), pch=19, cex=0.2))

with(data, cor.test(log2(dist+1), log2(amtgro+1))) # r = -0.36
with(data, cor.test(log2(dist+1), V5)) # r = -0.24
with(data, cor.test(V5, log2(1+amtgro), pch=19, cex=0.2)) # r = 0.65

a = with(data, lm(log2(dist+1) ~ log2(amtgro+1) + V5))
summary(a)  # the V5 variable does not improve the linear model



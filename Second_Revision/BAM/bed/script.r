etwd("/Users/adamjr/Dropbox/Apps/MyUploader_AJR_Khavarila/7SK/ChIP_edgeR/")

##############################
## Adam
# Baf
# counts = read.table("Annotate_mES_allEnh_centered_pm200bp_baf.bed.txt")
# counts = read.table("Annotate_mES_prom_startseq_centered_pm200bp_baf.bed.txt")

# # Baf - large region
# counts = read.table("Annotate_mES_allEnh_centered_pm1kb_baf.bed.txt")
# counts = read.table("Annotate_mES_prom_startseq_centered_pm1kb_baf.bed.txt")

# # Input A (Baf, Oct4)
# counts = read.table("Annotate_mES_allEnh_centered_pm200bp_inputA.bed.txt")
# counts = read.table("Annotate_mES_prom_startseq_centered_pm200bp_inputA.bed.txt")

# # Oct4
# counts = read.table("Annotate_mES_allEnh_centered_pm200bp_oct4.bed.txt")
# counts = read.table("Annotate_mES_prom_startseq_centered_pm200bp_oct4.bed.txt")

# # gH2Ax
# counts = read.table("Annotate_mES_allEnh_centered_pm500bp_gH2Ax.bed.txt")
# counts = read.table("Annotate_mES_prom_startseq_centered_p0_p500bp_gH2Ax.bed.txt")

# # Input B (gH2Ax)
# counts = read.table("Annotate_mES_allEnh_centered_pm500bp_inputB.bed.txt")
# counts = read.table("Annotate_mES_prom_startseq_centered_p0_p500bp_inputB.bed.txt")

# Brian

#baf_enh_1kb_counts = read.table("BAF_allEnh_1kb.bed")
baf_enh_1kb_counts = read.table("../../AJR/ChIP_counts/Annotate_mES_allEnh_centered_pm1kb_baf.bed.txt")
#baf_enh_200bp_counts = read.table("BAF_allEnh_200bp.bed")
baf_enh_200bp_counts = read.table("../../AJR/ChIP_counts/Annotate_mES_allEnh_centered_pm200bp_baf.bed.txt")
#baf_prom_1kb_counts = read.table("BAF_prom_1kb.bed")
baf_prom_1kb_counts = read.table("../../AJR/ChIP_counts/Annotate_mES_prom_startseq_centered_pm1kb_baf.bed.txt")
#baf_prom_prox_counts = read.table("BAF_prom_prox.bed")
baf_prom_prox_counts = read.table("../../AJR/ChIP_counts/Annotate_mES_prom_startseq_centered_pm200bp_baf.bed.txt")
gh2ax_enh_1kb_counts = read.table("gH2AX_allEnh_1kb.bed")
gh2ax_enh_200bp_counts = read.table("gH2AX_allEnh_200bp.bed")
gh2ax_prom_1kb_counts = read.table("gH2AX_prom_1kb.bed")
gh2ax_prom_prox_counts = read.table("gH2AX_prom_prox.bed")
gro_enh_1kb_counts = read.table("GRO_allEnh_1kb.bed")
gro_enh_200bp_counts = read.table("GRO_allEnh_200bp.bed")
gro_prom_1kb_counts = read.table("GRO_prom_1kb.bed")
gro_prom_prox_counts = read.table("GRO_prom_prox.bed")
gro_prom_dist_counts = read.table("GRO_prom_dist.bed")
oct4_enh_1kb_counts = read.table("Oct4_allEnh_1kb.bed")
oct4_enh_200bp_counts = read.table("Oct4_allEnh_200bp.bed")
oct4_prom_1kb_counts = read.table("Oct4_prom_1kb.bed")
oct4_prom_prox_counts = read.table("Oct4_prom_prox.bed")

##############################

# Library sizes from Homer tagInfo.txt
# arranged ctl1, ctl2, aso1, aso2

lib.size_n = c(24156541, 21252957, 29303059, 20847033)  # Baf
# lib.size_n = c(5609875, 162, 8416414, 14322431) # input A 
lib.size_n = c(21690336, 37247718, 16717755, 27399659) # Oct4
lib.size_n = c(39648219, 36060025, 41679902, 35611134) # gH2Ax
# lib.size_n = c(31852036, 27440705, 25213007, 20900353) # input B 
lib.size_n = c(35844103,22886901,46072000,40614974) # GROseq; ctl3-ctl4, aso3-aso4

lib.size_n = c(1e6,1e6,1e6,1e6)

# Select counts columns in arrangement ctl1, ctl2, aso1, aso2
counts = gro_prom_dist_counts
prefix = "gro_prom_dist"
countsn = counts[,c(9,10,13,14)] # for Brian's GRO
countsn = counts[,c(7,8,9,10)] # for all Brian's files except GRO
countsn = counts[,8:11] # if using Adam's files
# ##############################

# # # Preliminary check of count distributions
# s1=4
# s2=1
# plot(log2((countsn[,s1]/lib.size_n[s1])*1e6+0.1),log2((countsn[,s2]/lib.size_n[s2])*1e6+0.1),pch=19,cex=0.2,main="Oct4: C2 vs A2, lib size norm")
# abline(0,1)

# s1=3
# s2=4
# plot(log2(countsn[,s1]+1),log2(countsn[,s2]+1),pch=19,cex=0.2,main="Oct4: C1 vs A1, r")
# abline(0,1)

# plot(density(asinh(countsn[,1])))
# lines(density(asinh(countsn[,2])),col="darkgreen")
# lines(density(asinh(countsn[,3])),col="darkred")
# lines(density(asinh(countsn[,4])),col="blue")

# plot(density(asinh(countsn[,1]/lib.size_n[1])),xlim=c(0,2e-5))
# lines(density(asinh(countsn[,2]/lib.size_n[2])),col="darkgreen")
# lines(density(asinh(countsn[,3]/lib.size_n[3])),col="darkred")
# lines(density(asinh(countsn[,4]/lib.size_n[4])),col="blue")

# # ########################

# # # More QC - clutering, histograms
# library(gplots)
# library(RColorBrewer)


# # hist(table[,8],n=5000,xlim=c(0,100))
# library(preprocessCore)

# counts_qNorm = normalize.quantiles(as.matrix(countsn))
# plot(log2(counts_qNorm[,4]+1),log2(counts_qNorm[,1]+1),pch=19,cex=0.2)
# abline(0,1)

# cor_raw = cor(countsn,method="pearson")
# cor_qNorm = cor(counts_qNorm)

# colors = c(seq(0.75,1,length=15))
# my_palette <- rev(colorRampPalette(brewer.pal(6,"RdYlBu"))(n = 14)) # zero-centered
# my_palette <- colorRampPalette(c("white","blue","darkblue"))(n = 14)

# heatmap.2(as.matrix(cor_qNorm), col=my_palette,breaks=colors,density.info="none",trace="none",Rowv=T,Colv=T,labRow=c("A2","C2","A1","C1"),labCol=c("A2","C2","A1","C1"),symm=T,margins=c(4,4),colsep=c(0:10),rowsep=c(0:8),sepcolor="black",sepwidth=c(0.01,0.01))


########################

# edgeR differential calls
library(edgeR)

#Brian
y = DGEList(countsn, group=c(2,1,2,1), lib.size=lib.size_n)  # 1 = WT, 2 = ASO. BAF (if using Adam's files)
y = DGEList(countsn, group=c(1,1,2,2), lib.size=lib.size_n)  # 1 = WT, 2 = ASO

#Adam
# y = DGEList(countsn, group=c(1,1,2,2))  # 1 = WT, 2 = ASO
# y = DGEList(countsn, group=c(2,1,2,1))  # 1 = WT, 2 = ASO ** BAF updated ordering
# y = DGEList(countsn[,c(2,3,4)], group=c(1,2,1),lib.size=lib.size_n[c(2,3,4)])  # 1 = WT, 2 = ASO ** BAF updated ordering # drop ASO1 (sample 3) based on clustering

# y = DGEList(countsn, group=c(2,1,2,1),lib.size=lib.size_n)  # 1 = WT, 2 = ASO ** BAF updated ordering # drop ASO1 (sample 3) based on clustering


y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y, pair=c("1", "2"))
edge <- as.data.frame(topTags(et, n= dim(countsn)[1]))

# Select FDR thresholded regions and write to file
edge_padj <- edge[edge$FDR < 0.1, ]
edge_padj_up <- edge_padj[edge_padj$logFC >= 0,]
edge_padj_down <- edge_padj[edge_padj$logFC <= -0,]
head(edge)

# MA plot
alpha=0.8
# xlim=c(-1,3)
# ylim=c(-4,4)
# cex=0.3
# #pdf("enh_ma.pdf", width=6,height=4)
# plot(edge$logCPM,edge$logFC, pch=19,cex=cex,col=rgb(0,0,0,alpha),xlim=xlim,ylim=ylim,xlab="Average Counts (log2)",ylab="ASO/wt (log2)")
# points(edge_padj_up$logCPM, edge_padj_up $logFC, pch=19, cex=cex,col=rgb(1,0,0,alpha))
# points(edge_padj_down$logCPM, edge_padj_down$logFC, pch=19, cex=cex,col=rgb(0,0,1,alpha))
# legend("topright",c("All Enhancers","ASO high","wt high"), col=c("black","red","blue"),inset=0.06,pch=19,cex=0.7)
# title("gH2Ax ChIP Signal at Promoters \n FDR 0.05, log2(FC) > 0.5",cex.main=0.8)
# abline(h=0)
# #dev.off()

alpha=0.5
cex=0.1
pdf(paste(prefix, "_volcano.pdf", sep=''), width=4,height=4)
plot(edge$logFC, -log10(edge$FDR), pch=19,cex=cex,col=rgb(0,0,0,alpha),ylim=c(0,30), xlab="ASO/CTL (log2)",ylab="-log10 FDR")
points(edge_padj_up$logFC, -log10(edge_padj_up$FDR),  pch=19, cex=cex,col=rgb(1,0,0,alpha))
points(edge_padj_down$logFC, -log10(edge_padj_down$FDR),  pch=19, cex=cex,col=rgb(0,0,1,alpha))
legend("topright",c("All Enh","ASO high","CTL high"), col=c("black","red","blue"),inset=0.06,pch=19,cex=0.7)
title(paste(prefix, "\n FDR 0.1, log2(FC) > 0 \n lib.size norm", sep=''),cex.main=0.8)
dev.off()


# Print regions
x = c(2,3,4,1,6,5)  # for Adam's things
x = 1:6
regions_up = cbind(counts[as.numeric(rownames(edge_padj_up)),x], edge_padj_up)
regions_down = cbind(counts[as.numeric(rownames(edge_padj_down)),x], edge_padj_down)
regions_all_annotation = counts[as.numeric(rownames(edge)),x]
regions_all_edge = cbind(regions_all_annotation,edge)

write.table(regions_up, paste(prefix, "_ASO_high_FDR10_FC0.bed", sep=''), quote=F,row.names=F,col.names=T,sep="\t")
write.table(regions_down, paste(prefix, "_CTL_high_FDR10_FC0.bed", sep=''), quote=F,row.names=F,col.names=T,sep="\t")
write.table(regions_all_edge, paste(prefix, "_edgeR_stats.bed", sep=''), quote=F,row.names=F,col.names=F,sep="\t")

## CT ##

/arrayAhome/raflynn/7SK/GROseq/Flynn/metagenes/regionLists

################ Regression #################
# here, get regions_all_edge for each file and set it to the right names

gh2ax_enh_1kb_all = regions_all_edge
baf_enh_1kb_all = regions_all_edge
gro_enh_1kb_all = regions_all_edge
gro_enh_1kb_raw = gro_enh_1kb_counts
baf_enh_1kb_raw = baf_enh_1kb_counts

rownames(gh2ax_enh_1kb_all) = gh2ax_enh_1kb_all$V4
rownames(gro_enh_1kb_all) = gro_enh_1kb_all$V4
rownames(baf_enh_1kb_all) = baf_enh_1kb_all$V1
rownames(gro_enh_1kb_raw) = gro_enh_1kb_raw[,4]
rownames(baf_enh_1kb_raw) = baf_enh_1kb_raw$V1

baf_enh_1kb_all = baf_enh_1kb_all[rownames(gh2ax_enh_1kb_all),]
gro_enh_1kb_all = gro_enh_1kb_all[rownames(gh2ax_enh_1kb_all),]
gro_enh_1kb_raw = gro_enh_1kb_raw[rownames(gh2ax_enh_1kb_all),]
baf_enh_1kb_raw = baf_enh_1kb_raw[rownames(gh2ax_enh_1kb_all),]

gh2ax_fc = gh2ax_enh_1kb_all$logFC
gro_fc = gro_enh_1kb_all$logFC
baf_fc = baf_enh_1kb_all$logFC
gro_raw = log2(apply(gro_enh_1kb_raw[,9:10],1,sum))
gro_raw[!is.finite(gro_raw)] = 0
gro_abs_chg = log2(apply(gro_enh_1kb_raw[,13:14],1,sum) - apply(gro_enh_1kb_raw[,9:10],1,sum))
gro_abs_chg[!is.finite(gro_abs_chg)] = 0
baf_raw = log2(apply(baf_enh_1kb_raw[,c(8,10)],1,sum))
baf_raw[!is.finite(baf_raw)] = 0

ct_re = read.delim("/arrayAhome/raflynn/7SK/GROseq/Flynn/metagenes/regionLists/RY_enh_centered_list_ConvT.txt", header=TRUE, row.names=1)
x = strsplit(rownames(ct_re), '__')
rownames(ct_re) = sapply(x, function(y) paste(y[4:length(y)], collapse='__'))
ct_se = read.delim("/arrayAhome/raflynn/7SK/GROseq/Flynn/metagenes/regionLists/SE_indiv_list_ConvT.txt", header=TRUE, row.names=1)
x = strsplit(rownames(ct_se), '__')
rownames(ct_se) = sapply(x, function(y) paste(y[4:length(y)], collapse='__'))
ct = rbind(ct_re, ct_se)
ct = ct[rownames(gh2ax_enh_1kb_all),]
ct_fc = log2(ct$X5ASOcomb_ConvT/ct$Scrcomb_ConvT)
ct_fc[!is.finite(ct_fc)] = 0

a = lm(gh2ax_fc ~ gro_fc + baf_fc + gro_raw + baf_raw + gro_abs_chg)












library(grDevices)
df_enh = data.frame(gh2ax_fc, gro_fc, baf_fc, gro_raw, ct_fc)
is_te = sapply(rownames(df_enh), function(x) grepl("chr", x))
is_te[is_te == TRUE] = adjustcolor("darkorange", alpha.f=0.4)
is_te[is_te == FALSE] = adjustcolor("deepskyblue2", alpha.f=0.4)
df_enh$is_te = is_te
x = median(gro_raw)

# BAF vs GH2AX
pdf("baf_vs_gh2ax.pdf", width=10, height=5)
par(mfrow=c(1,2))
plot(baf_fc[gro_raw>x], gh2ax_fc[gro_raw>x], pch=19, cex=0.5,
col=df_enh$is_te[gro_raw>x], xlab="log2 BAF FC", ylab="log2 gH2AX FC", main="BAF FC vs GH2AX FC,\nGROseq > median") # -0.22, <2.2e-16, m=-0.266
abline(v=0, lwd=2, lty=2)
abline(h=0, lwd=2, lty=2)
reg1 <- lm(gh2ax_fc[gro_raw>x]~baf_fc[gro_raw>x])
abline(reg1, lwd=2, col="darkred")
legend("topright", legend=c("TE", "SE"), pch=19, col=c("darkorange", "deepskyblue2"))
plot(baf_fc[gro_raw<x], gh2ax_fc[gro_raw<x], pch=19, cex=0.5, ylim=c(-1,2), xlim=c(-1.5,1),
col=df_enh$is_te[gro_raw<x], xlab="log2 BAF FC", ylab="log2 gH2AX FC", main="BAF FC vs GH2AX FC,\nGROseq < median") # -0.053, 2.8e-3, m = -0.0483
abline(v=0, lwd=2, lty=2)
abline(h=0, lwd=2, lty=2)
reg2 <- lm(gh2ax_fc[gro_raw<x]~baf_fc[gro_raw<x])
abline(reg2, lwd=2, col="darkred")
legend("topright", legend=c("TE", "SE"), pch = 19, col=c("darkorange", "deepskyblue2"))
dev.off()

# CT_FC vs GH2AX
pdf("ct_vs_gh2ax.pdf", width=10, height=5)
par(mfrow=c(1,2))
with(df_enh, plot(ct_fc[gro_raw>x], gh2ax_fc[gro_raw>x], pch=19, cex=0.5, 
col=is_te[gro_raw>x], xlab="log2 CT FC", ylab="log2 gH2AX FC", main="CT FC vs GH2AX FC,\nGROseq > median")) # -0.118, 2.8e-10; m = -0.053
abline(v=0, lwd=2, lty=2)
abline(h=0, lwd=2, lty=2)
reg1 <- lm(gh2ax_fc[gro_raw>x]~ct_fc[gro_raw>x])
abline(reg1, lwd=2, col="darkred")
legend("topright", legend=c("TE", "SE"), pch=19, col=c("darkorange", "deepskyblue2"))
with(df_enh, plot(ct_fc[gro_raw<x], gh2ax_fc[gro_raw<x], pch=19, cex=0.5, 
col=is_te[gro_raw<x], xlab="log2 CT FC", ylab="log2 gH2AX FC", main="CT FC vs GH2AX FC,\nGROseq > median")) # 0.242, <2.2e-16, m = 0.04
abline(v=0, lwd=2, lty=2)
abline(h=0, lwd=2, lty=2)
reg2 <- lm(gh2ax_fc[gro_raw<x]~ct_fc[gro_raw<x])
abline(reg2, lwd=2, col="darkred")
legend("topright", legend=c("TE", "SE"), pch = 19, col=c("darkorange", "deepskyblue2"))
dev.off()


# Gro_FC vs GH2AX
pdf("gro_vs_gh2ax.pdf", width=10, height=5)
par(mfrow=c(1,2))
with(df_enh, plot(gro_fc[gro_raw>x], gh2ax_fc[gro_raw>x], pch=19, cex=0.5, 
col=is_te[gro_raw>x], xlab="log2 GROseq FC", ylab="log2 gH2AX FC", main="GROseq FC vs GH2AX FC,\nGROseq > median")) # 0.225, <2.2e-16, m = 0.24
abline(v=0, lwd=2, lty=2)
abline(h=0, lwd=2, lty=2)
reg1 <- lm(gh2ax_fc[gro_raw>x]~gro_fc[gro_raw>x])
abline(reg1, lwd=2, col="darkred")
legend("topright", legend=c("TE", "SE"), pch=19, col=c("darkorange", "deepskyblue2"))
with(df_enh, plot(gro_fc[gro_raw<x], gh2ax_fc[gro_raw<x], pch=19, cex=0.5, 
col=is_te[gro_raw<x], xlab="log2 GROseq FC", ylab="log2 gH2AX FC", main="GROseq FC vs GH2AX FC,\nGROseq < median")) # 0.115, 7.55e-10, m = 0.041
abline(v=0, lwd=2, lty=2)
abline(h=0, lwd=2, lty=2)
reg2 <- lm(gh2ax_fc[gro_raw<x]~gro_fc[gro_raw<x])
abline(reg2, lwd=2, col="darkred")
legend("topright", legend=c("TE", "SE"), pch = 19, col=c("darkorange", "deepskyblue2"))
dev.off()

# Gro_FC vs BAF_FC
pdf("gro_vs_baf.pdf", width=10, height=5)
par(mfrow=c(1,2))
with(df_enh, plot(gro_fc[gro_raw>x], baf_fc[gro_raw>x], pch=19, cex=0.5, 
col=is_te[gro_raw>x], xlab="log2 GROseq FC", ylab="log2 BAF FC", main="GROseq FC vs BAF FC,\nGROseq > median")) # -0.196, <2.2e-16, m=-0.175
abline(v=0, lwd=2, lty=2)
abline(h=0, lwd=2, lty=2)
reg1 <- lm(baf_fc[gro_raw>x]~gro_fc[gro_raw>x])
abline(reg1, lwd=2, col="darkred")
legend("topright", legend=c("TE", "SE"), pch=19, col=c("darkorange", "deepskyblue2"))
with(df_enh, plot(gro_fc[gro_raw<x], baf_fc[gro_raw<x], pch=19, cex=0.5, 
col=is_te[gro_raw<x], xlab="log2 GROseq FC", ylab="log2 BAF FC", main="GROseq FC vs BAF FC,\nGROseq < median")) # -0.022, 0.245, m = -0.0090
abline(v=0, lwd=2, lty=2)
abline(h=0, lwd=2, lty=2)
reg2 <- lm(baf_fc[gro_raw<x]~gro_fc[gro_raw<x])
abline(reg2, lwd=2, col="darkred")
legend("topright", legend=c("TE", "SE"), pch = 19, col=c("darkorange", "deepskyblue2"))
dev.off()

# BAF vs CT
pdf("baf_vs_ct.pdf", width=10, height=5)
par(mfrow=c(1,2))
plot(baf_fc[gro_raw>x], ct_fc[gro_raw>x], pch=19, cex=0.5,
col=df_enh$is_te[gro_raw>x], xlab="log2 BAF FC", ylab="log2 CT FC", main="BAF FC vs CT FC,\nGROseq > median") # -0.0092, 0.625, m = -0.025
abline(v=0, lwd=2, lty=2)
abline(h=0, lwd=2, lty=2)
reg1 <- lm(ct_fc[gro_raw>x]~baf_fc[gro_raw>x])
abline(reg1, lwd=2, col="darkred")
legend("topright", legend=c("TE", "SE"), pch=19, col=c("darkorange", "deepskyblue2"))
plot(baf_fc[gro_raw<x], ct_fc[gro_raw<x], pch=19, cex=0.5, xlim=c(-1.5,1),
col=df_enh$is_te[gro_raw<x], xlab="log2 BAF FC", ylab="log2 CT FC", main="BAF FC vs CT FC,\nGROseq < median") # -0.116, 5.26e-10, m = -0.538
abline(v=0, lwd=2, lty=2)
abline(h=0, lwd=2, lty=2)
reg2 <- lm(ct_fc[gro_raw<x]~baf_fc[gro_raw<x])
abline(reg2, lwd=2, col="darkred")
legend("topright", legend=c("TE", "SE"), pch = 19, col=c("darkorange", "deepskyblue2"))
dev.off()


# Gro_FC vs CT_FC
pdf("gro_vs_ct.pdf", width=10, height=5)
par(mfrow=c(1,2))
with(df_enh, plot(gro_fc[gro_raw>x], ct_fc[gro_raw>x], pch=19, cex=0.5, 
col=is_te[gro_raw>x], xlab="log2 GROseq FC", ylab="log2 CT FC", main="GROseq FC vs CT FC,\nGROseq > median")) # 0.364, < 2.2e-16, m=0.868
abline(v=0, lwd=2, lty=2)
abline(h=0, lwd=2, lty=2)
reg1 <- lm(ct_fc[gro_raw>x]~gro_fc[gro_raw>x])
abline(reg1, lwd=2, col="darkred")
legend("topright", legend=c("TE", "SE"), pch=19, col=c("darkorange", "deepskyblue2"))
with(df_enh, plot(gro_fc[gro_raw<x], ct_fc[gro_raw<x], pch=19, cex=0.5, 
col=is_te[gro_raw<x], xlab="log2 GROseq FC", ylab="log2 CT FC", main="GROseq FC vs CT FC,\nGROseq < median")) # 0.242, <2.2e-16, m=0.463
abline(v=0, lwd=2, lty=2)
abline(h=0, lwd=2, lty=2)
reg2 <- lm(ct_fc[gro_raw<x]~gro_fc[gro_raw<x])
abline(reg2, lwd=2, col="darkred")
legend("topright", legend=c("TE", "SE"), pch = 19, col=c("darkorange", "deepskyblue2"))
dev.off()

# Prediction
a = lm(gh2ax_fc ~ baf_fc + gro_fc + ct_fc, data=df_enh[df_enh$gro_raw>x,]) # 0.35, <2.2e-16; all 3 factors contribute significantly; only ct_fc contributes sig for <x
pdf("prediction_gro_above_median.pdf", width=5, height=5)
pred = predict(a)
actual = df_enh[df_enh$gro_raw>x,]$gh2ax_fc
plot(pred, actual, xlim=c(0,1), col=df_enh$is_te, 
pch=19, cex=0.5, xlab="Predicted", ylab="Actual", main="Predicted vs Actual gH2AX FC,\nGROseq > median")
b = lm(actual ~ pred)
abline(b, lwd=2, col="darkred")
legend("topright", legend=c("TE", "SE"), pch=19, col=c("darkorange", "deepskyblue2"))
dev.off()

pdf("prediction_sequential_gh2ax.pdf", width=9, height=3)
par(mfrow=c(1,3))
actual = df_enh[df_enh$gro_raw>x,]$gh2ax_fc

a = lm(gh2ax_fc ~ gro_fc, data=df_enh[df_enh$gro_raw>x,])
pred = predict(a)
plot(pred, actual, xlim=c(0,1), col=df_enh$is_te, 
pch=19, cex=0.5, xlab="Predicted gH2AX FC", ylab="Actual gH2AX FC", main="GRO fit,\nR = 0.225")
aa = lm(actual ~ pred) # 0.225
#abline(aa, lwd=2, col="darkred")
legend("topright", legend=c("TE", "SE"), pch=19, col=c("darkorange", "deepskyblue2"))

b = lm(gh2ax_fc ~ gro_fc + baf_fc, data=df_enh[df_enh$gro_raw>x,])
pred = predict(b) # 0.288
plot(pred, actual, xlim=c(0,1), col=df_enh$is_te, 
pch=19, cex=0.5, xlab="Predicted gH2AX FC", ylab="Actual gH2AX FC", main="GRO + BAF fit,\nR = 0.288")
bb = lm(actual ~ pred)
#abline(bb, lwd=2, col="darkred")
legend("topright", legend=c("TE", "SE"), pch=19, col=c("darkorange", "deepskyblue2"))

c = lm(gh2ax_fc ~ gro_fc + baf_fc + ct_fc, data=df_enh[df_enh$gro_raw>x,])
pred = predict(c)
plot(pred, actual, xlim=c(0,1), col=df_enh$is_te, 
pch=19, cex=0.5, xlab="Predicted gH2AX FC", ylab="Actual gH2AX FC", main="GRO + BAF + CT fit,\nR = 0.352")
cc = lm(actual ~ pred) # 0.352
#abline(cc, lwd=2, col="darkred")
legend("topright", legend=c("TE", "SE"), pch=19, col=c("darkorange", "deepskyblue2"))

dev.off()

k36 = read.table("Enh_H3K36me3_edgeR_stats.bed")
plot(k36$V12, -log(K36$V15))
rownames(k36) = k36[,1]
k36 = k36[rownames(df_enh),]
df_enh$k36_fc = k36$V12
df_enh$k36_fc2 = 2^k36$V12


d = lm(gh2ax_fc ~ gro_fc + baf_fc + ct_fc + k36_fc, data=df_enh[df_enh$gro_raw>x,])
pred = predict(c)
plot(pred, actual, xlim=c(0,1), col=df_enh$is_te, 
pch=19, cex=0.5, xlab="Predicted gH2AX FC", ylab="Actual gH2AX FC", main="GRO + BAF + CT fit,\nR = 0.352")
cc = lm(actual ~ pred) # 0.352
#abline(cc, lwd=2, col="darkred")
legend("topright", legend=c("TE", "SE"), pch=19, col=c("darkorange", "deepskyblue2"))






# write to file
all_enhancer_bed = read.table("~/7SK/ChIRPseq/genes/mES_All_enhancers_centered_2kb.bed", stringsAsFactors=FALSE)
df_enh_reord = df_enh[all_enhancer_bed$V4,]
combined_df_enh = cbind(all_enhancer_bed, df_enh_reord)
write.table(combined_df_enh, "all_enhancer_data.txt", quote=FALSE, col.names=TRUE, row.names=FALSE, sep='\t')










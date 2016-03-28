# list of files
# files=GROseq_7SK_ASO_rep1_genome_sorted.bam GROseq_7SK_ASO_rep2_genome_sorted.bam GROseq_7SK_ASO_rep3_genome_sorted.bam GROseq_7SK_ASO_rep4_genome_sorted.bam GROseq_7SK_Scr_rep1_genome_sorted.bam GROseq_7SK_Scr_rep2_genome_sorted.bam GROseq_7SK_Scr_rep3_genome_sorted.bam GROseq_7SK_Scr_rep4_genome_sorted.bam
# bedtools multicov -bams $files -bed mm9_tss_startseq_centered_TR_prox30.bed > bed/tss_prox_gro.bed

# Process union of CHiCAGO call data for binomial DI calls
counts = read.table("tss_prox_gro.bed")
counts = read.table("enh_gro.bed")
# # Preliminary check of count distributions
# plot(log2(counts[,1]+1),log2(counts[,3]+1),pch=19,cex=0.2)
# abline(0,1)

# hist(table[,8],n=5000,xlim=c(0,100))
# library(preprocessCore)

# counts_qNorm = normalize.quantiles(as.matrix(counts))
# plot(log2(counts_qNorm[,1]+1),log2(counts_qNorm[,6]+1),pch=19,cex=0.2)
# abline(0,1)

# cor_raw = cor(counts,method="spearman")
# cor_qNorm = cor(counts_qNorm)
# library(gplots)
# library(RColorBrewer)
# colors = c(seq(0.9,1,length=15))
# my_palette <- rev(colorRampPalette(brewer.pal(6,"RdYlBu"))(n = 14)) # zero-centered
# my_palette <- colorRampPalette(c("white","blue","darkblue"))(n = 14) 

# heatmap.2(as.matrix(cor_raw),col=my_palette,breaks=colors,density.info="none",trace="none",Rowv=T,Colv=T,labRow=c(1:8),labCol=c(1:8),symm=T,margins=c(4,4))


# edgeR differential calls
library(edgeR)
countsn = counts[,c(13,14,9,10)]
y = DGEList(countsn, group=c(1,1,2,2), lib.size=c(35844103,22886901,46072000,40614974))  # 1 = WT, 2 = ASO

y <- estimateCommonDisp(y)  
y <- estimateTagwiseDisp(y)
et <- exactTest(y, pair=c("1", "2"))
edge <- as.data.frame(topTags(et, n= dim(counts)[1]))

# Select FDR thresholded regions and write to file
edge_padj <- edge[edge$FDR < 0.05, ]
edge_padj_up <- edge_padj[edge_padj$logFC >= 0.5,]
edge_padj_down <- edge_padj[edge_padj$logFC <= -0.5,]

# MA plot
alpha=0.8
xlim=c(0,12)
ylim=c(-3,3)
cex=0.3
pdf("enh_ma.pdf", width=6,height=4)
plot(edge$logCPM,edge$logFC, pch=19,cex=cex,col=rgb(0,0,0,alpha),xlim=xlim,ylim=ylim,xlab="Average Counts (log2)",ylab="ASO/wt (log2)")
points(edge_padj_up$logCPM, edge_padj_up $logFC, pch=19, cex=cex,col=rgb(1,0,0,alpha))
points(edge_padj_down$logCPM, edge_padj_down$logFC, pch=19, cex=cex,col=rgb(0,0,1,alpha))
legend("topright",c("All Enhancers","ASO high","wt high"), col=c("black","red","blue"),inset=0.06,pch=19,cex=0.7)
title("GRO-seq Signal at Enhancers",cex.main=0.8)
abline(h=0)
dev.off()

pdf("enh_volcano.pdf", width=5,height=5)
plot(edge$logFC, -log10(edge$FDR), pch=19,cex=cex,col=rgb(0,0,0,alpha),ylim=c(0,8), xlab="ASO/WT (log2)",ylab="-log10 FDR")
points(edge_padj_up$logFC, -log10(edge_padj_up$FDR),  pch=19, cex=cex,col=rgb(1,0,0,alpha))
points(edge_padj_down$logFC, -log10(edge_padj_down$FDR),  pch=19, cex=cex,col=rgb(0,0,1,alpha))
legend("topright",c("All TSS","ASO high","wt high"), col=c("black","red","blue"),inset=0.06,pch=19,cex=0.7)
title("GRO-seq Signal at Enhancers",cex.main=0.8)
abline(h=0)
dev.off()



# Print regions
regions_up = counts[as.numeric(rownames(edge_padj_up)),1:6]
regions_down = counts[as.numeric(rownames(edge_padj_down)),1:6]

write.table(regions_up,"Enh_ASO_high_FDR5_FC0.5.bed",quote=F,row.names=F,col.names=F,sep="\t")
write.table(regions_down,"Enh_WT_high_FDR5_FC0.5.bed",quote=F,row.names=F,col.names=F,sep="\t")






process = function(fn) {
	a=read.delim(fn, header=FALSE,colClasses=rep("character",4))[,4]
	a[a=='.']=0
	a=as.numeric(a)
	return(a)
}

h3k27ac=process("ENCODE_Bruce4_H3K27ac_genome_shifted_norm_1kb.bed")
h3k27me3=process("ENCODE_Bruce4_H3K27me3_genome_shifted_norm_1kb.bed")
h3k4me1=process("ENCODE_E14_H3K4me1_genome_shifted_norm_1kb.bed")
h3k4me3=process("ENCODE_E14_H3K4me3_genome_shifted_norm_1kb.bed")
h3k9me3=process("ENCODE_E14_H3K9me3_genome_shifted_norm_1kb.bed")
ctr9=process("Rahl_CTR9_genome_shifted_norm_1kb.bed")
nelfa=process("Rahl_NELFA_genome_shifted_norm_1kb.bed")
pol2_ser2p=process("Rahl_Pol2_Ser2p_genome_shifted_norm_1kb.bed")
pol2_ser5p=process("Rahl_Pol2_Ser5p_genome_shifted_norm_1kb.bed")
spt5=process("Rahl_SPT5_genome_shifted_norm_1kb.bed")
wt_7sk=process("WT_genome_merged_norm_1kb.bed")
brd4=process("rahl_brd4_1kb.bed")
tbp=process("rahl_tbp_1kb.bed")
atacseq_12c=process("ATACseq_12C_merge_norm_1kb.bed")
gro_pos=process("GRO_12Ccomb_positive_1kb.bed")
gro_neg=process("GRO_12Ccomb_negative_1kb.bed")
gro_12c=(gro_pos+gro_neg)/2

data=data.frame(cbind(h3k27ac,h3k27me3,h3k4me1,h3k4me3,h3k9me3,ctr9,nelfa,pol2_ser2p,pol2_ser5p,spt5,wt_7sk,gro_12c,brd4,tbp,atacseq_12c))
c=cor(data)
c_spe=cor(data,method="spe")

require(pheatmap)
require(RColorBrewer)
colorRamp = colorRampPalette(c("white","red"))
pdf("pearson_correlation_heatmap.pdf",width=5,height=5)
pheatmap(c,col=colorRamp(50))
dev.off()
pdf("spearman_correlation_heatmap.pdf",width=5,height=5)
pheatmap(c_spe,col=colorRamp(50))
dev.off()

save.image("data.RData")



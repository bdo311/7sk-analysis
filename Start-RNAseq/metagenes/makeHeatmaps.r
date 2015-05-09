# Regular enhancers
start_sense = read.delim("/arrayAhome/raflynn/7SK/Start-RNAseq/Start-Seq_sense/bins/RY_enh_startseq_centered/allchr_sorted.txt",header=FALSE,colClasses=c(rep("character",7),rep("numeric",400)))
start_as = read.delim("/arrayAhome/raflynn/7SK/Start-RNAseq/Start-Seq_antisense/bins/RY_enh_startseq_centered/allchr_sorted.txt",header=FALSE,colClasses=c(rep("character",7),rep("numeric",400)))

start_sense = start_sense[order(start_sense[,4]),]
start_as = start_as[order(start_as[,4]),]
startseq = start_sense[,8:407] + start_as[,8:407]
rownames(startseq) = start_sense[,4]
startseq$dist = as.numeric(start_sense[,5])
startseq = startseq[order(startseq$dist),]
startseq_vals = log2(startseq[,1:400]+1)

require(pheatmap)
require(RColorBrewer)
colorRamp = colorRampPalette(c("white","blue"))
png("ry_enh_startseq_ctr.png")
pheatmap(startseq_vals,scale="none",cluster_rows=FALSE,cluster_cols=FALSE,col=colorRamp(50),show_rownames=FALSE,show_colnames=FALSE,width=6,height=12,breaks=c(seq(0,4,length.out=50),20))
dev.off()

# Super enhancers
start_sense = read.delim("/arrayAhome/raflynn/7SK/Start-RNAseq/Start-Seq_sense/bins/SE_indiv_startseq_centered/allchr_sorted.txt",header=FALSE,colClasses=c(rep("character",7),rep("numeric",400)))
start_as = read.delim("/arrayAhome/raflynn/7SK/Start-RNAseq/Start-Seq_antisense/bins/SE_indiv_startseq_centered/allchr_sorted.txt",header=FALSE,colClasses=c(rep("character",7),rep("numeric",400)))

start_sense = start_sense[order(start_sense[,4]),]
start_as = start_as[order(start_as[,4]),]
startseq = start_sense[,8:407] + start_as[,8:407]
rownames(startseq) = start_sense[,4]
startseq$dist = as.numeric(start_sense[,5])
startseq = startseq[order(startseq$dist),]
startseq_vals = log2(startseq[,1:400]+1)

require(pheatmap)
require(RColorBrewer)
colorRamp = colorRampPalette(c("white","blue"))
png("se_indiv_startseq_ctr.png")
pheatmap(startseq_vals,scale="none",cluster_rows=FALSE,cluster_cols=FALSE,col=colorRamp(50),show_rownames=FALSE,show_colnames=FALSE,width=6,height=12,breaks=c(seq(0,10,length.out=50),20))
dev.off()


# TSS
start_sense = read.delim("/arrayAhome/raflynn/7SK/Start-RNAseq/metagenes/Start-Seq_sense/bins/TSS_startseq_centered/allchr_sorted.txt",header=FALSE,colClasses=c(rep("character",7),rep("numeric",400)))
start_as = read.delim("/arrayAhome/raflynn/7SK/Start-RNAseq/metagenes/Start-Seq_antisense/bins/TSS_startseq_centered/allchr_sorted.txt",header=FALSE,colClasses=c(rep("character",7),rep("numeric",400)))

start_sense = start_sense[order(start_sense[,4]),]
start_as = start_as[order(start_as[,4]),]
startseq = start_sense[,8:407] + start_as[,8:407]
rownames(startseq) = start_sense[,4]
startseq$dist = as.numeric(start_sense[,5])
startseq = startseq[order(startseq$dist),]
startseq_vals = log2(startseq[,1:400]+1)

require(pheatmap)
require(RColorBrewer)
colorRamp = colorRampPalette(c("white","blue"))
png("tss_startseq_ctr.png")
pheatmap(startseq_vals,scale="none",cluster_rows=FALSE,cluster_cols=FALSE,col=colorRamp(50),show_rownames=FALSE,show_colnames=FALSE,width=6,height=12,breaks=c(seq(0,10,length.out=50),20))
dev.off()


# Random files
data = read.delim("/arrayAhome/raflynn/7SK/ChIRPseq/metagenes/7SK_mES_WT/bins/TSS_centered/allchr_sorted.txt",header=FALSE,colClasses=c(rep("character",7),rep("numeric",400)))
rownames(data) = data[,4]
data$dist = as.numeric(data[,5])
data = data[order(data$dist),]
data_vals = log2(data[,8:407]+1)

require(pheatmap)
require(RColorBrewer)
colorRamp = colorRampPalette(c("white","blue"))
png("tss_chirpseq.png")
pheatmap(data_vals,scale="none",cluster_rows=TRUE,cluster_cols=FALSE,col=colorRamp(50),show_rownames=FALSE,show_colnames=FALSE,width=6,height=12,breaks=c(seq(0,4,length.out=50),20))
dev.off()
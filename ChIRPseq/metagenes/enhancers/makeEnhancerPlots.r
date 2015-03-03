# makeEnhancerPlots.r
# 2/27/15
# combine super enhancer and enhancer plots

# name="7SK_HeLa"
# se_name="HeLa_SE_indiv"
# ry_name="HeLa_enhancers"
name="7SK_H1"
se_name="H1_SE_indiv"
ry_name="H1_enhancers"
# name="7SK_mES_WT"
# se_name="SE_indiv"
# ry_name="RY_enh"
se = read.delim(paste("~/7SK/ChIRPseq/metagenes/",name,"/bins/",se_name,"/allchr_sorted.txt",sep=''),colClasses=c(rep("character",6),rep("numeric",401)),header=FALSE)
ry = read.delim(paste("~/7SK/ChIRPseq/metagenes/",name,"/bins/",ry_name,"/allchr_sorted.txt",sep=''),colClasses=c(rep("character",6),rep("numeric",401)),header=FALSE)

#elbow plot
se.len=se[,7]
ry.len=ry[,7]
comb=c(se.len,ry.len)
o=order(comb,decreasing=TRUE)
color=c(rep(2,length(se.len)),rep(1,length(ry.len)))
comb=comb[o]
color=color[o]
# plot(1:length(comb),comb,col=color,pch=19)

#metagene
allenh=rbind(se,ry)
allenh=allenh[o,]
metagene=apply(allenh[,8:407],2,mean)
pdf(paste(name,"_enhancer_metagene.pdf",sep=''),width=6,height=6)
plot(metagene,type='l',lwd=2,ylim=c(0,max(metagene)*1.3))
dev.off()

#heatmap
allenh.log=log2(allenh[,8:407]+1)
require(pheatmap)
require(RColorBrewer)
colorRamp = colorRampPalette(c("white","red"))
png(paste(name,"_enhancer_heatmap.png",sep=''),width=800,height=800)
pheatmap(allenh.log,col=colorRamp(50),cluster_col=FALSE,cluster_row=FALSE,scale='none',show_rownames=FALSE,show_colnames=FALSE)
dev.off()
png(paste(name,"_enhancer_heatmap_legend.png",sep=''),width=100,height=800)
pheatmap(as.matrix(color),col=colorRamp(2),cluster_col=FALSE,cluster_row=FALSE,scale='none',show_rownames=FALSE,show_colnames=FALSE)
dev.off()
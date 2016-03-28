tss.baf.wt = read.delim("/home/raflynn/7SK/ChIPseq/Flynn/metagenes/BAF155_Scr_merge/bins/TSS_centered/allchr_sorted.txt", colClasses=c(rep("character",6),rep("numeric",401)),header=FALSE)
tss.baf.aso = read.delim("/home/raflynn/7SK/ChIPseq/Flynn/metagenes/BAF155_5p_merge/bins/TSS_centered/allchr_sorted.txt", colClasses=c(rep("character",6),rep("numeric",401)),header=FALSE)
re.baf.wt = read.delim("/home/raflynn/7SK/ChIPseq/Flynn/metagenes/BAF155_Scr_merge/bins/RY_enh_centered/allchr_sorted.txt", colClasses=c(rep("character",6),rep("numeric",401)),header=FALSE)
re.baf.aso = read.delim("/home/raflynn/7SK/ChIPseq/Flynn/metagenes/BAF155_5p_merge/bins/RY_enh_centered/allchr_sorted.txt", colClasses=c(rep("character",6),rep("numeric",401)),header=FALSE)
se.baf.wt = read.delim("/home/raflynn/7SK/ChIPseq/Flynn/metagenes/BAF155_Scr_merge/bins/SE_indiv/allchr_sorted.txt", colClasses=c(rep("character",6),rep("numeric",401)),header=FALSE)
se.baf.aso = read.delim("/home/raflynn/7SK/ChIPseq/Flynn/metagenes/BAF155_5p_merge/bins/SE_indiv/allchr_sorted.txt", colClasses=c(rep("character",6),rep("numeric",401)),header=FALSE)
rownames(tss.baf.wt) = tss.baf.wt[,4]
rownames(tss.baf.aso) = tss.baf.aso[,4]
rownames(re.baf.wt) = re.baf.wt[,4]
rownames(re.baf.aso) = re.baf.aso[,4]
rownames(se.baf.wt) = se.baf.wt[,4]
rownames(se.baf.aso) = se.baf.aso[,4]

enh.baf.wt = rbind(re.baf.wt[,8:407]/0.255, se.baf.wt[,8:407]/0.374)
enh.baf.aso = rbind(re.baf.aso[,8:407]/0.242, se.baf.aso[,8:407]/0.325)
tss.baf.wt = tss.baf.wt[,8:407]/0.249
tss.baf.aso = tss.baf.aso[,8:407]/0.244

pdf("baf.pdf", width=9, height=4)
par(mfrow=c(1,2))
plot(apply(enh.baf.wt,2,mean),type='l',ylim=c(1,3), main="Enhancers")
lines(apply(enh.baf.aso,2,mean),col="red")
plot(apply(tss.baf.wt,2,mean),type='l',ylim=c(1,3), main="TSS")
lines(apply(tss.baf.aso,2,mean),col="red")
dev.off()

require(ICSNP)
HotellingsT2(enh.baf.wt, enh.baf.aso) #2.9e-10, -500 to +500: <2.2e-16
HotellingsT2(tss.baf.wt, tss.baf.aso) #2.2e-16




tss.oct4.wt = read.delim("/home/raflynn/7SK/ChIPseq/Flynn/metagenes/Pou5f1_Scr_mergeAD/bins/TSS_centered/allchr_sorted.txt", colClasses=c(rep("character",6),rep("numeric",401)),header=FALSE)
tss.oct4.aso = read.delim("/home/raflynn/7SK/ChIPseq/Flynn/metagenes/Pou5f1_5p_mergeAD/bins/TSS_centered/allchr_sorted.txt", colClasses=c(rep("character",6),rep("numeric",401)),header=FALSE)
re.oct4.wt = read.delim("/home/raflynn/7SK/ChIPseq/Flynn/metagenes/Pou5f1_Scr_mergeAD/bins/RY_enh_centered/allchr_sorted.txt", colClasses=c(rep("character",6),rep("numeric",401)),header=FALSE)
re.oct4.aso = read.delim("/home/raflynn/7SK/ChIPseq/Flynn/metagenes/Pou5f1_5p_mergeAD/bins/RY_enh_centered/allchr_sorted.txt", colClasses=c(rep("character",6),rep("numeric",401)),header=FALSE)
se.oct4.wt = read.delim("/home/raflynn/7SK/ChIPseq/Flynn/metagenes/Pou5f1_Scr_mergeAD/bins/SE_indiv/allchr_sorted.txt", colClasses=c(rep("character",6),rep("numeric",401)),header=FALSE)
se.oct4.aso = read.delim("/home/raflynn/7SK/ChIPseq/Flynn/metagenes/Pou5f1_5p_mergeAD/bins/SE_indiv/allchr_sorted.txt", colClasses=c(rep("character",6),rep("numeric",401)),header=FALSE)
rownames(tss.oct4.wt) = tss.oct4.wt[,4]
rownames(tss.oct4.aso) = tss.oct4.aso[,4]
rownames(re.oct4.wt) = re.oct4.wt[,4]
rownames(re.oct4.aso) = re.oct4.aso[,4]
rownames(se.oct4.wt) = se.oct4.wt[,4]
rownames(se.oct4.aso) = se.oct4.aso[,4]

enh.oct4.wt = rbind(re.oct4.wt[,8:407]/0.239, se.oct4.wt[,8:407]/0.334)
enh.oct4.aso = rbind(re.oct4.aso[,8:407]/0.235, se.oct4.aso[,8:407]/0.313)
tss.oct4.wt = tss.oct4.wt[,8:407]/0.235
tss.oct4.aso = tss.oct4.aso[,8:407]/0.234

pdf("oct4.pdf", width=9, height=4)
par(mfrow=c(1,2))
plot(apply(enh.oct4.wt,2,mean),type='l',ylim=c(1,5), main="Enhancers")
lines(apply(enh.oct4.aso,2,mean),col="red")
plot(apply(tss.oct4.wt,2,mean),type='l',ylim=c(1,5), main="TSS")
lines(apply(tss.oct4.aso,2,mean),col="red")
dev.off()

require(ICSNP)
HotellingsT2(enh.oct4.wt, enh.oct4.aso) #7.59e-12, -500 to +500: <2.2e-16
HotellingsT2(tss.oct4.wt, tss.oct4.aso) #2.2e-16





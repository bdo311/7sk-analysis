### DATASETS ###

# gro

wt.sense = read.delim("/home/raflynn/7SK/GROseq/Flynn/metagenes/GRO_Scr_comb_norm_sense/bins/RY_enh_centered_4k/allchr_sorted.txt", header=FALSE, colClasses=c(rep("character",6),rep("numeric",401)))
wt.antisense = read.delim("/home/raflynn/7SK/GROseq/Flynn/metagenes/GRO_Scr_comb_norm_antisense/bins/RY_enh_centered_4k/allchr_sorted.txt", header=FALSE, colClasses=c(rep("character",6),rep("numeric",401)))
rownames(wt.sense) = wt.sense[,4]
rownames(wt.antisense) = wt.antisense[,4]
wt.antisense = wt.antisense[rownames(wt.sense),]
wt = wt.sense[,8:407] + wt.antisense[,8:407]

aso.sense = read.delim("/home/raflynn/7SK/GROseq/Flynn/metagenes/GRO_5ASO_comb_norm_sense/bins/RY_enh_centered_4k/allchr_sorted.txt", header=FALSE, colClasses=c(rep("character",6),rep("numeric",401)))
aso.antisense = read.delim("/home/raflynn/7SK/GROseq/Flynn/metagenes/GRO_5ASO_comb_norm_antisense/bins/RY_enh_centered_4k/allchr_sorted.txt", header=FALSE, colClasses=c(rep("character",6),rep("numeric",401)))
rownames(aso.sense) = aso.sense[,4]
rownames(aso.antisense) = aso.antisense[,4]
aso.antisense = aso.antisense[rownames(aso.sense),]
aso = aso.sense[,8:407] + aso.antisense[,8:407]
aso = aso[rownames(wt),]

# chip
baf.wt = read.delim("/home/raflynn/7SK/ChIPseq/Flynn/metagenes/BAF155_Scr_merge/bins/RY_enh_centered/allchr_sorted.txt", header=FALSE, colClasses=c(rep("character",6),rep("numeric",401)))
baf.aso = read.delim("/home/raflynn/7SK/ChIPseq/Flynn/metagenes/BAF155_5p_merge/bins/RY_enh_centered/allchr_sorted.txt", header=FALSE, colClasses=c(rep("character",6),rep("numeric",401)))
rownames(baf.wt) = baf.wt[,4]
rownames(baf.aso) = baf.aso[,4]
baf.wt = baf.wt[,8:407]
baf.aso = baf.aso[,8:407]
baf.wt = baf.wt[rownames(wt),]
baf.aso = baf.aso[rownames(wt),]

h2ax.wt = read.delim("/home/raflynn/7SK/ChIPseq/Flynn/metagenes/mES_7SKaso_Scr1_gH2AX/bins/RY_enh_centered/allchr_sorted.txt", header=FALSE, colClasses=c(rep("character",6),rep("numeric",401)))
h2ax.aso = read.delim("/home/raflynn/7SK/ChIPseq/Flynn/metagenes/mES_7SKaso_5ASO1_gH2AX/bins/RY_enh_centered/allchr_sorted.txt", header=FALSE, colClasses=c(rep("character",6),rep("numeric",401)))
rownames(h2ax.wt) = h2ax.wt[,4]
rownames(h2ax.aso) = h2ax.aso[,4]
h2ax.wt = h2ax.wt[,8:407]
h2ax.aso = h2ax.aso[,8:407]
h2ax.wt = h2ax.wt[rownames(wt),]
h2ax.aso = h2ax.aso[rownames(wt),]


chg = log2(apply(aso[,150:250],1,sum)/apply(wt[,150:250],1,sum))
chg = chg[order(chg, decreasing=TRUE)]
chg = chg[is.finite(chg)]  # 5347


a = log2(apply(aso[,150:250],1,sum))
b = log2(apply(aso[,150:250],1,sum)/apply(wt[,150:250],1,sum))
plot(a,b)

# top 500
wt_top500 = wt[names(chg[1:500]),]
aso_top500 = aso[names(chg[1:500]),]
plot(apply(aso_top500,2,mean), type='l', ylim=c(0,4), col="red")
lines(apply(wt_top500,2,mean))

baf.wt_top500 = baf.wt[names(chg[1:500]),]
baf.aso_top500 = baf.aso[names(chg[1:500]),]
plot(apply(baf.aso_top500,2,mean), type='l', ylim=c(0,1), col="red")
lines(apply(baf.wt_top500,2,mean))

chg = log2(apply(aso,1,sum)/apply(wt,1,sum))
baf.chg = log2(apply(baf.aso[,158:257],1,sum)/apply(baf.wt[,158:257],1,sum))

df = data.frame(chg, baf.chg)
df = df[is.finite(apply(df,1,sum)),]

baf.wt_top500 = baf.wt[names(chg[1:5347]),]
baf.aso_top500 = baf.aso[names(chg[1:5347]),]
plot(apply(baf.aso_top500,2,mean), type='l', ylim=c(0,1), col="red")
lines(apply(baf.wt_top500,2,mean))

x = apply(wt[,150:250], 1, sum)
wt_top500 = wt[order(x, decreasing=TRUE),]

baf.wt_top500 = baf.wt[rownames(wt_top500[1:500,]),]
baf.aso_top500 = baf.aso[rownames(wt_top500[1:500,]),]
plot(apply(baf.aso_top500,2,mean), type='l', ylim=c(0,1), col="red")
lines(apply(baf.wt_top500,2,mean))

h2ax.wt_top500 = h2ax.wt[names(chg[4847:5347]),]
h2ax.aso_top500 = h2ax.aso[names(chg[4847:5347]),]
plot(apply(h2ax.aso_top500,2,mean), type='l', ylim=c(0,1000), col="red")
lines(apply(h2ax.wt_top500,2,mean))

h2ax.wt_top500 = h2ax.wt[names(chg[1:500]),]
h2ax.aso_top500 = h2ax.aso[names(chg[1:500]),]
plot(apply(h2ax.aso_top500,2,mean), type='l', ylim=c(0,1000), col="red")
lines(apply(h2ax.wt_top500,2,mean))


# gro_chip_chirp.r
# 5/24/14
# to explore the contribution of divtx/bidir to the upstream peak

require(amap)
load("data.RData")

data.s = data.extrap[,order(apply(data.extrap, 2, sum), decreasing=TRUE)]
data.s = data.extrap[100:300,order(apply(data.extrap, 2, sum), decreasing=TRUE)]
# data.s=data.s[,1:5000]
data.ss = data.s[, apply(data.s, 2, max)!=0]
data.sss = apply(data.ss, 2, function(x) x/max(x))

bigenes = read.delim("/home/raflynn/7SK/ChIRPseq/genes/mm9_bidir_genes.txt", sep='\t')
bigenes.200 = bigenes[bigenes[,6]-bigenes[,3]<200,]
bidir = c(as.character(bigenes.200[,1]), as.character(bigenes.200[,4]))

annot = read.delim("/home/raflynn/7SK/ChIRPseq/genes/annotated_genes.txt", sep='\t', header=TRUE)
div = as.character(annot[annot$udRNA.lfc!=0,]$gene)
ft = as.character(annot[annot$termination.defect!=0,]$gene)

`%ni%` = Negate(`%in%`)
data.nodb = data.ss[,which(colnames(data.ss) %ni% c(bidir, div))] #no bidir, div
data.nodbs = apply(data.nodb, 2, function(x) x/max(x)) #scaled
collapseTo = 50
data.nodbs.coll = apply(data.nodbs, 2, function(x) approx(1:length(x), x, n=collapseTo)$y) #collapsed to 50 bins

# k-means
d = Dist(t(data.nodbs.coll), method="pearson", nbproc=10)
km = kmeans(t(data.nodbs), centers=3, nstart=3)

par(mfrow=c(1,3))
plot(type='l', apply(data.nodb[,km$cluster==1], 1, mean))
plot(type='l', apply(data.nodb[,km$cluster==2], 1, mean))
plot(type='l', apply(data.nodb[,km$cluster==3], 1, mean))

### ChIRP

# chirp heatmap
cl1 = data.nodb[,km$cluster==1]
cl1 = cl1[,order(apply(cl1, 2, sum), decreasing=TRUE)]
cl2 = data.nodb[,km$cluster==2]
cl2 = cl2[,order(apply(cl2, 2, sum), decreasing=TRUE)]
cl3 = data.nodb[,km$cluster==3]
cl3 = cl3[,order(apply(cl3, 2, sum), decreasing=TRUE)]

chirp.tss = t(cbind(cl2, cl1, cl3))
chirp.tss.sc = log10(chirp.tss + 1)

require(RColorBrewer)
colorRamp = colorRampPalette(c("white", "blue"))
require(pheatmap)
png("chirp_tss_all.png", width=1000, height=2000)
pheatmap(chirp.tss.sc, cluster_rows = FALSE, cluster_cols = FALSE, scale="none", show_rownames=FALSE, show_colnames=FALSE,col=colorRamp(50))
dev.off()

### Pol II

# polII heatmap
pol2.raw = read.table("/arrayAhome/raflynn/mmchipseq/RNAPII_all/bins/tss/allchr.txt", header=FALSE, nrows=5)
classes <- sapply(pol2.raw, class)
pol2.raw = read.table("/arrayAhome/raflynn/mmchipseq/RNAPII_all/bins/tss/allchr.txt", header=FALSE, colClasses=classes)
genenames = pol2.raw[,1]
pol2.raw = t(pol2.raw[,108:308])
colnames(pol2.raw) = genenames
pol2.raw = t(pol2.raw)
pol2.sc = log10(pol2.raw[rownames(chirp.tss),]+1)

colorRamp = colorRampPalette(c("white", "darkgreen"))
png("pol2_all_tss.png", width=1000, height=2000)
pheatmap(pol2.sc, cluster_rows = FALSE, cluster_cols = FALSE, scale="none", show_rownames=FALSE, show_colnames=FALSE,col=colorRamp(50))
dev.off()

### GRO Sense

# gro heatmap
gro12c.sense.raw = read.table("/arrayAhome/raflynn/7SK/GROseq/combined/GRO_12C_tss_sense.txt", header=FALSE, nrows=5)
classes <- sapply(gro12c.sense.raw, class)
gro12c.sense.raw = read.table("/arrayAhome/raflynn/7SK/GROseq/combined/GRO_12C_tss_sense.txt", header=FALSE, colClasses=classes)
genenames = gro12c.sense.raw[,1]
gro12c.sense.raw = t(gro12c.sense.raw[,108:308])
colnames(gro12c.sense.raw) = genenames
gro12c.sense.raw = t(gro12c.sense.raw)
gro12c.sense.sc = log10(gro12c.sense.raw[rownames(chirp.tss),]+1)

colorRamp = colorRampPalette(c("white", "firebrick4"))
png("gro12c_tss_sense.png", width=1000, height=2000)
pheatmap(gro12c.sense.sc, cluster_rows = FALSE, cluster_cols = FALSE, scale="none", show_rownames=FALSE, show_colnames=FALSE,col=colorRamp(50))
dev.off()

# gro ASO heatmap
gro123.sense.raw = read.table("/arrayAhome/raflynn/7SK/GROseq/combined/GRO_123_tss_sense.txt", header=FALSE, nrows=5)
classes <- sapply(gro123.sense.raw, class)
gro123.sense.raw = read.table("/arrayAhome/raflynn/7SK/GROseq/combined/GRO_123_tss_sense.txt", header=FALSE, colClasses=classes)
genenames = gro123.sense.raw[,1]
gro123.sense.raw = t(gro123.sense.raw[,108:308])
colnames(gro123.sense.raw) = genenames
gro123.sense.raw = t(gro123.sense.raw)
gro123.sense.sc = log10(gro123.sense.raw[rownames(chirp.tss),]+1)

plot(apply(gro12c.sense.raw[colnames(cl3),],2,mean),type='l')
lines(apply(gro123.sense.raw[colnames(cl3),],2,mean),type='l',col="red")

# gro changes
gro12.chg = gro123.sense.sc - gro12c.sense.sc
breaks=c(seq(-2,-0.01,length.out=25),seq(0.01,2,length.out=25))
colorRamp2 = colorRampPalette(c("green4","green3","white","firebrick3","firebrick4"))
png("gro_12_chg.png", width=1000, height=2000)
pheatmap(gro12.chg, cluster_rows = FALSE, cluster_cols = FALSE, scale="none", show_rownames=FALSE, show_colnames=FALSE,col=colorRamp2(50), breaks=breaks)
dev.off()


### PolII ChIP
pol2.raw = read.table("/arrayAhome/raflynn/mmchipseq/mes_all/RNAPII_S2/tss/allchr.txt", header=FALSE, nrows=5)
classes <- sapply(pol2.raw, class)
pol2.raw = read.table("/arrayAhome/raflynn/mmchipseq/mes_all/RNAPII_S2/tss/allchr.txt", header=FALSE, colClasses=classes)
genenames = pol2.raw[,1]
pol2.raw = t(pol2.raw[,108:308])
colnames(pol2.raw) = genenames
pol2.raw = t(pol2.raw)
pol2.raw = pol2.raw[which(rownames(pol2.raw) %in% rownames(chirp.tss)),]
pol2.sc = log10(pol2.raw[which(rownames(chirp.tss) %in% rownames(pol2.raw)),]+1)

png("pol2.png", width=400, height=800)
pheatmap(pol2.sc, cluster_rows = FALSE, cluster_cols = FALSE, scale="none", show_rownames=FALSE, show_colnames=FALSE,col=colorRamp2(50), breaks=breaks)
dev.off()

### GRO AS

# gro heatmap
gro12c.as.raw = read.table("/arrayAhome/raflynn/7SK/GROseq/combined/GRO_12C_tss_antisense.txt", header=FALSE, nrows=5)
classes <- sapply(gro12c.as.raw, class)
gro12c.as.raw = read.table("/arrayAhome/raflynn/7SK/GROseq/combined/GRO_12C_tss_antisense.txt", header=FALSE, colClasses=classes)
genenames = gro12c.as.raw[,1]
gro12c.as.raw = t(gro12c.as.raw[,108:308])
colnames(gro12c.as.raw) = genenames
gro12c.as.raw = t(gro12c.as.raw)
gro12c.as.sc = log10(gro12c.as.raw[rownames(chirp.tss),]+1)

png("gro12c_tss_as.png", width=400, height=800)
pheatmap(gro12c.as.sc, cluster_rows = FALSE, cluster_cols = FALSE, scale="none", show_rownames=FALSE, show_colnames=FALSE,col=colorRamp(50))
dev.off()



# gro ASO heatmap
gro123.as.raw = read.table("/arrayAhome/raflynn/7SK/GROseq/combined/GRO_123_tss_antisense.txt", header=FALSE, nrows=5)
classes <- sapply(gro123.as.raw, class)
gro123.as.raw = read.table("/arrayAhome/raflynn/7SK/GROseq/combined/GRO_123_tss_antisense.txt", header=FALSE, colClasses=classes)
genenames = gro123.as.raw[,1]
gro123.as.raw = t(gro123.as.raw[,108:308])
colnames(gro123.as.raw) = genenames
gro123.as.raw = t(gro123.as.raw)
gro123.as.sc = log10(gro123.as.raw[rownames(chirp.tss),]+1)

# gro changes
gro12.as.chg = gro123.as.sc - gro12c.as.sc
breaks=c(seq(-2,-0.01,length.out=25),seq(0.01,2,length.out=25))
colorRamp2 = colorRampPalette(c("green4","green2","white","firebrick2","firebrick"))
png("gro_12_chg_as.png", width=400, height=800)
pheatmap(gro12.as.chg, cluster_rows = FALSE, cluster_cols = FALSE, scale="none", show_rownames=FALSE, show_colnames=FALSE,col=colorRamp2(50), breaks=breaks)
dev.off()
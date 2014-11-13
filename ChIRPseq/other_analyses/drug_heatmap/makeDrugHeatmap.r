alltss = as.matrix(read.delim("tss_avg7sk_drugs_matrix.txt", colClasses=c("character", rep("numeric", 4)), row.names=1, header=TRUE))


pdf("top400_scatterplot_withlines.pdf", width=6, height=6)
tss = alltss[1:4000,]
flavo = log2(tss[,3]/tss[,1])
jq1 = log2(tss[,4]/tss[,1])
actd = log2(tss[,2]/tss[,1])
plot(flavo, pch=19, cex=0.3, ylim=c(-4,4), col=rgb(0,1,0,alpha=0.2),  xlab="Drugs ranked by WT2 ChIRP density", ylab="log2 Treatment/WT 7SK read density", main="Treatmeant/WT 7SK read density for top 4000 genes")
points(jq1, pch=19, cex=0.3, col=rgb(0,0,1,alpha=0.2))
points(actd, pch=19, cex=0.3,  col=rgb(1,0,0,alpha=0.2))

segments(0, median(actd[1:200]), 200, median(actd[1:200]), col="red", lwd=5)
segments(0, median(flavo[1:200]), 200, median(flavo[1:200]), col="green", lwd=5)
segments(0, median(jq1[1:200]), 200, median(jq1[1:200]), col="blue", lwd=5)

segments(200, median(actd[200:1000]), 1000, median(actd[200:1000]), col="red", lwd=5)
segments(200, median(flavo[200:1000]), 1000, median(flavo[200:1000]), col="green", lwd=5)
segments(200, median(jq1[200:1000]), 1000, median(jq1[200:1000]), col="blue", lwd=5)


segments(1000, median(actd[1000:2000]), 2000, median(actd[1000:2000]), col="red", lwd=5)
segments(1000, median(flavo[1000:2000]), 2000, median(flavo[1000:2000]), col="green", lwd=5)
segments(1000, median(jq1[1000:2000]), 2000, median(jq1[1000:2000]), col="blue", lwd=5)


segments(2000, median(actd[2000:4000]), 4000, median(actd[2000:4000]), col="red", lwd=5)
segments(2000, median(flavo[2000:4000]), 4000, median(flavo[2000:4000]), col="green", lwd=5)
segments(2000, median(jq1[2000:4000]), 4000, median(jq1[2000:4000]), col="blue", lwd=5)

abline(h=0)
legend("topright", legend=c("ActD", "Flavo", "JQ1"), pch=19, col=c("red", "green", "blue"))
dev.off()

require(gplots)
require(RColorBrewer)

colors = c(seq(0, 2, length=50),seq(2, 6, length=50))
my_palette <- colorRampPalette(c("blue", "green", "yellow", "orange", "red"))(n = 99)

pdf("top1000_heatmap.pdf", width=8, height=8)
tss = alltss[1:1000,]
heatmap.2(log2(tss), scale="none", dendrogram="none", breaks=colors, col=my_palette, trace='none', Rowv=FALSE, Colv=FALSE)
dev.off()

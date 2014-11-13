# lookAtClusters.r
# 4/19/14
# Let's look at heatmaps and things in the three clusters

# Load data and process gene information
load("/home/raflynn/ChIRPseq/WT2_new/bins/tss/data.RData")
tss = t(data.extrap)

load("/home/raflynn/ChIRPseq/WT2_new/bins/tes/data.RData")
tes = t(data.extrap)


#load("tes_and_tss.RData") #has tss and tes data tables: 21674x400, with names
geneInfo = read.table("annotated_genes.txt", sep='\t', header=TRUE)
rownames(geneInfo) = geneInfo[,1]
geneInfo = geneInfo[,-1]
geneInfo$udRNA.lfc = as.numeric(as.character(geneInfo$udRNA.lfc))

# list of genes with failed termination
failedTerm = rownames(geneInfo[geneInfo[,2] ==1 ,])

# list of genes with expression changes, in order
orderedByExp = geneInfo[order(geneInfo[,1], decreasing=TRUE),]
upExp = rownames(orderedByExp[orderedByExp[,1] > 0,]) #sorted from most upreg to least up
downExp = rev(rownames(orderedByExp[orderedByExp[,1] < 0,])) #sorted from most downreg to least down

# list of genes with divergent transcription changes, in order
orderedByDiv = geneInfo[order(geneInfo[,3], decreasing=TRUE),]
upDiv = rownames(orderedByDiv[orderedByDiv[,3] > 0,]) #sorted from most upreg to least up
downDiv = rev(rownames(orderedByDiv[orderedByDiv[,3] < 0,])) #sorted from most downreg to least down

# all genes ordered by ChIRP level
tssSums = order(apply(tss, 1, sum), decreasing=TRUE)
tesSums = order(apply(tes, 1, sum), decreasing=TRUE)

# Draw mega boxplots
pdf("140707_7skrank_gbclusters_tss_boxplot.pdf", width=8, height=6)
par(mfrow=c(1,1))
boxplot(
21674-match(upDiv, rownames(tss[tssSums,])), 
21674-match(failedTerm, rownames(tss[tssSums,])),
21674-match(upDiv[which(upDiv %in% failedTerm)], rownames(tss[tssSums,])),
21674-match(upExp, rownames(tss[tssSums,])), 
21674-match(downExp, rownames(tss[tssSums,])),
names=c("Up Div", "Failed term", "Div+Term", "Up Expr", "Down Expr"), xlab="Class of genes", ylab="Rank of ChIRP signal", col=c("blue", "green", "seagreen", "red", "orange"), main="At the TSS")
dev.off()

# Draw mega boxplots
pdf("140707_7skrank_gbclusters_tes_boxplot.pdf", width=8, height=6)
par(mfrow=c(1,1))
boxplot(
21674-match(upDiv, rownames(tes[tesSums,])), 
21674-match(failedTerm, rownames(tes[tesSums,])),
21674-match(upDiv[which(upDiv %in% failedTerm)], rownames(tes[tesSums,])),
21674-match(upExp, rownames(tes[tesSums,])), 
21674-match(downExp, rownames(tes[tesSums,])),
names=c("Up Div", "Failed term", "Div+Term", "Up Expr", "Down Expr"), xlab="Class of genes", ylab="Rank of ChIRP signal", col=c("blue", "green", "seagreen", "red", "orange"), main="At the TES")
dev.off()

# old
# par(mfrow=c(1,1))
# boxplot(
# 21674-match(upDiv, rownames(tss[tssSums,])), 
# 21674-match(downDiv, rownames(tss[tssSums,])), 
# 21674-match(failedTerm, rownames(tss[tssSums,])),
# 21674-match(upExp, rownames(tss[tssSums,])), 
# 21674-match(downExp, rownames(tss[tssSums,])),
# names=c("Up Div", "Down Div", "Failed term", "Up Expr", "Down Expr"), xlab="Class of genes", ylab="Rank of ChIRP signal", col=c("blue", "blue", "green", "red", "red"), main="At the TSS")



# # Draw TSS and TES for divergent transcription
# par(mfrow=c(1,2))
# boxplot(21674-match(upDiv, rownames(tss[tssSums,])), 21674-match(downDiv, rownames(tss[tssSums,])), names=c("Up", "Down"), xlab="Direction of divergent transcription", ylab="Rank of ChIRP signal", main="TSS")
# boxplot(21674-match(upDiv, rownames(tes[tesSums,])), 21674-match(downDiv, rownames(tes[tesSums,])), names=c("Up", "Down"), xlab="Direction of divergent transcription", ylab="Rank of ChIRP signal", main="TES")

# ks.test(match(upDiv, rownames(tss[tssSums,])), sample(1:21674, length(upDiv))) #2.2e-16
# ks.test(match(downDiv, rownames(tss[tssSums,])), sample(1:21674, length(downDiv))) #1.2e-15
# ks.test(match(upDiv, rownames(tes[tesSums,])), sample(1:21674, length(upDiv))) #2.2e-16
# ks.test(match(downDiv, rownames(tes[tesSums,])), sample(1:21674, length(downDiv))) #9e-3

# # Draw TSS and TES for expression
# par(mfrow=c(1,2))
# boxplot(match(upExp, rownames(tss[tssSums,])), match(downExp, rownames(tss[tssSums,])), names=c("Up", "Down"), xlab="Direction of expression change", ylab="Rank of ChIRP signal", main="TSS")
# boxplot(match(upExp, rownames(tes[tesSums,])), match(downExp, rownames(tes[tesSums,])), names=c("Up", "Down"), xlab="Direction of expression change", ylab="Rank of ChIRP signal", main="TES")

# ks.test(match(upExp, rownames(tss[tssSums,])), 1:21674) #2.2e-16
# ks.test(match(downExp, rownames(tss[tssSums,])), 1:21674) #1.9e-2
# ks.test(match(upExp, rownames(tes[tesSums,])), 1:21674) #2.2e-16
# ks.test(match(downExp, rownames(tes[tesSums,])), 1:21674) #3.6e-4

# # Draw TSS and TES for failed termination
# par(mfrow=c(1,1))
# boxplot(match(failedTerm, rownames(tss[tssSums,])), match(failedTerm, rownames(tes[tesSums,])), names=c("TSS", "TES"), xlab="Region", ylab="Rank of ChIRP signal", main="Failed termination")
# # boxplot(match(failedTerm[which(failedTerm %in% upDiv)], rownames(tss[tssSums,])), match(failedTerm[which(failedTerm %in% upDiv)], rownames(tes[tesSums,])), names=c("TSS", "TES"), xlab="Region", ylab="Rank of ChIRP signal", main="Failed termination and div tx")

# ks.test(match(failedTerm, rownames(tss[tssSums,])), 1:21674) #2.2e-16
# ks.test(match(failedTerm, rownames(tes[tesSums,])), 1:21674) #2.2e-16

# Plotting metagenes for divergent transcription
pdf("140707_7skmetagene_divtx.pdf", width=10, height=5.5)
par(mfrow=c(1,2))
plot(apply(tss[which(rownames(tss) %in% upDiv),], 2, function(x) mean(x, na.rm=TRUE)), type='l', lwd=2, col="red",ylim=c(0,3), main="TSS for divergent transcription", xlab="Distance relative to TSS", ylab="Reads per nucleotide", xaxt='n')
#lines(apply(tss[which(rownames(tss) %in% downDiv),], 2, function(x) mean(x, na.rm=TRUE)), type='l', lwd=2, col="lightblue", ylim=c(0,3))
lines(apply(tss, 2, function(x) mean(x, na.rm=TRUE)), type='l', lwd=2, col="black", ylim=c(0,3))
#legend("bottomright", legend=c("All genes", "More divergent transcription", "Less divergent transcription"), col=c("black", "red", "lightblue"), lty=rep(1,3), lwd=rep(2,3))
legend("bottomright", legend=c("All genes", "Divergently transcribed genes"), col=c("black", "red"), lty=rep(1,2), lwd=rep(2,2))
axis(1, at = seq(0, 400, 100), labels=seq(-1000, 1000, 500))

plot(apply(tes[which(rownames(tes) %in% upDiv),], 2, function(x) mean(x, na.rm=TRUE)), type='l', lwd=2, col="red",ylim=c(0,3), main="TES for divergent transcription", xlab="Distance relative to TES", ylab="Reads per nucleotide", xaxt='n')
#lines(apply(tes[which(rownames(tes) %in% downDiv),], 2, function(x) mean(x, na.rm=TRUE)), type='l', lwd=2, col="lightblue", ylim=c(0,3))
lines(apply(tes, 2, function(x) mean(x, na.rm=TRUE)), type='l', lwd=2, col="black", ylim=c(0,3))
legend("bottomright", legend=c("All genes", "Divergently transcribed genes"), col=c("black", "red"), lty=rep(1,2), lwd=rep(2,2))
# legend("bottomright", legend=c("All genes", "More divergent transcription", "Less divergent transcription"), col=c("black", "red", "lightblue"), lty=rep(1,3), lwd=rep(2,3))
axis(1, at = seq(0, 400, 100), labels=seq(-1000, 1000, 500))
dev.off()

# plot(apply(tes[which(rownames(tes) %in% upDiv[which(upDiv %in% failedTerm)]),], 2, function(x) mean(x, na.rm=TRUE)), type='l', lwd=2, col="red",ylim=c(0,3), main="TES for divergent transcription and failed termination", xlab="Position", ylab="Reads per nucleotide")

# Plotting metagenes for expression changes
pdf("140707_7skmetagene_exprchanges.pdf", width=10, height=5.5)
par(mfrow=c(1,2))
plot(apply(tss[which(rownames(tss) %in% upExp),], 2, function(x) mean(x, na.rm=TRUE)), type='l', lwd=2, col="red",ylim=c(0,3), main="TSS for expression", xlab="Position relative to TSS", ylab="Reads per nucleotide", xaxt='n')
lines(apply(tss[which(rownames(tss) %in% downExp),], 2, function(x) mean(x, na.rm=TRUE)), type='l', lwd=2, col="lightblue", ylim=c(0,3))
lines(apply(tss, 2, function(x) mean(x, na.rm=TRUE)), type='l', lwd=2, col="black", ylim=c(0,3))
legend("topleft", legend=c("All genes", "Higher expression", "Lower expression"), col=c("black", "red", "lightblue"), lty=rep(1,3), lwd=rep(2,3))
axis(1, at = seq(0, 400, 100), labels=seq(-1000, 1000, 500))


plot(apply(tes[which(rownames(tes) %in% upExp),], 2, function(x) mean(x, na.rm=TRUE)), type='l', lwd=2, col="red",ylim=c(0,3), main="TES for expression", xlab="Position relative to TES", ylab="Reads per nucleotide", xaxt='n')
lines(apply(tes[which(rownames(tes) %in% downExp),], 2, function(x) mean(x, na.rm=TRUE)), type='l', lwd=2, col="lightblue", ylim=c(0,3))
lines(apply(tes, 2, function(x) mean(x, na.rm=TRUE)), type='l', lwd=2, col="black", ylim=c(0,3))
legend("topright", legend=c("All genes", "Higher expression", "Lower expression"), col=c("black", "red", "lightblue"), lty=rep(1,3), lwd=rep(2,3))
axis(1, at = seq(0, 400, 100), labels=seq(-1000, 1000, 500))

dev.off()

# Plotting metagenes for failed termination
pdf("140707_7skmetagene_failedterm.pdf", width=10, height=5.5)

par(mfrow=c(1,2))
plot(apply(tss[which(rownames(tss) %in% failedTerm),], 2, function(x) mean(x, na.rm=TRUE)), type='l', lwd=2, col="red",ylim=c(0,4), main="TSS for failed termination", xlab="Position relative to TSS", ylab="Reads per nucleotide", xaxt='n')
lines(apply(tss, 2, function(x) mean(x, na.rm=TRUE)), type='l', lwd=2, col="black", ylim=c(0,4))
legend("bottomright", legend=c("All genes", "Failed termination"), col=c("black", "red"), lty=rep(1,2), lwd=rep(2,2))
axis(1, at = seq(0, 400, 100), labels=seq(-1000, 1000, 500))

plot(apply(tes[which(rownames(tes) %in% failedTerm),], 2, function(x) mean(x, na.rm=TRUE)), type='l', lwd=2, col="red",ylim=c(0,3), main="TES for failed termination", xlab="Position relative to TES", ylab="Reads per nucleotide", xaxt='n')
lines(apply(tes, 2, function(x) mean(x, na.rm=TRUE)), type='l', lwd=2, col="black", ylim=c(0,3))
legend("bottomright", legend=c("All genes", "Failed termination"), col=c("black", "red"), lty=rep(1,2), lwd=rep(2,2))
axis(1, at = seq(0, 400, 100), labels=seq(-1000, 1000, 500))

dev.off()

# Plotting heatmaps

# Make heatmap
drawHeatmap = function(df, fnName, ...) {
	require(RColorBrewer)
	colors=colorRampPalette(c("white", "blue", "darkblue"))(100)  
	
	jpeg(fnName, ...)
	layout(matrix(c(1), 1, 1), widths=c(5), heights = c(10), respect=T)
	par(mar=c(1,1,1,1))

	df = df[order(apply(df, 1, sum), decreasing=TRUE),]
	data = log10(df + 1)
	heatmap(data, Rowv=NA, Colv="Rowv", labRow = NA, labCol = NA, col=colors, scale="none")
	dev.off()
}

drawHeatmap(tss, "140707_7skheatmap_alltss.jpg", width=1000, height=2000)
drawHeatmap(tes, "140707_7sk_alltes_heatmap.jpg", width=1000, height=2000)
# drawHeatmap(tss[which(rownames(tss) %in% upDiv),])
# drawHeatmap(tss[match(upDiv, rownames(tss)),])


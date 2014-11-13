# makeRateHeatmap.r
# 5/25/14
# make heatmap of gene bodies based on rate

load("data.RData")
rates = read.delim("/home/raflynn/ChIRPseq/groseq/gene_rates_with_7sk.txt", sep='\t', header=TRUE)
rates = rates[rates[,2]>0,]
rownames(rates) = rates[,1]
ordGenes = rownames(rates)[order(rates[,2], decreasing=TRUE)]

data.extrap.log = log10(data.extrap + 1)

colors=colorRampPalette(c("white", "blue", "darkblue"))(100)  
heatmap(t(data.extrap.log[,ordGenes[1:5]]), Rowv=NA, Colv="Rowv", labCol = NA, col=colors, scale="none")

data.coll = apply(data.extrap, 2, function(x) approx(1:length(x), x, n=50)$y)
data.coll.log = log10(data.coll + 1)

heatmap(t(data.coll.log[,ordGenes]), Rowv=NA, Colv="Rowv", labRow = NA, labCol = NA, col=colors, scale="none")

# for tss
heatmap(t(data.extrap.log[150:250,ordGenes]), Rowv=NA, Colv="Rowv", labRow = NA, labCol = NA, col=colors, scale="none")

left = apply(data.extrap.log[185:200,], 2, mean)
leftDensity = left[ordGenes]
leftOrder = order(leftDensity)

right = apply(data.extrap.log[200:215,], 2, mean)
rightDensity = right[ordGenes]
rightOrder = order(rightDensity)
plot(rightDensity - leftDensity, rates[ordGenes,2]/1000)
plot(rates[ordGenes,2], log(middleDensity))
cor.test(rates[ordGenes,2], middleDensity, method="spearman")

par(mfrow=c(1,2))
heatmap(t(data.extrap.log[175:225, ordGenes]), Rowv=NA, Colv="Rowv", labRow = NA, labCol = NA, col=colors, scale="none")
barplot(abs(rates[rev(ordGenes),2]), horiz=TRUE)


topRates = ordGenes[1:100]
bottomRates = ordGenes[length(ordGenes)-100,length(ordGenes)]
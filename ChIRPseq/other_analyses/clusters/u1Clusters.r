# u1Clusters.r
# May 2014
# different ways to understand U1/PAS and its relationship to the set of genes
# that undergo failed termination after 7SK knockdown in mESCs

data = read.table("gene_classes_new_50.txt", sep='\t', header=TRUE)
rownames(data) = data[,1]
data.short = data[data$length<30000,]

data = data[,-1]
data.prob = data[(data$class!="N" & data$class !="B"),]

# require(MASS)
# data.prob$up_ratio = log(data.prob$up_ratio)
# train = sample(1:nrow(data.prob), 1500, replace=FALSE)

# model = glm(class ~ length, data=data.prob, subset=train, family="binomial")
# data.pred = predict(model, data.prob[-train,], type="response")
# data.pred=ifelse(data.pred>0.5, "D", "T")
# table(data.pred, "obs" = data.prob[-train,]$class)

# mean(data.prob[-train,]$class == data.pred)


# u1 and pas files; can use gene_up_u1.txt and gene_up_pas.txt for upstream things
require(parallel)
no_col = count.fields("gene_u1.txt", sep = "\t") 
u1.data = read.table("gene_u1.txt", sep='\t', fill=TRUE, col.names=1:max(no_col))

no_col = count.fields("gene_pas.txt", sep = "\t") 
pas.data = read.table("gene_pas.txt", sep='\t', fill=TRUE, col.names=1:max(no_col))

rownames(u1.data) = u1.data[,1]
u1.data = u1.data[,-1]

rownames(pas.data) = pas.data[,1]
pas.data = pas.data[,-1]

# u1 and pas fill in as 'na'
makeBins = function(x) {
	x=na.omit(as.numeric(x))
	bins = rep(0, 20)
	numSites = length(x)
	if (numSites == 0) {return(bins)}
	
	fakeHist = hist(x, breaks=c((0:20)/20))
	return(fakeHist$counts/numSites)	
}

# bin making
cl = makeCluster(10)
u1.bins = t(parApply(cl, u1.data[which(rownames(u1.data) %in% rownames(data.prob)),], 1, makeBins))
u1.bins.all = t(parApply(cl, u1.data, 1, makeBins))
u1.bins.cumul = t(parApply( cl, u1.bins, 1, cumsum))
u1.bins.all.cumul = t(parApply(cl, u1.bins.all, 1, cumsum))

pas.bins = t(parApply(cl, pas.data[which(rownames(pas.data) %in% rownames(data.prob)),], 1, makeBins))
pas.bins.all = t(parApply(cl, pas.data, 1, makeBins))
pas.bins.cumul = t(parApply(cl, pas.bins, 1, cumsum))
pas.bins.all.cumul = t(parApply(cl, pas.bins.all, 1, cumsum))

ratio.bins.cumul = u1.bins.cumul - pas.bins.cumul
ratio.bins.all.cumul = u1.bins.all.cumul - pas.bins.all.cumul





# divtx has longer length
boxplot(data.prob[data.prob$class=="D",]$length, data.prob[data.prob$class=="T",]$length, ylim=c(0,100000))

# upstream things don't seem to matter
boxplot(data.prob[data.prob$class=="D",]$up_ratio, data.prob[data.prob$class=="T",]$up_ratio, ylim=c(0,5))

# gene body things do though
boxplot(data.prob[data.prob$class=="D",]$gene_pasrate, data.prob[data.prob$class=="T",]$gene_pasrate, ylim=c(0,5), ylab="U1/PAS ratio", main="U1/PAS ratio")
axis(1, at=c(1,2), labels=c("Divergent\ntranscription", c("Termination\ndefect")))

boxplot(data[data$length>20000,]$gene_ratio, data[data$length<20000,]$gene_ratio, ylim=c(0,5), ylab="U1/PAS ratio", main="U1/PAS ratio", xlab="Length of gene")
axis(1, at=c(1,2), labels=c("> 20000", c("< 20000")))

# and so does length
boxplot(log(data.prob[data.prob$class=="D",]$length), log(data.prob[data.prob$class=="T",]$length), ylim=c(0,15), ylab="Log of length", main="Gene length")
axis(1, at=c(1,2), labels=c("Divergent\ntranscription", c("Termination\ndefect")))

# 5/14: looking at all classes
boxplot(log(data.short[data.short$class=="T",]$length), log(data.short[data.short$class=="D",]$length), log(data.short[data.short$class=="N",]$length), ylim=c(0,15), ylab="Log of length", main="Gene length")

boxplot(data.short[data.short$class=="T",]$length, data.short[data.short$class=="D",]$length, data.short[data.short$class=="N",]$length, ylim=c(0,30000), ylab="Log of length", main="Gene length")

boxplot(data.short[data.short$class=="T",]$gene_u1rate, data.short[data.short$class=="D",]$gene_u1rate, data.short[data.short$class=="N",]$gene_u1rate, ylim=c(0,5), ylab="Log of length", main="Gene length")

boxplot(data[data$class=="T",]$gene_ratio, data[data$class=="D",]$gene_ratio, data[data$class=="N",]$gene_ratio, ylim=c(0,5), ylab="Log of length", main="Gene length")

# by length
boxplot(data[data$length<30000,]$gene_pasrate, data[data$length>30000,]$gene_pasrate, ylim=c(0,5), ylab="Log of length", main="Gene length")

boxplot(data[data$length<30000,]$gene_ratio, data[data$length>30000,]$gene_ratio, ylim=c(0,5), ylab="Log of length", main="Gene length")

# looking only at short genes
boxplot(data.short[data.short$class=="T",]$gene_u1rate, data.short[data.short$class!="T",]$gene_u1rate, ylim=c(0,5), ylab="U1/PAS ratio", main="U1/PAS ratio")
axis(1, at=c(1,2), labels=c("Divergent\ntranscription", c("Termination\ndefect")))







rocCurve = function(data.byRatio) {
	tp = cumsum(rev(data.byRatio)/sum(data.byRatio))
	n = nrow(data)
	fp = cumsum((n/10 - rev(data.byRatio))/(n-sum(data.byRatio)))
	roc = cbind(fp, tp)
	return(rbind(c(0,0),roc))
	
}

# pdf("u1_pas_ratio_termdefects.pdf", width=12, height=5)
# par(mfrow=c(1,3))
# prediction: u1/pas --> 'T'
data.byRatio = data[order(data$gene_ratio),]
numPerBin = floor(nrow(data.byRatio)/10)

nums = rep(0, 10)
for(i in 0:9) {
	nums[i+1] = sum(data.byRatio[(i * numPerBin):((i+1)*numPerBin),]$class=="T")
}
allROC = rocCurve(nums)
plot(nums/sum(nums), type='b', pch=19, lwd=2, ylim=c(0,0.2),xlab="Decile of U1/PAS ratio", ylab="Percentage of failed term genes in bin", main="Proportion of failed termination genes \nas a function of U1/PAS ratio",xaxt='n')
axis(1, at=c(1:10),labels=c(1:10))





# barplot(nums/numPerBin*100, ylim=c(0, 10), names.arg=1:10, main="Gene-wide U1/PAS ratio", xlab="Bins of genes with increasing gene-wide U1/PAS ratio", ylab="Percentage of genes with termination defects")
# prediction: 5'
data.byRatio = data[order(data$X5_ratio),]
numPerBin = floor(nrow(data.byRatio)/10)

nums = rep(0, 10)
for(i in 0:9) {
	nums[i+1] = sum(data.byRatio[(i * numPerBin):((i+1)*numPerBin),]$class=="T")
}
fiveROC = rocCurve(nums)
plot(nums/sum(nums), type='b', pch=19, lwd=2, ylim=c(0,0.2),col='blue',xlab="Decile of U1/PAS ratio", ylab="Percentage of failed term genes in bin", main="Proportion of failed termination genes \nas a function of U1/PAS ratio",xaxt='n')

# barplot(nums/numPerBin*100, xlab="Bins of genes with increasing 5' U1/PAS ratio", main=c("5' U1/PAS ratio"), names.arg=1:10, ylim=c(0,10))

# prediction: 3'
data.byRatio = data[order(data$X3_ratio),]
numPerBin = floor(nrow(data.byRatio)/10)

nums = rep(0, 10)
for(i in 0:9) {
	nums[i+1] = sum(data.byRatio[(i * numPerBin):((i+1)*numPerBin),]$class=="T")
}
threeROC = rocCurve(nums)
lines(nums/sum(nums), type='b', pch=19, lwd=2, col='red')
legend("topleft", legend=c("First 25% of gene (5' end)", "Last 25% of gene (3' end)"), col=c("blue", "red"), lty=c(1,1), lwd=c(2,2))

axis(1, at=c(1:10),labels=c(1:10))



# barplot(nums/numPerBin*100, xlab="Bins of genes with increasing 3' U1/PAS ratio", main=c("3' U1/PAS ratio"), names.arg=1:10, ylim=c(0,10))
#dev.off()

plot(allROC, type='b', pch=19, col="black", xlim=c(0,1), ylim=c(0,1))
points(fiveROC, type='b', pch=19, col="blue")
points(threeROC, type='b', pch=19, col="red")



# plotting bins
plot(data$length, pas.bins.all[which(rownames(pas.bins.all) %in% rownames(data)),20])

# U1
div = parApply(cl, u1.bins[which(rownames(u1.bins) %in% rownames(data[data.prob$class=="D",])),], 2, mean)
plot(div/sum(div), type='l', pch=19, lwd=2, ylim=c(0,0.08), xaxt='n', xlab=NA, ylab="Average U1 site density", main="Average distribution of U1 sites along gene")

term = parApply(cl, u1.bins[which(rownames(u1.bins) %in% rownames(data[data$class=="T",])),], 2, mean) 
lines(term/sum(term), type='l', pch=19, lwd=2, col="red")

none = parApply(cl, u1.bins.all[which(rownames(u1.bins.all) %in% rownames(data[data$class=="N",])),], 2, mean)
lines(none/sum(none), type='l', pch=19, lwd=2, col="blue")

legend("topright", legend=c("Div tx", "Term defects", "Unaffected"), col=c("black", "red", "blue"), lty=c(1,1), lwd=c(2,2))
axis(1, at=c(1,20), labels=c("TSS", "TES"))

# PAS
div = parApply(cl, pas.bins[which(rownames(pas.bins) %in% rownames(data[data$class=="D",])),], 2, mean)
plot(div/sum(div), type='l', pch=19, lwd=2, ylim=c(0,0.3), xaxt='n', xlab=NA, ylab="Average PAS site density", main="Average distribution of PAS sites along gene")

term = parApply(cl, pas.bins[which(rownames(pas.bins) %in% rownames(data[data$class=="T",])),], 2, mean)
lines(term/sum(term), type='l', lwd=2, pch=19, col="red")

none = parApply(cl, pas.bins.all[which(rownames(pas.bins.all) %in% rownames(data[data$class=="N",])),], 2, mean)
lines(none/sum(none), type='l', lwd=2, pch=19, col="blue")

legend("topleft", legend=c("Div tx", "Term defects", "Unaffected"), col=c("black", "red", "blue"), lty=c(1,1), lwd=c(2,2))
axis(1, at=c(1,20), labels=c("TSS", "TES"))

# Ratios
is.na(ratio.bins.cumul) = sapply(ratio.bins.cumul, is.infinite)
is.na(ratio.bins.all.cumul) = sapply(ratio.bins.all.cumul, is.infinite)

div = parApply(cl, ratio.bins.cumul[which(rownames(ratio.bins.cumul) %in% rownames(data.short[data.short$class=="D",])),], 2, function(x) mean(x, na.rm=TRUE))
plot(div, type='l', pch=19, lwd=2, ylim=c(0,2), xaxt='n', xlab=NA, ylab="Average PAS site density", main="Average distribution of PAS sites along gene")

term = parApply(cl, ratio.bins.cumul[which(rownames(ratio.bins.cumul) %in% rownames(data.short[data.short$class=="T",])),], 2, function(x) mean(x, na.rm=TRUE))
lines(term, type='l', lwd=2, pch=19, col="red")

none = parApply(cl, ratio.bins.all.cumul[which(rownames(ratio.bins.all.cumul) %in% rownames(data.short[data.short$class=="N",])),], 2, function(x) mean(x, na.rm=TRUE))
lines(none, type='l', lwd=2, pch=19, col="blue")

legend("topleft", legend=c("Div tx", "Term defects", "Unaffected"), col=c("black", "red", "blue"), lty=c(1,1), lwd=c(2,2))
axis(1, at=c(1,20), labels=c("TSS", "TES"))

#scratch  
divless = parApply(cl, pas.bins.all[which(rownames(pas.bins.all) %in% rownames(data[data$length<20000,])),], 2, mean)
plot(divless/sum(divless), type='l', pch=19, lwd=2, ylim=c(0,0.3), xaxt='n', xlab=NA, ylab="Average PAS site density", main="Average density of PAS sites across gene")
divmore = parApply(cl, pas.bins.all[which(rownames(pas.bins.all) %in% rownames(data[data$length>20000,])),], 2, mean)
lines(divmore/sum(divmore), type='l', lwd=2, ylim=c(0,0.3), col="red")
legend("topleft", legend=c("< 20000 nt", "> 20000 nt"), col=c("black", "red"), lty=c(1,1), lwd=c(2,2))
axis(1, at=c(1,20), labels=c("start", "end"))


stopCluster(cl)

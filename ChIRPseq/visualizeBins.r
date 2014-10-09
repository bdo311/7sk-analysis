# File: visualizeBins.r
# 2/28/14
# Takes in bin files and visualizes average for top/bottom

library(data.table)

args = commandArgs(TRUE)
exprFile = args[1]




# for tss/tes
data = t(fread(exprFile, sep='\t'))
colnames(data) = data[1,]
data = data[8:nrow(data),]

data = apply(data, c(1,2), as.numeric)
data = apply(data, 2, function(x) x/max(x)) #scales

plot(apply(data[,1:1000], 1, mean), type='l', lwd=2, ylab = "Average reads per nt", main="Top 1000", col="red")

# for genebody
extrap = function(x) {
  num = sum(!is.na(x))
  approx(1:num, x[1:num], n=1000)$y
}

no_col = count.fields(exprFile, sep = "\t") 
data = t(read.table(exprFile, sep='\t', fill=TRUE, col.names=1:max(no_col))) #can't use fread

colnames(data) = data[1,]
data = data[8:nrow(data),]

data.extrap = apply(data, 2, extrap)
data.extrap = apply(data.extrap, 2, function(x) x/max(x))

plot(apply(data[,1:1000], 1, mean), type='l', lwd=2, ylab = "Average reads per nt", main="Top 1000", col="red")

# shrinking list so that instead of 1000 bins it's 100
shrinkList = function(x) {
	lenPerInt = length(x)/100
	newList = rep(0, 100)

	for(i in 1:100) {
		newList[i] = mean(x[i:i+lenPerInt-1])
	}
	return(newList)
}

library(parallel)
clus = makeCluster(getOption("cl.cores", 10))

data.new = parApply(cl=clus, X=data.extrap, MARGIN=2, FUN=shrinkList)
colnames(data.new) = colnames(data.extrap)

data.new.t = t(data.new)
kobj = kmeans(data.new.t[1:20000,], centers=3)
plot(kobj$centers[1,], type="l", lwd=2)
table(kobj$cluster)

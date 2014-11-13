# File: visualizeBins_new.r
# 3/5/14
# In a given folder, produces plot files and saves data into an RData

library(data.table)

args = commandArgs(TRUE)
typeOfFile = args[1]
folderName = args[2]

# load("data.RData")

# Extrapolate bins to exactly 1000
extrap = function(x) {
  num = sum(!is.na(x))
  approx(1:num, x[1:num], n=1000)$y
}

# Gene file requires extrapolation and fread cannot be used. TSS/TES files can
if(typeOfFile == "gene") {
	no_col = count.fields("allchr_sorted.txt", sep = "\t") 
	data = t(read.table("allchr_sorted.txt", sep='\t', fill=TRUE,
	 col.names=1:max(no_col))) #can't use fread

	colnames(data) = data[1,]
	data = data[8:nrow(data),]

	data.extrap = apply(data, 2, extrap) # coerce to exactly 1000 bins
} else {
	data = t(fread("allchr_sorted.txt", sep='\t'))
	colnames(data) = data[1,]
	data = data[8:nrow(data),]

	data.extrap = apply(data, c(1,2), as.numeric) # doesn't need extrapolation, but needs conversion to numeric
}

# Scale and collapse data, and print averages to files
# data.scaled = apply(data.extrap, 2, function(x) x/max(x))
# data.collapsed = apply(data.scaled[,1:19000], 1, function(x) mean(x, na.rm=TRUE))
# write.table(data.collapsed, paste('avgscaled_', typeOfFile, '_', folderName, ".txt", sep=''), sep='\t')
# data.extrap = data.extrap[,-(which(colnames(data.extrap) %in% c("C920006O11Rik")))]
data.avg = apply(data.extrap[,50:ncol(data.extrap)], 1, function(x) mean(x, na.rm=TRUE))
write.table(data.avg, paste('avgraw_', typeOfFile, '_', folderName, ".txt", sep=''), sep='\t')

# Make plots of the scaled average data
pdf("plot.pdf")
plot(data.avg, type='l', lwd=2, 
	ylab = "Average reads per nt", main="All genes", col="red")
dev.off()

save(data, data.extrap, data.avg, file="data.RData")

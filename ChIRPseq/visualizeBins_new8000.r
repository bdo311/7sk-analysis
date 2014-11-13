# File: visualizeBins_new8000.r
# 4/15/14
# In a given folder, produces plot files and saves data into an RData for the top 8000 genes only

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
	no_col = count.fields("allchr_sorted_byWT2.txt", sep = "\t") 
	data = t(read.table("allchr_sorted_byWT2.txt", sep='\t', fill=TRUE,
	 col.names=1:max(no_col))) #can't use fread

	colnames(data) = data[1,]
	data = data[8:nrow(data),]

	data.extrap = apply(data, 2, extrap) # coerce to exactly 1000 bins
} else {
	data = t(fread("allchr_sorted_byWT2.txt", sep='\t'))
	colnames(data) = data[1,]
	data = data[8:nrow(data),]

	data.extrap = apply(data, c(1,2), as.numeric) # doesn't need extrapolation, but needs conversion to numeric
}

# Scale and collapse data, and print averages to files
data.scaled = apply(data.extrap, 2, function(x) x/max(x))
data.collapsed = apply(data.scaled[,1:7999], 1, function(x) mean(x, na.rm=TRUE))
write.table(data.collapsed, paste('avgscaled_8000_', typeOfFile, '_', folderName, ".txt", sep=''), sep='\t')
write.table(apply(data.extrap[,1:7999], 1, function(x) mean(x, na.rm=TRUE)), paste('avgraw_8000', typeOfFile, '_', folderName, ".txt", sep=''), sep='\t')

# Make plots of the scaled average data
pdf("plot_8000.pdf")
plot(data.collapsed, type='l', lwd=2, 
	ylab = "Average reads per nt", main="All genes", col="red")
dev.off()

save(data, data.extrap, data.scaled, file="data_8000.RData")

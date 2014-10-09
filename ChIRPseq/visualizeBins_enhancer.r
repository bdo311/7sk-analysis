# File: visualizeBins_enhancer.r
# 4/8/14, adapted from visualizeBins_new.r
# In a given folder, produces plot files and saves data into an RData

library(data.table)

args = commandArgs(TRUE)
typeOfFile = args[1]
folderName = args[2]
start = strtoi(args[3])


# load("data.RData")

# Extrapolate bins to exactly 400
extrap = function(x) {
  num = sum(!is.na(x))
  approx(1:num, x[1:num], n=400)$y
}

# Super file requires extrapolation
if(typeOfFile == "super") {
	no_col = count.fields("allchr_sorted.txt", sep = "\t") 
	data = t(read.table("allchr_sorted.txt", sep='\t', fill=TRUE,
	 col.names=1:max(no_col))) #can't use fread

	colnames(data) = data[1,]
	data = data[start:nrow(data),]

	data.extrap = apply(data, 2, extrap) # coerce to exactly 400 bins
} else {
	data = t(fread("allchr_sorted.txt", sep='\t'))
	colnames(data) = data[1,]
	data = data[start:nrow(data),]

	data.extrap = apply(data, c(1,2), as.numeric) # doesn't need extrapolation, but needs conversion to numeric
}

# Scale and collapse data, and print averages to files
if(typeOfFile == "reg") {
	enh = "non super enhancers"
} else {
	enh = "super enhancers"
}
# data.scaled = apply(data.extrap, 2, function(x) x/max(x))
# data.collapsed = apply(data.scaled, 1, function(x) mean(x, na.rm=TRUE))
# write.table(data.collapsed, paste('avgscaled_', typeOfFile, '_', folderName, ".txt", sep=''), sep='\t')
# `%ni%` = Negate(`%in%`)
# data.extrap = data.extrap[,which(colnames(data.extrap) %ni% c("INT_STITCHED_4655"))]
data.avg = apply(data.extrap, 1, function(x) mean(x, na.rm=TRUE))
write.table(data.avg, paste('avgraw_', folderName, ".txt", sep=''), sep='\t')

# Make plots of the scaled average data
pdf("plot.pdf")
plot(data.avg, type='l', lwd=2, 
	ylab = "Average reads per nt", main=paste("All ", enh, sep=''), col="red")
dev.off()

save(data, data.extrap, data.avg, file="data.RData")

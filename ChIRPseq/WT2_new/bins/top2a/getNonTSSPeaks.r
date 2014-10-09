# getNonTSSPeaks.r
# 8/9/14
# compares nonTSS and nongenic Top2A peak metagene to peak metagene across all peaks

nontss = read.delim("/home/raflynn/mmchipseq/top2/Top2a_provided_peaks_notss.bed", colClasses=rep("character",5), header=FALSE)
nongenic = read.delim("/home/raflynn/mmchipseq/top2/Top2a_provided_peaks_nogenic.bed", colClasses=rep("character",5), header=FALSE)
nontss = nontss$V4
nongenic = nongenic$V4

# Extrapolate bins to exactly 1000
extrap = function(x) {
  num = sum(!is.na(x))
  approx(1:num, x[1:num], n=1000)$y
}

no_col = count.fields("allchr_sorted.txt", sep = "\t") 
data = t(read.table("allchr_sorted.txt", sep='\t', fill=TRUE,
 col.names=1:max(no_col))) #can't use fread

colnames(data) = data[4,]
data = data[8:nrow(data),]

data.extrap = apply(data, 2, extrap) # coerce to exactly 1000 bins

data.nontss = data.extrap[,which(colnames(data.extrap) %in% nontss)]
data.nongenic = data.extrap[,which(colnames(data.extrap) %in% nongenic)]

data.avg.all = apply(data.extrap, 1, function(x) mean(x, na.rm=TRUE))
data.avg.nt = apply(data.nontss, 1, function(x) mean(x, na.rm=TRUE))
data.avg.ng = apply(data.nongenic, 1, function(x) mean(x, na.rm=TRUE))

df = data.frame("All"=data.avg.all, "Non-TSS" = data.avg.nt, "Non-genic" = data.avg.ng)
write.table(df, 'top2a_peak_densities', row.names=FALSE, sep='\t')

# Make plots of the scaled average data
pdf("plot.pdf")
plot(data.avg, type='l', lwd=2, 
	ylab = "Average reads per nt", main="All genes", col="red")
dev.off()

pdf("plot_nontss_nongenic.pdf", width=6, height=6)
plot(data.avg.all, type='l', ylim=c(0,1.2), lwd=2, ylab = "Average reads per nt", main="7SK read density at TOP2A peaks", col="black", xaxt='n', xlab=NA)
lines(data.avg.nt, lwd=2, col="blue")
lines(data.avg.ng, lwd=2, col="cyan")
legend("topright", legend=c("All peaks", "Non-TSS peaks", "Non-genic peaks"), col=c("black", "blue", "cyan"), lty=rep(1,3), lwd=rep(2,3))
axis(1, at=c(333, 666), labels=c("Peak start", "Peak end"))
dev.off()




dte = apply(data.extrap[333:666,], 2, mean)
dtg = apply(data.nongenic[333:666,], 2, mean)
dts = apply(data.nontss[333:666,], 2, mean)
pdf("plot_nontss_nongenic_density.pdf", width=6, height=6)

plot(density(dtg), xlim=c(0,3), col="cyan", lwd=2, xlab="Mean 7SK density across Top2A peak", ylab="Density", main="Distribution of 7SK density at Top2A peaks")
lines(density(dte), col="black", lwd=2)
lines(density(dts), col="blue", lwd=2)
legend("topright", legend=c("All peaks", "Non-TSS peaks", "Non-genic peaks"), col=c("black", "blue", "cyan"), lty=rep(1,3), lwd=rep(2,3))

abline(v=0.3, lty=2)
dev.off()


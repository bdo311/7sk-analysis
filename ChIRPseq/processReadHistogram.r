# processReadHistogram.r
# 3/22/14
# Makes a cumulative histogram of reads per bin to total reads up to that point

library(data.table)

data = fread("/home/raflynn/7SK_ChIRPseq/WT2/odd_hist.txt")

cumuldata = cbind(data$V1, rep(0, nrow(data)))

cumul = 0
for (i in 1:nrow(data)) {
	cumul = cumul + data$V1[i] * data$V2[i]
	cumuldata[i,2] = cumul
}

pdf("/home/raflynn/7SK_ChIRPseq/WT2/odd_cumul_100.pdf")
write.table(cumuldata, '/home/raflynn/7SK_ChIRPseq/WT2/odd_hist_cumul_100.txt', sep='\t')
plot(cumuldata, type='l',xlim=c(0,100), xlab="Reads per base", ylab="Cumulative number of reads")
dev.off()

	
normalize = function(x) {
z = x - min(x)
z/max(z)
}

data = read.delim("/home/raflynn/Mitosis/NETseq/metagenes/cNET_HeLa_Rep1_sense/bins/Churchman_CTCF_all/allchr_sorted.txt", header=FALSE, colClasses=c(rep("character",6),rep("numeric",401)))
data2 = data[,128:287]
data3 = apply(data2, 1, normalize)
data3[is.nan(data3)] = 0
data3 = t(data3)

data4 = read.delim("/home/raflynn/Mitosis/NETseq/metagenes/cNET_HeLa_Rep1_antisense/bins/Churchman_CTCF_all/allchr_sorted.txt", header=FALSE, colClasses=c(rep("character",6),rep("numeric",401)))
data5 = data4[,128:287]
data6 = apply(data5, 1, normalize)
data6[is.nan(data6)] = 0
data6 = t(data6)


ma <- function(arr, n){
  res = arr
  for(i in 1:length(arr)){
	res[i] = ifelse(i < n, mean(arr[1:i]), mean(arr[(i-n):i]))
  }
  res
}

# if I want to remove 0s
data60 = data6[apply(data6,1,sum) > 0,]
data30 = data3[apply(data3,1,sum) > 0,]
#

pdf("replicate_churchman_ctcf_fig6a.pdf", width=5, height=4)
plot(ma(apply(data30,2,mean), 5), type='l', lwd=2, col="red", xaxt='n',
xlab="Position", ylab="Mean normalized density")
title("Replicating Figure 6a: NETseq\n around 16,339 CTCF sites")
axis(1, at=seq(0, 160, 40), labels=seq(-400,400,200))
lines(ma(apply(data60,2,mean), 5), lwd=2, col="green")
legend("topleft", legend=c("Sense","Antisense"), lty=1, lwd=2, col=c("red","green"))
abline(v=80)
dev.off()

# only the top 500
data7 = data2[order(apply(data2,1,sum),decreasing=TRUE)[1:5000],]
data8 = apply(data7,1,normalize)
data8[is.nan(data8)] = 0
data8 = t(data8)
plot(apply(data8,2,mean),type='l')

data9 = data5[order(apply(data5,1,sum),decreasing=TRUE)[1:5000],]
data10 = apply(data9,1,normalize)
data10[is.nan(data10)] = 0
data10 = t(data10)
plot(apply(data10,2,mean),type='l')
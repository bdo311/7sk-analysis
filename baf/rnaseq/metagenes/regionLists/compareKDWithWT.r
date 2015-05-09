se = read.delim("mES_SEindiv_noTSS_1kb_list.txt",header=TRUE,row.names=1,colClasses=c("character",rep("numeric",4)))
re = read.delim("RY_enh_list.txt",header=TRUE,row.names=1,colClasses=c("character",rep("numeric",4)))

se_chg = log2(se$SRR1265785_rmdup+se$SRR1265787_rmdup)-log2(se$SRR1265786_rmdup+se$SRR1265788_rmdup)
se_chg_f = se_chg[is.finite(se_chg)] #5107, p = 0.00834 (T-test), median = 0.134

re_chg = log2(re$SRR1265785_rmdup+re$SRR1265787_rmdup)-log2(re$SRR1265786_rmdup+re$SRR1265788_rmdup)
re_chg_f = re_chg[is.finite(re_chg)] #349, p < 2.2e-16 (T-test), median = 0.661

pdf("erna_change.pdf",width=6,height=6)
plot(density(re_chg_f),xlim=c(-5,5),ylim=c(0,0.6),lwd=2,xlab="log2 Read Density Change",col="red",main="eRNA density change")
lines(density(se_chg_f),col="green",lwd=2)
abline(v=0,col="lightblue",lty=2)
legend("topleft",legend=c("RE","SE"),lwd=2,col=c("red","green"))
dev.off()
re = read.delim("re_matrix.txt",header=TRUE,row.names=1)
se = read.delim("se_matrix.txt",header=TRUE,row.names=1)
tss = read.delim("tss_matrix.txt",header=TRUE,row.names=1)

re.oct4.aso = with(re, log2(Pou5f1_5p_mergeAD/Input_5p_merge))
re.oct4.scr = with(re, log2(Pou5f1_Scr_mergeAD/Input_Scr_merge))

se.oct4.aso = with(se, log2(Pou5f1_5p_mergeAD/Input_5p_merge))
se.oct4.scr = with(se, log2(Pou5f1_Scr_mergeAD/Input_Scr_merge))

tss.oct4.aso = with(tss, log2(Pou5f1_5p_mergeAD/Input_5p_merge))
tss.oct4.scr = with(tss, log2(Pou5f1_Scr_mergeAD/Input_Scr_merge))


# re.oct4.aso = with(re, log2(BAF155_5p_merge/Input_5p_merge))
# re.oct4.scr = with(re, log2(BAF155_Scr_merge/Input_Scr_merge))
                            
# se.oct4.aso = with(se, log2(BAF155_5p_merge/Input_5p_merge))
# se.oct4.scr = with(se, log2(BAF155_Scr_merge/Input_Scr_merge))

# tss.oct4.aso = with(tss, log2(BAF155_5p_merge/Input_5p_merge))
# tss.oct4.scr = with(tss, log2(BAF155_Scr_merge/Input_Scr_merge))


re.chg = re.oct4.aso - re.oct4.scr
se.chg = se.oct4.aso - se.oct4.scr
tss.chg = tss.oct4.aso - tss.oct4.scr

pdf("region_comparison_oct4.pdf",width=10,height=10)
par(mfrow=c(2,2))

plot(tss.oct4.scr, tss.chg, pch=19, cex=0.1, xlab="WT log2 signal", ylab="log2 ASO/WT change", main="Oct4 WT binding vs ASO chg")
points(re.oct4.scr, re.chg, pch=19, cex=0.1, col="red")
points(se.oct4.scr, se.chg, pch=19, cex=0.2, col="green")
abline(h=0, col="lightblue",lty=2,lwd=2)
legend("topright", legend=c("TSS","RE","SE"), col=c("black","red","green"))

boxplot(tss.chg, re.chg, se.chg, ylim=c(-2,2), ylab="log2 ASO/WT change", main="Oct4 binding change")
axis(1, at=1:3, labels=c("TSS","RE","SE"))
abline(h=0, col="lightblue",lty=2,lwd=2)

# density plot of Oct4 change
plot(density(re.chg),col="red", xlab="log2 ASO/WT", main="Change in Oct4 binding after ASO")
lines(density(tss.chg[is.finite(tss.chg)]))
lines(density(se.chg[is.finite(se.chg)]),col="green")
legend("topright", legend=c("TSS","RE","SE"), col=c("black","red","green"), lwd=2, lty=1)
abline(v=0,col="lightblue",lty=2,lwd=2)

dev.off()

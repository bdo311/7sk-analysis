# makeCDF.r
# 3/4/15

#re = read.delim("re_matrix.txt",header=TRUE,row.names=1)
#re = read.delim("se_matrix.txt",header=TRUE,row.names=1)
re = read.delim("tss_matrix.txt",header=TRUE,row.names=1)

BAF155_3p = with(re, log2(BAF155_3p_merge/Input_3p_merge))
BAF155_5p = with(re, log2(BAF155_5p_merge/Input_5p_merge))
BAF155_Scr = with(re, log2(BAF155_Scr_merge/Input_Scr_merge))

# Pou5f1_3p_mergeAB = with(re, log2(Pou5f1_3p_mergeAB/Input_3p_merge))
# Pou5f1_5p_mergeAB = with(re, log2(Pou5f1_5p_mergeAB/Input_5p_merge))
# Pou5f1_Scr_mergeAB = with(re, log2(Pou5f1_Scr_mergeAB/Input_Scr_merge))

Pou5f1_3p_mergeAD = with(re, log2(Pou5f1_3p_mergeAD/Input_3p_merge))
Pou5f1_5p_mergeAD = with(re, log2(Pou5f1_5p_mergeAD/Input_5p_merge))
Pou5f1_Scr_mergeAD = with(re, log2(Pou5f1_Scr_mergeAD/Input_Scr_merge))

baf=BAF155_5p-BAF155_Scr
# pouAB=Pou5f1_5p_mergeAB-Pou5f1_Scr_mergeAB
pouAD=Pou5f1_5p_mergeAD-Pou5f1_Scr_mergeAD
Fn1=ecdf(baf)
# Fn2=ecdf(pouAB)
Fn3=ecdf(pouAD)
pdf("tss_chg.pdf",width=5,height=5)
plot(Fn1,lwd=2,xlim=c(-3,3),main="5p vs ASO Enrichment",xlab="log2 Enrichment 5p/ASO")
# lines(Fn2,col="red",lwd=2)
lines(Fn3,col="red",lwd=2)
abline(v=0,lty=2)
legend("topleft",legend=c("BAF155","Pou5f1_mergeAD"),lwd=2,col=c("black","red"))
dev.off()

pdf("tss_cdfs.pdf",width=5,height=8)
par(mfrow=c(2,1))
Fn1=ecdf(BAF155_5p)
Fn2=ecdf(BAF155_Scr)
plot(Fn2,lwd=2,xlab="log2 IP/Input",xlim=c(-1.5,2.5))	
lines(Fn1,col="red",lwd=2)
abline(v=0,lty=2,col="lightblue")
ks.test(BAF155_5p,BAF155_Scr) 

# Fn1=ecdf(Pou5f1_5p_mergeAB)
# Fn2=ecdf(Pou5f1_Scr_mergeAB)
# plot(Fn2,lwd=2,xlab="log2 IP/Input",xlim=c(-1.5,3.5))
# lines(Fn1,col="red",lwd=2)
# abline(v=0,lty=2,col="lightblue")
# ks.test(Pou5f1_5p_mergeAB,Pou5f1_Scr_mergeAB)

Fn1=ecdf(Pou5f1_5p_mergeAD)
Fn2=ecdf(Pou5f1_Scr_mergeAD)
plot(Fn2,lwd=2,xlab="log2 IP/Input",xlim=c(-1.5,3.5))
lines(Fn1,col="red",lwd=2)
abline(v=0,lty=2,col="lightblue")
ks.test(Pou5f1_5p_mergeAD,Pou5f1_Scr_mergeAD)

dev.off()

region="TSS_centered"
region="RY_enh_centered"
region="SE_indiv"
chirp = read.delim(paste("/home/raflynn/7SK/ChIRPseq/metagenes/7SK_mES_WT/bins/",region,"/allchr_sorted.txt",sep=''),header=FALSE,colClasses=c(rep("character",6),rep("numeric",401)))

if (region=="TSS_centered") {
	gro_wt = 
	read.delim(paste("/home/raflynn/7SK/GROseq/Flynn/metagenes/GRO_Scr_comb_norm_sense/bins/",region,"/allchr_sorted.txt",sep=''),header=FALSE,colClasses=c(rep("character",6),rep("numeric",401)))
	rownames(gro_wt) = gro_wt[,4]
	gro_aso = 
	read.delim(paste("/home/raflynn/7SK/GROseq/Flynn/metagenes/GRO_5ASO_comb_norm_sense/bins/",region,"/allchr_sorted.txt",sep=''),header=FALSE,colClasses=c(rep("character",6),rep("numeric",401)))
	rownames(gro_aso) = gro_aso[,4]
} else {
	wt_sense = 
	read.delim(paste("/home/raflynn/7SK/GROseq/Flynn/metagenes/GRO_Scr_comb_norm_sense/bins/",region,"/allchr_sorted.txt",sep=''),header=FALSE,colClasses=c(rep("character",6),rep("numeric",401)))
	rownames(wt_sense) = wt_sense[,4]
	wt_antisense = 
	read.delim(paste("/home/raflynn/7SK/GROseq/Flynn/metagenes/GRO_Scr_comb_norm_antisense/bins/",region,"/allchr_sorted.txt",sep=''),header=FALSE,colClasses=c(rep("character",6),rep("numeric",401)))
	rownames(wt_antisense) = wt_antisense[,4]
	wt_antisense = wt_antisense[rownames(wt_sense),]
	aso_sense = 
	read.delim(paste("/home/raflynn/7SK/GROseq/Flynn/metagenes/GRO_5ASO_comb_norm_sense/bins/",region,"/allchr_sorted.txt",sep=''),header=FALSE,colClasses=c(rep("character",6),rep("numeric",401)))
	rownames(aso_sense) = aso_sense[,4]
	aso_sense = aso_sense[rownames(wt_sense),]
	aso_antisense = 
	read.delim(paste("/home/raflynn/7SK/GROseq/Flynn/metagenes/GRO_5ASO_comb_norm_antisense/bins/",region,"/allchr_sorted.txt",sep=''),header=FALSE,colClasses=c(rep("character",6),rep("numeric",401)))
	rownames(aso_antisense) = aso_antisense[,4]
	aso_antisense = aso_antisense[rownames(wt_sense),]

	wt_vals = wt_sense[,8:407] + wt_antisense[,8:407]
	aso_vals = aso_sense[,8:407] + aso_antisense[,8:407]
	gro_wt = cbind(wt_sense[,1:7], wt_vals)
	gro_aso = cbind(aso_sense[,1:7], aso_vals)	
}

# ordering by ChIRP signal
rownames(chirp) = chirp[,4]
signal = apply(chirp[,158:257],1,sum)
chirp = chirp[order(signal, decreasing=TRUE),]
gro_wt = gro_wt[rownames(chirp),]
gro_aso = gro_aso[rownames(chirp),]

grochg = log2(apply(gro_aso[,108:307],1,sum)/apply(gro_wt[,108:307],1,sum))
grosig = log2(apply(gro_wt[,108:307],1,sum))
grosig_aso = log2(apply(gro_aso[,108:307],1,sum))
chirpsig = log2(apply(chirp[,108:307],1,sum))

chirpsig[chirpsig<3] = 3
chirpsig[chirpsig>11] = 11
grochg_binary = (grochg < -0.5)
df = data.frame(cbind(grochg, grosig, chirpsig, grochg_binary))
df = df[(is.finite(df$chirpsig) & is.finite(df$grochg) & is.finite(grosig)),]

# cool scatterplot

require(ggplot2)
df$grochg_lim = df$grochg
df$grochg_lim[df$grochg_lim < -1.2] = -1.2
df$grochg_lim[df$grochg_lim > 1.2] = 1.2
pdf("tss_rainbow.pdf", width=7, height=4)
ggplot(df, aes(x=grosig, y=chirpsig)) + 
geom_point(aes(color=grochg_lim), size=1) +
scale_colour_gradientn(colours=c("orange","black","cyan"), limits=c(-1.2,1.2)) +
xlim(-1,15) +
theme_bw()
dev.off()

df.hii=df[df$grochg>1,]
df.hi = df[df$grochg > 0 & df$grochg < 1,]
df.lo = df[df$grochg < 0 & df$grochg > -1,]
df.loo=df[df$grochg< -1,]

hii=loess(chirpsig~grosig,data=df.hii)
hi = loess(chirpsig~grosig, data=df.hi)
lo = loess(chirpsig~grosig, data=df.lo)
loo=loess(chirpsig~grosig,data=df.loo)



pdf(paste(region, "chirp_vs_gro_byGroChg.pdf", sep='_'), width=5, height=5)
plot(df$grosig,df$chirpsig,cex=0.2,pch=19, main="GROseq vs ChIRP signal,\nstratified by GROseq ASO/WT change",
xlab="log2 WT GROseq signal", ylab="log2 WT ChIRPseq signal")
o1 = order(df.hii$grosig)
lines(df.hii$grosig[o1], predict(hii)[o1], lwd=3, col="red")
o1 = order(df.hi$grosig)
lines(df.hi$grosig[o1], predict(hi)[o1], lwd=3, col="orange")
o1 = order(df.lo$grosig)
lines(df.lo$grosig[o1], predict(lo)[o1], lwd=3, col="green")
o1 = order(df.loo$grosig)
lines(df.loo$grosig[o1], predict(loo)[o1], lwd=3, col="blue")
legend("topleft", legend=c("> 1", "0-1", "-1-0", "< -1"), col=c("red", "orange", "green", "blue"), lwd=3)
dev.off()


  





pdf(paste(region,"_grochg_vs_chirpsig_boxplot.pdf",sep=''),width=18,height=5)
par(mfrow=c(1,7))

for (i in seq(3,9)) {
	df.num=df[df$grosig > i & df$grosig < (i+2),]
	boxplot(df.num[df.num$grochg< -1,]$chirpsig,
	df.num[df.num$grochg> -1 & df.num$grochg < 0,]$chirpsig,
	df.num[df.num$grochg>0 & df.num$grochg < 1,]$chirpsig,
	df.num[df.num$grochg>1,]$chirpsig,
	ylim=c(3,10),
	main=paste("GROseq ",i,'-',i+1,sep=''),
	xlab="log2 GROseq\nchange", ylab="log2 ChIRPseq signal",
	xaxt='n')
	axis(1, at=1:4,labels=c("< -1", "-1 to 0","0 to 1", "> 1"))
}
dev.off()

# all significant by KS test
ks.test(df3[df3$grochg<0,]$chirpsig,df3[df3$grochg>0,]$chirpsig)$p.value
ks.test(df4[df4$grochg<0,]$chirpsig,df4[df4$grochg>0,]$chirpsig)$p.value
ks.test(df5[df5$grochg<0,]$chirpsig,df5[df5$grochg>0,]$chirpsig)$p.value
ks.test(df6[df6$grochg<0,]$chirpsig,df6[df6$grochg>0,]$chirpsig)$p.value
ks.test(df7[df7$grochg<0,]$chirpsig,df7[df7$grochg>0,]$chirpsig)$p.value
ks.test(df8[df8$grochg<0,]$chirpsig,df8[df8$grochg>0,]$chirpsig)$p.value
ks.test(df9[df9$grochg<0,]$chirpsig,df9[df9$grochg>0,]$chirpsig)$p.value



# Density plots
plot(density(2^grosig),main="Reads at SE", ylim=c(0,2e-3),xlab="Reads -500 to +500 around SE",ylab="Density")
lines(density(2^grosig_aso),col="red")
plot(density(grosig),main="Reads at SE", xlab="log2 Reads -500 to +500 around SE",ylab="Density")
lines(density(grosig_aso),col="red")









### Scatterplots for change ###

pdf(paste(region, "_chg.pdf", sep=''), width=13, height=4)
par(mfrow=c(1,3))
plot(grosig, grochg, pch=19, cex=0.3,
main="GROseq change vs GRO signal", xlab="log2 GROseq signal WT", ylab="log2 GROseq change ASO/WT")
abline(h=0,lty=2,col="lightblue")
plot(chirpsig, grochg, pch=19, cex=0.3, xlim=c(2,12),
main="GROseq change vs ChIRP signal", xlab="log2 ChIRPseq signal WT", ylab="log2 GROseq change ASO/WT")
abline(h=0,lty=2,col="lightblue")
plot(chirpsig, grosig, pch=19, cex=0.3, xlim=c(2,12),
main="GROseq signal vs ChIRP signal", xlab="log2 ChIRPseq signal WT", ylab="log2 GROseq signal WT")
dev.off()

require(ggplot2)
require(RColorBrewer)

df2 = df[abs(df$grochg) > 0.5,]
pdf(paste(region, "_chg_color.pdf", sep=''), width=5, height=4)
ggplot(df, aes(grosig, grochg)) + geom_point(aes(color=chirpsig), size=1) + 
theme_bw() + scale_colour_gradient(low="red",high="darkblue",limits=c(3,11)) +
ggtitle("GROseq change vs GROseq signal")
ggplot(df, aes(chirpsig, grosig)) + geom_point(aes(color=grochg), size=1) + 
theme_bw() + scale_colour_gradient(low="red",high="darkblue") +
ggtitle("GROseq signal vs ChIRPseq signal")
ggplot(df2, aes(grosig, chirpsig)) + geom_point(aes(color=grochg_binary), size=1) + 
theme_bw() + ggtitle("ChIRPseq signal vs GROseq signal")
dev.off()
  
  
  
# ranked everything
dfo = as.data.frame(cbind(rank(df$chirpsig),rank(df$grosig),df$grochg))
names(dfo) = c("chirpsig","grosig","grochg")
dfo$cminusg = dfo$chirpsig - dfo$grosig
plot(rank(dfo$cminusg), dfo$grochg)

# z-scored everything
newchirpsig = (df$chirpsig-mean(df$chirpsig))/mean(df$chirpsig)
newgrosig = (df$grosig-mean(df$grosig))/mean(df$grosig)

dfm = as.data.frame(cbind(newchirpsig,newgrosig,df$grochg))
names(dfm) = c("chirpsig","grosig","grochg")
dfm$cminusg = dfm$chirpsig - dfm$grosig
plot(rank(dfm$cminusg), dfm$grochg)



pdf(paste(region,"_grochg_vs_chirpsig_scatter.pdf",sep=''),width=10,height=5)
par(mfrow=c(1,2))
plot(log2(chirpsig/grosig),grochg,cex=0.3,pch=19,xlab="log2 ChIRP signal / GRO signal", ylab="log2 GROseq change")
dfx=df[order(log2(df$chirpsig/df$grosig)),]
plot(dfx$grochg,cex=0.3,pch=19,xlab="Rank of log2 ChIRP signal/GRO signal", ylab="log2 GROseq change")
dev.off()

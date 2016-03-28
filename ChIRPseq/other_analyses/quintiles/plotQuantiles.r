region="TSS_centered"
region="RY_enh_centered"
region="SE_indiv"

chirp = 
read.delim(paste("/home/raflynn/7SK/ChIRPseq/metagenes/7SK_mES_WT/bins/",region,"/allchr_sorted.txt",sep=''),header=FALSE,colClasses=c(rep("character",6),rep("numeric",401)))
rownames(chirp) = chirp[,4]
signal = apply(chirp[,158:257],1,sum)
chirp = chirp[order(signal, decreasing=TRUE),]

if (region=="TSS_centered") {
	gro_wt = 
	read.delim(paste("/home/raflynn/7SK/GROseq/Flynn/metagenes/GRO_Scr_comb_norm_sense/bins/",region,"/allchr_sorted.txt",sep=''),header=FALSE,colClasses=c(rep("character",6),rep("numeric",401)))
	gro_aso = 
	read.delim(paste("/home/raflynn/7SK/GROseq/Flynn/metagenes/GRO_5ASO_comb_norm_sense/bins/",region,"/allchr_sorted.txt",sep=''),header=FALSE,colClasses=c(rep("character",6),rep("numeric",401)))
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

rownames(gro_wt) = gro_wt[,4]
rownames(gro_aso) = gro_aso[,4]
gro_wt_ord = gro_wt[rownames(chirp),]
gro_aso_ord = gro_aso[rownames(chirp),]

pdf(paste(region, 'chirp_and_gro.pdf', sep='_'), width=12, height=4)
par(mfrow=c(1,3))

# Plot 1: quintiles for ChIRP

n = 5
x = ceiling(nrow(chirp)/n)
plot(apply(chirp[1:x,8:407],2,mean),type='l',ylim=c(0,3), col="firebrick", lwd=2,
main="ChIRPseq signal per quintile", ylab="ChIRPseq signal")
lines(apply(chirp[(x+1):(2*x),8:407],2,mean), col="orange", lwd=2)
lines(apply(chirp[(2*x+1):(3*x),8:407],2,mean), col="goldenrod", lwd=2)
lines(apply(chirp[(3*x+1):(4*x-10),8:407],2,mean), col="seagreen", lwd=2)
lines(apply(chirp[(4*x+1):(5*x-10),8:407],2,mean), col="lightblue", lwd=2)

# Plot 2: ChIRPseq vs GROseq
chirp_vals = log2(apply(chirp[,8:407], 1, mean))
gro_wt_vals = log2(apply(gro_wt_ord[,8:407], 1, mean))
plot(chirp_vals, gro_wt_vals, pch=19, cex=0.3, main="ChIRPseq vs GROseq signal",
xlab="log2 ChIRPseq signal", ylab="log2 GROseq signal", xlim=c(-6,4))

wt1 = apply(gro_wt_ord[1:x,8:407],2,mean)
aso1 = apply(gro_aso_ord[1:x,8:407],2,mean)
wt2 = apply(gro_wt[(x+1):(2*x),8:407],2,mean)
aso2 = apply(gro_aso[(x+1):(2*x),8:407],2,mean)
wt3 = apply(gro_wt[(2*x+1):(3*x),8:407],2,mean)
aso3 = apply(gro_aso[(2*x+1):(3*x),8:407],2,mean)
wt4 = apply(gro_wt[(3*x+1):(4*x),8:407],2,mean)
aso4 = apply(gro_aso[(3*x+1):(4*x),8:407],2,mean)
wt5 = apply(gro_wt[(4*x+1):(5*x-10),8:407],2,mean)
aso5 = apply(gro_aso[(4*x+1):(5*x-10),8:407],2,mean)

# Plot 3: GROseq change per ChIRP quintile
plot(log2(aso1/wt1), type='l', lwd=2, col="firebrick", 
ylim=c(-0.6, 0.6), ylab="log2 GROseq ASO/WT",
main="GROseq ASO/WT change per quintile")
lines(log2(aso2/wt2), lwd=2, col="orange")
lines(log2(aso3/wt3), lwd=2, col="goldenrod")
lines(log2(aso4/wt4), lwd=2, col="seagreen")
lines(log2(aso5/wt5), lwd=2, col="lightblue")
abline(h=0, lty=2, lwd=2, col="darkblue")

dev.off()

pdf(paste(region, "gro_by_quintile.pdf", sep='_'),width=10,height=3)
par(mfrow=c(1,5))
plot(wt1,col="black",type='l',ylim=c(0,1.3*max(wt1)), main="Quintile 1")
lines(aso1,col="red")
plot(wt2,col="black",type='l',ylim=c(0,1.3*max(wt2)), main="Quintile 2")
lines(aso2,col="red")
plot(wt3,col="black",type='l',ylim=c(0,1.3*max(wt3)), main="Quintile 3")
lines(aso3,col="red")
plot(wt4,col="black",type='l',ylim=c(0,1.3*max(wt4)), main="Quintile 4")
lines(aso4,col="red")
plot(wt5,col="black",type='l',ylim=c(0,1.3*max(wt5)), main="Quintile 5")
lines(aso5,col="red")

dev.off()



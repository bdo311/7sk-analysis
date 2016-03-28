require(parallel)
require(ggplot2)

	
# bootstrap
getMean = function(x, i) {
	nums = sample(1:nrow(x), size=nrow(x),replace=TRUE)
	return(apply(x[nums,8:ncol(x)],2,mean))
}
	
region="TSS_centered"
region="RY_enh_centered"
region="SE_indiv"

wt_sense = 
read.delim(paste("/home/raflynn/7SK/GROseq/Flynn/metagenes/GRO_Scr_comb_norm_sense/bins/",region,"/allchr_sorted.txt",sep=''),header=FALSE,colClasses=c(rep("character",6),rep("numeric",401)))
rownames(wt_sense) = wt_sense[,4]
if (region == "TSS_centered") {wt_sense = wt_sense[-c(1),]}
if (region == "RY_enh_centered") {wt_sense = wt_sense[-c(3,4),]}

wt_antisense = 
read.delim(paste("/home/raflynn/7SK/GROseq/Flynn/metagenes/GRO_Scr_comb_norm_antisense/bins/",region,"/allchr_sorted.txt",sep=''),header=FALSE,colClasses=c(rep("character",6),rep("numeric",401)))
rownames(wt_antisense) = wt_antisense[,4]
if (region == "TSS_centered") {wt_antisense = wt_antisense[-c(1),]}
if (region == "RY_enh_centered") {wt_antisense = wt_antisense[-c(3,4),]}

aso_sense = 
read.delim(paste("/home/raflynn/7SK/GROseq/Flynn/metagenes/GRO_5ASO_comb_norm_sense/bins/",region,"/allchr_sorted.txt",sep=''),header=FALSE,colClasses=c(rep("character",6),rep("numeric",401)))
rownames(aso_sense) = aso_sense[,4]
if (region == "TSS_centered") {aso_sense = aso_sense[-c(1),]}
if (region == "RY_enh_centered") {aso_sense = aso_sense[-c(3,4),]}

aso_antisense = 
read.delim(paste("/home/raflynn/7SK/GROseq/Flynn/metagenes/GRO_5ASO_comb_norm_antisense/bins/",region,"/allchr_sorted.txt",sep=''),header=FALSE,colClasses=c(rep("character",6),rep("numeric",401)))
rownames(aso_antisense) = aso_antisense[,4]
if (region == "TSS_centered") {aso_antisense = aso_antisense[-c(1),]}
if (region == "RY_enh_centered") {aso_antisense = aso_antisense[-c(3,4),]}


wt_sense_means = t(as.data.frame(mclapply(1:100,function(i) getMean(wt_sense, i),mc.cores=10)))
wt_antisense_means = t(as.data.frame(mclapply(1:100,function(i) getMean(wt_antisense, i),mc.cores=10)))
aso_sense_means = t(as.data.frame(mclapply(1:100,function(i) getMean(aso_sense, i),mc.cores=10)))
aso_antisense_means = t(as.data.frame(mclapply(1:100,function(i) getMean(aso_antisense, i),mc.cores=10)))

wt_sense_mean = apply(wt_sense_means,2,mean)
wt_sense_975 = apply(wt_sense_means,2,function(x) quantile(x,0.975)[[1]])
wt_sense_025 = apply(wt_sense_means,2,function(x) quantile(x,0.025)[[1]])

aso_sense_mean = apply(aso_sense_means,2,mean)
aso_sense_975 = apply(aso_sense_means,2,function(x) quantile(x,0.975)[[1]])
aso_sense_025 = apply(aso_sense_means,2,function(x) quantile(x,0.025)[[1]])

wt_antisense_mean = apply(wt_antisense_means,2,mean)
wt_antisense_975 = apply(wt_antisense_means,2,function(x) quantile(x,0.975)[[1]])
wt_antisense_025 = apply(wt_antisense_means,2,function(x) quantile(x,0.025)[[1]])

aso_antisense_mean = apply(aso_antisense_means,2,mean)
aso_antisense_975 = apply(aso_antisense_means,2,function(x) quantile(x,0.975)[[1]])
aso_antisense_025 = apply(aso_antisense_means,2,function(x) quantile(x,0.025)[[1]])

pdf(paste(region,"_gro_metagene.pdf",sep=''),width=5,height=4)
means = data.frame(bins=seq(-1000,995,5),wt_sense_mean,aso_sense_mean,wt_antisense_mean,aso_antisense_mean)
(plt = ggplot(means,aes(bins)) + 

geom_line(aes(y=wt_sense_mean,colour="WT"),color="black") + 
geom_line(aes(y=-wt_antisense_mean,colour="WT"),color="black") + 
geom_line(aes(y=aso_sense_mean,color="ASO"),color="red") + 
geom_line(aes(y=-aso_antisense_mean,color="ASO"),color="red") + 

geom_ribbon(aes(ymin=wt_sense_025,ymax=wt_sense_975),fill="black",alpha=0.3) + 
geom_ribbon(aes(ymin=-wt_antisense_975,ymax=-wt_antisense_025),fill="black",alpha=0.3) + 
geom_ribbon(aes(ymin=aso_sense_025,ymax=aso_sense_975),fill="red",alpha=0.3) + 
geom_ribbon(aes(ymin=-aso_antisense_975,ymax=-aso_antisense_025),fill="red",alpha=0.3) + 
theme_bw() + 
ylim(-max(means[,2:ncol(means)])*1.3,max(means[,2:ncol(means)])*1.3) +
#ylim(-6, 17) +
ggtitle("Metagenes around RE") + 
ylab("Value") + xlab("Distance from center")) +
geom_hline(yintercept=0)
dev.off()




# Hotelling's
HotellingsT2(wt_sense[,358:407], aso_sense[,358:407])
HotellingsT2(wt_sense[,308:407], aso_sense[,308:407])
HotellingsT2(wt_antisense[,8:57], aso_antisense[,8:57])
HotellingsT2(wt_antisense[,8:107], aso_antisense[,8:107])
HotellingsT2(rbind(as.matrix(wt_sense[,407:308]),as.matrix(wt_antisense[,8:107])), rbind(as.matrix(aso_sense[,407:308]),as.matrix(aso_antisense[,8:107])))
HotellingsT2(rbind(as.matrix(wt_sense[,407:358]),as.matrix(wt_antisense[,8:57])), rbind(as.matrix(aso_sense[,407:358]),as.matrix(aso_antisense[,8:57])))



# TSS
data=read.delim("startRNA_NELFwt_tss_peaks.txt",na.strings='.',header=FALSE)

# Basic stats
nrow(data) #14116
# bidirectional : 2058
# divergent     :10672
# unidirectional: 1386

# Distance between peaks
distBetweenPeaks = data$V12-data$V8
distBetweenPeaks = abs(distBetweenPeaks[is.finite(distBetweenPeaks)])
plot(density(distBetweenPeaks))
median(distBetweenPeaks) #184 bp

# Size of sense peak
senseSize = log2(as.numeric(apply(data,1,function(x) {ifelse(x[6]=="+",x[14],x[10])})))
plot(density(senseSize))
median(senseSize) #6.64 --> 99.7

# Ratio of sense peak size over antisense peak
ratio = log2(apply(data,1,function(x) {ifelse(x[6]=="+",as.numeric(x[14])/as.numeric(x[10]),as.numeric(x[10])/as.numeric(x[14]))}))
ratio.all = ratio[is.finite(ratio.all)]
plot(density(ratio.all))
abline(v=0)
median(ratio.all) #2.092 --> 4.26

ratio.bidir = ratio[data[,15]=="bidirectional"]
ratio.bidir = ratio.bidir[is.finite(ratio.bidir)]
median(ratio.bidir) #0.025

ratio.div = ratio[data[,15]=="divergent"]
ratio.div = ratio.div[is.finite(ratio.div)]
median(ratio.div) #2.458 --> 5.49

# Distance from annotated TSS
tssDist = apply(data,1,function(x) {ifelse(x[6]=="+",as.numeric(x[12])-(as.numeric(x[2])+1000),as.numeric(x[2])+1000-as.numeric(x[8]))})
plot(density(tssDist))
nums = seq(0,1000,by=50)
numGreater = sapply(nums,function(x)sum(abs(tssDist)>=x))
plot(nums,numGreater,type='l',xlab="Distance from TSS",ylab="Number of genes with equal or greater distance")

## RE ##
data=read.delim("startRNA_NELFwt_RYenh_peaks.txt",na.strings='.',header=FALSE)
nrow(data) #5356
distBetweenPeaks = data$V15
plot(density(distBetweenPeaks))
median(distBetweenPeaks) #273 bp

# Size of peaks
senseSize = log2(c(data$V14,data$V10))
median(senseSize) # 2.807 --> 6.99

## SE ##
data=read.delim("startRNA_NELFwt_SEindiv_peaks.txt",na.strings='.',header=FALSE)
nrow(data) #360
distBetweenPeaks = data$V18
plot(density(distBetweenPeaks))
median(distBetweenPeaks) #314 bp

# Size of peaks
senseSize = log2(c(data$V17,data$V13))
median(senseSize) # 5.156 --> 35.7





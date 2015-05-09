gro_se_sense = read.delim("/home/raflynn/7SK/GROseq/Flynn/metagenes/GRO_12Ccomb_sense/bins/mES_SEindiv_noTSS_1kb/allchr_sorted.txt",header=FALSE,colClasses=c(rep("character",7),rep("numeric",400)))
gro_se_antisense = read.delim("/home/raflynn/7SK/GROseq/Flynn/metagenes/GRO_12Ccomb_antisense/bins/mES_SEindiv_noTSS_1kb/allchr_sorted.txt",header=FALSE,colClasses=c(rep("character",7),rep("numeric",400)))

gro_enh_sense = read.delim("/home/raflynn/7SK/GROseq/Flynn/metagenes/GRO_12Ccomb_sense/bins/RY_enh/allchr_sorted.txt",header=FALSE,colClasses=c(rep("character",7),rep("numeric",400)))
gro_enh_antisense = read.delim("/home/raflynn/7SK/GROseq/Flynn/metagenes/GRO_12Ccomb_antisense/bins/RY_enh/allchr_sorted.txt",header=FALSE,colClasses=c(rep("character",7),rep("numeric",400)))

chirp_se = read.delim("/home/raflynn/7SK/ChIRPseq/metagenes/7SK_mES_WT/bins/mES_SEindiv_noTSS_1kb/allchr_sorted.txt",header=FALSE,colClasses=c(rep("character",7),rep("numeric",400)))
chirp_enh = read.delim("/home/raflynn/7SK/ChIRPseq/metagenes/7SK_mES_WT/bins/RY_enh/allchr_sorted.txt",header=FALSE,colClasses=c(rep("character",7),rep("numeric",400)))



gro_se_sense_vals = apply(gro_se_sense[,8:ncol(gro_se_sense)],1,mean)
names(gro_se_sense_vals) = apply(gro_se_sense[,1:6],1,function(x) paste(x,collapse='__'))
gro_se_sense_vals = gro_se_sense_vals[order(names(gro_se_sense_vals))]

gro_se_antisense_vals = apply(gro_se_antisense[,8:ncol(gro_se_antisense)],1,mean)
names(gro_se_antisense_vals) = apply(gro_se_antisense[,1:6],1,function(x) paste(x,collapse='__'))
gro_se_antisense_vals = gro_se_antisense_vals[order(names(gro_se_antisense_vals))]

gro_enh_sense_vals = apply(gro_enh_sense[,8:ncol(gro_enh_sense)],1,mean)
names(gro_enh_sense_vals) = apply(gro_enh_sense[,1:6],1,function(x) paste(x,collapse='__'))
gro_enh_sense_vals = gro_enh_sense_vals[order(names(gro_enh_sense_vals))]

gro_enh_antisense_vals = apply(gro_enh_antisense[,8:ncol(gro_enh_antisense)],1,mean)
names(gro_enh_antisense_vals) = apply(gro_enh_antisense[,1:6],1,function(x) paste(x,collapse='__'))
gro_enh_antisense_vals = gro_enh_antisense_vals[order(names(gro_enh_antisense_vals))]

gro_se_vals = gro_se_sense_vals + gro_se_antisense_vals
gro_enh_vals = gro_enh_sense_vals + gro_enh_antisense_vals

# ---

chirp_se_vals = apply(chirp_se[,8:ncol(chirp_se)],1,mean)
names(chirp_se_vals) = apply(chirp_se[,1:6],1,function(x) paste(x,collapse='__'))
chirp_se_vals = chirp_se_vals[order(names(chirp_se_vals))]

chirp_enh_vals = apply(chirp_enh[,8:ncol(chirp_enh)],1,mean)
names(chirp_enh_vals) = apply(chirp_enh[,1:6],1,function(x) paste(x,collapse='__'))
chirp_enh_vals = chirp_enh_vals[order(names(chirp_enh_vals))]

# boxplot
pdf("gro_chirp_se_enh_comparison.pdf",width=7,height=6)
par(mfrow=c(1,2))
boxplot(log10(chirp_se_vals),log10(chirp_enh_vals),ylab="log10 Mean 7SK ChIRPseq Read Density/nt",xaxt='n',main="7SK ChIRPseq")
axis(1,at=1:2,labels=c("Super\nEnhancer","Regular\nEnhancer"), tck=0)
ks.test(chirp_se_vals, chirp_enh_vals)
boxplot(log10(gro_se_vals),log10(gro_enh_vals),ylab="log10 Mean GROseq Read Density/nt",xaxt='n',main="GROseq")
axis(1,at=1:2,labels=c("Super\nEnhancer","Regular\nEnhancer"), tck=0)
ks.test(gro_se_vals, gro_enh_vals)
dev.off()

# hockey plot
pdf("gro_chirp_se_enh_hockey.pdf",width=7,height=6)
par(mfrow=c(1,2))
chirp.df = rbind(cbind(chirp_enh_vals,1),cbind(chirp_se_vals,2))
chirp.df = chirp.df[order(chirp.df[,1],decreasing=FALSE),]
plot(chirp.df[,1],pch=19,col=chirp.df[,2],main="7SK ChIRPseq",xlab="Ranked enhancers",ylab="ChIRPseq read density")
gro.df = rbind(cbind(gro_enh_vals,1),cbind(gro_se_vals,2))
gro.df = gro.df[order(gro.df[,1],decreasing=FALSE),]
plot(gro.df[,1],pch=19,col=gro.df[,2],main="GROseq",xlab="Ranked enhancers",ylab="GROseq read density")
dev.off()

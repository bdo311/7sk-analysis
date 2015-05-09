data=read.delim("TSS_centered_list.txt",header=TRUE,row.names=1)

data$aso1 = log2(with(data,GRO_1251_sense+GRO_1251_antisense)+1)
data$aso2 = log2(with(data,GRO_1252_sense+GRO_1252_antisense)+1)
data$ctrl1 = log2(with(data,GRO_12C1_sense+GRO_12C1_antisense)+1)
data$ctrl2 = log2(with(data,GRO_12C2_sense+GRO_12C2_antisense)+1)

pdf("groseq_replicates.pdf",width=8,height=5)
par(mfrow=c(1,2))
plot(data$aso1,data$aso2,pch=19,cex=0.2,xlab="Replicate 1, log2 signal",ylab="Replicate 2, log2 signal",main="7SK ASO GROseq TSS signal")
plot(data$ctrl1,data$ctrl2,pch=19,cex=0.2,xlab="Replicate 1, log2 signal",ylab="Replicate 2, log2 signal",main="7SK WT GROseq TSS signal")
dev.off()


cor.test(log2(data$aso1+1),log2(data$aso2+1),method="spe")
cor.test(log2(data$ctrl1+1),log2(data$ctrl2+1),method="spe")
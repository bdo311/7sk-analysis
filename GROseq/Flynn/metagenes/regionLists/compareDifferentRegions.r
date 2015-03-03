# compareDifferentRegions.r
# to compare GRO-seq signal at different region types during 7SK ASO 

aso = c("GRO_12Ccomb","GRO_123comb","GRO_125comb","GRO_6Ccomb","GRO_63comb","GRO_65comb")
#aso = c("GRO_12C1","GRO_12C2","GRO_1231","GRO_1232","GRO_1251","GRO_1252","GRO_6C1","GRO_6C2","GRO_631","GRO_632","GRO_651","GRO_652")

### make bargraphs

means = list()
ks = list()
#for (regionType in c("TSS_trans_paused","TSS_transcribed","TSS_paused","GeneBody","TES")) {
for (regionType in c("TSS_trans_paused","GeneBody_trans_paused","TES_trans_paused")) {
	fn = paste(regionType, "_list.txt",sep='')
	data = read.delim(fn,header=TRUE,row.names=1,colClasses=c("character",rep("numeric",12)))
	
	sampleNames = as.character(sapply(aso, function(x) paste(x,"_sense",sep='')))
	data=data[,sampleNames]
	names(data) = aso

	meansForSample = apply(data, 2, mean)
	
	ksForSample = c(ks.test(data$GRO_12Ccomb,data$GRO_123comb)$p.value,
					ks.test(data$GRO_12Ccomb,data$GRO_125comb)$p.value,
					ks.test(data$GRO_6Ccomb,data$GRO_63comb)$p.value,
					ks.test(data$GRO_6Ccomb,data$GRO_65comb)$p.value)
	means[[regionType]] = meansForSample
	ks[[regionType]] = ksForSample	
}

for (regionType in c("RY_enh","SE_indiv")) {
	fn = paste(regionType, "_list.txt",sep='')
	data = read.delim(fn,header=TRUE,row.names=1,colClasses=c("character",rep("numeric",12)))
	
	sense = as.character(sapply(aso, function(x) paste(x,"_sense",sep='')))
	antisense = as.character(sapply(aso, function(x) paste(x,"_antisense",sep='')))
	
	combined = data[,sense] + data[,antisense]
	names(combined) = aso
	meansForSample = apply(combined, 2, mean)
	
	ksForSample = c(ks.test(combined$GRO_12Ccomb,combined$GRO_123comb)$p.value,
					ks.test(combined$GRO_12Ccomb,combined$GRO_125comb)$p.value,
					ks.test(combined$GRO_6Ccomb,combined$GRO_63comb)$p.value,
					ks.test(combined$GRO_6Ccomb,combined$GRO_65comb)$p.value)
	means[[regionType]] = meansForSample
	ks[[regionType]] = ksForSample	
}


# table for bedgraphs

means=do.call(rbind, means)
ks=do.call(rbind, ks)
write.table(means,"means_comb.txt",sep='\t',quote=FALSE)
write.table(ks,"ks_comb.txt",sep='\t',quote=FALSE)


### make boxplots of change
aso = c("GRO_12Ccomb","GRO_123comb","GRO_125comb","GRO_6Ccomb","GRO_63comb","GRO_65comb")
ttest = list()

list.125 = list()
list.65 = list()

pdf("change_comb.pdf",width=5,10)
par(mfrow=c(3,2))
par(mar=c(3,4,2,2))
#for (regionType in c("TSS_trans_paused","TSS_transcribed","TSS_paused","GeneBody","TES")) {
for (regionType in c("TSS_trans_paused","GeneBody_trans_paused","TES_trans_paused")) {
	fn = paste(regionType, "_list.txt",sep='')
	data = read.delim(fn,header=TRUE,row.names=1,colClasses=c("character",rep("numeric",12)))
	
	sampleNames = as.character(sapply(aso, function(x) paste(x,"_sense",sep='')))
	data=data[,sampleNames]
	names(data) = aso
	
	gro123 = log2(data$GRO_123comb/data$GRO_12Ccomb)
	gro125 = log2(data$GRO_125comb/data$GRO_12Ccomb)
	gro63 = log2(data$GRO_63comb/data$GRO_6Ccomb)
	gro65 = log2(data$GRO_65comb/data$GRO_6Ccomb)
	
	gro123=gro123[is.finite(gro123)]
	gro125=gro125[is.finite(gro125)]
	gro63=gro63[is.finite(gro63)]
	gro65=gro63[is.finite(gro65)]
	
	list.125[[regionType]] = gro125
	list.65[[regionType]] = gro65
	
	boxplot(gro123,gro125,gro63,gro65,
			main=regionType, xaxt='n',ylab="log2 Change versus WT",ylim=c(-2,2))
	axis(1,at=1:4,labels=c("123","125","63","65"))
	abline(h=0,col="lightblue",lty=2)
	
	tResults = c(t.test(gro123,mu=0)$p.value,
				 t.test(gro125,mu=0)$p.value,
				 t.test(gro63,mu=0)$p.value,
				 t.test(gro65,mu=0)$p.value)
	ttest[[regionType]] = tResults
}

for (regionType in c("RY_enh","SE_indiv")) {
	fn = paste(regionType, "_list.txt",sep='')
	data = read.delim(fn,header=TRUE,row.names=1,colClasses=c("character",rep("numeric",12)))
	
	sense = as.character(sapply(aso, function(x) paste(x,"_sense",sep='')))
	antisense = as.character(sapply(aso, function(x) paste(x,"_antisense",sep='')))
	
	data = data[,sense] + data[,antisense]
	names(data) = aso

	gro123 = log2(data$GRO_123comb/data$GRO_12Ccomb)
	gro125 = log2(data$GRO_125comb/data$GRO_12Ccomb)
	gro63 = log2(data$GRO_63comb/data$GRO_6Ccomb)
	gro65 = log2(data$GRO_65comb/data$GRO_6Ccomb)
	
	gro123=gro123[is.finite(gro123)]
	gro125=gro125[is.finite(gro125)]
	gro63=gro63[is.finite(gro63)]
	gro65=gro63[is.finite(gro65)]

	list.125[[regionType]] = gro125
	list.65[[regionType]] = gro65
		
	boxplot(gro123,gro125,gro63,gro65,
			main=regionType, xaxt='n',ylab="log2 Change versus WT",ylim=c(-2,2))
	axis(1,at=1:4,labels=c("123","125","63","65"))
	abline(h=0,col="lightblue",lty=2)
	
	tResults = c(t.test(gro123,mu=0)$p.value,
				 t.test(gro125,mu=0)$p.value,
				 t.test(gro63,mu=0)$p.value,
				 t.test(gro65,mu=0)$p.value)
	ttest[[regionType]] = tResults
}

dev.off()

ttest=do.call(rbind, ttest)
write.table(ttest,"ttest_comb.txt",sep='\t',quote=FALSE)

# make combined boxplots
pdf("5paso.pdf",width=6,height=5)
par(mfrow=c(1,2))
with(list.65,boxplot(TSS_trans_paused,RY_enh,SE_indiv,main="6hr ASO",ylab='log2 Change in GROseq Read Density versus WT',ylim=c(-3,3),col=c("gray","red","green")))
axis(1,at=1:3,labels=c("TSS","RE","SE"))
abline(h=0,col="lightblue",lty=2)
with(list.125,boxplot(TSS_trans_paused,RY_enh,SE_indiv,main="12hr ASO",ylab='log2 Change in GROseq Read Density versus WT',ylim=c(-3,3),col=c("gray","red","green")))
axis(1,at=1:3,labels=c("TSS","RE","SE"))
abline(h=0,col="lightblue",lty=2)
dev.off()
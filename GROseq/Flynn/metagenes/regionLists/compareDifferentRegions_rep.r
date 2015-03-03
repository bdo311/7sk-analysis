# compareDifferentRegions.r
# to compare GRO-seq signal at different region types during 7SK ASO 

#aso = c("GRO_12Ccomb","GRO_123comb","GRO_125comb","GRO_6Ccomb","GRO_63comb","GRO_65comb")
aso = c("GRO_12C1","GRO_12C2","GRO_1231","GRO_1232","GRO_1251","GRO_1252","GRO_6C1","GRO_6C2","GRO_631","GRO_632","GRO_651","GRO_652")

### make bargraphs

# means = list()
# ks = list()
# for (regionType in c("TSS_trans_paused","TSS_transcribed","TSS_paused","GeneBody","TES")) {
	# fn = paste(regionType, "_list.txt",sep='')
	# data = read.delim(fn,header=TRUE,row.names=1,colClasses=c("character",rep("numeric",12)))
	
	# sampleNames = as.character(sapply(aso, function(x) paste(x,"_sense",sep='')))
	# data=data[,sampleNames]
	# names(data) = aso

	# meansForSample = apply(data, 2, mean)
	
	# ksForSample = c(ks.test(data$GRO_12Ccomb,data$GRO_123comb)$p.value,
					# ks.test(data$GRO_12Ccomb,data$GRO_125comb)$p.value,
					# ks.test(data$GRO_6Ccomb,data$GRO_63comb)$p.value,
					# ks.test(data$GRO_6Ccomb,data$GRO_65comb)$p.value)
	# means[[regionType]] = meansForSample
	# ks[[regionType]] = ksForSample	
# }

# for (regionType in c("RY_enh","SE_indiv")) {
	# fn = paste(regionType, "_list.txt",sep='')
	# data = read.delim(fn,header=TRUE,row.names=1,colClasses=c("character",rep("numeric",12)))
	
	# sense = as.character(sapply(aso, function(x) paste(x,"_sense",sep='')))
	# antisense = as.character(sapply(aso, function(x) paste(x,"_antisense",sep='')))
	
	# combined = data[,sense] + data[,antisense]
	# names(combined) = aso
	# meansForSample = apply(combined, 2, mean)
	
	# ksForSample = c(ks.test(combined$GRO_12Ccomb,combined$GRO_123comb)$p.value,
					# ks.test(combined$GRO_12Ccomb,combined$GRO_125comb)$p.value,
					# ks.test(combined$GRO_6Ccomb,combined$GRO_63comb)$p.value,
					# ks.test(combined$GRO_6Ccomb,combined$GRO_65comb)$p.value)
	# means[[regionType]] = meansForSample
	# ks[[regionType]] = ksForSample	
# }


# # table for bedgraphs

# means=do.call(rbind, means)
# ks=do.call(rbind, ks)
# write.table(means,"means.txt",sep='\t',quote=FALSE)
# write.table(ks,"ks.txt",sep='\t',quote=FALSE)


### make boxplots of change
#aso = c("GRO_12Ccomb","GRO_123comb","GRO_125comb","GRO_6Ccomb","GRO_63comb","GRO_65comb")
aso = c("GRO_12C1","GRO_12C2","GRO_1231","GRO_1232","GRO_1251","GRO_1252","GRO_6C1","GRO_6C2","GRO_631","GRO_632","GRO_651","GRO_652")
ttest = list()

pdf("change_reps.pdf",width=12,10)
par(mfrow=c(3,3))
par(mar=c(3,4,2,2))
for (regionType in c("TSS_trans_paused","TSS_transcribed","TSS_paused","GeneBody","TES")) {
	fn = paste(regionType, "_list.txt",sep='')
	data = read.delim(fn,header=TRUE,row.names=1,colClasses=c("character",rep("numeric",36)))
	
	sampleNames = as.character(sapply(aso, function(x) paste(x,"_sense",sep='')))
	data=data[,sampleNames]
	names(data) = aso
	
	gro1231 = log2(data$GRO_1231/data$GRO_12C1)
	gro1251 = log2(data$GRO_1251/data$GRO_12C1)
	gro631 = log2(data$GRO_631/data$GRO_6C1)
	gro651 = log2(data$GRO_651/data$GRO_6C1)	
	gro1232 = log2(data$GRO_1232/data$GRO_12C2)
	gro1252 = log2(data$GRO_1252/data$GRO_12C2)
	gro632 = log2(data$GRO_632/data$GRO_6C2)
	gro652 = log2(data$GRO_652/data$GRO_6C2)
	
	gro1231=gro1231[is.finite(gro1231)]
	gro1251=gro1251[is.finite(gro1251)]
	gro631=gro631[is.finite(gro631)]
	gro651=gro651[is.finite(gro651)]
	gro1232=gro1232[is.finite(gro1232)]
	gro1252=gro1252[is.finite(gro1252)]
	gro632=gro632[is.finite(gro632)]
	gro652=gro652[is.finite(gro652)]
	
	boxplot(gro1231,gro1232,gro1251,gro1252,gro631,gro632,gro651,gro652,
			main=regionType, xaxt='n',ylab="log2 Change versus WT",ylim=c(-3,3))
	axis(1,at=1:8,labels=c("1231","1232","1251","1252","631","632","651","652"))
	abline(h=0,col="lightblue",lty=2)
	
	tResults = c(t.test(gro1231,mu=0)$p.value,
				 t.test(gro1232,mu=0)$p.value,
				 t.test(gro1251,mu=0)$p.value,
				 t.test(gro1252,mu=0)$p.value,
				 t.test(gro631,mu=0)$p.value,
				 t.test(gro632,mu=0)$p.value,
				 t.test(gro651,mu=0)$p.value,
				 t.test(gro652,mu=0)$p.value)
	ttest[[regionType]] = tResults
}

for (regionType in c("RY_enh","SE_indiv")) {
	fn = paste(regionType, "_list.txt",sep='')
	data = read.delim(fn,header=TRUE,row.names=1,colClasses=c("character",rep("numeric",36)))
	
	sense = as.character(sapply(aso, function(x) paste(x,"_sense",sep='')))
	antisense = as.character(sapply(aso, function(x) paste(x,"_antisense",sep='')))
	
	data = data[,sense] + data[,antisense]
	names(data) = aso
	
	gro1231 = log2(data$GRO_1231/data$GRO_12C1)
	gro1251 = log2(data$GRO_1251/data$GRO_12C1)
	gro631 = log2(data$GRO_631/data$GRO_6C1)
	gro651 = log2(data$GRO_651/data$GRO_6C1)	
	gro1232 = log2(data$GRO_1232/data$GRO_12C2)
	gro1252 = log2(data$GRO_1252/data$GRO_12C2)
	gro632 = log2(data$GRO_632/data$GRO_6C2)
	gro652 = log2(data$GRO_652/data$GRO_6C2)
	
	gro1231=gro1231[is.finite(gro1231)]
	gro1251=gro1251[is.finite(gro1251)]
	gro631=gro631[is.finite(gro631)]
	gro651=gro651[is.finite(gro651)]
	gro1232=gro1232[is.finite(gro1232)]
	gro1252=gro1252[is.finite(gro1252)]
	gro632=gro632[is.finite(gro632)]
	gro652=gro652[is.finite(gro652)]
	
	boxplot(gro1231,gro1232,gro1251,gro1252,gro631,gro632,gro651,gro652,
			main=regionType, xaxt='n',ylab="log2 Change versus WT",ylim=c(-3,3))
	axis(1,at=1:8,labels=c("1231","1232","1251","1252","631","632","651","652"))
	abline(h=0,col="lightblue",lty=2)
	
	tResults = c(t.test(gro1231,mu=0)$p.value,
				 t.test(gro1232,mu=0)$p.value,
				 t.test(gro1251,mu=0)$p.value,
				 t.test(gro1252,mu=0)$p.value,
				 t.test(gro631,mu=0)$p.value,
				 t.test(gro632,mu=0)$p.value,
				 t.test(gro651,mu=0)$p.value,
				 t.test(gro652,mu=0)$p.value)
	ttest[[regionType]] = tResults

}

dev.off()

ttest=do.call(rbind, ttest)
write.table(ttest,"ttest_reps.txt",sep='\t',quote=FALSE)

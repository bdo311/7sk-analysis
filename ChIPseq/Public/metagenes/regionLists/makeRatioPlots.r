comp = read.table("comparisons_with7sk_allchip.txt",header=FALSE)
#comp = read.table("comparisons_with7sk_selectchip.txt",header=FALSE)
denom = read.table("denominators_allchip.txt",header=FALSE)

# large boxplots
regionTypes = c("TSS_trans_paused","SE_indiv","RY_enh")
comparisons = list()

for (regionType in regionTypes) {
	print(regionType)
	fn = paste(regionType, "_list_with7sk_atac.txt",sep='')
	data = read.delim(fn,header=TRUE,row.names=1,colClasses=c("character",rep("numeric",74)))
	
	newDataList = list()
	for (i in 1:nrow(comp)) {
		top = as.character(comp[i,1])
		bottom = as.character(comp[i,2])
		if (substring(top,1,1)=='7') top=paste("X",top,sep='')
		if (substring(bottom,1,1)=='7') bottom=paste("X",bottom,sep='')
		if (bottom=='none') {
			x = data[,top]*5
		} else {
			x = data[,top]/data[,bottom]
		}
		x[!is.finite(x)] = 0
		newDataList[[top]] = x	
	}
	
	# pdf(paste(regionType,"_boxplot.pdf",sep=''),width=20,height=6)
	# par(mar=c(7,4,4,2))
	for (i in 1:nrow(denom)){
		bottom = as.character(denom[i,1])
		if (!(bottom %in% names(comparisons))) comparisons[[bottom]] = list()
		divided = lapply(newDataList,function(x) log2(x/newDataList[[bottom]]))
		divided_zero = lapply(divided,function(x) x[is.finite(x)])
		medians = lapply(divided_zero, median)
		for(x in names(medians)) {
			#print(c(bottom, x, medians[[x]]))
			if (!(x %in% names(comparisons[[bottom]]))) comparisons[[bottom]][[x]] = vector()

			comparisons[[bottom]][[x]] = c(comparisons[[bottom]][[x]],medians[[x]])	
			
		}
		# boxplot(divided_zero,main=bottom,xaxt='n')			
		# axis(1,1:length(divided_zero),labels=names(newDataList),las=2,cex.axis=0.7)
		# abline(h=0,col="blue",lty=2)
	}	

	# dev.off()
}

makeDataFrame = function(x) {
	a=t(as.data.frame(x))
	colnames(a) = regionTypes
	return(a)
}

# heatmap: uses boxplot info above
comparisons.df = lapply(comparisons, makeDataFrame)
#comparisons.df = lapply(comparisons.df, function(x) colnames(x) = regionTypes)
require(pheatmap)
require(RColorBrewer)
colorRamp = colorRampPalette(c("blue","white","red"))

#pdf("ratio_heatmaps_with7sk_atac_selectchip.pdf",width=4,height=5)
pdf("ratio_heatmaps_with7sk_atac_allchip.pdf",width=4,height=10)

for (i in 1:nrow(denom)) {
	d = as.character(denom[i,1])
	pheatmap(comparisons.df[[d]],main=d,scale="none",col=colorRamp(50),breaks=seq(-4,4,8/50))
	#for allchip (bigger)
	#pheatmap(comparisons.df[[d]],cluster_rows=FALSE,cluster_cols=FALSE,main=d,scale="none",col=colorRamp(50),breaks=seq(-2,2,4/50)) #for selectchip (smaller)

}
dev.off()


### PCA, based on allchip
regionTypes = c("TSS_trans_paused","RY_enh","SE_indiv")

newDataList = list()
dataTypes = list()
for (regionType in regionTypes) {
	print(regionType)
	fn = paste(regionType, "_list_with7sk_atac.txt",sep='')
	data = read.delim(fn,header=TRUE,row.names=1,colClasses=c("character",rep("numeric",74)))
	dataTypes[[regionType]] = nrow(data)
	
	for (i in 1:nrow(comp)) {
		top = as.character(comp[i,1])
		bottom = as.character(comp[i,2])
		if (substring(top,1,1)=='7') top=paste("X",top,sep='')
		if (substring(bottom,1,1)=='7') bottom=paste("X",bottom,sep='')
		if (bottom=='none') {
			x = data[,top]*5
		} else {
			x = data[,top]/data[,bottom]
		}
		x[!is.finite(x)] = 0
		newDataList[[top]] = c(newDataList[[top]],x)
	}
}

newData = as.data.frame(newDataList)

#pca
addTrans <- function(color,trans)
{
  # This function adds transparency to a color.
  # Define transparency with an integer between 0 and 255
  # 0 being fully transparent and 255 being fully visible
  # Works with either color and trans a vector of equal length,
  # or one of the two of length 1.

  if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
  if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
  if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))

  num2hex <- function(x)
  {
    hex <- unlist(strsplit("0123456789ABCDEF",split=""))
    return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
  }
  rgb <- rbind(col2rgb(color),trans)
  res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
  return(res)
}

#TR for super enhancers
data=read.delim("/home/raflynn/7SK/GROseq/Flynn/grochg/se_tr_withnames.txt",header=TRUE,row.names=1)
tr=log2(data$GRO_125/data$GRO_12C)
names(tr)=rownames(data)

#reorder TR list for our super enhancers
fn = "SE_indiv_list_with7sk_atac.txt"
data = read.delim(fn,header=TRUE,row.names=1,colClasses=c("character",rep("numeric",74)))
tr=tr[rownames(data)]

# tr_bin=0
# tr_bin[abs(tr)<0.5]=3 #3 for >0
# tr_bin[abs(tr)>0.5]=4 #4 for <0

#whiten the data by making all factors mean 0
newData.log = apply(newData,c(1,2),function(x) ifelse(x!=0,log2(x),0))
newData.log=as.data.frame(newData.log)
newData.means = apply(newData.log,2,mean)

newData.logm = as.data.frame(t(apply(newData.log,1,function(x) x - newData.means)))
#dataList = c(rep(1,dataTypes[["TSS_trans_paused"]]),rep(2,dataTypes[["RY_enh"]]),rep(3,dataTypes[["SE_indiv"]]))
dataList = c(rep(1,dataTypes[["TSS_trans_paused"]]),rep(2,dataTypes[["RY_enh"]]),tr_bin) #use tr to mark SE TRs

pr = prcomp(newData.logm) #pca with RAW VALUES
pdf("first2PCs_raw_allchip_tr.pdf",width=7,height=7)
plot(pr$x[,1:2],col=addTrans(dataList,200),pch=19,cex=0.7,main="First 2 PCs")
legend("topleft",legend=c("TSS","Super enhancer","Regular enhancer"),col=c("black","green","red"),pch=19)
dev.off()

# pca with RATIOS
newData.log.scannell = apply(newData.log,2,function(x)x-newData.log$Scannell_TBP)
newData.means = apply(newData.log.scannell,2,mean)
newData.logm.scannell = as.data.frame(t(apply(newData.log.scannell,1,function(x) x - newData.means)))
pr = prcomp(newData.logm.scannell) 
pdf("first2PCs_scannellTBP_allchip_tr.pdf",width=7,height=7)
plot(pr$x[,1:2],col=addTrans(dataList,200),pch=19,cex=0.7,main="First 2 PCs")
legend("topleft",legend=c("TSS","Super enhancer","Regular enhancer"),col=c("black","green","red"),pch=19)
dev.off()

newData.log.rahl = apply(newData.log,2,function(x)x-newData.log$Rahl_TBP)
newData.means = apply(newData.log.rahl,2,mean)
newData.logm.rahl = as.data.frame(t(apply(newData.log.rahl,1,function(x) x - newData.means)))
pr = prcomp(newData.logm.rahl) 
pdf("first2PCs_rahlTBP_allchip_tr.pdf",width=7,height=7)
plot(pr$x[,1:2],col=addTrans(dataList,200),pch=19,cex=0.7,main="First 2 PCs")
legend("topleft",legend=c("TSS","Super enhancer","Regular enhancer"),col=c("black","green","red"),pch=19)
dev.off()


library(ggplot2)
library(RColorBrewer)
se.pca = pr$x[(nrow(pr$x)-559):nrow(pr$x),1:2]
se.pca=cbind(se.pca,tr)
colnames(se.pca)=c("pc1","pc2","tr")
write.table(se.pca,"se_pca.txt",sep='\t')
qplot(se.pca,color=tr)+scale_color_brewer()

# get loadings
pdf("loadings.pdf",width=20,height=10)
par(mar=c(13,4,2,1),las=2)
barplot(t(pr$rotation[,1:2]),beside=TRUE)
#axis(1,at=1:45,labels=rownames(pr$rotation),las=2)
dev.off()

# make boxplot
pdf("boxplot_firstPC.pdf",width=5,height=5)
boxplot(pr$x[1:dataTypes[["TSS_trans_paused"]],1],
pr$x[(dataTypes[["TSS_trans_paused"]]+1+dataTypes[["RY_enh"]]+1):nrow(pr$x),1],
pr$x[(dataTypes[["TSS_trans_paused"]]+1):(dataTypes[["TSS_trans_paused"]]+1+dataTypes[["RY_enh"]]),1],
ylim=c(-20,20),
xaxt='n')
axis(1,at=1:3,labels=c("TSS","Super\nEnhancer","Regular\nEnhancer"))
dev.off()




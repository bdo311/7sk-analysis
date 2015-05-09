# makeScatter.r

se=read.delim("se_info.txt",header=TRUE)
se$conv_chg = log2(se$conv_aso/se$conv_wt)

tss=read.delim("tss_info.txt",header=TRUE)
tss$conv_chg = log2(tss$conv_aso/tss$conv_wt)

re=read.delim("re_info.txt",header=TRUE)
re$conv_chg = log2(re$conv_aso/re$conv_wt)

# ctcf=read.delim("ctcf_info.txt",header=TRUE)
# ctcf$conv_chg = log2(ctcf$conv_aso/ctcf$conv_wt)

# rand3=read.delim("rand3_info.txt",header=TRUE)
# rand3$conv_chg = log2(rand3$conv_aso/rand3$conv_wt)

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

pdf("gro_vs_ct_chg.pdf",width=6,height=6)
logt = log2(tss$gro)
logr = log2(re$gro)
logs = log2(se$gro)
plot(logt,tss$conv_chg,pch=19,cex=0.5,col=addTrans("black",100),xlim=c(0,7),ylim=c(-4,4),xlab="log2 GRO read density",ylab="CT log2 fold change",main="GROseq density (WT) vs CT change")
points(logr,re$conv_chg,pch=19,cex=0.5,col=addTrans("red",100))
#points(log2(ctcf$gro),ctcf$conv_chg,pch=19,cex=0.4,col=addTrans("blue",40))
#points(log2(rand3$gro),rand3$conv_chg,pch=19,cex=0.5,col=addTrans("blue",100))
points(logs,se$conv_chg,pch=19,cex=0.5,col=addTrans("green",250))
abline(h=0,lty=6,lwd=2, col="cyan")
points(median(logs[is.finite(logs)]),median(se$conv_chg[is.finite(se$conv_chg)]),pch=21,cex=4,col="yellow",bg="darkgreen")
points(median(logr[is.finite(logr)]),median(re$conv_chg[is.finite(re$conv_chg)]),pch=21,cex=4,col="yellow",bg="darkred")
points(median(logt[is.finite(logt)]),median(tss$conv_chg[is.finite(tss$conv_chg)]),pch=21,cex=4,col="yellow",bg="gray")
#points(median(is.finite(log2(rand3$gro[rand3$gro]))),median(rand3$conv_chg[is.finite(rand3$conv_chg)]),pch=21,cex=4,col="yellow",bg="blue")
abline(h=0,lty=6,lwd=2, col="cyan")
legend("topright",legend=c("TSS","Super enhancer","Regular enhancer"),col=c("black","green","red"),pch=19)
dev.off()

pdf("ct_vs_ct_chg.pdf",width=6,height=6)
logt = log2(tss$conv_wt)
logr = log2(re$conv_wt)
logs = log2(se$conv_wt)
plot(logt,tss$conv_chg,pch=19,cex=0.5,col=addTrans("black",100),xlim=c(-8,3),ylim=c(-4,4),xlab="log2 CT",ylab="CT log2 fold change",main="CT (WT) vs CT change")
points(logr,re$conv_chg,pch=19,cex=0.5,col=addTrans("red",100))
#points(log2(ctcf$conv_wt),ctcf$conv_chg,pch=19,cex=0.4,col=addTrans("blue",40))
#points(log2(rand3$conv_wt),rand3$conv_chg,pch=19,cex=0.5,col=addTrans("blue",100))
points(logs,se$conv_chg,pch=19,cex=0.5,col=addTrans("green",250))
abline(h=0,lty=6,lwd=2, col="cyan")
points(median(logs[is.finite(logs)]),median(se$conv_chg[is.finite(se$conv_chg)]),pch=21,cex=4,col="yellow",bg="darkgreen")
points(median(logr[is.finite(logr)]),median(re$conv_chg[is.finite(re$conv_chg)]),pch=21,cex=4,col="yellow",bg="darkred")
points(median(logt[is.finite(logt)]),median(tss$conv_chg[is.finite(tss$conv_chg)]),pch=21,cex=4,col="yellow",bg="gray")
#points(median(is.finite(log2(rand3$conv_wt[rand3$conv_wt]))),median(rand3$conv_chg[is.finite(rand3$conv_chg)]),pch=21,cex=4,col="yellow",bg="blue")
abline(h=0,lty=6,lwd=2, col="cyan")
legend("topright",legend=c("TSS","Super enhancer","Regular enhancer"),col=c("black","green","red"),pch=19)
dev.off()
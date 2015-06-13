re=read.delim("re_matrix.txt",header=TRUE,row.names=1)
se=read.delim("se_matrix.txt",header=TRUE,row.names=1)

re$rna = log2(with(re,(SRR1265785_rmdup+SRR1265787_rmdup)/(SRR1265786_rmdup+SRR1265788_rmdup)))
re$gro = log2(with(re,(GRO_125comb_sense+GRO_125comb_antisense)/(GRO_12Ccomb_sense+GRO_12Ccomb_antisense)))

se$rna = log2(with(se,(SRR1265785_rmdup+SRR1265787_rmdup)/(SRR1265786_rmdup+SRR1265788_rmdup)))
se$gro = log2(with(se,(GRO_125comb_sense+GRO_125comb_antisense)/(GRO_12Ccomb_sense+GRO_12Ccomb_antisense)))

re.stats = re[is.finite(re$rna) & is.finite(re$gro),][,c("rna","gro")]
se.stats = se[is.finite(se$rna) & is.finite(se$gro),][,c("rna","gro")]

cor.test(re.stats$rna,re.stats$gro)
cor.test(se.stats$rna,se.stats$gro)

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

pdf("rnaseq_vs_groseq.pdf",width=5,height=5)
plot(re.stats$rna,re.stats$gro,cex=0.4,pch=19,col=addTrans("red",100),xlab="Brg1 KD RNAseq log2 change",ylab="7SK ASO GROseq log2 change",xlim=c(-4,4),ylim=c(-4,4))
points(se.stats$rna,se.stats$gro,cex=0.4,pch=19,col=addTrans("green",150))
abline(v=0,h=0,lty=2,lwd=4,col="lightblue")

# Regression lines
re.reg = lm(re.stats$gro ~ re.stats$rna)
se.reg = lm(se.stats$gro ~ se.stats$rna)
abline(re.reg,col="firebrick",lwd=4)
abline(se.reg,col="darkgreen",lwd=4)
dev.off()

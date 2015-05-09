# fn="promoter_tr.txt"
# ofn="promoter.pdf"

fn="se_tr.txt"
ofn="se.pdf"

data=read.delim(fn,header=FALSE,row.names=1)
gro123=as.numeric(log2(data["GRO_123",]))
gro125=as.numeric(log2(data["GRO_125",]))
gro12C=as.numeric(log2(data["GRO_12C",]))
gro63=as.numeric(log2(data["GRO_63",]))
gro65=as.numeric(log2(data["GRO_65",]))
gro6C=as.numeric(log2(data["GRO_6C",]))


gro12C=gro12C[gro12C!=0 & is.finite(gro12C)]
gro125=gro125[gro125!=0 & is.finite(gro125)]
gro6C=gro6C[gro6C!=0 & is.finite(gro6C)]
gro65=gro65[gro65!=0  & is.finite(gro65)]

pdf(ofn,width=7,height=10)
par(mfrow=c(2,1))
Fn1=ecdf(gro12C)
Fn2=ecdf(gro125)
plot(Fn2,col="red",xlim=c(-4,9),xlab="log2 TR",main="Gene TR 12hr ASO")
plot(Fn1,col="black",add=TRUE)
legend("topleft",legend=c("12hr Ctrl","12hr 5'ASO"),col=c("black","red"),lty=1)
Fn1=ecdf(gro6C)
Fn2=ecdf(gro65)
plot(Fn2,col="red",xlim=c(-4,9),xlab="log2 TR",main="Gene TR 6hr ASO")
plot(Fn1,col="black",add=TRUE)
legend("topleft",legend=c("6hr Ctrl","6hr 5'ASO"),col=c("black","red"),lty=1)
dev.off()





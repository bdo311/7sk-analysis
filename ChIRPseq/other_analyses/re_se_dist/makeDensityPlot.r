data=read.delim("re_se_dist.txt")
re=log10(data[data[,7]=="RE",9])
se=log10(data[data[,7]=="SE",9])
re=re[is.finite(re)]
se=se[is.finite(se)]

pdf("peak_distance.pdf",width=5,height=5)
plot(density(re),col="red",xlim=c(2,7),xlab="Peak distance",xaxt='n',lwd=2,main="Peak-to-peak distance")
lines(density(se),col="green",lwd=2)
legend("topleft",legend=c("RE","SE"),col=c("red","green"),lwd=2)
axis(1,at=2:7,labels=c(100,1000,10000,100000,1000000,10000000))
dev.off()
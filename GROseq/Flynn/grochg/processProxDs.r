# processProxDs.r

promProx = read.delim("promoter_prox_withnames.txt",row.names=1,header=TRUE)
promDs = read.delim("promoter_ds_withnames.txt",row.names=1,header=TRUE)

reProx = read.delim("re_prox_withnames.txt",row.names=1,header=TRUE)
reDs = read.delim("re_ds_withnames.txt",row.names=1,header=TRUE)

seProx = read.delim("se_prox_withnames.txt",row.names=1,header=TRUE)
seDs = read.delim("se_ds_withnames.txt",row.names=1,header=TRUE)

# ## Boxplots

# pdf("proximal_distal_density.pdf",width=5,height=8)
# par(mfrow=c(2,1))
# boxplot(log2(promProx$GRO_125),log2(promProx$GRO_12C),log2(reProx$GRO_125),log2(reProx$GRO_12C),log2(seProx$GRO_125),log2(seProx$GRO_12C),xaxt='n',ylab="log2 density",main="Proximal read density",col=c("gray","gray","red","red","green","green"))
# axis(1,at=1:6,labels=c("Prom\nASO","Prom\nWT","RE\nASO","RE\nWT","SE\nASO","SE\nWT"))
# boxplot(log2(promDs$GRO_125),log2(promDs$GRO_12C),log2(reDs$GRO_125),log2(reDs$GRO_12C),log2(seDs$GRO_125),log2(seDs$GRO_12C),xaxt='n',ylab="log2 density",main="Downstream read density",col=c("gray","gray","red","red","green","green"))
# axis(1,at=1:6,labels=c("Prom\nASO","Prom\nWT","RE\nASO","RE\nWT","SE\nASO","SE\nWT"))
# dev.off()


## New numbers from Start-seq centering

median(promProx$GRO_125)	#4.68
median(promProx$GRO_12C)	#5.42
median(reProx$GRO_125)	#0.55
median(reProx$GRO_12C)	#0.36
median(seProx$GRO_125)	#3.74
median(seProx$GRO_12C)	#2.76

median(promDs$GRO_125)	#0.75
median(promDs$GRO_12C)	#0.75
median(reDs$GRO_125)	#0.25
median(reDs$GRO_12C)	#0.12
median(seDs$GRO_125)	#3.33
median(seDs$GRO_12C)	#1.89


mean(promProx$GRO_125)	#8.33
mean(promProx$GRO_12C)	#9.28
mean(reProx$GRO_125)	#1.85
mean(reProx$GRO_12C)	#1.52
mean(seProx$GRO_125)	#5.38
mean(seProx$GRO_12C)	#4.34

mean(promDs$GRO_125)	#1.92
mean(promDs$GRO_12C)	#1.94
mean(reDs$GRO_125)	#0.92
mean(reDs$GRO_12C)	#0.68
mean(seDs$GRO_125)	#4.78
mean(seDs$GRO_12C)	#2.99

ks.test(promProx$GRO_125,promProx$GRO_12C)	#2.1e-13
ks.test(reProx$GRO_125,reProx$GRO_12C)	#<2.2e-16
ks.test(seProx$GRO_125,seProx$GRO_12C)	#1.1e-4

ks.test(promDs$GRO_125,promDs$GRO_12C)	#0.016
ks.test(reDs$GRO_125,reDs$GRO_12C)	#<2.2e-16
ks.test(seDs$GRO_125,seDs$GRO_12C)	#5.9e-14

# Scatterplot/Distributions
tss_prox_chg = log2(promProx$GRO_125/promProx$GRO_12C)
tss_dist_chg = log2(promDs$GRO_125/promDs$GRO_12C)
re_prox_chg = log2(reProx$GRO_125/reProx$GRO_12C)
re_dist_chg = log2(reDs$GRO_125/reDs$GRO_12C)
se_prox_chg = log2(seProx$GRO_125/seProx$GRO_12C)
se_dist_chg = log2(seDs$GRO_125/seDs$GRO_12C)

tss_prox_chg = tss_prox_chg[is.finite(tss_prox_chg)]
re_prox_chg = re_prox_chg[is.finite(re_prox_chg)]
se_prox_chg = se_prox_chg[is.finite(se_prox_chg)]
tss_dist_chg = tss_dist_chg[is.finite(tss_dist_chg)]
re_dist_chg = re_dist_chg[is.finite(re_dist_chg)]
se_dist_chg = se_dist_chg[is.finite(se_dist_chg)]

pdf("prox_ds_changes.pdf",width=5,height=8)
par(mfrow=c(2,1))
plot(density(tss_prox_chg),ylim=c(0,0.4),xlim=c(-7,7),main="Proximal",lwd=2,xlab="log2 change in proximal signal")
lines(density(re_prox_chg),col="red",lwd=2)
lines(density(se_prox_chg),col="green",lwd=2)
legend("topright",legend=c("TSS","RE","SE"),col=c("black","red","green"),lwd=2,lty=1)
abline(v=0)

plot(density(tss_dist_chg),ylim=c(0,0.4),xlim=c(-7,7),main="Distal",lwd=2,xlab="log2 change in distal signal")
lines(density(re_dist_chg),col="red",lwd=2)
lines(density(se_dist_chg),col="green",lwd=2)
legend("topright",legend=c("TSS","RE","SE"),col=c("black","red","green"),lwd=2,lty=1)
abline(v=0)
dev.off()

t.test(tss_prox_chg)	#1.09e-15
t.test(re_prox_chg)	#2.7e-12
t.test(se_prox_chg)	#1.25e-8
t.test(tss_dist_chg)	#0.014
t.test(re_dist_chg)	#<2.2e-16
t.test(se_dist_chg)	#<2.2e-16

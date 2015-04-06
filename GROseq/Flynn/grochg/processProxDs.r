# processProxDs.r

promProx = read.delim("promoter_prox_withnames.txt",row.names=1,header=TRUE)
promDs = read.delim("promoter_ds_withnames.txt",row.names=1,header=TRUE)

reProx = read.delim("re_prox_withnames.txt",row.names=1,header=TRUE)
reDs = read.delim("re_ds_withnames.txt",row.names=1,header=TRUE)

seProx = read.delim("se_prox_withnames.txt",row.names=1,header=TRUE)
seDs = read.delim("se_ds_withnames.txt",row.names=1,header=TRUE)

pdf("proximal_distal_density.pdf",width=5,height=8)
par(mfrow=c(2,1))
boxplot(log2(promProx$GRO_125),log2(promProx$GRO_12C),log2(reProx$GRO_125),log2(reProx$GRO_12C),log2(seProx$GRO_125),log2(seProx$GRO_12C),xaxt='n',ylab="log2 density",main="Proximal read density",col=c("gray","gray","red","red","green","green"))
axis(1,at=1:6,labels=c("Prom\nASO","Prom\nWT","RE\nASO","RE\nWT","SE\nASO","SE\nWT"))
boxplot(log2(promDs$GRO_125),log2(promDs$GRO_12C),log2(reDs$GRO_125),log2(reDs$GRO_12C),log2(seDs$GRO_125),log2(seDs$GRO_12C),xaxt='n',ylab="log2 density",main="Downstream read density",col=c("gray","gray","red","red","green","green"))
axis(1,at=1:6,labels=c("Prom\nASO","Prom\nWT","RE\nASO","RE\nWT","SE\nASO","SE\nWT"))
dev.off()

median(promProx$GRO_125) #7.03
median(promProx$GRO_12C) #8.09
median(reProx$GRO_125) #0.49
median(reProx$GRO_12C) #0.33
median(seProx$GRO_125) #6.78
median(seProx$GRO_12C) #5.29

median(promDs$GRO_125) #1.20
median(promDs$GRO_12C) #1.21
median(reDs$GRO_125) #0.051
median(reDs$GRO_12C) #0.026
median(seDs$GRO_125) #1.28
median(seDs$GRO_12C) #0.68


mean(promProx$GRO_125) #10.49
mean(promProx$GRO_12C) #11.83
mean(reProx$GRO_125) #1.58
mean(reProx$GRO_12C) #1.28
mean(seProx$GRO_125) #9.68
mean(seProx$GRO_12C) #8.10

mean(promDs$GRO_125) #2.53
mean(promDs$GRO_12C) #2.48
mean(reDs$GRO_125) #0.31
mean(reDs$GRO_12C) #0.24
mean(seDs$GRO_125) #2.23
mean(seDs$GRO_12C) #1.44

ks.test(promProx$GRO_125,promProx$GRO_12C) #1.2e-13
ks.test(reProx$GRO_125,reProx$GRO_12C) #<2.2e-16
ks.test(seProx$GRO_125,seProx$GRO_12C) #3.2e-4

ks.test(promDs$GRO_125,promDs$GRO_12C) #0.46
ks.test(reDs$GRO_125,reDs$GRO_12C) #<2.2e-16
ks.test(seDs$GRO_125,seDs$GRO_12C) #2.4e-8

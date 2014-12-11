if (Sys.glob("gro_tr_comb.RData")) {
	load("gro_tr_comb.RData")
} else {
	fn="tr_comb.txt"
	no_col = count.fields(fn, sep = "\t") 
	data = t(read.table(fn, sep='\t', fill=TRUE, row.names=1)) 
}

fn="enh_tr_comb_10000.txt"
no_col = count.fields(fn, sep = "\t") 
data = t(read.table(fn, sep='\t', fill=TRUE, row.names=1)) 

gro_12C = data[,"GRO_12Ccomb"]
gro_12C = gro_12C[!is.na(gro_12C)]
gro_123 = data[,"GRO_123comb"]
gro_123 = gro_123[!is.na(gro_123)]
gro_125 = data[,"GRO_125comb"]
gro_125 = gro_125[!is.na(gro_125)]
gro_6C = data[,"GRO_6Ccomb"]
gro_6C = gro_6C[!is.na(gro_6C)]
gro_63 = data[,"GRO_63comb"]
gro_63 = gro_63[!is.na(gro_63)]
gro_65 = data[,"GRO_65comb"]
gro_65 = gro_65[!is.na(gro_65)]

pdf("gro_tr_comb_enh_10000_boxplot.pdf", width=8, height=5)
par(mfrow=c(1,2))
boxplot(log2(gro_12C), log2(gro_123), log2(gro_125), xaxt='n', ylab="log2 TR", main="12hr ASO", ylim=c(-4,8))
axis(1, at=1:3, labels=c("Control", "3' ASO", "5' ASO"))
abline(h=median(log2(gro_12C)), lty=2, col="red")

boxplot(log2(gro_6C), log2(gro_63), log2(gro_65), xaxt='n', ylab="log2 TR", main="6hr ASO", ylim=c(-4,8))
axis(1, at=1:3, labels=c("Control", "3' ASO", "5' ASO"))
abline(h=median(log2(gro_6C)), lty=2, col="red")
dev.off()

pdf("gro_tr_comb_enh_10000_ecdf.pdf", width=8, height=8)
par(mfrow=c(2, 1))
Fn1 = ecdf(log2(gro_12C))
Fn2 = ecdf(log2(gro_123))
Fn3 = ecdf(log2(gro_125))
plot(Fn1, xlab="log2 Traveling Ratio", ylab="Cumulative Proportion of genes", main="12hr ASO", xlim=c(-8,8))
lines(Fn2, col="red")
lines(Fn3, col="blue")
legend("topleft", legend=c("Control", "3' ASO", "5' ASO"), col=c("black", "red", "blue"), lty=1)

Fn1 = ecdf(log2(gro_6C))
Fn2 = ecdf(log2(gro_63))
Fn3 = ecdf(log2(gro_65))
plot(Fn1, xlab="log2 Traveling Ratio", ylab="Cumulative Proportion of genes", main="6hr ASO", xlim=c(-8,8))
lines(Fn2, col="red")
lines(Fn3, col="blue")
legend("topleft", legend=c("Control", "3' ASO", "5' ASO"), col=c("black", "red", "blue"),lty=1)
dev.off()

save.image("gro_tr_comb_enh_10000.RData")

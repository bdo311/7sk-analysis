# reading in files
sampleList = list()
samples = c("123", "125", "12C", "63", "65", "6C")
for (ns in 1:length(samples)) {
	print(ns)
	fn = paste("GRO_", samples[ns], "_repeat_norm.bedGraph", sep='')
	data = read.delim(fn, colClasses=c("character", "numeric", "numeric", "numeric"), header=FALSE)
	vals = list()
	for (row in 1:nrow(data)) {
		start = data[row,2]
		end = data[row,3]
		val = data[row,4]
		while(end > start) {
			vals[[start + 1]] = val
			start = start + 1
		}
	}
	vals = sapply(vals, function(x)ifelse(is.null(x), 0, x)) #bedgraphs omit positions where read density is 0, so I need to put them back in
	sampleList[[ns]] = vals
}

df <- data.frame(matrix(unlist(sampleList), ncol=length(samples), byrow=FALSE))
names(df) = samples


## repeats
repeatPos = read.delim("Mm_repeatIndex_positions.txt", colClasses=c("character", "numeric", "numeric", "numeric"), header=FALSE, row.names=1)
names(repeatPos) = c("len", "start", "stop")
pdf("12hr_repeats.pdf", width=15, height=10)
par(mfrow=c(4,5))
for (i in 1:nrow(repeatPos)) {
	start = repeatPos$start[i] + 5
	stop = repeatPos$stop[i] - 5
	maxY = max(c(df$"12C"[start:stop], df$"123"[start:stop], df$"125"[start:stop])) * 1.3
	plot(df$"12C"[start:stop], type='l', main=rownames(repeatPos)[i], ylim=c(0, maxY), xlab="Position", ylab="GRO-seq density")
	lines(df$"123"[start:stop], col="red")
	lines(df$"125"[start:stop], col="blue")
	legend("topright", legend=c("Control", "5' ASO", "3' ASO"), lty=rep(1, 3), col=c("black", "blue", "red"))
}
dev.off()

pdf("6hr_repeats.pdf", width=15, height=10)
par(mfrow=c(4,5))
for (i in 1:nrow(repeatPos)) {
	start = repeatPos$start[i] + 5
	stop = repeatPos$stop[i] - 5
	maxY = max(c(df$"6C"[start:stop], df$"63"[start:stop], df$"65"[start:stop])) * 1.3
	plot(df$"6C"[start:stop], type='l', main=rownames(repeatPos)[i], ylim=c(0, maxY), xlab="Position", ylab="GRO-seq density")
	lines(df$"63"[start:stop], col="red")
	lines(df$"65"[start:stop], col="blue")
	legend("topright", legend=c("Control", "5' ASO", "3' ASO"), lty=rep(1, 3), col=c("black", "blue", "red"))
}
dev.off()

## rRNA
rRNA_start = repeatPos$start[17]
start18s=rRNA_start +4007
end18s=rRNA_start +5876
start5s=rRNA_start +6877
end5s=rRNA_start +7033
start28s=rRNA_start +8123
end28s=rRNA_start +12836
rRNAend=rRNA_start +13401

# whole rRNA
pdf("12hr_rrna.pdf", width=8, height=3)
maxY = max(c(df$"12C"[rRNA_start:rRNAend], df$"123"[rRNA_start:rRNAend], df$"125"[rRNA_start:rRNAend])) * 1.3
plot(df$"12C"[rRNA_start:rRNAend], type='l', main="rRNA", ylim=c(0, maxY), xlab="Position", ylab="GRO-seq density")
lines(df$"123"[rRNA_start:rRNAend], col="red")
lines(df$"125"[rRNA_start:rRNAend], col="blue")
abline(v=c(start18s, end18s, start5s, end5s, start28s, end28s)-rRNA_start, col="green")
legend("topleft", legend=c("Control", "5' ASO", "3' ASO"), lty=rep(1, 3), col=c("black", "blue", "red"))
dev.off()

pdf("6hr_rrna.pdf", width=8, height=3)
maxY = max(c(df$"6C"[rRNA_start:rRNAend], df$"63"[rRNA_start:rRNAend], df$"65"[rRNA_start:rRNAend])) * 1.3
plot(df$"6C"[rRNA_start:rRNAend], type='l', main="rRNA", ylim=c(0, maxY), xlab="Position", ylab="GRO-seq density")
lines(df$"63"[rRNA_start:rRNAend], col="red")
lines(df$"65"[rRNA_start:rRNAend], col="blue")
abline(v=c(start18s, end18s, start5s, end5s, start28s, end28s)-rRNA_start, col="green")
legend("topleft", legend=c("Control", "5' ASO", "3' ASO"), lty=rep(1, 3), col=c("black", "blue", "red"))
dev.off()


# regions of rRNA
pdf("12hr_rrna_regions.pdf", width=8, height=3)
par(mfrow=c(1, 3))

maxY = max(c(df$"12C"[start18s:end18s], df$"123"[start18s:end18s], df$"125"[start18s:end18s])) * 1.3
plot(df$"12C"[start18s:end18s], type='l', main="18S rRNA", ylim=c(0, maxY), xlab="Position", ylab="GRO-seq density")
lines(df$"123"[start18s:end18s], col="red")
lines(df$"125"[start18s:end18s], col="blue")
legend("topright", legend=c("Control", "5' ASO", "3' ASO"), lty=rep(1, 3), col=c("black", "blue", "red"))

maxY = max(c(df$"12C"[start5s:end5s], df$"123"[start5s:end5s], df$"125"[start5s:end5s])) * 1.3
plot(df$"12C"[start5s:end5s], type='l', main="5S rRNA", ylim=c(0, maxY), xlab="Position", ylab="GRO-seq density")
lines(df$"123"[start5s:end5s], col="red")
lines(df$"125"[start5s:end5s], col="blue")
legend("topright", legend=c("Control", "5' ASO", "3' ASO"), lty=rep(1, 3), col=c("black", "blue", "red"))

maxY = max(c(df$"12C"[start28s:end28s], df$"123"[start28s:end28s], df$"125"[start28s:end28s])) * 1.3
plot(df$"12C"[start28s:end28s], type='l', main="28s rRNA", ylim=c(0, maxY), xlab="Position", ylab="GRO-seq density")
lines(df$"123"[start28s:end28s], col="red")
lines(df$"125"[start28s:end28s], col="blue")
legend("topright", legend=c("Control", "5' ASO", "3' ASO"), lty=rep(1, 3), col=c("black", "blue", "red"))

dev.off()


pdf("6hr_rrna_regions.pdf", width=8, height=3)
par(mfrow=c(1, 3))

maxY = max(c(df$"6C"[start18s:end18s], df$"63"[start18s:end18s], df$"65"[start18s:end18s])) * 1.3
plot(df$"6C"[start18s:end18s], type='l', main="18S rRNA", ylim=c(0, maxY), xlab="Position", ylab="GRO-seq density")
lines(df$"63"[start18s:end18s], col="red")
lines(df$"65"[start18s:end18s], col="blue")
legend("topright", legend=c("Control", "5' ASO", "3' ASO"), lty=rep(1, 3), col=c("black", "blue", "red"))

maxY = max(c(df$"6C"[start5s:end5s], df$"63"[start5s:end5s], df$"65"[start5s:end5s])) * 1.3
plot(df$"6C"[start5s:end5s], type='l', main="5S rRNA", ylim=c(0, maxY), xlab="Position", ylab="GRO-seq density")
lines(df$"63"[start5s:end5s], col="red")
lines(df$"65"[start5s:end5s], col="blue")
legend("topright", legend=c("Control", "5' ASO", "3' ASO"), lty=rep(1, 3), col=c("black", "blue", "red"))

maxY = max(c(df$"6C"[start28s:end28s], df$"63"[start28s:end28s], df$"65"[start28s:end28s])) * 1.3
plot(df$"6C"[start28s:end28s], type='l', main="28s rRNA", ylim=c(0, maxY), xlab="Position", ylab="GRO-seq density")
lines(df$"63"[start28s:end28s], col="red")
lines(df$"65"[start28s:end28s], col="blue")
legend("topright", legend=c("Control", "5' ASO", "3' ASO"), lty=rep(1, 3), col=c("black", "blue", "red"))

dev.off()

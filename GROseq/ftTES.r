# ftTES.r
# 10/5/14
# looks at GRO data in ASO/not in ASO to see whether there are changes at 3' to FT genes

annot = read.delim("/home/raflynn/ChIRPseq/genes/annotated_genes.txt", sep='\t', header=TRUE)
div = as.character(annot[annot$udRNA.lfc!=0,]$gene)
ft = as.character(annot[annot$termination.defect!=0,]$gene)

setwd("combined")

process = function(treatment) {
	fn = paste(treatment, '_tes_sense.txt', sep='')
	data = t(read.table(fn))
	colnames(data) = data[1,]
	data = data[8:nrow(data),]
	data = apply(data, c(1,2), as.numeric)
	
	ftcols = data[,which(colnames(data) %in% ft)]
	restcols = data[,which(!(colnames(data) %in% ft))]

	return(list(ftcols, restcols))
}

twelvec = process("GRO_12C")
ft.12c = twelvec[[1]]
rest.12c = twelvec[[2]]
twelve3 = process("GRO_123")
ft.123 = twelve3[[1]]
rest.123 = twelve3[[2]]
twelve5 = process("GRO_125")
ft.125 = twelve5[[1]]
rest.125 = twelve5[[2]]

sixc = process("GRO_6C")
ft.6c = sixc[[1]]
rest.6c = sixc[[2]]
six3 = process("GRO_63")
ft.63 = six3[[1]]
rest.63 = six3[[2]]
six5 = process("GRO_65")
ft.65 = six5[[1]]
rest.65 = six5[[2]]

pdf("ft.pdf",width=10, height=6)
par(mfrow=c(2,2))
plot(apply(ft.12c, 1, mean), type='l', ylim=c(0,8), main="12hr GRO FT genes", ylab=c("GROseq density"))
lines(apply(ft.123, 1, mean), col="blue")
lines(apply(ft.125, 1, mean), col="red")
legend("topleft", legend=c("Ctrl", "3'ASO", "5'ASO"),lty=1, col=c("black","blue","red"))
plot(apply(ft.6c, 1, mean), type='l', ylim=c(0,8), main="6hr GRO FT genes", ylab=c("GROseq density"))
lines(apply(ft.63, 1, mean), col="blue")
lines(apply(ft.65, 1, mean), col="red")
legend("topleft", legend=c("Ctrl", "3'ASO", "5'ASO"),lty=1, col=c("black","blue","red"))


plot(apply(rest.12c, 1, mean), type='l', ylim=c(0,3), main="12hr GRO Non-FT genes", ylab=c("GROseq density"))
lines(apply(rest.123, 1, mean), col="blue")
lines(apply(rest.125, 1, mean), col="red")
legend("topleft", legend=c("Ctrl", "3'ASO", "5'ASO"),lty=1, col=c("black","blue","red"))

plot(apply(rest.6c, 1, mean), type='l', ylim=c(0,3), main="6hr GRO Non-FT genes", ylab=c("GROseq density"))
lines(apply(rest.63, 1, mean), col="blue")
lines(apply(rest.65, 1, mean), col="red")
legend("topleft", legend=c("Ctrl", "3'ASO", "5'ASO"),lty=1, col=c("black","blue","red"))
dev.off()

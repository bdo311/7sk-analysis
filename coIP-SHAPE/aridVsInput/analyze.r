data = read.delim("gene_exp.diff", header=TRUE)
data.tr = data[data$value_1>0.1 & data$value_2>0.1,]
data.tr = data.tr[order(data.tr$log2.fold_change., decreasing=TRUE),]
plot(log2(data.tr$value_1), log2(data.tr$value_2))

arid = read.delim("aridVsInput_allreads/gene_exp.diff", header=TRUE, row.names=1)
igg = read.delim("iggVsInput_allreads/gene_exp.diff", header=TRUE, row.names=1)

df = data.frame(arid$gene, arid_arid = arid$value_1, arid_input = arid$value_2, arid_lfc = -arid$log2.fold_change., igg_igg = igg$value_1, igg_input = igg$value_2, igg_lfc = -igg$log2.fold_change.)
df = df[df[,2]>0.1 & df[,3]>0.1 & df[,5]>0.1  & df[,6]>0.1 ,] #15745 genes left
norm.lfc = df$arid_lfc - df$igg_lfc
df = cbind(df, norm.lfc )
df = df[order(df$norm.lfc, decreasing=TRUE),]
write.table(df, "arid_norm.txt", row.names=FALSE, sep='\t')
# arid.tr = arid[arid$value_1>0.1 & arid$value_2>0.1,]
# igg.tr = igg[igg$value_1>0.1 & igg$value_2>0.1,]

# arid.tr = arid.tr[order(arid.tr$log2.fold_change., decreasing=TRUE),]
# igg.tr = igg.tr[order(igg.tr$log2.fold_change., decreasing=TRUE),]

# plot(log2(data.tr$value_1), log2(data.tr$value_2))


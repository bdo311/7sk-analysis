library(edgeR)

get_region_edge = function(libsize, countsn, grp, x) {
	y = DGEList(countsn, group=grp, lib.size=libsize)  # 1 = WT, 2 = ASO

	y <- estimateCommonDisp(y)
	y <- estimateTagwiseDisp(y)
	et <- exactTest(y, pair=c("1", "2"))
	edge <- as.data.frame(topTags(et, n= dim(countsn)[1]))

	# Select FDR thresholded regions and write to file
	edge_padj <- edge[edge$FDR < 0.1, ]
	edge_padj_up <- edge_padj[edge_padj$logFC >= 0,]
	edge_padj_down <- edge_padj[edge_padj$logFC <= -0,]

	regions_up = cbind(counts[as.numeric(rownames(edge_padj_up)),x], edge_padj_up)
	regions_down = cbind(counts[as.numeric(rownames(edge_padj_down)),x], edge_padj_down)
	regions_all_annotation = counts[as.numeric(rownames(edge)),x]
	return(cbind(regions_all_annotation,edge))
}

lib.size_n = c(24156541, 21252957, 29303059, 20847033)  # Baf
baf_prom_prox_counts = read.table("../../AJR/ChIP_counts/Annotate_mES_prom_startseq_centered_pm200bp_baf.bed.txt")
baf_prom_prox_counts = baf_prom_prox_counts[,8:11] # if using Adam's files
grp = c(2,1,2,1)
x = c(2,3,4,1,6,5)  # for Adam's things
baf_prom_prox = get_region_edge(lib.size_n, baf_prom_prox_counts, grp, x)

lib.size_n = c(39648219, 36060025, 41679902, 35611134) # gH2Ax
gh2ax_prom_prox_counts = read.table("gH2AX_prom_prox.bed")
gh2ax_prom_prox_counts = gh2ax_prom_prox_counts[,c(7,8,9,10)] # for all Brian's files except GRO
grp = c(1,1,2,2)
x = 1:6
gh2ax_prom_prox = get_region_edge(lib.size_n, gh2ax_prom_prox_counts, grp, x)

lib.size_n = c(35844103,22886901,46072000,40614974) # GROseq; ctl3-ctl4, aso3-aso4
gro_prom_prox_raw = read.table("GRO_prom_prox.bed")
gro_prom_prox_counts = gro_prom_prox_raw[,c(9,10,13,14)] # for Brian's GRO
grp = c(1,1,2,2)
x = 1:6
gro_prom_prox = get_region_edge(lib.size_n, gro_prom_prox_counts, grp, x)

rownames(gh2ax_prom_prox) = gh2ax_prom_prox$V4
rownames(gro_prom_prox) = gro_prom_prox$V4
rownames(baf_prom_prox) = baf_prom_prox$V4
rownames(gro_prom_prox_raw) = gro_prom_prox_raw[,4]

baf_prom_prox = baf_prom_prox[rownames(gh2ax_prom_prox),]
gro_prom_prox = gro_prom_prox[rownames(gh2ax_prom_prox),]
gro_prom_prox_raw = gro_prom_prox_raw[rownames(gh2ax_prom_prox),]

gh2ax_fc = gh2ax_prom_prox$logFC
gro_fc = gro_prom_prox$logFC
baf_fc = baf_prom_prox$logFC
gro_raw = log2(apply(gro_prom_prox_raw[,9:10],1,sum))
gro_raw[!is.finite(gro_raw)] = 0

ct_promprox = read.delim("/home/raflynn/7SK/GROseq/Flynn/metagenes/regionLists/TR_prox30_TSS_centered_list_CT_promprox.txt", header=TRUE, row.names=1)
x = strsplit(rownames(ct_promprox), '__')
rownames(ct_promprox) = sapply(x, function(y) y[4])
ct_promprox = ct_promprox[rownames(gh2ax_prom_prox),]
ct_fc = log2(ct_promprox$X5ASOcomb_ConvT/ct_promprox$Scrcomb_ConvT)
ct_fc[!is.finite(ct_fc)] = 0

df_promprox = data.frame(gh2ax_fc, gro_fc, baf_fc, gro_raw, ct_fc)

# write to file
all_promprox_bed = read.table("~/7SK/ChIRPseq/genes/mm9_tss_startseq_centered_1kb.bed", stringsAsFactors=FALSE)
df_promprox_reord = df_promprox[all_promprox_bed$V4,]
combined_df_promprox = cbind(all_promprox_bed, df_promprox_reord)
write.table(combined_df_promprox, "all_promprox_data.txt", quote=FALSE, col.names=TRUE, row.names=FALSE, sep='\t')



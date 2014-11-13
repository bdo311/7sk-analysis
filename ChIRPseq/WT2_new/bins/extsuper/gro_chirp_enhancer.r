# gro_chirp_enhancer.r
# 10/9/14
# to make heatmaps of enhancers in gro aso

require(pheatmap)
require(RColorBrewer)
load("data.RData")

getMatrix = function(path) {
	load(paste(path, "data.RData", sep=''))
	return(data.extrap)
}

### ChIRP

# chirp heatmap
chirp_reg = getMatrix("/home/raflynn/7SK/ChIRPseq/WT2_new/bins/extreg/")
chirp_super = getMatrix("/home/raflynn/7SK/ChIRPseq/WT2_new/bins/extsuper/")

regSortOrder = order(apply(chirp_reg, 2, sum), decreasing=TRUE)
superSortOrder = order(apply(chirp_super, 2, sum), decreasing=TRUE)

chirp_reg = t(chirp_reg[,regSortOrder])
chirp_super = t(chirp_super[,superSortOrder])

chirp_reg_sc = log10(chirp_reg + 1)

colorRamp = colorRampPalette(c("white", "blue"))
png("chirp_reg.png", width=1000, height=2000)
pheatmap(chirp_reg_sc, cluster_rows = FALSE, cluster_cols = FALSE, scale="none", show_rownames=FALSE, show_colnames=FALSE,col=colorRamp(50))
dev.off()

### GRO Sense/Antisense

# gro heatmap
gro12c_reg_sense = getMatrix("/arrayAhome/raflynn/7SK/GROseq/GRO_12C_plus/bins/extreg/")
gro12c_reg_as = getMatrix("/arrayAhome/raflynn/7SK/GROseq/GRO_12C_minus/bins/extreg/")
gro12c_reg = gro12c_reg_as + gro12c_reg_sense
gro12c_reg = t(gro12c_reg[,regSortOrder])

gro12c_reg_sc = log10(gro12c_reg +1)

colorRamp = colorRampPalette(c("white", "firebrick4"))
png("gro12c_reg.png", width=1000, height=2000)
pheatmap(gro12c_reg_sc, cluster_rows = FALSE, cluster_cols = FALSE, scale="none", show_rownames=FALSE, show_colnames=FALSE,col=colorRamp(50))
dev.off()

# gro ASO heatmap
gro125_reg_sense = getMatrix("/arrayAhome/raflynn/7SK/GROseq/GRO_125_plus/bins/extreg/")
gro125_reg_as = getMatrix("/arrayAhome/raflynn/7SK/GROseq/GRO_125_minus/bins/extreg/")
gro125_reg = gro125_reg_as + gro125_reg_sense
gro125_reg = t(gro125_reg[,regSortOrder])

gro123_reg_sense = getMatrix("/arrayAhome/raflynn/7SK/GROseq/GRO_123_plus/bins/extreg/")
gro123_reg_as = getMatrix("/arrayAhome/raflynn/7SK/GROseq/GRO_123_minus/bins/extreg/")
gro123_reg = gro123_reg_as + gro123_reg_sense
gro123_reg = t(gro123_reg[,regSortOrder])

gro12a_reg = (gro125_reg + gro123_reg)/2
gro12a_reg_sc = log10(gro12a_reg +1)

colorRamp = colorRampPalette(c("white", "firebrick4"))
png("gro12a_reg.png", width=1000, height=2000)
pheatmap(gro12a_reg_sc, cluster_rows = FALSE, cluster_cols = FALSE, scale="none", show_rownames=FALSE, show_colnames=FALSE,col=colorRamp(50))
dev.off()

# gro changes
gro12_reg_chg = gro12a_reg_sc - gro12c_reg_sc
breaks=c(seq(-2,-0.01,length.out=25),seq(0.01,2,length.out=25))
colorRamp2 = colorRampPalette(c("dodgerblue4","skyblue3","white","yellow1","yellow4"))
png("gro_12_reg_chg.png", width=1000, height=2000)
pheatmap(gro12_reg_chg, cluster_rows = FALSE, cluster_cols = FALSE, scale="none", show_rownames=FALSE, show_colnames=FALSE,col=colorRamp2(50), breaks=breaks)
dev.off()





# SUPER
chirp_super_sc = log10(chirp_super + 1)

colorRamp = colorRampPalette(c("white", "blue"))
png("chirp_super.png", width=1000, height=200)
pheatmap(chirp_super_sc, cluster_rows = FALSE, cluster_cols = FALSE, scale="none", show_rownames=FALSE, show_colnames=FALSE,col=colorRamp(50))
dev.off()

### GRO Sense/Antisense

# gro heatmap
gro12c_super_sense = getMatrix("/arrayAhome/raflynn/7SK/GROseq/GRO_12C_plus/bins/extsuper/")
gro12c_super_as = getMatrix("/arrayAhome/raflynn/7SK/GROseq/GRO_12C_minus/bins/extsuper/")
gro12c_super = gro12c_super_as + gro12c_super_sense
gro12c_super = t(gro12c_super[,superSortOrder])

gro12c_super_sc = log10(gro12c_super +1)

colorRamp = colorRampPalette(c("white", "firebrick4"))
png("gro12c_super.png", width=1000, height=200)
pheatmap(gro12c_super_sc, cluster_rows = FALSE, cluster_cols = FALSE, scale="none", show_rownames=FALSE, show_colnames=FALSE,col=colorRamp(50))
dev.off()

# gro ASO heatmap
gro125_super_sense = getMatrix("/arrayAhome/raflynn/7SK/GROseq/GRO_125_plus/bins/extsuper/")
gro125_super_as = getMatrix("/arrayAhome/raflynn/7SK/GROseq/GRO_125_minus/bins/extsuper/")
gro125_super = gro125_super_as + gro125_super_sense
gro125_super = t(gro125_super[,superSortOrder])

gro123_super_sense = getMatrix("/arrayAhome/raflynn/7SK/GROseq/GRO_123_plus/bins/extsuper/")
gro123_super_as = getMatrix("/arrayAhome/raflynn/7SK/GROseq/GRO_123_minus/bins/extsuper/")
gro123_super = gro123_super_as + gro123_super_sense
gro123_super = t(gro123_super[,superSortOrder])

gro12a_super = (gro125_super + gro123_super)/2
gro12a_super_sc = log10(gro12a_super +1)

colorRamp = colorRampPalette(c("white", "firebrick4"))
png("gro12a_super.png", width=1000, height=200)
pheatmap(gro12a_super_sc, cluster_rows = FALSE, cluster_cols = FALSE, scale="none", show_rownames=FALSE, show_colnames=FALSE,col=colorRamp(50))
dev.off()

# gro changes
gro12_super_chg = gro12a_super_sc - gro12c_super_sc
breaks=c(seq(-2,-0.01,length.out=25),seq(0.01,2,length.out=25))
colorRamp2 = colorRampPalette(c("dodgerblue4","skyblue3","white","yellow1","yellow4"))
png("gro_12_super_chg.png", width=1000, height=200)
pheatmap(gro12_super_chg, cluster_rows = FALSE, cluster_cols = FALSE, scale="none", show_rownames=FALSE, show_colnames=FALSE,col=colorRamp2(50), breaks=breaks)
dev.off()


# compareFazzioToGRO.r
# 3/21/15
#d
gro125_plus = read.delim("/arrayAhome/raflynn/7SK/GROseq/Flynn/metagenes/GRO_125comb_sense/bins/Fazzio_DHS_ext/allchr_sorted.txt",header=FALSE,row.names=4)
gro125_minus = read.delim("/arrayAhome/raflynn/7SK/GROseq/Flynn/metagenes/GRO_125comb_antisense/bins/Fazzio_DHS_ext/allchr_sorted.txt",header=FALSE,row.names=4)
gro12C_plus = read.delim("/arrayAhome/raflynn/7SK/GROseq/Flynn/metagenes/GRO_12Ccomb_sense/bins/Fazzio_DHS_ext/allchr_sorted.txt",header=FALSE,row.names=4)
gro12C_minus = read.delim("/arrayAhome/raflynn/7SK/GROseq/Flynn/metagenes/GRO_12Ccomb_antisense/bins/Fazzio_DHS_ext/allchr_sorted.txt",header=FALSE,row.names=4)

rnaseq = read.delim("GFP and Brg1_DNaseI_averages_matched_log2_sortBrg1_onlysig.txt",header=TRUE)
names(rnaseq) = c("name","gfp","brg1","brg1_kd_chg")
newname = paste("chr",rnaseq$name,sep='')
newname = gsub(":","_",newname)
newname = gsub("-",'_',newname)
rnaseq$name = newname

gro125 = gro125_plus[newname,7] + gro125_minus[newname,7]
gro12C = gro12C_plus[newname,7] + gro12C_minus[newname,7]

rnaseq = cbind(rnaseq, gro125,gro12C,"chg"=log2(gro125/gro12C))

png("brg1_vs_groseq.png",width=1000,height=1000)
plot(rnaseq$brg1_kd_chg,rnaseq$chg,cex=0.7,pch=19)
dev.off()
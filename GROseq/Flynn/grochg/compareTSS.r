# compareTSS.r
# 9/13/14
# compares TSS from GRO and GRO ASO

data = read.delim("tss_compared.txt", header=TRUE, row.names=1, colClasses=c("character", rep("numeric", 12)))
data[data<0.1]=0.1

#[1] "senseGRO_123" "senseGRO_6C"  "senseGRO_125" "senseGRO_12C" "senseGRO_63"
#[6] "senseGRO_65"  "divGRO_123"   "divGRO_6C"    "divGRO_125"   "divGRO_12C"
#[11] "divGRO_63"    "divGRO_65"

#aso.chg.12 = data$senseGRO_12C/((data$senseGRO_123 +data$senseGRO_125)/2)


# genes whose TSS GRO increases or decreases most upon 7SK ASO
aso.chg.12 = apply(data, 1, function(x) ((x[1]+x[3])/2)/x[4])
aso.chg.12 = log2(aso.chg.12[aso.chg.12 !=1])
aso.chg.12 = aso.chg.12[order(aso.chg.12, decreasing=TRUE)]
aso.chg.12 = cbind(aso.chg.12, gro = data[names(aso.chg.12),]$senseGRO_12C)
boxplot(log2(aso.chg.12[1:2000,2]), log2(aso.chg.12[15000:17059,2])) # the genes that increase upon ASO often have low ctrl GRO, and vv

plot(aso.chg.12[,1], log2(aso.chg.12[,2]), cex=0.2, pch=19, xlab="log2 ASO/Ctrl", ylab="log2 Ctrl GRO")

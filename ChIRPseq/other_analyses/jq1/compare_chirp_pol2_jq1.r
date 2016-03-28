


tss_chirp = read.delim("/home/raflynn/7SK/ChIRPseq/metagenes/regionLists/TSS_centered_list.txt", header=TRUE, row.names=1)
tss_chip = read.delim("/home/raflynn/7SK/ChIPseq/Flynn/metagenes/regionLists/TSS_centered_1k_list.txt", header=TRUE, row.names=1)
 
re_chirp = read.delim("/home/raflynn/7SK/ChIRPseq/metagenes/regionLists/RY_enh_centered_list.txt", header=TRUE, row.names=1)
re_chip = read.delim("/home/raflynn/7SK/ChIPseq/Flynn/metagenes/regionLists/RY_enh_centered_1k_list.txt", header=TRUE, row.names=1)
 
se_chirp = read.delim("/home/raflynn/7SK/ChIRPseq/metagenes/regionLists/SE_indiv_list.txt", header=TRUE, row.names=1)
se_chip = read.delim("/home/raflynn/7SK/ChIPseq/Flynn/metagenes/regionLists/SE_indiv_1k_list.txt", header=TRUE, row.names=1)

tss_wt_ratio = 
log2(with(tss_chirp, X7SK_mES_WT/1.42) / 
with(tss_chip, (mES_PolII_DMSO1 + mES_PolII_DMSO2)/(2*1363)))

tss_jq1_ratio = 
log2(with(tss_chirp, X7SK_mES_JQ11/1.12) / 
with(tss_chip, (mES_PolII_JQ11_1 + mES_PolII_JQ11_2)/(2*1383)))

tss = data.frame(wt=tss_wt_ratio, jq1=tss_jq1_ratio)
tss$sig = (abs(tss$wt-tss$jq1)>1)

pdf("plot.pdf", width=5, height=5)
with(tss, plot(wt, jq1, col=factor(sig), pch=19, cex=0.5))
abline(0,1)
dev.off()

re_wt_ratio = 
log2(with(re_chirp, X7SK_mES_WT/1.42) / 
with(re_chip, (mES_PolII_DMSO1 + mES_PolII_DMSO2)/(2*1258)))

re_jq1_ratio = 
log2(with(re_chirp, X7SK_mES_JQ11/1.08) / 
with(re_chip, (mES_PolII_JQ11_1 + mES_PolII_JQ11_2)/(2*1241)))

re = data.frame(wt=re_wt_ratio, jq1=re_jq1_ratio)
re$sig = (abs(re$wt-re$jq1)>1)

pdf("plot.pdf", width=5, height=5)
with(re, plot(wt, jq1, col=factor(sig), pch=19, cex=0.5))
abline(0,1)
dev.off()


se_wt_ratio = 
log2(with(se_chirp, X7SK_mES_WT/2.18) / 
with(se_chip, (mES_PolII_DMSO1 + mES_PolII_DMSO2)/(2*1833)))

se_jq1_ratio = 
log2(with(se_chirp, X7SK_mES_JQ11/1.21) / 
with(se_chip, (mES_PolII_JQ11_1 + mES_PolII_JQ11_2)/(2*1422)))

se = data.frame(wt=se_wt_ratio, jq1=se_jq1_ratio)
se$sig = (abs(se$wt-se$jq1)>1)

pdf("plot.pdf", width=5, height=5)
with(se, plot(wt, jq1, col=factor(sig), pch=19, cex=0.5))
abline(0,1)
dev.off()

pdf("jq1_vs_wt.pdf", width=5, height=9)
par(mfrow=c(2,1))
boxplot(tss$wt, tss$jq1, re$wt, re$jq1, se$wt, se$jq1, xaxt='n', ylab="log2 7SK/Pol2 ratio", ylim=c(-6,6))
axis(1, at=1:6, labels=c("TSS WT", "TSS JQ1", "RE WT", "RE JQ1", "SE WT", "SE JQ1"))

tss.chg = tss$jq1 - tss$wt
tss.chg = tss.chg[is.finite(tss.chg)]
re.chg = re$jq1 - re$wt
re.chg = re.chg[is.finite(re.chg)]
se.chg = se$jq1 - se$wt
se.chg = se.chg[is.finite(se.chg)]
boxplot(tss.chg, re.chg, se.chg, xaxt='n', ylab="log2 JQ1/WT CC ratio", ylim=c(-5,4))
abline(h=0, lty=2)
axis(1, at=1:3, labels=c("TSS", "RE", "SE"))
dev.off()

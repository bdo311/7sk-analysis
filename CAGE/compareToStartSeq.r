cage = read.delim("mm9_tss_cage_centered_1kb.bed",header=FALSE)
startseq = read.delim("/arrayAhome/raflynn/7SK/Start-RNAseq/mm9_tss_startseq_centered_1kb.bed",header=FALSE)

# Sense TSS peak
cage$ctr = cage[,2]+1000
startseq$ctr = startseq[,2]+1000
rownames(cage) = as.character(cage[,4])
rownames(startseq) = as.character(startseq[,4])

genenames = cage[,4][which(cage[,4] %in% startseq[,4])]
cage.ctr = cage[genenames,]$ctr
startseq.ctr = startseq[genenames,]$ctr
hist(abs(cage.ctr-startseq.ctr),breaks=seq(from=0,to=2000,by=5),xlim=c(0,100))

# Antisense TSS peak
cage=cage[cage[,5]>0,]
startseq=startseq[startseq[,5]>0,]
genenames = cage[,4][which(cage[,4] %in% startseq[,4])]
cage.div = cage[genenames,5]
startseq.div = startseq[genenames,5]
hist(abs(cage.div-startseq.div),breaks=seq(from=0,to=2000,by=5),xlim=c(0,300))


# h1_hockeyplot.r

re=read.delim("/arrayAhome/raflynn/7SK/ChIRPseq/metagenes/7SK_H1/bins/H1_regenh_all/allchr_sorted.txt",header=FALSE)
se=read.delim("/arrayAhome/raflynn/7SK/ChIRPseq/metagenes/7SK_H1/bins/H1_superenh_all/allchr_sorted.txt",header=FALSE)

re_count=re[,7]*(re[,3]-re[,2])
se_count=se[,7]*(se[,3]-se[,2])
nameList = c(as.character(re[,4]),as.character(se[,4]))
re_df = data.frame(re[,1],re[,2],re[3],re[,4],re_count)
se_df = data.frame(se[,1],se[,2],se[3],se[,4],se_count)

counts = as.data.frame(t(cbind(t(re_df),t(se_df))))
colnames(counts)=c("chr","start","stop","name","total_count")
counts$total_count=as.numeric(as.character(counts$total_count))
counts$type = c(rep(1,length(re_count)),rep(2,length(se_count)))
counts=counts[order(counts$total_count),]
write.table(counts, "h1_hockeyplot_list.txt",row.names=TRUE,sep='\t')

pdf("h1_hockeyplot.pdf",width=5,height=6)
plot(1:nrow(counts),counts$total_count,pch=19,col=counts$type,xlab="Enhancers ranked by read count",ylab="7SK read count",main="7SK ChIRP read count per enhancer")
legend("topleft",legend=c("Regular enhancer","Super enhancer"),col=1:2,pch=19)
dev.off()
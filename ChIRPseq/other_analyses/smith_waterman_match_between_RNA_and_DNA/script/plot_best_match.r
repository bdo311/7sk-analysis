args <- commandArgs(TRUE)
inputfolder=args[1]
labels = c('RNA-DNA','RNA-DNA*','RNA*-DNA','RNA*-DNA*')
hit_th = 0 # only count hits with score at least this
nbin = 100

# combining two strands
for (rnafile in c('RNA','RNA.shuffled')){
    for (dnafile in c('DNA','DNA.shuffled')){
	pair = paste(rnafile,dnafile,sep='-')
        print(pair)
	minus = read.table(paste(inputfolder,'/',pair,'.alignment.bed',sep=''),header=F,stringsAsFactors=F)
	plus = read.table(paste(inputfolder,'/',pair,'.alignment.rc.bed',sep=''),header=F,stringsAsFactors=F)
	minus = data.frame(minus,stringsAsFactors=F)
	plus = data.frame(plus,stringsAsFactors=F)
	best = minus[,5] > plus[,5]
	plus[best,] = minus[best,]
	write.table(plus,file=paste(inputfolder,'/',pair,'.best',sep=''),quote=F,sep='\t',row.names=F,col.names=F)
    }
}


boxplot
pdf(paste(inputfolder,'/plot.pdf',sep=''))
par(cex=1.5)
scores = numeric(0)
for (rnafile in c('RNA','RNA.shuffled')){
    for (dnafile in c('DNA','DNA.shuffled')){
	x = read.table(paste(inputfolder,'/',rnafile,'-',dnafile,'.best',sep=''),header=F,stringsAsFactors=F)
	scores = cbind(scores,x[,5])
    }
}
scores[scores>20]=20
par(mar=c(6,4,4,2))
boxplot(scores,names=labels,las=2,ylab="alignment score")
dev.off()

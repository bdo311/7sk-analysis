library(boot)

# args <- commandArgs(TRUE)
# inputfolder=args[1]
# labels = c('RNA-DNA','RNA-DNA*','RNA*-DNA','RNA*-DNA*')
# hit_th = 0 # only count hits with score at least this
# nbin = 100


# # combining two strands
# for (rnafile in c('RNA')){
    # for (dnafile in c('DNA')){
	# pair = paste(rnafile,dnafile,sep='-')
        # print(pair)
	# minus = read.table(paste(inputfolder,'/',pair,'.alignment_new.bed',sep=''),header=F,stringsAsFactors=F)
	# plus = read.table(paste(inputfolder,'/',pair,'.alignment.rc_new.bed',sep=''),header=F,stringsAsFactors=F)
	# minus = data.frame(minus,stringsAsFactors=F)
	# plus = data.frame(plus,stringsAsFactors=F)
	# best = minus[,5] > plus[,5]
	# plus[best,] = minus[best,]
	# write.table(plus,file=paste(inputfolder,'/',pair,'.best',sep=''),quote=F,sep='\t',row.names=F,col.names=F)
    # }
# }

data = read.table("RNA-DNA.best", sep='\t')
sw = data[,5]
mean(sw)

# getting a bootstrap)
getMean = function(data, indices) {
	d = data[indices]
	mean(d)
}
	
bootobj = boot(data=sw, statistic=getMean, R=1000)
bootobj
plot(bootobj)
boot.ci(bootobj, conf=0.95)

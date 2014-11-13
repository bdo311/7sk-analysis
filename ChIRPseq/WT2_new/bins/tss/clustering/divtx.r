# divtx.r
# 5/24/14
# to explore the contribution of divtx/bidir to the upstream peak

require(amap)
load("data.RData")

data.s = data.extrap[,order(apply(data.extrap, 2, sum), decreasing=TRUE)]
# data.s = data.extrap[100:300,order(apply(data.extrap, 2, sum), decreasing=TRUE)]
# data.s=data.s[,1:5000]
data.ss = data.s[, apply(data.s, 2, max)!=0]
data.sss = apply(data.ss, 2, function(x) x/max(x))

bigenes = read.delim("/home/raflynn/ChIRPseq/genes/mm9_bidir_genes.txt", sep='\t')
bigenes.200 = bigenes[bigenes[,6]-bigenes[,3]<200,]
bidir = c(as.character(bigenes.200[,1]), as.character(bigenes.200[,4]))

annot = read.delim("/home/raflynn/ChIRPseq/genes/annotated_genes.txt", sep='\t', header=TRUE)
div = as.character(annot[annot$udRNA.lfc!=0,]$gene)
ft = as.character(annot[annot$termination.defect!=0,]$gene)

`%ni%` = Negate(`%in%`)
data.nodb = data.ss[,which(colnames(data.ss) %ni% c(bidir, div))] #no bidir, div
data.nodbs = apply(data.nodb, 2, function(x) x/max(x)) #scaled
collapseTo = 50
data.nodbs.coll = apply(data.nodbs, 2, function(x) approx(1:length(x), x, n=collapseTo)$y) #collapsed to 50 bins

# # plot decreasing upstream peak
# plot(type='l', apply(data.ss, 1, mean), lwd=2, ylim=c(1,3.5), xlab="Position", ylab="Density", main="7SK at TSS, -500 to +500")
# points(type='l', apply(data.ss[,which(colnames(data.ss) %ni% bidir)], 1, mean), lwd=2, col="blue")
# points(type='l', apply(data.ss[,which(colnames(data.ss) %ni% div)], 1, mean), lwd=2, col="red")
# points(type='l', apply(data.ss[,which(colnames(data.ss) %ni% c(bidir, div))], 1, mean), lwd=2, col="green")
# legend("topright", legend=c("all", "no bidir", "no divtx", "no bidir/divtx"), col=c("black", "blue", "red", "green"), lty=rep(1,4), lwd=rep(2,4))

# get input line
input = read.delim("/home/raflynn/ChIRPseq/7SK_Input/bins/tss/avgraw_7SK_Input_tss.txt", header=TRUE)


# plot three things on one plot
pdf("all_div_ft.pdf", width=7, height=7)
plot(type='l', apply(data.ss, 1, mean), lwd=2, main="TSS at genes", ylab="7SK read density/nt", xlab="Distance from TSS", xaxt='n', ylim=c(0,3))
lines(apply(data.ss[,which(colnames(data.ss)%in% div)],1,mean), lwd=2, col="green")
lines(apply(data.ss[,which(colnames(data.ss)%in% ft)],1,mean), lwd=2, col="red")
lines(input, lwd=2, col="blue")
axis(1, at=c(0,100,200,300,400),labels=c(-1000,-500,0,"+500","+1000"))
dev.off()


dev.off()

# plot individual classes
pdf("all.pdf", width=7, height=7)
par(mfrow=c(2,2))
n = ncol(data.ss)
plot(type='l', apply(data.ss, 1, mean), lwd=2, col="blue", main=paste("TSS at top 10K 7SK-bound genes\nn = ", n, sep=''), ylab="7SK read density/nt", xlab="Distance from TSS", xaxt='n', ylim=c(0,4))
axis(1, at=c(0,100,200),labels=c(-500,0,"+500"))
#pdf("div.pdf", width=4, height=4)
n=ncol(data.ss[,which(colnames(data.ss) %in% div)])
plot(type='l', apply(data.ss[,which(colnames(data.ss) %in% div)], 1, mean), lwd=2, col="blue", xlab="Distance from TSS", xaxt='n', ylim=c(0,4), ylab="7SK read density/nt", main=paste("TSS at divergent genes,\nn = ", n, sep=''))
axis(1, at=c(0,100,200),labels=c(-500,0,"+500"))
#pdf("bidir.pdf", width=4, height=4)
n=ncol(data.ss[,which(colnames(data.ss) %in% bidir)])
plot(type='l', apply(data.ss[,which(colnames(data.ss) %in% bidir)], 1, mean), lwd=2, col="blue", xlab="Distance from TSS", xaxt='n', ylim=c(0,4), ylab="7SK read density/nt", main=paste("TSS at bidirectional genes,\nn = ", n, sep=''))
axis(1, at=c(0,100,200),labels=c(-500,0,"+500"))
#pdf("rest.pdf", width=4, height=4)
n=ncol(data.nodb)
plot(type='l', apply(data.nodb, 1, mean), lwd=2, col="blue", ylim=c(0,4), ylab="7SK read density/nt", xaxt='n', xlab="Distance from TSS", main=paste("TSS at non-divergent, \nnon-bidirectional genes,\nn = ", n, sep=''))
axis(1, at=c(0,100,200),labels=c(-500,0,"+500"))
dev.off()

pdf("all.pdf", width=5, height=5)
plot(type='l', apply(data.nodb, 1, mean), lwd=2, col="blue", ylim=c(0,3), ylab="7SK read density/nt", xaxt='n', xlab="Distance from TSS", main="TSS at genes")
lines(apply(data.ss[,which(colnames(data.ss) %in% div)], 1, mean), lwd=2, col="green")
lines(apply(data.ss[,which(colnames(data.ss) %in% bidir)], 1, mean), lwd=2, col="red")
axis(1, at=c(0, 200, 400),labels=c(-1000,0,"+1000"))
legend("topleft", legend=c("Divergent txn", "Bidirectional", "All others"), lwd=2, lty=1, col=c("green", "red", "blue"))
dev.off()

# k-means
d = Dist(t(data.nodbs.coll), method="pearson", nbproc=10)


km = kmeans(t(data.nodbs), centers=3, nstart=3)

par(mfrow=c(1,3))
plot(type='l', apply(data.nodb[,km$cluster==1], 1, mean))
plot(type='l', apply(data.nodb[,km$cluster==2], 1, mean))
plot(type='l', apply(data.nodb[,km$cluster==3], 1, mean))

# graphic of all the clusters.
pdf("kmeans.pdf", width=12, height=10)
par(mfrow=c(5,6),mar=c(2,2,2,2))
kmeans.objects = list()
#kmeans.objects[[1]] = kmeans(t(data.nodbs.coll), centers=1, nstart=3)
kmeans.objects[[1]] = Kmeans(t(data.nodbs.coll), centers=1, nstart=5, iter.max=50, method="pearson") #

for(i in 2:6) {
	# km = kmeans(t(data.nodbs.coll), centers=i, nstart=3)
	km=Kmeans(t(data.nodbs.coll), centers=i, nstart=5, iter.max=50, method="pearson")
	kmeans.objects[[i]] = km
	
	# for(j in 1:i) {
		# plot(type='l', apply(data.nodb[,km$cluster==j], 1, mean), ylab='n', xlab='n',xaxt='n',ylim=c(0,6),lwd=2)
		# axis(1, at=c(0,100,200), labels=c("-500", 0, "+500"))
		# text(x=100,y=0.5, labels=as.character(sum(km$cluster==j)), col="red", cex=2)
	# }
	# if (i < 6) {
		# for(k in (i+1):6) {	frame()}	
	# }
}
dev.off()

# writing gene names
threemeans = kmeans.objects[[3]]
allgenes = colnames(data.nodb)
center = colnames(data.nodb[,threemeans$cluster==1])
left = colnames(data.nodb[,threemeans$cluster==2])
right = colnames(data.nodb[,threemeans$cluster==3])

write.table(allgenes, "allgenes.txt", row.names=F,col.names=F,sep="\t")
write.table(center, "center.txt", row.names=F,col.names=F,sep="\t")
write.table(left, "left.txt", row.names=F,col.names=F,sep="\t")
write.table(right, "right.txt", row.names=F,col.names=F,sep="\t")

# gap statistic
require(parallel)
cl=makeCluster(10)
unif.Wk = matrix(0, nrow=6, ncol=100)
for(i in 1:ncol(unif.Wk)) {
	print(i)
	set.seed(i)
	rand = t(parSapply(cl, 1:ncol(data.nodbs), function(x) runif(10, 0, 1)))
	# rand = matrix(runif(collapseTo*ncol(data.nodbs), 0, 1), nrow=ncol(data.nodbs), ncol=collapseTo)
	for(k in 1:6) {
		km = kmeans(rand, centers=k, nstart=3)
		unif.Wk[k,i] = km$tot.withinss
	}
}

e.log.wk = apply(log(unif.Wk), 1, mean)
sd.log.wk = apply(log(unif.Wk), 1, sd)
log.obs.wk = sapply(1:6, function(x) log(kmeans.objects[[x]]$tot.withinss))

gapStat = e.log.wk - log.obs.wk
upSD = gapStat + sd.log.wk
downSD = gapStat - sd.log.wk

# plot gap statistic with error bars
pdf("gap.pdf", width=4, height=4)
par(mfrow=c(1,1))
plot(gapStat, type='b', pch=19, xlab="Number of clusters", ylab="Gap statistic", main="Gap statistic as a function of k")
arrows(1:6, upSD, 1:6, downSD,code=3,length=0.05,angle=90,col='red')
dev.off()




# another metric: PDN method
# obs.wk = sapply(1:6, function(x) kmeans.objects[[x]]$tot.withinss)
obs.wk = sapply(1:6, function(x) sum(kmeans.objects[[x]]$withinss))

n_d = 50
getAK = function(largestK) {
	if(largestK == 2) { 
		return(1 - 3/(4*n_d))
	}
	a_kminus1 = getAK(largestK - 1)
	return(a_kminus1 + (1 - a_kminus1)/6)
}

getFK = function(largestK) {
	if(largestK == 1) { 
		return(1)
	}
	fk = obs.wk[largestK]/(a_k[largestK]*obs.wk[largestK - 1])
	return(fk)
}

a_k = c(0, sapply(2:6, getAK))
f_k = sapply(1:6, getFK)

pdf("pdn.pdf", width=4, height=4)
plot(f_k, type='b',pch=19, main="PDN statistic as \nfunction of cluster number", xlab="Cluster number", ylab="PDN statistic")
dev.off()









# CH
chList = rep(0,6)
for (i in 2:6) {
  between = kmeans.objects[[i]]$betweenss
  within = kmeans.objects[[i]]$tot.withinss
  ch = (between/(i - 1))/(within/(7291 - i)) #CH = (b/(k-1))/(w/(n-k)) and n = 20
  chList[i] = ch
}

# plot of CH statistic
par(mfrow=c(1,1))
plot(2:7, chList[2:7], type='b', pch=19, xlab='Number of clusters', ylab = "CH statistic", main="CH statistic as a function of k")

# gap statistic
require(cluster)
clusObj = clusGap(t(data.nodbs.coll), FUN=kmeans, nstart=3, K.max = 7)

# get gap statistic and s[k]
gapStat = clusObj$Tab[,3]
upSD = gapStat + clusObj$Tab[,4]
downSD = gapStat - clusObj$Tab[,4]

# plot gap statistic with error bars
plot(gapStat, type='b', pch=19, xlab="Number of clusters", ylab="Gap statistic", main="Gap statistic as a function of k", ylim=c(-0.3, 0.5))
arrows(1:7, upSD, 1:7, downSD,code=3,length=0.05,angle=90,col='red')

# hierarchical clustering
d = dist(t(data.nodbs))
attributes(d)$Size
x=hclust(d, method="single", )
par(mfrow=c(1,1))
plot(x, labels=FALSE) #kinda crappy

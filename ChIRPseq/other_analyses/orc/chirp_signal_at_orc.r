pm2 = read.delim("/home/raflynn/7SK/ChIRPseq/metagenes/7SK_mES_WT/bins/Mouse_ORC_pm2kb/allchr_sorted.txt", colClasses=c(rep("character",4),"numeric","character",rep("numeric",401)), header=FALSE)
pm2 = pm2[order(pm2[,5], decreasing=TRUE),]

pdf("orc_dist_to_7sk_signal.pdf", width=6, height=6)
plot(log10(pm2[,5]+1), log10(apply(pm2[,108:307],1,sum)+1), pch=19, cex=0.1,
xlab="log10(Dist + 1)", ylab="log10(7SK signal + 1)")
dev.off()

pdf("orc_7sk_metagene.pdf", width=6, height=6)
plot(apply(pm2[1:10000,8:407],2,mean),type='l', col="red", ylim=c(0,1))
lines(apply(pm2[10000:20000,8:407],2,mean),type='l', col="orange")
lines(apply(pm2[20000:30000,8:407],2,mean),type='l', col="gold")
lines(apply(pm2[30000:40000,8:407],2,mean),type='l', col="green")
lines(apply(pm2[40000:50000,8:407],2,mean),type='l', col="cyan")
lines(apply(pm2[50000:60000,8:407],2,mean),type='l', col="blue")
lines(apply(pm2[60000:65000,8:407],2,mean),type='l', col="purple")
title("Red = farthest, blue = closest")
dev.off()
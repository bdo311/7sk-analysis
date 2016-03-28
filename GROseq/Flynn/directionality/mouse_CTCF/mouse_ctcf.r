## set up data

chirp_7sk = read.delim("/home/raflynn/7SK/GROseq/Flynn/directionality/mouse_CTCF/metagenes/7SK/bins/CTCF_all/allchr_sorted.txt", header=FALSE, colClasses=c(rep("character",6),rep("numeric",401)))
rownames(chirp_7sk) = apply(chirp_7sk, 1, function(x) paste(x[1], x[2], x[3], sep='_'))
chirp_input = read.delim("/home/raflynn/7SK/GROseq/Flynn/directionality/mouse_CTCF/metagenes/7SK_Input/bins/CTCF_all/allchr_sorted.txt", header=FALSE, colClasses=c(rep("character",6),rep("numeric",401)))
rownames(chirp_input) = apply(chirp_input, 1, function(x) paste(x[1], x[2], x[3], sep='_'))
chirp_input = chirp_input[rownames(chirp_7sk),]


hela_7sk = read.delim("/home/raflynn/Mitosis/NETseq/metagenes/7SK_HeLa/bins/Churchman_CTCF_all/allchr_sorted.txt", header=FALSE, colClasses=c(rep("character",6),rep("numeric",401)))
rownames(hela_7sk) = apply(hela_7sk, 1, function(x) paste(x[1], x[2], x[3], sep='_'))
hela_input = read.delim("/home/raflynn/Mitosis/NETseq/metagenes/7SK_HeLa_Input/bins/Churchman_CTCF_all/allchr_sorted.txt", header=FALSE, colClasses=c(rep("character",6),rep("numeric",401)))
rownames(hela_input) = apply(hela_input, 1, function(x) paste(x[1], x[2], x[3], sep='_'))
hela_input = hela_input[rownames(hela_7sk),]



## Inputs 
plot(apply(hela_input[1:1000,8:407],2,mean),type='l',ylim=c(0.35,0.6))
lines(apply(chirp_input[1:1000,8:407],2,mean),type='l', col="red")
abline(v=200, col="lightblue")
legend("topright", legend=c("mES","HeLa"), col=c("red","black"),lty=1,lwd=2)
title("ChIRP inputs")

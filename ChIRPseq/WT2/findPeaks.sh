perl /seq/chirpseq/bedGraph2sam.pl ../genes/mm9.sizes combined.bedGraph combined.sam
macs14 -t combined.sam -c input.sam -f SAM -n combined --bw 200 -m 10,50

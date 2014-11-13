#1. bowtie
#parallel -j 3 "bowtie2 -p 4 -x /home/raflynn/Scripts/repeat_index/mm9/rep -k 1 -U {} --un {.}_unm.fastq -S {.}.sam" ::: *2_trimmed.fastq
#parallel -j 3 "bowtie2 -p 4 -x /seq/bowtie2-2.1.0/indexes/mm9 -k 1 -U {} -S {.}.sam" ::: *2_trimmed_unm.fastq

#2. samtools to sort the reads mapped to genome
# parallel -j 6 "samtools view -Suo - {} | samtools sort - {.}_sorted.bam" ::: *_trimmed_unm.sam
#parallel -j 6 "samtools view -Suo - {} | samtools sort - {.}_sorted.bam" ::: *trimmed.sam
#parallel -j 8 "samtools view -F 0x10 -Suo - {} | samtools sort - {.}_positive_sorted" ::: *trimmed_unm.sam
#parallel -j 8 "samtools view -f 0x10 -Suo - {} | samtools sort - {.}_negative_sorted" ::: *trimmed_unm.sam

#3. bedgraph making
# parallel -j 8 "bedtools genomecov -ibam {} -bg > {.}.bedGraph; norm_bedGraph.pl {.}.bedGraph {.}_norm.bedGraph" ::: *2*sorted.bam
parallel -j 8 "bedGraphToBigWig {} /seq/chromosome/mm9/mm9.sizes {.}.bw" ::: *2*unm*_norm.bedGraph
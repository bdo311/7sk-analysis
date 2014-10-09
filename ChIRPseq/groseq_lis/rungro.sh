# bowtie2 -p 4 -x /seq/bowtie2-2.1.0/indexes/mm9 -k 1 -U SRR935093.fastq -S SRR935093.sam
# samtools view -uS SRR935093.sam | samtools sort - SRR935093sorted
# samtools view -F 0x10 -bo SRR935093_positive.bam SRR935093sorted.bam
# samtools view -f 0x10 -bo SRR935093_negative.bam SRR935093sorted.bam
parallel "bedtools genomecov -ibam {} -bg > {.}.bedGraph" ::: SRR935093_*.bam
parallel "norm_bedGraph.pl {} {.}_norm.bedGraph" ::: SRR935093*.bedGraph


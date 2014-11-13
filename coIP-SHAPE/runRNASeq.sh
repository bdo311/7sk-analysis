fn=$1
x=${fn%.*}
basename=${x##*/}

index_7sk="~/Scripts/repeat_index/mm9_7sk/mm_7sk"
repeat_index="~/Scripts/repeat_index/mm9/rep"
genome_index="/seq/bowtie2-2.1.0/indexes/mm9"
sizes="/seq/chromosome/mm9/mm9.sizes"

#tophat -p 4 -o ${basename}_allreads mm9 ${basename}.fastq

echo $fn
#bowtie2 -p 4 -x $index_7sk -k 1 -U $fn -S ${basename}_7sk.sam

# #1. bowtie
# bowtie2 -p 4 -x $repeat_index -k 1 -U $fn --un ${basename}_genome.fastq -S ${basename}_repeat.sam
# tophat -p 4 -o $basename mm9 ${basename}_genome.fastq

# #6. samtools to make bedgraph for repeat index
samtools view -Suo - ${basename}.sam | samtools sort - ${basename}_sorted
samtools index ${basename}_sorted.bam
samtools flagstat ${basename}_sorted.bam

# bedtools genomecov -ibam ${basename}_repeat_sorted.bam -bg > ${basename}_repeat_sorted.bedGraph
# norm_bedGraph.pl ${basename}_repeat_sorted.bedGraph ${basename}_repeat_sorted_norm.bedGraph

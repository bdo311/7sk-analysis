#!/bin/sh
# runchirp.sh

# parameters
if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <organism> <name> <fastq>"
  exit 1
fi

org=$1
name=$2
fastq=$3

if [ $org == "mouse" ]; then
	repeat_pos="~/Scripts/repeat_index/mm9/Mm_repeatIndex_spaced_positions.txt"
	repeat_index="~/Scripts/repeat_index/mm9/rep_spaced"
	genome_index="/seq/bowtie2-2.1.0/indexes/mm9"
	sizes="/seq/chromosome/mm9/mm9.sizes"
elif [ $org == "human" ]; then
	repeat_pos="~/Scripts/repeat_index/hg19/Hs_repeatIndex_spaced_positions.txt"
	repeat_index="~/Scripts/repeat_index/hg19/rep_spaced"
	genome_index="/seq/bowtie2-2.1.0/indexes/hg19"
	sizes="/seq/chromosome/hg19/hg19.sizes"
else
	echo "Organism must be human or mouse."
	exit 1
fi

#1. bowtie
# parallel "bowtie2 -p 4 -x $repeat_index -q {} --un {.}_genome.fastq -S {.}_repeat.sam" ::: $fastq
# parallel "bowtie2 -p 4 -x $genome_index -q {.}_genome.fastq -S {.}_genome.sam" ::: $fastq

# #2. samtools - sorting and removing duplicates
# parallel "cat {.}_genome.sam | samtools view -Suo - - | samtools sort - {.}_genome_sorted" ::: $fastq
# parallel "samtools rmdup -s {.}_genome_sorted.bam {.}_genome_sorted_rmdup.bam" ::: $fastq

#3. make bed files and get 1nt
parallel "bedtools bamtobed -i {.}_genome_sorted.bam | awk -F '\t' 'BEGIN {OFS=\"\t\"} {if (\$6==\"+\") print \$1,\$2,\$2+1,\$4,\$5,\$6;}' > {.}_genome_sorted_1nt_plus.bed" ::: $fastq
parallel "bedtools bamtobed -i {.}_genome_sorted.bam | awk -F '\t' 'BEGIN {OFS=\"\t\"} {if (\$6==\"-\") print \$1,\$3-1,\$3,\$4,\$5,\$6;}' > {.}_genome_sorted_1nt_minus.bed" ::: $fastq
parallel "bedtools bamtobed -i {.}_genome_sorted_rmdup.bam | awk -F '\t' 'BEGIN {OFS=\"\t\"} {if (\$6==\"+\") print \$1,\$2,\$2+1,\$4,\$5,\$6;}' > {.}_genome_sorted_rmdup_1nt_plus.bed" ::: $fastq
parallel "bedtools bamtobed -i {.}_genome_sorted_rmdup.bam | awk -F '\t' 'BEGIN {OFS=\"\t\"} {if (\$6==\"-\") print \$1,\$3-1,\$3,\$4,\$5,\$6;}' > {.}_genome_sorted_rmdup_1nt_minus.bed" ::: $fastq

#4. make bedgraphs for everything
parallel "bedtools genomecov -i {.}_genome_sorted_1nt_plus.bed -g $sizes -bg > {.}_genome_plus.bedGraph; norm_bedGraph.pl {.}_genome_plus.bedGraph {.}_genome_plus_norm.bedGraph" ::: $fastq
parallel "bedtools genomecov -i {.}_genome_sorted_1nt_minus.bed -g $sizes -bg > {.}_genome_minus.bedGraph; norm_bedGraph.pl {.}_genome_minus.bedGraph {.}_genome_minus_norm.bedGraph" ::: $fastq
parallel "bedtools genomecov -i {.}_genome_sorted_rmdup_1nt_plus.bed -g $sizes -bg > {.}_genome_rmdup_plus.bedGraph; norm_bedGraph.pl {.}_genome_rmdup_plus.bedGraph {.}_genome_rmdup_plus_norm.bedGraph" ::: $fastq
parallel "bedtools genomecov -i {.}_genome_sorted_rmdup_1nt_minus.bed -g $sizes -bg > {.}_genome_rmdup_minus.bedGraph; norm_bedGraph.pl {.}_genome_rmdup_minus.bedGraph {.}_genome_rmdup_minus_norm.bedGraph" ::: $fastq

#6. make bedgraphs for repeat index; not merging repeat bedgraph
#parallel "bedtools genomecov -ibam {.}_genome_sorted_rmdup.bam -bg > {.}_genome.bedGraph; norm_bedGraph.pl {.}_genome.bedGraph {.}_genome_norm.bedGraph" ::: $fastq
# parallel "bedtools genomecov -ibam {.}_genome_sorted.bam -bg > {.}_genome.bedGraph; norm_bedGraph.pl {.}_genome.bedGraph {.}_genome_norm.bedGraph" ::: $fastq

#7. get stats for everything
parallel "samtools flagstat {.}_genome_sorted.bam > {.}_genome_sorted_stats.txt" ::: $fastq
parallel "samtools flagstat {.}_genome_sorted_rmdup.bam > {.}_genome_sorted_rmdup_stats.txt" ::: $fastq


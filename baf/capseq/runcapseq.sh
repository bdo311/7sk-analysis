#!/bin/sh
# runchirp.sh

# parameters
if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <organism> <name> <fasta>"
  exit 1
fi

org=$1
name=$2
fasta=$3

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
parallel "bowtie2 -p 4 -x $repeat_index -k 1 -f {} --un {.}_genome.fastq -S {.}_repeat.sam" ::: $fasta
parallel "bowtie2 -p 4 -x $genome_index -k 1 -f {.}_genome.fastq -S {.}_genome.sam" ::: $fasta

#2. samtools - sorting and removing duplicates
parallel "cat {.}_genome.sam | samtools view -Suo - - | samtools sort - {.}_genome_sorted" ::: $fasta
parallel "samtools rmdup -s {.}_genome_sorted.bam {.}_genome_sorted_rmdup.bam" ::: $fasta

#6. make bedgraphs for repeat index; not merging repeat bedgraph
#parallel "bedtools genomecov -ibam {.}_genome_sorted_rmdup.bam -bg > {.}_genome.bedGraph; norm_bedGraph.pl {.}_genome.bedGraph {.}_genome_norm.bedGraph" ::: $fasta
parallel "bedtools genomecov -ibam {.}_genome_sorted.bam -bg > {.}_genome.bedGraph; norm_bedGraph.pl {.}_genome.bedGraph {.}_genome_norm.bedGraph" ::: $fasta

#7. get stats for everything
parallel "samtools flagstat {.}_genome_sorted.bam > {.}_genome_sorted_stats.txt" ::: $even $odd
parallel "samtools flagstat {.}_genome_sorted_rmdup.bam > {.}_genome_sorted_rmdup_stats.txt" ::: $even $odd


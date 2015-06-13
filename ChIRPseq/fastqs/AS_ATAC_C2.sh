# for mapping ATAC-seq data against the Cast/129 SNP-masked genome
# after mapping to this custom genome, use Split_snp.pl to divide into 2 genomes
##parameters#######################################
file1="/home/raflynn/7SK/ChIRPseq/fastq/C2_7SK_Even_R1.fastq.gz"
file2="/home/raflynn/7SK/ChIRPseq/fastq/C2_7SK_Even_R2.fastq.gz"
output="C2_Even"
thread=4 # you can use more if only one sample
ref="/home/jinxu/DB/mouse/mmu9/129S1_CASTEiJ_mm9/129S1_CASTEiJ_mm9"
ref_size="/home/jinxu/DB/mouse/mmu9/129S1_CASTEiJ_mm9/129S1_CASTEiJ_mm9.chrsize"
##Mapping step###########################################
echo $file1
echo $file2
#echo $file1.test
echo "trimming adaptor"
time=`date`
echo $time
# the path of cutadapt should be change into your own path
#/home/raflynn/Scripts/cutadapt-1.7.1/bin/cutadapt -b CTGTCTCTTATACACATCTCCGAGCCCACGAGA  -m 25   -o $file1.tmp.fastq  -p $file2.tmp.fastq  $file1 $file2
#/home/raflynn/Scripts/cutadapt-1.7.1/bin/cutadapt  -b CTGTCTCTTATACACATCTGACGCTGCCGACGA  -m 25   -o $file2.trimmed.fastq.gz  -p $file1.trimmed.fastq.gz   $file2.tmp.fastq $file1.tmp.fastq
#rm $file1.tmp.fastq
#rm $file2.tmp.fastq
echo "trimming adaptor finished "
time=`date`
echo $time
echo "Mapping by bowtie2"
time=`date`
echo $time
#echo $file1.trimmed.fastq.gz
#echo $file2.trimmed.fastq.gz
#bowtie2   -p $thread   --very-sensitive   -x $ref -1 $file1.trimmed.fastq.gz -2 $file2.trimmed.fastq.gz  -S  $output.sam 
# edit a little bit to avoid ATACseq specific analysis 
#awk '$3!="chrM"' $output.sam |samtools view -S -b -f 0x2 -q 10 - |samtools sort -  $output.pe.q10.sort
#samtools rmdup -S $output.pe.q10.sort.bam $output.pe.q10.rmdup.bam 
#awk '$3=="chrM" || NF<10 ' $output.sam |samtools view -S -b -  > $output.chrM.bam
##rm $output.sam
######new command line ######
#samtools view -S -b -f 0x2 -q 10 $output.sam |samtools sort -  $output.pe.q10.sort
#samtools rmdup -S $output.pe.q10.sort.bam $output.pe.q10.rmdup.bam
# this script is ending for mapping step. split SNP step is not including. 
#echo "Mapping by bowtie2 finished"

perl shift_sam_bases_edit.pl /home/raflynn/7SK/ChIRPseq/genes/mm9.sizes ATACseq_12C1.trim.sort.nuc.uniq.rmdup.bam.sam ATACseq_12C1.sam
cat header.sam ATACseq_12C1.sam > tmp.12C1.shifted.sam
mv tmp.12C1.shifted.sam ATACseq_12C1.sam
samtools view -b -S -o ATACseq_12C1.bam ATACseq_12C1.sam


#perl shift_sam_bases_edit.pl /home/raflynn/7SK/ChIRPseq/genes/mm9.sizes ATACseq_1252.trim.sort.nuc.uniq.rmdup.bam.sam ATACseq_1252.sam
cat header.sam ATACseq_1252.sam > tmp.1252.shifted.sam
mv tmp.1252.shifted.sam ATACseq_1252.sam
samtools view -b -S -o ATACseq_1252.bam ATACseq_1252.sam


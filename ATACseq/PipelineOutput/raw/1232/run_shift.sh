#perl shift_sam_bases_edit.pl /home/raflynn/7SK/ChIRPseq/genes/mm9.sizes ATACseq_1232.trim.sort.nuc.uniq.rmdup.bam.sam ATACseq_1232.sam
cat header.sam ATACseq_1232.sam > tmp.1232.shifted.sam
mv tmp.1232.shifted.sam ATACseq_1232.sam
samtools view -b -S -o ATACseq_1232.bam ATACseq_1232.sam


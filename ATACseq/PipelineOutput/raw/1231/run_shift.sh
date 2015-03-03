#perl shift_sam_bases_edit.pl /home/raflynn/7SK/ChIRPseq/genes/mm9.sizes ATACseq_1231.trim.sort.nuc.uniq.rmdup.bam.sam ATACseq_1231.sam
cat header.sam ATACseq_1231.sam > tmp.1231.shifted.sam
mv tmp.1231.shifted.sam ATACseq_1231.sam
samtools view -b -S -o ATACseq_1231.bam ATACseq_1231.sam


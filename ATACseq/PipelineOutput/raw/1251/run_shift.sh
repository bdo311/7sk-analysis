#perl shift_sam_bases_edit.pl /home/raflynn/7SK/ChIRPseq/genes/mm9.sizes ATACseq_1251.trim.sort.nuc.uniq.rmdup.bam.sam ATACseq_1251.sam
cat header.sam ATACseq_1251.sam > tmp.1251.shifted.sam
mv tmp.1251.shifted.sam ATACseq_1251.sam
samtools view -b -S -o ATACseq_1251.bam ATACseq_1251.sam


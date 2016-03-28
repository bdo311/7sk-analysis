N {OFS="\t"}; {print $1,$2+800,$3-800,$4,$5,$6}' mES_allEnh_centered_1kb.bed > mES_allEnh_centered_200bp.bed

# BAF - done
nohup bedtools multicov -bams ChIP1_Baf155_12hrControlASO_rep1_genome_sorted.bam ChIP1_Baf155_12hrControlASO_rep2_genome_sorted.bam ChIP1_Baf155_12hr7SKASO_rep1_genome_sorted.bam ChIP1_Baf155_12hr7SKASO_rep2_genome_sorted.bam -bed ../mES_allEnh_centered_1kb.bed > bed/BAF_allEnh_1kb.bed &

nohup bedtools multicov -bams ChIP1_Baf155_12hrControlASO_rep1_genome_sorted.bam ChIP1_Baf155_12hrControlASO_rep2_genome_sorted.bam ChIP1_Baf155_12hr7SKASO_rep1_genome_sorted.bam ChIP1_Baf155_12hr7SKASO_rep2_genome_sorted.bam -bed ../mES_allEnh_centered_200bp.bed > bed/BAF_allEnh_200bp.bed &

nohup bedtools multicov -bams ChIP1_Baf155_12hrControlASO_rep1_genome_sorted.bam ChIP1_Baf155_12hrControlASO_rep2_genome_sorted.bam ChIP1_Baf155_12hr7SKASO_rep1_genome_sorted.bam ChIP1_Baf155_12hr7SKASO_rep2_genome_sorted.bam -bed ../mm9_tss_startseq_centered_TR_prox30.bed > bed/BAF_prom_prox.bed &

nohup bedtools multicov -bams ChIP1_Baf155_12hrControlASO_rep1_genome_sorted.bam ChIP1_Baf155_12hrControlASO_rep2_genome_sorted.bam ChIP1_Baf155_12hr7SKASO_rep1_genome_sorted.bam ChIP1_Baf155_12hr7SKASO_rep2_genome_sorted.bam -bed ../mES_prom_startseq_centered_1kb.bed > bed/BAF_prom_1kb.bed &

# Pou5f1
nohup bedtools multicov -bams ChIP1_Pou5f1_12hrControlASO_rep1_genome_sorted.bam ChIP1_Pou5f1_12hrControlASO_rep2_genome_sorted.bam ChIP1_Pou5f1_12hr7SKASO_rep1_genome_sorted.bam ChIP1_Pou5f1_12hr7SKASO_rep2_genome_sorted.bam -bed ../mES_allEnh_centered_1kb.bed > bed/Oct4_allEnh_1kb.bed &

nohup bedtools multicov -bams ChIP1_Pou5f1_12hrControlASO_rep1_genome_sorted.bam ChIP1_Pou5f1_12hrControlASO_rep2_genome_sorted.bam ChIP1_Pou5f1_12hr7SKASO_rep1_genome_sorted.bam ChIP1_Pou5f1_12hr7SKASO_rep2_genome_sorted.bam -bed ../mES_allEnh_centered_200bp.bed > bed/Oct4_allEnh_200bp.bed &

nohup bedtools multicov -bams ChIP1_Pou5f1_12hrControlASO_rep1_genome_sorted.bam ChIP1_Pou5f1_12hrControlASO_rep2_genome_sorted.bam ChIP1_Pou5f1_12hr7SKASO_rep1_genome_sorted.bam ChIP1_Pou5f1_12hr7SKASO_rep2_genome_sorted.bam -bed ../mm9_tss_startseq_centered_TR_prox30.bed > bed/Oct4_prom_prox.bed &

nohup bedtools multicov -bams ChIP1_Pou5f1_12hrControlASO_rep1_genome_sorted.bam ChIP1_Pou5f1_12hrControlASO_rep2_genome_sorted.bam ChIP1_Pou5f1_12hr7SKASO_rep1_genome_sorted.bam ChIP1_Pou5f1_12hr7SKASO_rep2_genome_sorted.bam -bed ../mES_prom_startseq_centered_1kb.bed > bed/Oct4_prom_1kb.bed &

# gH2AX
nohup bedtools multicov -bams ChIP3_gH2AX_12hrControlASO_rep1_genome_sorted.bam ChIP3_gH2AX_12hrControlASO_rep2_genome_sorted.bam ChIP3_gH2AX_12hr7SKASO_rep1_genome_sorted.bam ChIP3_gH2AX_12hr7SKASO_rep2_genome_sorted.bam -bed ../mES_allEnh_centered_1kb.bed > bed/gH2AX_allEnh_1kb.bed &

nohup bedtools multicov -bams ChIP3_gH2AX_12hrControlASO_rep1_genome_sorted.bam ChIP3_gH2AX_12hrControlASO_rep2_genome_sorted.bam ChIP3_gH2AX_12hr7SKASO_rep1_genome_sorted.bam ChIP3_gH2AX_12hr7SKASO_rep2_genome_sorted.bam -bed ../mES_allEnh_centered_200bp.bed > bed/gH2AX_allEnh_200bp.bed &

nohup bedtools multicov -bams ChIP3_gH2AX_12hrControlASO_rep1_genome_sorted.bam ChIP3_gH2AX_12hrControlASO_rep2_genome_sorted.bam ChIP3_gH2AX_12hr7SKASO_rep1_genome_sorted.bam ChIP3_gH2AX_12hr7SKASO_rep2_genome_sorted.bam -bed ../mm9_tss_startseq_centered_TR_prox30.bed > bed/gH2AX_prom_prox.bed &

nohup bedtools multicov -bams ChIP3_gH2AX_12hrControlASO_rep1_genome_sorted.bam ChIP3_gH2AX_12hrControlASO_rep2_genome_sorted.bam ChIP3_gH2AX_12hr7SKASO_rep1_genome_sorted.bam ChIP3_gH2AX_12hr7SKASO_rep2_genome_sorted.bam -bed ../mES_prom_startseq_centered_1kb.bed > bed/gH2AX_prom_1kb.bed &

# GRO - done
nohup bedtools multicov -bams GROseq_7SK_Scr_rep1_genome_sorted.bam GROseq_7SK_Scr_rep2_genome_sorted.bam GROseq_7SK_Scr_rep3_genome_sorted.bam GROseq_7SK_Scr_rep4_genome_sorted.bam GROseq_7SK_ASO_rep1_genome_sorted.bam GROseq_7SK_ASO_rep2_genome_sorted.bam GROseq_7SK_ASO_rep3_genome_sorted.bam GROseq_7SK_ASO_rep4_genome_sorted.bam -bed ../mES_allEnh_centered_1kb.bed > bed/GRO_allEnh_1kb.bed &

nohup bedtools multicov -bams GROseq_7SK_Scr_rep1_genome_sorted.bam GROseq_7SK_Scr_rep2_genome_sorted.bam GROseq_7SK_Scr_rep3_genome_sorted.bam GROseq_7SK_Scr_rep4_genome_sorted.bam GROseq_7SK_ASO_rep1_genome_sorted.bam GROseq_7SK_ASO_rep2_genome_sorted.bam GROseq_7SK_ASO_rep3_genome_sorted.bam GROseq_7SK_ASO_rep4_genome_sorted.bam -bed ../mES_allEnh_centered_200bp.bed > bed/GRO_allEnh_200bp.bed &

nohup bedtools multicov -bams GROseq_7SK_Scr_rep1_genome_sorted.bam GROseq_7SK_Scr_rep2_genome_sorted.bam GROseq_7SK_Scr_rep3_genome_sorted.bam GROseq_7SK_Scr_rep4_genome_sorted.bam GROseq_7SK_ASO_rep1_genome_sorted.bam GROseq_7SK_ASO_rep2_genome_sorted.bam GROseq_7SK_ASO_rep3_genome_sorted.bam GROseq_7SK_ASO_rep4_genome_sorted.bam -bed ../mm9_tss_startseq_centered_TR_dist30.bed > bed/GRO_prom_dist.bed &

nohup bedtools multicov -bams GROseq_7SK_Scr_rep1_genome_sorted.bam GROseq_7SK_Scr_rep2_genome_sorted.bam GROseq_7SK_Scr_rep3_genome_sorted.bam GROseq_7SK_Scr_rep4_genome_sorted.bam GROseq_7SK_ASO_rep1_genome_sorted.bam GROseq_7SK_ASO_rep2_genome_sorted.bam GROseq_7SK_ASO_rep3_genome_sorted.bam GROseq_7SK_ASO_rep4_genome_sorted.bam -bed ../mES_prom_startseq_centered_1kb.bed > bed/GRO_prom_1kb.bed &

nohup bedtools multicov -bams GROseq_7SK_Scr_rep1_genome_sorted.bam GROseq_7SK_Scr_rep2_genome_sorted.bam GROseq_7SK_Scr_rep3_genome_sorted.bam GROseq_7SK_Scr_rep4_genome_sorted.bam GROseq_7SK_ASO_rep1_genome_sorted.bam GROseq_7SK_ASO_rep2_genome_sorted.bam GROseq_7SK_ASO_rep3_genome_sorted.bam GROseq_7SK_ASO_rep4_genome_sorted.bam -bed ../mES_prom_startseq_centered_1kb.bed > bed/GRO_prom_1kb.bed &


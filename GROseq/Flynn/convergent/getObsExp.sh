#convergent_wt: 963013*100=96,301,300
#convergent_aso: 1349785*100=134,978,500

### SUPER ENHANCERS ###
bedtools intersect -a convergent_wt.bedGraph -b mES_SE_individual_1kb_BED6_sorted.bed | wc -l 
#obs: 5251*100=525,100
#total nt: 361*2000=722,000
#exp in genome: 96,301,000
#total in genome: 2,800,000,000
#exp: 24,832
#enrichment: 21.15
bedtools intersect -a convergent_aso.bedGraph -b mES_SE_individual_1kb_BED6_sorted.bed | wc -l 
#obs: 6418*100=641,800
#total nt: 361*2000=722,000
#exp in genome: 134,978,500
#total in genome: 2,800,000,000
#exp: 34,805
#enrichment: 18.44

### REGULAR ENHANCERS ###
#Old: bedtools intersect -a convergent_wt.bedGraph -b ~/7SK/ChIRPseq/genes/mES_reg_enhancers_RY_1kb_BED6.bed | wc -l
bedtools intersect -a convergent_wt.bedGraph -b mES_reg_enhancers_RY_startseq_centered_1kb_sorted.bed | wc -l
#obs: 25795*100 = 2,579,500
#total nt: 5356*2000=10,712,000
#exp in genome: 96,301,000
#total in genome: 2,800,000,000
#exp: 368420
#enrichment: 7.00
bedtools intersect -a convergent_aso.bedGraph -b mES_reg_enhancers_RY_startseq_centered_1kb_sorted.bed | wc -l 
#obs: 41034*100 = 4,103,400
#total nt: 5356*2000=10,712,000
#exp in genome: 134,978,500
#total in genome: 2,800,000,000
#exp: 516389
#enrichment: 7.95

### TSS ###
# Old: python makeTransPausedBed.py ~/7SK/ChIRPseq/genes/mm9_refseq_tss_BED6.bed ~/7SK/ChIRPseq/genes/subsets/mm9_refseq_transcribedOrPausedGenes.txt mm9_refseq_tss_transPaused.bed #7867
bedtools intersect -a convergent_wt.bedGraph -b mm9_tss_startseq_centered_1kb_sorted.bed | wc -l 
#obs: 124465*100 = 12,446,500
#total nt: 14234*2000=28,468,000
#exp in genome: 96,301,000
#total in genome: 2,800,000,000
#exp: 979,106
#enrichment: 12.71
bedtools intersect -a convergent_aso.bedGraph -b mm9_tss_startseq_centered_1kb_sorted.bed | wc -l 
#obs: 147047*100 = 14,704,700
#total nt: 7867*2000=28,468,000
#exp in genome: 134,978,500
#total in genome: 2,800,000,000
#exp: 1372346
#enrichment: 10.72

### CTCF SITES ###
bedtools intersect -a convergent_wt.bedGraph -b /home/raflynn/7SK/ATACseq/Nucleosome_signal_v3/Regions/CTCF_motif-peaks_noTSS_win2k.bed | wc -l 
#obs: 51344*100 = 5,134,400
#total nt: 32710*2000 = 65,420,000
#exp in genome: 96,301,000
#total in genome: 2,800,000,000
#exp: 2250004
#enrichment: 2.28
bedtools intersect -a convergent_aso.bedGraph -b /home/raflynn/7SK/ATACseq/Nucleosome_signal_v3/Regions/CTCF_motif-peaks_noTSS_win2k.bed | wc -l 
#obs: 74479*100 = 7,447,900
#total nt: 32710*2000=65,420,000
#exp in genome: 134,978,500
#total in genome: 2,800,000,000
#exp: 3153676
#enrichment: 2.36

### RANDOM INTERGENIC REGIONS ###
bedtools random -l 2000 -n 10000 -g /seq/chromosome/mm9/mm9.sizes | sort -k1,1 -k2,2n > mm9_random_loci.bed
bedtools shuffle -i mm9_random_loci.bed -g /seq/chromosome/mm9/mm9.sizes -seed 1 -excl /home/raflynn/7SK/ChIRPseq/genes/mm9_refseq_coll_gene.bed | sort -k1,1 -k2,2n > mm9_random_loci_1.bed
bedtools shuffle -i mm9_random_loci.bed -g /seq/chromosome/mm9/mm9.sizes -seed 2 -excl /home/raflynn/7SK/ChIRPseq/genes/mm9_refseq_coll_gene.bed | sort -k1,1 -k2,2n > mm9_random_loci_2.bed
bedtools shuffle -i mm9_random_loci.bed -g /seq/chromosome/mm9/mm9.sizes -seed 3 -excl /home/raflynn/7SK/ChIRPseq/genes/mm9_refseq_coll_gene.bed | sort -k1,1 -k2,2n > mm9_random_loci_3.bed
bedtools intersect -a convergent_wt.bedGraph -b mm9_random_loci_1.bed | wc -l 
#obs: 457,800 419,000 427,600 
#total nt: 10000*2000 = 20,000,000
#exp in genome: 96,301,000
#total in genome: 2,800,000,000
#exp: 687864
#enrichment: 0.67 0.61 0.62
bedtools intersect -a convergent_aso.bedGraph -b mm9_random_loci_1.bed | wc -l 
#obs: 665300 590000 620500
#total nt: 10000*2000 = 20,000,000
#exp in genome: 134,978,500
#total in genome: 2,800,000,000
#exp: 964132
#enrichment: 0.69 0.61 0.64

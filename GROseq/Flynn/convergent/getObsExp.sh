#convergent_wt: 963013*100=96,301,300
#convergent_aso: 1349785*100=134,978,500

### SUPER ENHANCERS ###
bedtools intersect -a convergent_wt.bedGraph -b ~/7SK/ChIRPseq/genes/mES_SE_individual_1kb_BED6.bed | wc -l 
#obs: 7841*100=784,100
#total nt: 559*2000=1118000
#exp in genome: 96,301,000
#total in genome: 2,800,000,000
#exp: 38,452
#enrichment: 20.39
bedtools intersect -a convergent_aso.bedGraph -b ~/7SK/ChIRPseq/genes/mES_SE_individual_1kb_BED6.bed | wc -l #
#obs: 9571*100=957,100
#total nt: 559*2000=1118000
#exp in genome: 134,978,500
#total in genome: 2,800,000,000
#exp: 53,895
#enrichment: 17.76

### REGULAR ENHANCERS ###
bedtools intersect -a convergent_wt.bedGraph -b ~/7SK/ChIRPseq/genes/mES_reg_enhancers_RY_1kb_BED6.bed | wc -l 
#obs: 28149*100 = 2,814,900
#total nt: 8563*2000=17,126,000
#exp in genome: 96,301,000
#total in genome: 2,800,000,000
#exp: 589,018
#enrichment: 4.78
bedtools intersect -a convergent_aso.bedGraph -b ~/7SK/ChIRPseq/genes/mES_reg_enhancers_RY_1kb_BED6.bed | wc -l 
#obs: 46574*100 = 4,657,400
#total nt: 8563*2000=17,126,000
#exp in genome: 134,978,500
#total in genome: 2,800,000,000
#exp: 825586
#enrichment: 5.64

### TSS ###
python makeTransPausedBed.py ~/7SK/ChIRPseq/genes/mm9_refseq_tss_BED6.bed ~/7SK/ChIRPseq/genes/subsets/mm9_refseq_transcribedOrPausedGenes.txt mm9_refseq_tss_transPaused.bed #7867
bedtools intersect -a convergent_wt.bedGraph -b mm9_refseq_tss_transPaused.bed | wc -l 
#obs: 83159*100 = 8,315,900
#total nt: 7867*2000=15,734,000
#exp in genome: 96,301,000
#total in genome: 2,800,000,000
#exp: 541,143
#enrichment: 15.37
bedtools intersect -a convergent_aso.bedGraph -b mm9_refseq_tss_transPaused.bed | wc -l 
#obs: 96300*100 = 9,630,000
#total nt: 7867*2000=15,734,000
#exp in genome: 134,978,500
#total in genome: 2,800,000,000
#exp: 785483
#enrichment: 12.26

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

### TELOMERES ###




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

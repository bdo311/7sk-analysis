#############
### SETUP ###
#############

# these can be multiplied by 100 since all blocks are 100bp each
cat convergent_wt_min.bedGraph | awk '{ sum += $4 } END { print sum }' 
#80729.4*100=8,072,940
cat convergent_aso_min.bedGraph | awk '{ sum += $4 } END { print sum }' 
#100776*100=10,077,600

cat ../norm_bedGraphs/GRO_125comb_negative_norm.bedGraph | python ~/Scripts/chirpseq_analysis/getTotalCount.py
cat ../norm_bedGraphs/GRO_125comb_positive_norm.bedGraph | python ~/Scripts/chirpseq_analysis/getTotalCount.py
cat ../norm_bedGraphs/GRO_12Ccomb_negative_norm.bedGraph | python ~/Scripts/chirpseq_analysis/getTotalCount.py
cat ../norm_bedGraphs/GRO_12Ccomb_positive_norm.bedGraph | python ~/Scripts/chirpseq_analysis/getTotalCount.py
# all 500,000,000 so total = 1,000,000,000

#######################
### SUPER ENHANCERS ###
#######################

## WT
bedtools intersect -a convergent_wt_min.bedGraph -b mES_SE_individual_1kb_BED6_sorted.bed | awk '{ sum += $4 } END { print sum }' 
#obs: 1109.45*100=110,945

bedtools intersect -a ../norm_bedGraphs/GRO_12Ccomb_positive_norm.bedGraph -b mES_SE_individual_1kb_BED6_sorted.bed | python ~/Scripts/chirpseq_analysis/getTotalCount.py
#obs: 2073866.21

bedtools intersect -a ../norm_bedGraphs/GRO_12Ccomb_negative_norm.bedGraph -b mES_SE_individual_1kb_BED6_sorted.bed | python ~/Scripts/chirpseq_analysis/getTotalCount.py
#obs: 1925190.79

## Total for gro: 3999057
## Percent of gro: 0.400%
## Percent of conv: 1.37%
## Conv%/Gro% = 3.425


## ASO
bedtools intersect -a convergent_aso_min.bedGraph -b mES_SE_individual_1kb_BED6_sorted.bed | awk '{ sum += $4 } END { print sum }' 
#obs: 2853.91*100=285,391

bedtools intersect -a ../norm_bedGraphs/GRO_125comb_positive_norm.bedGraph -b mES_SE_individual_1kb_BED6_sorted.bed | python ~/Scripts/chirpseq_analysis/getTotalCount.py
#obs: 2917985.15

bedtools intersect -a ../norm_bedGraphs/GRO_125comb_negative_norm.bedGraph -b mES_SE_individual_1kb_BED6_sorted.bed | python ~/Scripts/chirpseq_analysis/getTotalCount.py
#obs: 2775337.40

## Total for gro: 5693322.55
## Percent of gro: 0.569%
## Percent of conv: 2.83%
## Conv%/Gro% = 4.974

#########################
### REGULAR ENHANCERS ###
#########################

## WT
bedtools intersect -a convergent_wt_min.bedGraph -b mES_reg_enhancers_RY_startseq_centered_1kb_sorted.bed |awk '{ sum += $4 } END { print sum }' 
#obs: 2379.1*100=237,910

bedtools intersect -a ../norm_bedGraphs/GRO_12Ccomb_positive_norm.bedGraph -b mES_reg_enhancers_RY_startseq_centered_1kb_sorted.bed | python ~/Scripts/chirpseq_analysis/getTotalCount.py
#obs: 8225590.68986

bedtools intersect -a ../norm_bedGraphs/GRO_12Ccomb_negative_norm.bedGraph -b mES_reg_enhancers_RY_startseq_centered_1kb_sorted.bed | python ~/Scripts/chirpseq_analysis/getTotalCount.py
#obs: 9074629.84756

## Total for gro: 17300221
## Percent of gro: 1.73%
## Percent of conv: 2.95%
## Conv%/Gro% = 1.705

## ASO
bedtools intersect -a convergent_aso_min.bedGraph -b mES_reg_enhancers_RY_startseq_centered_1kb_sorted.bed | awk '{ sum += $4 } END { print sum }' 
#obs: 5105.9*100=510,590

bedtools intersect -a ../norm_bedGraphs/GRO_125comb_positive_norm.bedGraph -b mES_reg_enhancers_RY_startseq_centered_1kb_sorted.bed | python ~/Scripts/chirpseq_analysis/getTotalCount.py
#obs: 10363379.8813

bedtools intersect -a ../norm_bedGraphs/GRO_125comb_negative_norm.bedGraph -b mES_reg_enhancers_RY_startseq_centered_1kb_sorted.bed | python ~/Scripts/chirpseq_analysis/getTotalCount.py
#obs: 11307656.9809

## Total for gro: 21671037
## Percent of gro: 2.17%
## Percent of conv: 5.07%
## Conv%/Gro% = 2.336

###########
### TSS ###
###########

## WT

bedtools intersect -a convergent_wt_min.bedGraph -b mm9_tss_startseq_centered_1kb_sorted.bed | awk '{ sum += $4 } END { print sum }'  
#obs: 20573.6*100=2,057,360

bedtools intersect -a ../norm_bedGraphs/GRO_12Ccomb_positive_norm.bedGraph -b mm9_tss_startseq_centered_1kb_sorted.bed | python ~/Scripts/chirpseq_analysis/getTotalCount.py
#obs: 66450874.7295

bedtools intersect -a ../norm_bedGraphs/GRO_12Ccomb_negative_norm.bedGraph -b mm9_tss_startseq_centered_1kb_sorted.bed | python ~/Scripts/chirpseq_analysis/getTotalCount.py
#obs: 67870228.5552

## Total for gro: 134321104
## Percent of gro: 13.43%
## Percent of conv: 25.48%
## Conv%/Gro% = 1.897

## ASO
bedtools intersect -a convergent_aso_min.bedGraph -b mm9_tss_startseq_centered_1kb_sorted.bed | awk '{ sum += $4 } END { print sum }'  
#obs: 22813.4*100=2,281,340

bedtools intersect -a ../norm_bedGraphs/GRO_125comb_positive_norm.bedGraph -b mm9_tss_startseq_centered_1kb_sorted.bed | python ~/Scripts/chirpseq_analysis/getTotalCount.py
#obs: 66513815.6778

bedtools intersect -a ../norm_bedGraphs/GRO_125comb_negative_norm.bedGraph -b mm9_tss_startseq_centered_1kb_sorted.bed | python ~/Scripts/chirpseq_analysis/getTotalCount.py
#obs: 67597006.1964

## Total for gro: 134110822
## Percent of gro: 13.41%
## Percent of conv: 22.64%
## Conv%/Gro% = 1.689

##################
### CTCF SITES ###
##################

## WT
bedtools intersect -a convergent_wt_min.bedGraph -b /home/raflynn/7SK/ATACseq/Nucleosome_signal_v3/Regions/CTCF_motif-peaks_noTSS_win2k.bed | awk '{ sum += $4 } END { print sum }'  
#obs: 2377.43*100=237,743

bedtools intersect -a ../norm_bedGraphs/GRO_12Ccomb_positive_norm.bedGraph -b /home/raflynn/7SK/ATACseq/Nucleosome_signal_v3/Regions/CTCF_motif-peaks_noTSS_win2k.bed | python ~/Scripts/chirpseq_analysis/getTotalCount.py
#obs: 23150713.8693

bedtools intersect -a ../norm_bedGraphs/GRO_12Ccomb_negative_norm.bedGraph -b /home/raflynn/7SK/ATACseq/Nucleosome_signal_v3/Regions/CTCF_motif-peaks_noTSS_win2k.bed | python ~/Scripts/chirpseq_analysis/getTotalCount.py
#obs: 23049690.1998

## Total for gro: 46200404.0691
## Percent of gro: 4.62%
## Percent of conv: 2.94%

## ASO
bedtools intersect -a convergent_aso_min.bedGraph -b /home/raflynn/7SK/ATACseq/Nucleosome_signal_v3/Regions/CTCF_motif-peaks_noTSS_win2k.bed | awk '{ sum += $4 } END { print sum }'  
#obs: 3389.02*100=338,902

bedtools intersect -a ../norm_bedGraphs/GRO_125comb_positive_norm.bedGraph -b /home/raflynn/7SK/ATACseq/Nucleosome_signal_v3/Regions/CTCF_motif-peaks_noTSS_win2k.bed | python ~/Scripts/chirpseq_analysis/getTotalCount.py 
#obs: 23346693.7563

bedtools intersect -a ../norm_bedGraphs/GRO_125comb_negative_norm.bedGraph -b /home/raflynn/7SK/ATACseq/Nucleosome_signal_v3/Regions/CTCF_motif-peaks_noTSS_win2k.bed | python ~/Scripts/chirpseq_analysis/getTotalCount.py
#obs: 23122459.3069

## Total for gro: 46469153.0632
## Percent of gro: 4.65%
## Percent of conv: 3.36%


#################################
### RANDOM INTERGENIC REGIONS ###
#################################

## WT
bedtools random -l 2000 -n 10000 -g /seq/chromosome/mm9/mm9.sizes | sort -k1,1 -k2,2n > mm9_random_loci.bed
bedtools shuffle -i mm9_random_loci.bed -g /seq/chromosome/mm9/mm9.sizes -seed 1 -excl /home/raflynn/7SK/ChIRPseq/genes/mm9_refseq_coll_gene.bed | sort -k1,1 -k2,2n > mm9_random_loci_1.bed
bedtools shuffle -i mm9_random_loci.bed -g /seq/chromosome/mm9/mm9.sizes -seed 2 -excl /home/raflynn/7SK/ChIRPseq/genes/mm9_refseq_coll_gene.bed | sort -k1,1 -k2,2n > mm9_random_loci_2.bed
bedtools shuffle -i mm9_random_loci.bed -g /seq/chromosome/mm9/mm9.sizes -seed 3 -excl /home/raflynn/7SK/ChIRPseq/genes/mm9_refseq_coll_gene.bed | sort -k1,1 -k2,2n > mm9_random_loci_3.bed
# Random1 is messed up
bedtools intersect -a convergent_wt_min.bedGraph -b mm9_random_loci_2.bed | awk '{ sum += $4 } END { print sum }'  
#obs: 173.14 160.221
bedtools intersect -a ../norm_bedGraphs/GRO_12Ccomb_positive_norm.bedGraph -b mm9_random_loci_2.bed | python ~/Scripts/chirpseq_analysis/getTotalCount.py
#obs: 1887048.79738 1547942.01301
bedtools intersect -a ../norm_bedGraphs/GRO_12Ccomb_negative_norm.bedGraph -b mm9_random_loci_2.bed | python ~/Scripts/chirpseq_analysis/getTotalCount.py
#obs: 1524587.32713 1699303.45516

## Rand2
## Total for gro: 3411636.12451
## Percent of gro: 0.34%
## Percent of conv: 0.21%

## Rand3
## Total for gro: 3247245.46817
## Percent of gro: 0.32%
## Percent of conv: 0.20%

## ASO
bedtools intersect -a convergent_aso_min.bedGraph -b mm9_random_loci_1.bed | awk '{ sum += $4 } END { print sum }'  
#obs: 289.759 238.688
bedtools intersect -a ../norm_bedGraphs/GRO_125comb_negative_norm.bedGraph -b mm9_random_loci_2.bed | python ~/Scripts/chirpseq_analysis/getTotalCount.py
#obs: 1587216.67451 1793585.03174
bedtools intersect -a ../norm_bedGraphs/GRO_125comb_positive_norm.bedGraph -b mm9_random_loci_2.bed | python ~/Scripts/chirpseq_analysis/getTotalCount.py 
#obs: 1891477.06693 1655390.23845

## Rand2
## Total for gro: 3478693.74144
## Percent of gro: 0.35%
## Percent of conv: 0.29%

## Rand3
## Total for gro: 3448975.27019
## Percent of gro: 0.35%
## Percent of conv: 0.24%
########
## SE ##
########

nohup bedtools map -a mES_SE_individual_1kb_BED6_sorted.bed -b convergent_wt_min_sorted.bedGraph -c 4 -o mean -null '0' > scatter/se_conv_wt.txt &
nohup bedtools map -a mES_SE_individual_1kb_BED6_sorted.bed -b convergent_aso_min_sorted.bedGraph -c 4 -o mean -null '0' > scatter/se_conv_aso.txt &
nohup bedtools map -a mES_SE_individual_1kb_BED6_sorted.bed -b old/GRO_12Ccomb_positive_norm_sorted.bedGraph -c 4 -o mean -null '0' > scatter/se_gro_wt_pos.txt &
nohup bedtools map -a mES_SE_individual_1kb_BED6_sorted.bed -b old/GRO_12Ccomb_negative_norm_sorted.bedGraph -c 4 -o mean -null '0' > scatter/se_gro_wt_neg.txt &


nohup bedtools map -a mES_reg_enhancers_RY_startseq_centered_1kb_sorted.bed  -b convergent_wt_min_sorted.bedGraph -c 4 -o mean -null '0' > scatter/re_conv_wt.txt &
nohup bedtools map -a mES_reg_enhancers_RY_startseq_centered_1kb_sorted.bed  -b convergent_aso_min_sorted.bedGraph -c 4 -o mean -null '0' > scatter/re_conv_aso.txt &
nohup bedtools map -a mES_reg_enhancers_RY_startseq_centered_1kb_sorted.bed  -b old/GRO_12Ccomb_positive_norm_sorted.bedGraph -c 4 -o mean -null '0' > scatter/re_gro_wt_pos.txt &
nohup bedtools map -a mES_reg_enhancers_RY_startseq_centered_1kb_sorted.bed  -b old/GRO_12Ccomb_negative_norm_sorted.bedGraph -c 4 -o mean -null '0' > scatter/re_gro_wt_neg.txt &


nohup bedtools map -a mm9_tss_startseq_centered_1kb_sorted.bed -b convergent_wt_min_sorted.bedGraph -c 4 -o mean -null '0' > scatter/tss_conv_wt.txt &
nohup bedtools map -a mm9_tss_startseq_centered_1kb_sorted.bed -b convergent_aso_min_sorted.bedGraph -c 4 -o mean -null '0' > scatter/tss_conv_aso.txt &
nohup bedtools map -a mm9_tss_startseq_centered_1kb_sorted.bed -b old/GRO_12Ccomb_positive_norm_sorted.bedGraph -c 4 -o mean -null '0' > scatter/tss_gro_wt_pos.txt &
nohup bedtools map -a mm9_tss_startseq_centered_1kb_sorted.bed -b old/GRO_12Ccomb_negative_norm_sorted.bedGraph -c 4 -o mean -null '0' > scatter/tss_gro_wt_neg.txt &

nohup bedtools map -a /home/raflynn/7SK/ATACseq/Nucleosome_signal_v3/Regions/CTCF_motif-peaks_noTSS_win2k.bed -b convergent_wt_min_sorted.bedGraph -c 4 -o mean -null '0' > scatter/ctcf_conv_wt.txt &
nohup bedtools map -a /home/raflynn/7SK/ATACseq/Nucleosome_signal_v3/Regions/CTCF_motif-peaks_noTSS_win2k.bed -b convergent_aso_min_sorted.bedGraph -c 4 -o mean -null '0' > scatter/ctcf_conv_aso.txt &
nohup bedtools map -a /home/raflynn/7SK/ATACseq/Nucleosome_signal_v3/Regions/CTCF_motif-peaks_noTSS_win2k.bed -b old/GRO_12Ccomb_positive_norm_sorted.bedGraph -c 4 -o mean -null '0' > scatter/ctcf_gro_wt_pos.txt &
nohup bedtools map -a /home/raflynn/7SK/ATACseq/Nucleosome_signal_v3/Regions/CTCF_motif-peaks_noTSS_win2k.bed -b old/GRO_12Ccomb_negative_norm_sorted.bedGraph -c 4 -o mean -null '0' > scatter/ctcf_gro_wt_neg.txt &


nohup bedtools map -a mm9_random_loci_2.bed -b convergent_wt_min_sorted.bedGraph -c 4 -o mean -null '0' > scatter/rand2_conv_wt.txt &
nohup bedtools map -a mm9_random_loci_2.bed -b convergent_aso_min_sorted.bedGraph -c 4 -o mean -null '0' > scatter/rand2_conv_aso.txt &
nohup bedtools map -a mm9_random_loci_2.bed -b old/GRO_12Ccomb_positive_norm_sorted.bedGraph -c 4 -o mean -null '0' > scatter/rand2_gro_wt_pos.txt &
nohup bedtools map -a mm9_random_loci_2.bed -b old/GRO_12Ccomb_negative_norm_sorted.bedGraph -c 4 -o mean -null '0' > scatter/rand2_gro_wt_neg.txt &

nohup bedtools map -a mm9_random_loci_3.bed -b convergent_wt_min_sorted.bedGraph -c 4 -o mean -null '0' > scatter/rand3_conv_wt.txt &
nohup bedtools map -a mm9_random_loci_3.bed -b convergent_aso_min_sorted.bedGraph -c 4 -o mean -null '0' > scatter/rand3_conv_aso.txt &
nohup bedtools map -a mm9_random_loci_3.bed -b old/GRO_12Ccomb_positive_norm_sorted.bedGraph -c 4 -o mean -null '0' > scatter/rand3_gro_wt_pos.txt &
nohup bedtools map -a mm9_random_loci_3.bed -b old/GRO_12Ccomb_negative_norm_sorted.bedGraph -c 4 -o mean -null '0' > scatter/rand3_gro_wt_neg.txt &


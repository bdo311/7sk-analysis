/home/raflynn/Scripts/MACS-master/bin/macs2 bdgbroadcall -i ~/7SK/ChIRPseq/bedgraphs/mES_Untreated_genome_merged_norm.bedGraph -o mES_Untreated_genome_merged_norm_broad.bed
/home/raflynn/Scripts/MACS-master/bin/macs2 bdgbroadcall -i ~/7SK/ChIPseq/Flynn/bedgraphs/BAF155_Scr_merged.bedGraph -o BAF155_Scr_merged_broad.bed
/home/raflynn/Scripts/MACS-master/bin/macs2 bdgbroadcall -i ~/7SK/ChIPseq/Flynn/bedgraphs/mES_WT_DMSO1_BAF155_genome_shifted_norm.bedGraph -o mES_WT_DMSO1_BAF155_genome_shifted_norm_broad.bed
/home/raflynn/Scripts/MACS-master/bin/macs2 bdgbroadcall -i ~/7SK/ChIPseq/Flynn/bedgraphs/mES_WT_DMSO2_BAF155_genome_shifted_norm.bedGraph -o mES_WT_DMSO2_BAF155_genome_shifted_norm_broad.bed
/home/raflynn/Scripts/MACS-master/bin/macs2 bdgbroadcall -i ~/7SK/ChIPseq/Flynn/bedgraphs/mES_Scr1_HEXIM1_genome_shifted_norm.bedGraph -o mES_Scr1_HEXIM1_genome_shifted_norm_broad.bed
/home/raflynn/Scripts/MACS-master/bin/macs2 bdgbroadcall -i ~/7SK/ChIPseq/Flynn/bedgraphs/mES_Scr2_HEXIM1_genome_shifted_norm.bedGraph -o mES_Scr2_HEXIM1_genome_shifted_norm_broad.bed
/home/raflynn/Scripts/MACS-master/bin/macs2 bdgbroadcall -i ~/7SK/ChIPseq/Flynn/bedgraphs/mES_HEXIM1_genome_shifted_norm.bedGraph -o mES_HEXIM1_genome_shifted_norm_broad.bed
/home/raflynn/Scripts/MACS-master/bin/macs2 bdgbroadcall -i ~/7SK/ChIPseq/Flynn/bedgraphs/mES_ChIPseq_DDX21_2i_genome_shifted_norm.bedGraph -o mES_ChIPseq_DDX21_2i_genome_shifted_norm_broad.bed


awk -F '\t' 'BEGIN {OFS="\t"} {print $1,$2,$3,"7sk__"NR"__"$5,$5,$6}' mES_Untreated_genome_merged_norm.bed > named_mES_Untreated_genome_merged_norm.bed
awk -F '\t' 'BEGIN {OFS="\t"} {print $1,$2,$3,"BAF155_Scr__"NR"__"$5,$5,$6}' BAF155_Scr_merged.bed > named_BAF155_Scr_merged_norm.bed
awk -F '\t' 'BEGIN {OFS="\t"} {print $1,$2,$3,"BAF155_DMSO1__"NR"__"$5,$5,$6}' mES_WT_DMSO1_BAF155_genome_shifted_norm.bed > named_mES_WT_DMSO1_BAF155_genome_shifted_norm.bed
awk -F '\t' 'BEGIN {OFS="\t"} {print $1,$2,$3,"BAF155_DMSO2__"NR"__"$5,$5,$6}' mES_WT_DMSO2_BAF155_genome_shifted_norm.bed > named_mES_WT_DMSO2_BAF155_genome_shifted_norm.bed
awk -F '\t' 'BEGIN {OFS="\t"} {print $1,$2,$3,"HEXIM1_Scr1__"NR"__"$5,$5,$6}' mES_Scr1_HEXIM1_genome_shifted_norm.bed > named_mES_Scr1_HEXIM1_genome_shifted_norm.bed
awk -F '\t' 'BEGIN {OFS="\t"} {print $1,$2,$3,"HEXIM1_Scr2__"NR"__"$5,$5,$6}' mES_Scr2_HEXIM1_genome_shifted_norm.bed > named_mES_Scr2_HEXIM1_genome_shifted_norm.bed
awk -F '\t' 'BEGIN {OFS="\t"} {print $1,$2,$3,"HEXIM1__"NR"__"$5,$5,$6}' mES_HEXIM1_genome_shifted_norm.bed > named_mES_HEXIM1_genome_shifted_norm.bed
awk -F '\t' 'BEGIN {OFS="\t"} {print $1,$2,$3,"DDX21__"NR"__"$5,$5,$6}' mES_ChIPseq_DDX21_2i_genome_shifted_norm.bed > named_mES_ChIPseq_DDX21_2i_genome_shifted_norm.bed

files=$(echo named*norm.bed)
bedops --everything $files > MasterList_narrow.bed
bedmap --echo --count --echo-map-id MasterList_narrow.bed > MasterList_narrow_overlap.bed
python collapseOverlapPeaks.py



/home/raflynn/Scripts/MACS-master/bin/macs2 bdgpeakcall -i ~/7SK/ChIRPseq/bedgraphs/mES_Untreated_genome_merged_norm.bedGraph -o mES_Untreated_genome_merged_norm.bed
/home/raflynn/Scripts/MACS-master/bin/macs2 bdgpeakcall -i ~/7SK/ChIPseq/Flynn/bedgraphs/BAF155_Scr_merged.bedGraph -o BAF155_Scr_merged.bed
/home/raflynn/Scripts/MACS-master/bin/macs2 bdgpeakcall -i ~/7SK/ChIPseq/Flynn/bedgraphs/mES_WT_DMSO1_BAF155_genome_shifted_norm.bedGraph -o mES_WT_DMSO1_BAF155_genome_shifted_norm.bed
/home/raflynn/Scripts/MACS-master/bin/macs2 bdgpeakcall -i ~/7SK/ChIPseq/Flynn/bedgraphs/mES_WT_DMSO2_BAF155_genome_shifted_norm.bedGraph -o mES_WT_DMSO2_BAF155_genome_shifted_norm.bed
/home/raflynn/Scripts/MACS-master/bin/macs2 bdgpeakcall -i ~/7SK/ChIPseq/Flynn/bedgraphs/mES_Scr1_HEXIM1_genome_shifted_norm.bedGraph -o mES_Scr1_HEXIM1_genome_shifted_norm.bed
/home/raflynn/Scripts/MACS-master/bin/macs2 bdgpeakcall -i ~/7SK/ChIPseq/Flynn/bedgraphs/mES_Scr2_HEXIM1_genome_shifted_norm.bedGraph -o mES_Scr2_HEXIM1_genome_shifted_norm.bed
/home/raflynn/Scripts/MACS-master/bin/macs2 bdgpeakcall -i ~/7SK/ChIPseq/Flynn/bedgraphs/mES_HEXIM1_genome_shifted_norm.bedGraph -o mES_HEXIM1_genome_shifted_norm.bed
/home/raflynn/Scripts/MACS-master/bin/macs2 bdgpeakcall -i ~/7SK/ChIPseq/Flynn/bedgraphs/mES_ChIPseq_DDX21_2i_genome_shifted_norm.bedGraph -o mES_ChIPseq_DDX21_2i_genome_shifted_norm.bed
/home/raflynn/Scripts/MACS-master/bin/macs2 bdgpeakcall -i ~/7SK/ChIPseq/Flynn/bedgraphs/mES_ChIPseq_DDX21_AFK_genome_shifted_norm.bedGraph -o mES_ChIPseq_DDX21_AFK_genome_shifted_norm.bed

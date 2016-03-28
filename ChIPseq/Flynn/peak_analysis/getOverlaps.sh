awk -F '\t' 'BEGIN {OFS="\t"} {print $1,$2,$3,"BAF155_rep1__"NR"__"$8,$8,"+"}' baf155_rep1.txt > baf155_rep1.bed
awk -F '\t' 'BEGIN {OFS="\t"} {print $1,$2,$3,"BAF155_rep2__"NR"__"$8,$8,"+"}' baf155_rep2.txt > baf155_rep2.bed
awk -F '\t' 'BEGIN {OFS="\t"} {print $1,$2,$3,"HEXIM1_rep1__"NR"__"$8,$8,"+"}' hexim1_rep1.txt > hexim1_rep1.bed
awk -F '\t' 'BEGIN {OFS="\t"} {print $1,$2,$3,"HEXIM1_rep2__"NR"__"$8,$8,"+"}' hexim1_rep2.txt > hexim1_rep2.bed
awk -F '\t' 'BEGIN {OFS="\t"} {print $1,$2,$3,"7SK_WT__"NR"__"10,10,"+"}' distinct_peaks_3.txt > 7SK_WT3.bed

# parallel "sed -i 's/^M//g' {}" ::: *.bed
# parallel "sed -i 's/+/+\n/g' {}" ::: *.bed

files=$(echo *.bed)
bedops --everything $files > MasterList.bed
bedmap --echo --count --echo-map-id MasterList.bed > MasterList_overlap.bed

# python ~/7SK/ChIRPseq/other_analyses/overlap_peaks/collapseOverlapPeaks.py
# python ~/7SK/ChIRPseq/other_analyses/overlap_peaks/jaccardPeaks.py
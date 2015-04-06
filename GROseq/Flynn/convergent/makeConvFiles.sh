# makeConvFiles.sh

# make windows
bedtools makewindows -g /seq/chromosome/mm9/mm9.sizes -w 100 > mm9_100.bed
bedtools makewindows -g /seq/chromosome/mm9/mm9.sizes -w 100 -s 50 > mm9_100_50.bed

# map gro to our 100bp windows
parallel "bedtools map -a mm9_100.bed -b {} -c 4 -null "0" -g /seq/chromosome/mm9/mm9.sizes > {.}_100.bed" ::: ../norm_bedGraphs/GRO_125comb_positive_norm.bedGraph  ../norm_bedGraphs/GRO_125comb_negative_norm.bedGraph ../norm_bedGraphs/GRO_12Ccomb_positive_norm.bedGraph ../norm_bedGraphs/GRO_12Ccomb_negative_norm.bedGraph
parallel "bedtools map -a mm9_100_50.bed -b {} -c 4 -null "0" -g /seq/chromosome/mm9/mm9.sizes > {.}_100_50.bed" ::: ../norm_bedGraphs/GRO_125comb_positive_norm.bedGraph  ../norm_bedGraphs/GRO_125comb_negative_norm.bedGraph ../norm_bedGraphs/GRO_12Ccomb_positive_norm.bedGraph ../norm_bedGraphs/GRO_12Ccomb_negative_norm.bedGraph
mv ../norm_bedGraphs/*.bed ./

# columns
nohup cut -f4 GRO_125comb_negative_norm_100_50.bed > 125n_50.txt &
nohup cut -f4 GRO_125comb_positive_norm_100_50.bed > 125p_50.txt &
nohup cut -f4 GRO_12Ccomb_negative_norm_100_50.bed > 12Cn_50.txt &
nohup cut -f4 GRO_12Ccomb_positive_norm_100_50.bed > 12Cp_50.txt &
nohup cut -f4 GRO_125comb_negative_norm_100.bed > 125n.txt &
nohup cut -f4 GRO_125comb_positive_norm_100.bed > 125p.txt &
nohup cut -f4 GRO_12Ccomb_negative_norm_100.bed > 12Cn.txt &
nohup cut -f4 GRO_12Ccomb_positive_norm_100.bed > 12Cp.txt &

# combine everything
paste mm9_100.bed 125p.txt 125n.txt 12Cp.txt 12Cn.txt > gro_100.bed 
paste mm9_100_50.bed 125p_50.txt 125n_50.txt 12Cp_50.txt 12Cn_50.txt > gro_100_50.bed

# run python script
python countChanges.py
# processHumanChIRP.sh
# 5/22/14
# not intended to be run all at once (because of the nohups), but this is what was run

# 1 tabbing the files
cd H1
#sed -i 's/ /\t/g' H1_even.bedGraph > even_tabbed.bedGraph
nohup sed -i 's/ /\t/g' H1_odd.bedGraph > odd_tabbed.bedGraph &
nohup sed -i 's/ /\t/g' H1_input.bedGraph > input_tabbed.bedGraph &

cd ../HeLa
nohup sed -i 's/ /\t/g' HeLa_even.bedGraph > even_tabbed.bedGraph &
nohup sed -i 's/ /\t/g' HeLa_odd.bedGraph > odd_tabbed.bedGraph &
nohup sed -i 's/ /\t/g' HeLa_input.bedGraph > input_tabbed.bedGraph &

# 2 combining even and odd
cd ../H1
nohup unionBedGraphs -i even_tabbed.bedGraph odd_tabbed.bedGraph -header -names Even Odd > combined_unnorm.bedGraph &
cd ../HeLa
nohup unionBedGraphs -i even_tabbed.bedGraph odd_tabbed.bedGraph -header -names Even Odd > combined_unnorm.bedGraph &

# 2.1 merging even and odd
cd ../H1
nohup python ../mergeCombined.py combined_unnorm.bedGraph combined.bedGraph &
cd ../HeLa
nohup python ../mergeCombined.py combined_unnorm.bedGraph combined.bedGraph &

# 3 normalize to 750 million total read count; make bedgraph and bigwig
cd ..
python normalizeBedGraph.py

# 4 bw and upload
cd ../H1
nohup ../aws put "x-amz-acl: public-read" changseq/bdo/H1_combined.bw combined_norm.bw &

cd ../HeLa
nohup ../aws put "x-amz-acl: public-read" changseq/bdo/HeLa_combined.bw combined_norm.bw &
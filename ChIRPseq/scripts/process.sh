# parallel "cd {}; macs14 -t even.sam -n even_macs -g mm -B -S" ::: ActD DRB Flavo JQ11 JQ14 
# parallel "cd {}; macs14 -t odd.sam -n odd_macs -g mm -B -S" ::: ActD DRB Flavo JQ11 JQ14 
# parallel "cd {}/even_macs_MACS_bedGraph/treat/; gunzip *; python ~/ChIRPseq/normalizeBedgraph_new.py *.bdg even_norm.bedGraph" ::: ActD DRB Flavo JQ11 JQ14
# parallel "cd {}/odd_macs_MACS_bedGraph/treat/; gunzip *; python ~/ChIRPseq/normalizeBedgraph_new.py *.bdg odd_norm.bedGraph" ::: ActD DRB Flavo JQ11 JQ14
# parallel "cd {}/; unionBedGraph -i even_macs_MACS_bedGraph/treat/even_norm.bedGraph odd_macs_MACS_bedGraph/treat/odd_norm.bedGraph -header > norm_union.bedGraph" ::: ActD DRB Flavo JQ11 JQ14
# parallel "cd {}; python ~/ChIRPseq/takeLower.py norm_union.bedGraph norm.bedGraph" ::: ActD DRB Flavo JQ11 JQ14 
# parallel -j 1 "cd {}; sort -k1,1 -k2,2n norm.bedGraph > norm_sorted.bedGraph" ::: ActD DRB Flavo JQ11 JQ14

bedtools unionbedg -i norm_sorted.bedGraph /home/raflynn/ChIRPseq/7sk_removal/mouse_7sk_withsides_all_sorted.bedGraph -header > norm_union_7sk.bedGraph
python ~/ChIRPseq/7sk_removal/remove7SK.py norm_union_7sk.bedGraph norm_no7sk.bedGraph
python ~/ChIRPseq/normalizeBedgraph_new.py norm_no7sk.bedGraph norm_combined.bedGraph
python ~/ChIRPseq/collapseBedgraph.py norm_combined.bedGraph norm_collapsed.bedGraph
bedGraphToBigWig norm_combined.bedGraph /seq/chromosome/mm9/mm9.sizes norm_combined.bw

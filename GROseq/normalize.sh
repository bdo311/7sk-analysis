python normalizeGroBedgraph.py GRO_123_mm9.bedGraph GRO_123_mm9_norm.bedGraph 24825211 &
python normalizeGroBedgraph.py GRO_125_mm9.bedGraph GRO_125_mm9_norm.bedGraph 22629580 &
python normalizeGroBedgraph.py GRO_12C_mm9.bedGraph GRO_12C_mm9_norm.bedGraph 19579653 &
python normalizeGroBedgraph.py GRO_63_mm9.bedGraph GRO_63_mm9_norm.bedGraph 28706630 &
python normalizeGroBedgraph.py GRO_65_mm9.bedGraph GRO_65_mm9_norm.bedGraph 25828161 &
python normalizeGroBedgraph.py GRO_6C_mm9.bedGraph GRO_6C_mm9_norm.bedGraph 29965093 &

mv GRO_123_mm9.bedGraph old
mv GRO_125_mm9.bedGraph old
mv GRO_12C_mm9.bedGraph old
mv GRO_65_mm9.bedGraph old
mv GRO_63_mm9.bedGraph old
mv GRO_6C_mm9.bedGraph old

parallel "bedGraphToBigWig {} /seq/chromosome/mm9/mm9.sizes {.}.bw" ::: *_mm9_norm.bedGraph
parallel "~/Scripts/aws put \"x-amz-acl: public-read\" changseq/rflynn/GROseq/{} {}" ::: *.bw

python normalizeGroBedgraph.py GRO_123_repeat.bedGraph GRO_123_repeat_norm.bedGraph 14466067 &
python normalizeGroBedgraph.py GRO_125_repeat.bedGraph GRO_125_repeat_norm.bedGraph 13214754 &
python normalizeGroBedgraph.py GRO_12C_repeat.bedGraph GRO_12C_repeat_norm.bedGraph 24988826 &
python normalizeGroBedgraph.py GRO_63_repeat.bedGraph GRO_63_repeat_norm.bedGraph 9341506 &
python normalizeGroBedgraph.py GRO_65_repeat.bedGraph GRO_65_repeat_norm.bedGraph 8377713 &
python normalizeGroBedgraph.py GRO_6C_repeat.bedGraph GRO_6C_repeat_norm.bedGraph 12425087 &

mv GRO_123_repeat.bedGraph old
mv GRO_125_repeat.bedGraph old
mv GRO_12C_repeat.bedGraph old
mv GRO_65_repeat.bedGraph old
mv GRO_63_repeat.bedGraph old
mv GRO_6C_repeat.bedGraph old
# here is how to use this pipeline
# 1. get the fasta files for ncRNA and target DNA
# 2. put the files in the folder 'sequences'
# 3. write down file name below
RNAfile=Mm_7SK.fa
DNAfile=WT2_combined_peaks.fasta

# 4. set output folder below
output=Mm_7SK_WT2_best_match

# 5. modify the parameters only if necessary
gap_penalty=100

script=./script
sequences=./sequences

$script/smith_waterman $sequences/$RNAfile.shuffled $sequences/$DNAfile.shuffled $script/identity.scoremat $gap_penalty  rc > $output/RNA.shuffled-DNA.shuffled.alignment.rc
python $script/sw2bed.py $output/RNA.shuffled-DNA.shuffled.alignment.rc rc

# plot
Rscript $script/plot_best_match.r $output  

# here is how to use this pipeline
# 1. get the fasta files for ncRNA and target DNA
# 2. put the files in the folder 'sequences'
# 3. write down file name below
RNAfile=Mm_Line1spa_2.fa
DNAfile=Mm_Line1ChIRP_odd.fasta

# 4. set output folder below
output=Mm_Line1ChIRP_odd

# 5. modify the parameters only if necessary
gap_penalty=100

# 6. save this file and run it by typing './pipeline.sh', will take a few minutes to run.

########## 

script=./script
sequences=./sequences

# make the output folder
mkdir $output

# dinucleotide shuffle the sequence
python $script/fasta-dinucleotide-shuffle.py -f $sequences/$RNAfile -t -shuffled > $sequences/$RNAfile.shuffled
python $script/fasta-dinucleotide-shuffle.py -f $sequences/$DNAfile -t -shuffled > $sequences/$DNAfile.shuffled

# fix shuffle
python $script/fixShuffled.py $sequences/$DNAfile $sequences/$DNAfile.shuffled

# find best matches
$script/smith_waterman $sequences/$RNAfile $sequences/$DNAfile $script/identity.scoremat $gap_penalty > $output/RNA-DNA.alignment
python $script/sw2bed.py $output/RNA-DNA.alignment norc
$script/smith_waterman $sequences/$RNAfile $sequences/$DNAfile $script/identity.scoremat $gap_penalty  rc > $output/RNA-DNA.alignment.rc
python $script/sw2bed.py $output/RNA-DNA.alignment.rc rc
$script/smith_waterman $sequences/$RNAfile $sequences/$DNAfile.shuffled $script/identity.scoremat $gap_penalty > $output/RNA-DNA.shuffled.alignment
python $script/sw2bed.py $output/RNA-DNA.shuffled.alignment norc
$script/smith_waterman $sequences/$RNAfile $sequences/$DNAfile.shuffled $script/identity.scoremat $gap_penalty  rc > $output/RNA-DNA.shuffled.alignment.rc
python $script/sw2bed.py $output/RNA-DNA.shuffled.alignment.rc rc
$script/smith_waterman $sequences/$RNAfile.shuffled $sequences/$DNAfile $script/identity.scoremat $gap_penalty > $output/RNA.shuffled-DNA.alignment
python $script/sw2bed.py $output/RNA.shuffled-DNA.alignment norc
$script/smith_waterman $sequences/$RNAfile.shuffled $sequences/$DNAfile $script/identity.scoremat $gap_penalty  rc > $output/RNA.shuffled-DNA.alignment.rc
python $script/sw2bed.py $output/RNA.shuffled-DNA.alignment.rc rc
$script/smith_waterman $sequences/$RNAfile.shuffled $sequences/$DNAfile.shuffled $script/identity.scoremat $gap_penalty > $output/RNA.shuffled-DNA.shuffled.alignment
python $script/sw2bed.py $output/RNA.shuffled-DNA.shuffled.alignment norc
$script/smith_waterman $sequences/$RNAfile.shuffled $sequences/$DNAfile.shuffled $script/identity.scoremat $gap_penalty  rc > $output/RNA.shuffled-DNA.shuffled.alignment.rc
python $script/sw2bed.py $output/RNA.shuffled-DNA.shuffled.alignment.rc rc

# plot
Rscript $script/plot_best_match.r $output  

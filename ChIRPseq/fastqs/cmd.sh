/home/raflynn/Scripts/cutadapt-1.7.1/bin/cutadapt -b CTGTCTCTTATACACATCTCCGAGCCCACGAGA -m 25  -o C2_7SK_Even_R1.tmp.fastq -p C2_7SK_Even_R2.tmp.fastq   C2_7SK_Even_R1.fastq.gz C2_7SK_Even_R2.fastq.gz  
/home/raflynn/Scripts/cutadapt-1.7.1/bin/cutadapt  -b CTGTCTCTTATACACATCTGACGCTGCCGACGA  -m 25   -o C2_7SK_Even_R2.trimmed.fastq.gz  -p C2_7SK_Even_R1.trimmed.fastq.gz   C2_7SK_Even_R2.tmp.fastq C2_7SK_Even_R1.tmp.fastq

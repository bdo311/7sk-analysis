### CONVERT TO TAGALIGN ###

# individual files
	# ChIP
nohup samtools view -b -F 1548 -q 30 ChIP2_BAF_DMSO1_genome_sorted.bam | bamToBed -i stdin | awk 'BEGIN{FS="\t";OFS="\t"}{$4="N"; print $0}' | gzip -c > BAM/peaks/ChIP2_BAF_DMSO1.tagAlign.gz &
nohup samtools view -b -F 1548 -q 30 ChIP2_BAF_DMSO2_genome_sorted.bam | bamToBed -i stdin | awk 'BEGIN{FS="\t";OFS="\t"}{$4="N"; print $0}' | gzip -c > BAM/peaks/ChIP2_BAF_DMSO2.tagAlign.gz &
nohup samtools view -b -F 1548 -q 30 ChIP2_Input_genome_sorted.bam | bamToBed -i stdin | awk 'BEGIN{FS="\t";OFS="\t"}{$4="N"; print $0}' | gzip -c > peaks/ChIP2_Input.tagAlign.gz &
nohup samtools view -b -F 1548 -q 30 ChIP3_HEXIM1_Scr1_genome_sorted.bam | bamToBed -i stdin | awk 'BEGIN{FS="\t";OFS="\t"}{$4="N"; print $0}' | gzip -c > BAM/peaks/ChIP3_HEXIM1_Scr1.tagAlign.gz &  # ChIP3_Scr1 == ChIP3.pr1
nohup samtools view -b -F 1548 -q 30 ChIP3_HEXIM1_Scr2_genome_sorted.bam | bamToBed -i stdin | awk 'BEGIN{FS="\t";OFS="\t"}{$4="N"; print $0}' | gzip -c > BAM/peaks/ChIP3_HEXIM1_Scr2.tagAlign.gz &  # ChIP3_Scr2 == ChIP3.pr2
nohup zcat ChIP3_HEXIM1_Scr1.tagAlign.gz ChIP3_HEXIM1_Scr2.tagAlign.gz | gzip -c > ChIP3_HEXIM1_Scr.tagAlign.gz &  # one replicate from these two
nohup samtools view -b -F 1548 -q 30 ChIP2_Hexim1_genome_sorted.bam | bamToBed -i stdin | awk 'BEGIN{FS="\t";OFS="\t"}{$4="N"; print $0}' | gzip -c > peaks/ChIP2_Hexim1.tagAlign.gz &  

	# ChIRP
nohup samtools view -b -F 1548 -q 30 Mouse_7SK_WT_even_genome_sorted.bam | bamToBed -i stdin | awk 'BEGIN{FS="\t";OFS="\t"}{$4="N"; print $0}' | gzip -c > peaks/Mouse_7SK_WT_even.tagAlign.gz &
nohup samtools view -b -F 1548 -q 30 Mouse_7SK_WT_odd_genome_sorted.bam | bamToBed -i stdin | awk 'BEGIN{FS="\t";OFS="\t"}{$4="N"; print $0}' | gzip -c > peaks/Mouse_7SK_WT_odd.tagAlign.gz &
nohup samtools view -b -F 1548 -q 30 Mouse_7SK_Input_even_genome_sorted.bam | bamToBed -i stdin | awk 'BEGIN{FS="\t";OFS="\t"}{$4="N"; print $0}' | gzip -c > peaks/Mouse_7SK_Input_even.tagAlign.gz &
nohup samtools view -b -F 1548 -q 30 Mouse_7SK_Input_odd_genome_sorted.bam | bamToBed -i stdin | awk 'BEGIN{FS="\t";OFS="\t"}{$4="N"; print $0}' | gzip -c > peaks/Mouse_7SK_Input_odd.tagAlign.gz &

# pooled
nohup zcat ChIP2_BAF_DMSO1.tagAlign.gz ChIP2_BAF_DMSO2.tagAlign.gz | gzip -c > ChIP2_BAF_DMSO.tagAlign.gz &
nohup zcat ChIP2_Hexim1.tagAlign.gz ChIP3_HEXIM1_Scr1.tagAlign.gz ChIP3_HEXIM1_Scr2.tagAlign.gz | gzip -c > ChIP23_HEXIM1.tagAlign.gz &  # ChIP23 is the pooled HEXIM
nohup zcat Mouse_7SK_WT_even.tagAlign.gz Mouse_7SK_WT_odd.tagAlign.gz | gzip -c > Mouse_7SK_WT.tagAlign.gz &
nohup zcat Mouse_7SK_Input_even.tagAlign.gz Mouse_7SK_Input_odd.tagAlign.gz | gzip -c > Mouse_7SK_Input.tagAlign.gz &

### PSEUDOREPLICATES FOR CHIP SAMPLES ### # ChIP3 already done
zcat ChIP2_BAF_DMSO1.tagAlign.gz | wc -l  # 33209465 --> 16604733
zcat ChIP2_BAF_DMSO1.tagAlign.gz | shuf | split -l 16604733
nohup cat xaa | gzip -c > ChIP2_BAF_DMSO1.pr1.tagAlign.gz &
nohup cat xab | gzip -c > ChIP2_BAF_DMSO1.pr2.tagAlign.gz &

zcat ChIP2_BAF_DMSO2.tagAlign.gz | wc -l  # 36364929 --> 18182500
zcat ChIP2_BAF_DMSO2.tagAlign.gz | shuf | split -l 18182500
nohup cat xaa | gzip -c > ChIP2_BAF_DMSO2.pr1.tagAlign.gz &
nohup cat xab | gzip -c > ChIP2_BAF_DMSO2.pr2.tagAlign.gz &

zcat ChIP2_Hexim1.tagAlign.gz | wc -l  # 30499504 --> 15250000
zcat ChIP2_Hexim1.tagAlign.gz | shuf | split -l 15250000
nohup cat xaa | gzip -c > ChIP2_Hexim1.pr1.tagAlign.gz &
nohup cat xab | gzip -c > ChIP2_Hexim1.pr2.tagAlign.gz &

zcat ChIP23_HEXIM1.tagAlign.gz | wc -l # 110094263 --> 55050000
zcat ChIP23_HEXIM1.tagAlign.gz | shuf | split -l 55050000
nohup cat xaa | gzip -c > ChIP23_HEXIM1.pr1.tagAlign.gz &
nohup cat xab | gzip -c > ChIP23_HEXIM1.pr2.tagAlign.gz &

zcat ChIP2_BAF_DMSO.tagAlign.gz | shuf | split -l 35000000
nohup cat xaa | gzip -c > ChIP2_BAF_DMSO.pr1.tagAlign.gz &
nohup cat xab | gzip -c > ChIP2_BAF_DMSO.pr2.tagAlign.gz &
rm -f xaa xab 

zcat Mouse_7SK_WT_even.tagAlign.gz | wc -l  # 45911357 --> 23000000
zcat Mouse_7SK_WT_even.tagAlign.gz | shuf | split -l 23000000
nohup cat xaa | gzip -c > Mouse_7SK_WT_even.pr1.tagAlign.gz & 
nohup cat xab | gzip -c > Mouse_7SK_WT_even.pr2.tagAlign.gz & 

zcat Mouse_7SK_WT_odd.tagAlign.gz | wc -l  # 74750268 --> 37380000
zcat Mouse_7SK_WT_odd.tagAlign.gz | shuf | split -l 37380000
nohup cat xaa | gzip -c > Mouse_7SK_WT_odd.pr1.tagAlign.gz &
nohup cat xab | gzip -c > Mouse_7SK_WT_odd.pr2.tagAlign.gz &

zcat Mouse_7SK_WT.tagAlign.gz | shuf | split -l 60350000
nohup cat xaa | gzip -c > Mouse_7SK_WT.pr1.tagAlign.gz & 
nohup cat xab | gzip -c > Mouse_7SK_WT.pr2.tagAlign.gz &

### MACS2 ###

# individual replicates vs pooled input
nohup macs2 callpeak -t ChIP2_BAF_DMSO1.tagAlign.gz -c ChIP2_Input.tagAlign.gz -f BED -n ChIP2_BAF_DMSO1_VS_ChIP2_Input -g mm -p 1e-3 --to-large & # done
nohup macs2 callpeak -t ChIP2_BAF_DMSO2.tagAlign.gz -c ChIP2_Input.tagAlign.gz -f BED -n ChIP2_BAF_DMSO2_VS_ChIP2_Input -g mm -p 1e-3 --to-large & # done
nohup macs2 callpeak -t ChIP3_HEXIM1_Scr.tagAlign.gz -c ChIP2_Input.tagAlign.gz -f BED -n ChIP3_HEXIM1_Scr_VS_ChIP2_Input -g mm -p 1e-3 --to-large &
nohup macs2 callpeak -t ChIP2_Hexim1.tagAlign.gz -c ChIP2_Input.tagAlign.gz -f BED -n ChIP2_Hexim1_VS_ChIP2_Input -g mm -p 1e-3 --to-large &
nohup macs2 callpeak -t Mouse_7SK_WT_even.tagAlign.gz -c Mouse_7SK_Input.tagAlign.gz -f BED -n Mouse_7SK_WT_even_VS_Mouse_7SK_Input -g mm -p 1e-3 --to-large & # done
nohup macs2 callpeak -t Mouse_7SK_WT_odd.tagAlign.gz -c Mouse_7SK_Input.tagAlign.gz -f BED -n Mouse_7SK_WT_odd_VS_Mouse_7SK_Input -g mm -p 1e-3 --to-large & # done

# pooled replicates vs pooled input
nohup macs2 callpeak -t ChIP2_BAF_DMSO.tagAlign.gz ll
-c ChIP2_Input.tagAlign.gz -f BED -n ChIP2_BAF_DMSO_VS_ChIP2_Input -g mm -p 1e-3 --to-large & # done
nohup macs2 callpeak -t ChIP23_HEXIM1.tagAlign.gz -c ChIP23_HEXIM1.tagAlign.gz -f BED -n ChIP23_HEXIM1_VS_ChIP2_Input -g mm -p 1e-3 --to-large &
nohup macs2 callpeak -t Mouse_7SK_WT.tagAlign.gz -c Mouse_7SK_Input.tagAlign.gz -f BED -n Mouse_7SK_WT_VS_Mouse_7SK_Input -g mm -p 1e-3 --to-large & # done

# individual pseudoreplicates vs pooled input
nohup macs2 callpeak -t ChIP2_BAF_DMSO1.pr1.tagAlign.gz -c ChIP2_Input.tagAlign.gz -f BED -n ChIP2_BAF_DMSO1_Pr1_VS_ChIP2_Input -g mm -p 1e-3 --to-large & # done
nohup macs2 callpeak -t ChIP2_BAF_DMSO1.pr2.tagAlign.gz -c ChIP2_Input.tagAlign.gz -f BED -n ChIP2_BAF_DMSO1_Pr2_VS_ChIP2_Input -g mm -p 1e-3 --to-large & # done
nohup macs2 callpeak -t ChIP2_BAF_DMSO2.pr1.tagAlign.gz -c ChIP2_Input.tagAlign.gz -f BED -n ChIP2_BAF_DMSO2_Pr1_VS_ChIP2_Input -g mm -p 1e-3 --to-large & # done
nohup macs2 callpeak -t ChIP2_BAF_DMSO2.pr2.tagAlign.gz -c ChIP2_Input.tagAlign.gz -f BED -n ChIP2_BAF_DMSO2_Pr2_VS_ChIP2_Input -g mm -p 1e-3 --to-large & # done
nohup macs2 callpeak -t ChIP3_HEXIM1_Scr1.tagAlign.gz -c ChIP2_Input.tagAlign.gz -f BED -n ChIP3_HEXIM1_Scr_Pr1_VS_ChIP2_Input -g mm -p 1e-3 --to-large &
nohup macs2 callpeak -t ChIP3_HEXIM1_Scr2.tagAlign.gz -c ChIP2_Input.tagAlign.gz -f BED -n ChIP3_HEXIM1_Scr_Pr2_VS_ChIP2_Input -g mm -p 1e-3 --to-large &
nohup macs2 callpeak -t ChIP2_Hexim1.pr1.tagAlign.gz -c ChIP2_Input.tagAlign.gz -f BED -n ChIP2_Hexim1_Pr1_VS_ChIP2_Input -g mm -p 1e-3 --to-large &
nohup macs2 callpeak -t ChIP2_Hexim1.pr2.tagAlign.gz -c ChIP2_Input.tagAlign.gz -f BED -n ChIP2_Hexim1_Pr2_VS_ChIP2_Input -g mm -p 1e-3 --to-large &
nohup macs2 callpeak -t Mouse_7SK_WT_even.pr1.tagAlign.gz -c Mouse_7SK_Input.tagAlign.gz -f BED -n Mouse_7SK_WT_even_Pr1_VS_Mouse_7SK_Input -g mm -p 1e-3 --to-large & # done
nohup macs2 callpeak -t Mouse_7SK_WT_even.pr2.tagAlign.gz -c Mouse_7SK_Input.tagAlign.gz -f BED -n Mouse_7SK_WT_even_Pr2_VS_Mouse_7SK_Input -g mm -p 1e-3 --to-large & # done
nohup macs2 callpeak -t Mouse_7SK_WT_odd.pr1.tagAlign.gz -c Mouse_7SK_Input.tagAlign.gz -f BED -n Mouse_7SK_WT_odd_Pr1_VS_Mouse_7SK_Input -g mm -p 1e-3 --to-large & # done
nohup macs2 callpeak -t Mouse_7SK_WT_odd.pr2.tagAlign.gz -c Mouse_7SK_Input.tagAlign.gz -f BED -n Mouse_7SK_WT_odd_Pr2_VS_Mouse_7SK_Input -g mm -p 1e-3 --to-large & # done

# pooled pseudoreplicates vs pooled input
nohup macs2 callpeak -t ChIP2_BAF_DMSO.pr1.tagAlign.gz -c ChIP2_Input.tagAlign.gz -f BED -n ChIP2_BAF_DMSO_Pr1_VS_ChIP2_Input -g mm -p 1e-3 --to-large & # done
nohup macs2 callpeak -t ChIP2_BAF_DMSO.pr2.tagAlign.gz -c ChIP2_Input.tagAlign.gz -f BED -n ChIP2_BAF_DMSO_Pr2_VS_ChIP2_Input -g mm -p 1e-3 --to-large & # done
nohup macs2 callpeak -t ChIP23_HEXIM1.pr1.tagAlign.gz -c ChIP2_Input.tagAlign.gz -f BED -n ChIP23_HEXIM1_Pr1_VS_ChIP2_Input -g mm -p 1e-3 --to-large &
nohup macs2 callpeak -t ChIP23_HEXIM1.pr2.tagAlign.gz -c ChIP2_Input.tagAlign.gz -f BED -n ChIP23_HEXIM1_Pr2_VS_ChIP2_Input -g mm -p 1e-3 --to-large &
nohup macs2 callpeak -t Mouse_7SK_WT.pr1.tagAlign.gz -c Mouse_7SK_Input.tagAlign.gz -f BED -n Mouse_7SK_WT_Pr1_VS_Mouse_7SK_Input -g mm -p 1e-3 --to-large & # done
nohup macs2 callpeak -t Mouse_7SK_WT.pr2.tagAlign.gz -c Mouse_7SK_Input.tagAlign.gz -f BED -n Mouse_7SK_WT_Pr2_VS_Mouse_7SK_Input -g mm -p 1e-3 --to-large & # done

### IDR ###

# sort peaks, take top 100k

# IDR

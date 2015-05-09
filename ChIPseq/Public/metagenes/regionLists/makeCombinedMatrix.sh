cut -f2- ~/7SK/ChIRPseq/metagenes/regionLists/TSS_centered_list.txt > 7sk_tss.txt
cut -f2- ~/7SK/ChIRPseq/metagenes/regionLists/RY_enh_centered_list.txt > 7sk_re.txt
cut -f2- ~/7SK/ChIRPseq/metagenes/regionLists/SE_indiv_list.txt > 7sk_se.txt

cut -f2- ~/7SK/ATACseq/PipelineOutput/metagenes/regionLists/TSS_centered_list.txt > atac_tss.txt
cut -f2- ~/7SK/ATACseq/PipelineOutput/metagenes/regionLists/RY_enh_centered_list.txt > atac_re.txt
cut -f2- ~/7SK/ATACseq/PipelineOutput/metagenes/regionLists/SE_indiv_list.txt > atac_se.txt

cut -f2- ~/7SK/ChIPseq/Flynn/Second_run/metagenes/regionLists/TSS_centered_list.txt > flynn_tss.txt
cut -f2- ~/7SK/ChIPseq/Flynn/Second_run/metagenes/regionLists/RY_enh_centered_list.txt > flynn_re.txt
cut -f2- ~/7SK/ChIPseq/Flynn/Second_run/metagenes/regionLists/SE_indiv_list.txt > flynn_se.txt

paste TSS_centered_list.txt *_tss.txt > TSS_centered_list_all.txt
paste RY_enh_centered_list.txt *_re.txt > RY_enh_centered_list_all.txt
paste SE_indiv_list.txt *_se.txt > SE_indiv_list_all.txt

rm -f 7sk* atac* flynn*
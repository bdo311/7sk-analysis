for file in `ls  *.pe.q10.sort.bam`
do
echo "start $file"

/home/jinxu/bin/SNPsplit    --snp_file /home/jinxu/DB/mouse/mmu9/129_Cast_SNP_inform/129_Cast_SNP_mm9.txt2.chr  --paired  --conflicting $file  1>$file.log 2>$file.err 

echo "$file done "

done


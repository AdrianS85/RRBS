export FASTQ_DEDUP_FOL=/tmp/Analysis/ &&

cd Bismark &&

cat nudup.txt >> nudup_old.txt; rm nudup.txt
cat fq2.txt >> fq2_old.txt; rm fq2.txt
cat pairs.txt >> pairs_old.txt; rm pairs.txt
cat pairs_to_do.txt >> pairs_to_do_old.txt; rm pairs_to_do.txt
cat jolog.txt >> jolog_old.txt; rm jolog.txt
cat nudup_raport.txt >> nudup_raport_old.txt; rm nudup_raport.txt

#FROM INSIDE THE CONTAINER?
ls *stripped.sam | sort > nudup.txt &&
ls $FASTQ_DEDUP_FOL*_R2_* | sort > fq2.txt &&
paste fq2.txt nudup.txt > pairs.txt &&

ls *stripped.sorted.dedup.bam | sed 's/_R1_.*//' > already_done.txt &&
grep -v -f already_done.txt pairs.txt > pairs_to_do.txt &&

parallel --verbose --joblog jolog.txt --tmpdir /tmp/ --jobs 5 --colsep "\t" "python /nudup-master/nudup.py --rmdup-only  -T /tmp/ --paired-end -f {1} -o {2.} {2}" :::: pairs_to_do.txt &> nudup_raport.txt &&

parallel "samtools sort -n -o {.}.sorted.bam {}" ::: *sorted.dedup.bam &> sam2_raport.txt &&

parallel bamqc ::: *dedup.sorted.bam &&

parallel "bismark_methylation_extractor --bedGraph --paired-end --comprehensive --merge_non_CpG" ::: *dedup.sorted.bam &> bedgraph.txt &&

mkdir Bedgraph_rest && mv *.sorted.txt Bedgraph_rest &&

bismark2report &&
bismark2summary &&
multiqc -n Bismark_multiqc.html . &&
rm -R ../work/
rm -R /tmp/nudup*

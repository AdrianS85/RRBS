### This adress is needed to 
export FASTQ_DEDUP_FOL=/tmp/

cd Bismark &&

#FROM INSIDE THE CONTAINER?
ls *stripped.bam | sort >> nudup.txt &&
ls $FASTQ_DEDUP_FOL*_R2_* | sort >> fq2.txt &&
# paste fq2.txt nudup.txt >> pairs
parallel --verbose --link --joblog jolog.txt --memfree 20 --retries 10 --delay 1 --tmpdir /tmp/tmp/ --jobs 10  "python /nugentechnologies-nudup-7a126eb/nudup.py --rmdup-only  -T /tmp/tmp/ --paired-end -f {1} -o {2.} {2}" :::: fq2.txt nudup.txt &> nudup_raport.txt
#FROM OUTSIDE THE CONTAINER?
#parallel --verbose --joblog jolog.txt --colsep '\t' "singularity exec --cleanenv --bind ./Bind:/tmp --pwd /tmp/Analysis/Bismark --writable RRBS_Singularity_New.simg python /nudup-master/nudup.py --rmdup-only  --paired-end -f {1} -o {2.} {2}" :::: Bind/Analysis/Bismark/pairs &> Bind/Analysis/Bismark/nudup_raport.txt

parallel "samtools sort -n -o {.}.sorted.bam {}" ::: *sorted.dedup.bam &> sam2_raport.txt &&

parallel bamqc ::: *dedup.sorted.bam &&

parallel "bismark_methylation_extractor --bedGraph --paired-end --comprehensive --merge_non_CpG" ::: *dedup.sorted.bam &> bedgraph.txt &&

mkdir Bedgraph_rest && mv *.sorted.txt Bedgraph_rest &&

bismark2report &&
bismark2summary &&
multiqc -n Bismark_multiqc.html . &&
rm -R ../work/

#

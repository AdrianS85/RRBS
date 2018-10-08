### This adress is needed to 
export FASTQ_DEDUP_FOL=/tmp/Analysis/ &&

cd Bismark &&

#FROM INSIDE THE CONTAINER?
ls *stripped.bam | sort >> nudup.txt &&
ls $FASTQ_DEDUP_FOL*_R2_* | sort >> fq2.txt &&
paste fq2.txt nudup.txt >> pairs

parallel --verbose --link --joblog jolog.txt--tmpdir /tmp/tmp/ --jobs 5  "python /nugentechnologies-nudup-7a126eb/nudup.py --rmdup-only  -T /tmp/tmp/ --paired-end -f {1} -o {2.} {2}" :::: pairs &> nudup_raport.txt

parallel "samtools sort -n -o {.}.sorted.bam {}" ::: *sorted.dedup.bam &> sam2_raport.txt &&

parallel bamqc ::: *dedup.sorted.bam &&

parallel "bismark_methylation_extractor --bedGraph --paired-end --comprehensive --merge_non_CpG" ::: *dedup.sorted.bam &> bedgraph.txt &&

mkdir Bedgraph_rest && mv *.sorted.txt Bedgraph_rest &&

bismark2report &&
bismark2summary &&
multiqc -n Bismark_multiqc.html . &&
rm -R ../work/

#singularity shell --cleanenv --bind ./Analysis:/tmp --pwd /tmp/ --writable RRBS_Singularity_New.simg

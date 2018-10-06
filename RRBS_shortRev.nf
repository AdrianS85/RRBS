RP2 = Channel.fromFilePairs("*_R2_001.fastq.gz", size: -1, flat: true) 
RP3 = Channel.fromFilePairs("Bismark/*_R1_001_val_1.fq_trimmed_bismark_bt2_pe.bam", size: -1, flat: true) 
dparallel = 5 //This option marks how many nudup instances can I run concurrently


params.gen_fol = "/tmp/Analysis/Genome/"
params.outputFR = "Fastqc_Raw"
params.outputFRM = "Fastqc_Raw/Multiqc"
params.outputTG = "Trim_Galore"
params.outputTGM = "Trim_Galore/Multiqc"
params.outputD = "Diversity"
params.outputFD = "Fastqc_Div"
params.outputFDM = "Fastqc_Div/Multiqc"
params.outputB = "Bismark"
params.outputDD = "Deduplication"
params.outputC = "Calling"



process DeDuplicationPrep {
         input:
         set val(ID), file(ddp1) from RP3

         output:
         set ID, file("${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_pe.sam_stripped.sam") into DDP_out
          
         """
         samtools view -h -o ${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_pe.sam ${ddp1} &&
         sleep 2 && 
         strip_bismark_sam.sh ${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_pe.sam &&
         sleep 2 
         """
}

DDP2_outDD = Channel.create()
DDP2_outDD = DDP_out.join(RP2)

process DeDuplication {
         publishDir path: params.outputDD, mode: 'copy', pattern: "*log.txt"
         
         maxForks = dparallel

         input:
         set val(ID), file(dd1), file(dd3) from DDP2_outDD

         output:
         set ID, file("${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_pe.sam_stripped.sorted.dedup.bam") into (DD_out1, DD_out2)
         file("${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_pe.sam_stripped_dup_log.txt")

         """
         python /nudup-master/nudup.py -T /tmp/ --rmdup-only  --paired-end -f ${dd3} -o ${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_pe.sam_stripped ${dd1} &&
         sleep 10
         """
}

/BB_Biopsy_62_S15_R1_001_val_1.fq_trimmed_bismark_bt2_pe.sam_stripped.sorted.dedup.bam

DD_out1.subscribe onComplete: {

process DeDuplicationPOst {
         publishDir path: params.outputDD, mode: 'copy' 

         input:
         set val(ID), file(ddpo1) from (DD_out2)

         output:
         set ID, file("${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_pe.sam_stripped.sorted.dedup.sorted.bam") into DDPO_out
         set file("${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_pe.sam_stripped.sorted.dedup.sorted_bamqc.zip"), file("${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_pe.sam_stripped.sorted.dedup.sorted_bamqc.html")

         """
         samtools sort -n -o ${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_pe.sam_stripped.sorted.dedup.sorted.bam ${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_pe.sam_stripped.sorted.dedup.bam &&
         sleep 2 &&
         bamqc ${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_pe.sam_stripped.sorted.dedup.sorted.bam 
         """
}


process Calling {
         publishDir path: params.outputC, mode: 'copy', pattern: "*gz"
         publishDir path: params.outputB, mode: 'copy', pattern: "*txt" 

         input:
         set val(ID), file(c1) from DDPO_out	

         output:
         set file("${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_pe.sam_stripped.sorted.dedup.sorted.bismark.cov.gz"), file("${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_pe.sam_stripped.sorted.dedup.sorted.bedGraph.gz"), file("${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_pe.sam_stripped.sorted.dedup.sorted.M-bias.txt"), file("${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_pe.sam_stripped.sorted.dedup.sorted_splitting_report.txt")

         """
         bismark_methylation_extractor --bedGraph --paired-end --comprehensive --merge_non_CpG ${c1}
         """
}}

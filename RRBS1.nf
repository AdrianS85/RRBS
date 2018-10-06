#!/usr/bin/env nextflow

/* 
-with-report -with-trace -with-timeline -with-dag
-resume
nextflow run RRBS.nf



WARNING!! 
USE export SINGULARITY_CACHEDIR=/mnt/iscsi/Adrian/BRAIN_RRBS/Analysis/tmp/
export SINGULARITY_LOCALCACHEDIR=/mnt/iscsi/Adrian/BRAIN_RRBS/Analysis/tmp/

I am not sure if that was neccesary once -T command was established in nudup.py. Need to check enviroment of Brain analysis once its finished.
*/




params.all = "*.fastq.gz"
ALL = file(params.all).flatten()
RP = Channel.fromFilePairs("*_R{1,3}_001.fastq.gz", size: -1, flat: true) //https://groups.google.com/forum/#!topic/nextflow/X4YyYmLTbTo
RP2 = Channel.fromFilePairs("*_R2_001.fastq.gz", size: -1, flat: true) 
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



process Fastqc_Raw {
         publishDir path: params.outputFR, mode: 'copy' //INFO: for some reason path dir can only be established by referencing parameter  

         input:
         file ALL

         output:
         file '*_fastqc.{zip,html}' into (FR_out1, FR_out2) //INFO: {} brackets include OR operator

         """
         fastqc $ALL;
         sleep 1
         """
}

FR_out2.subscribe onComplete: {
         FR_out3 = Channel.create()
         FR_out3 = FR_out1.collect()

         process Fastqc_Raw_Multi {
                  publishDir path: params.outputFRM, mode: 'copy'

                  input:
                  file FR_out3

                  output:
                  set file("Raw_Fastqc_multiqc.html"), file("Raw_Fastqc_multiqc_data.tar.gz")

                  """
                  multiqc -n Raw_Fastqc_multiqc.html ${FR_out3}; tar -czf Raw_Fastqc_multiqc_data.tar.gz Raw_Fastqc_multiqc_data;
                  sleep 1
                  """
}}




process Trim_Galore {
         publishDir path: params.outputTG, mode: 'copy'  

         input:
         set val(ID), file(tg1), file(tg3) from RP

         output:
         set ID, file("${ID}_R1_001_val_1.fq.gz"), file("${ID}_R3_001_val_2.fq.gz") into (TG_out) //INFO: we can transfer values from input to output by directly referencing value: "ID"
         file "${ID}_R{1,3}_001.fastq.gz_trimming_report.txt" into (TG_out1, TG_out2)
         """
         trim_galore --paired -a AGATCGGAAGAGC -a2 AAATCAAAAAAAC ${tg1} ${tg3};
         sleep 1
         """
}

////WARNING! This will work only if Calling will be performed during given run
TG_out2.subscribe onComplete: {
         TG_out3 = Channel.create()
         TG_out3 = TG_out1.collect()

         process Trim_Multi {
         publishDir path: params.outputTGM, mode: 'copy'

         input:
         file TG_out3

         output:
         set file("Trim_multiqc.html"), file("Trim_multiqc_data.tar.gz")

         """
         multiqc -n Trim_multiqc.html ${TG_out3}; tar -czf Trim_multiqc_data.tar.gz Trim_multiqc_data;
         sleep 1
         """
}}




process Diversity {
         publishDir path: params.outputD, mode: 'copy'  

         input:
         set val(ID), file(d1), file(d2) from TG_out	

         output:
         set ID, file("${ID}_R1_001_val_1.fq_trimmed.fq.gz"), file("${ID}_R3_001_val_2.fq_trimmed.fq.gz") into (D_out, D_outB)

         """
         python /trimRRBSdiversityAdaptCustomers.py -1 $d1 -2 $d2;
         sleep 1
         """
}




process Fastqc_Div {
         publishDir path: params.outputFD, mode: 'copy' 

         input:
         set val(ID), file(fd1), file(fd2) from D_out

         output:
         file '*_fastqc.{zip,html}' into (FD_out1, FD_out2)

         """
         fastqc $fd1; fastqc $fd2;
         sleep 1
         """
}

////WARNING! This will work only if Calling will be performed during given run
FD_out2.subscribe onComplete: {
         FD_out3 = Channel.create()
         FD_out3 = FD_out1.collect()

         process Fastqc_Div_Multi {
                  publishDir path: params.outputFDM, mode: 'copy'

                  input:
                  file FD_out3

                  output:
                  set file("Fastqc_Div_multiqc.html"), file("Fastqc_Div_multiqc_data.tar.gz")

                  """
                  multiqc -n Fastqc_Div_multiqc.html ${FD_out3}; tar -czf Fastqc_Div_multiqc_data.tar.gz Fastqc_Div_multiqc_data;
                  sleep 1
                  """
}}




process Bismark {
         publishDir path: params.outputB, mode: 'copy' 

         input:
         set val(ID), file(b1), file(b2) from D_outB	

         output:
         set ID, file("${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_pe.bam"), file("${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_PE_report.txt") into (B_out, B_outDD)

         """
         bismark --bowtie2 --genome_folder ${params.gen_fol}  -1 $b1 -2 $b2;
         sleep 10
         """ //gen_fol should be always aviable based on how the singularity is made and run
}




process Bismark_Report {
         publishDir path: params.outputB, mode: 'copy' 

         input:
         set val(ID), file(br1), file(br2) from B_out	

         output:
         set val(ID), file("${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_pe.nucleotide_stats.txt"), file("${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_pe_bamqc.html"), file("${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_pe_bamqc.zip")

         """
         bam2nuc --genome_folder ${params.gen_fol} $br1; bamqc $br1;
         sleep 10
         """
}




process DeDuplicationPrep {
         input:
         set val(ID), file(ddp1), file(ddp3) from B_outDD

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
         set ID, file("${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_pe.sam_stripped.sorted.dedup.bam") into (DD_out1, DD_out2, DD_out3)
         file("${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_pe.sam_stripped_dup_log.txt")

         """
         python /nudup-master/nudup.py -T /tmp/ --rmdup-only  --paired-end -f ${dd3} -o ${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_pe.sam_stripped ${dd1} &&
         sleep 10
         """
}



DD_out2.subscribe onComplete: {
         DD_out4 = Channel.create()
         DD_out4 = DD_out1.collect()

process DeDuplicationPOst {
         publishDir path: params.outputDD, mode: 'copy' 

         input:
         set val(ID), file(ddpo1) from (DD_out3)
         file from DD_out4

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


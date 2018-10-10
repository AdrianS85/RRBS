// nextflow run -with-report -with-trace -with-timeline -with-deg -name ##name## RRBS.nf ##-resume##
params.all = "*.fastq.gz"
ALL = file(params.all).flatten()
RP = Channel.fromFilePairs("*_R{1,3}_001.fastq.gz", size: -1, flat: true) //https://groups.google.com/forum/#!topic/nextflow/X4YyYmLTbTo
RP2 = Channel.fromFilePairs("*_R2_001.fastq.gz", size: -1, flat: true) 



params.gen_fol = "/tmp/Analysis/Genome/"
params.outputFR = "Fastqc_Raw"
params.outputFRM = "Fastqc_Raw/Multiqc"
params.outputTG = "Trim_Galore"
params.outputTGM = "Trim_Galore/Multiqc"
params.outputD = "Diversity"
params.outputFD = "Fastqc_Div"
params.outputFDM = "Fastqc_Div/Multiqc"
params.outputB = "Bismark"



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
         publishDir path: params.outputB, mode: 'copy'

         input:
         set val(ID), file(ddp1), file(ddp3) from B_outDD

         output:
         set ID, file("${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_pe.sam_stripped.sam")
          
         """
         samtools view -h -o ${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_pe.sam ${ddp1} &&
         sleep 2 && 
         strip_bismark_sam.sh ${ID}_R1_001_val_1.fq_trimmed_bismark_bt2_pe.sam &&
         sleep 2 
         """
}

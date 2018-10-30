#GET THE TOOLS AND ANNOTATIONS
library(RnBeads)
REGION_SET <- "tiling200bp"
ASSEMBLY <- "mm10"
rnb.load.annotation.from.db(REGION_SET, assembly=ASSEMBLY) ## This loads data from RnBeads site. This annotation should be already included in region.types option if it is set as NULL.
#GET THE TOOLS AND ANNOTATIONS


#SETUP WORKING ENVIROMENT
#setwd("D:/Zycie_zawodowe/Fede_seq/rnbeads")
setwd(getwd())
data_dir <- paste0(getwd(), "/data")
dataset_dir <- data_dir
sample_annotation <- paste0(data_dir, "/annot.csv")
analysis_dir <- paste0(getwd(), "/analysis_dir")
report_dir <- paste0(getwd(), "/reports")
#SETUP WORKING ENVIROMENT


#SETUP PARALLELISM
parallel.setup(24)
parallel.isEnabled()
#SETUP PARALLELISM


#SETUP RUN OPTIONS
rnb.options(
  analysis.name = "Fede_Brain1",
  assembly = "mm10",
  region.aggregation = "coverage.weighted",
  identifiers.column = "sampleID", ###Can behave wierdly. doesnt line "_"? integers?
  #import.bed.test.only = T,
  #IMPORT
  import.bed.style = "bismarkCov",
  #import.default.data.type = "bs.bed.dir", #This only sets default
  #import.gender.prediction = F,
  #import.bed.frame.shift
  #QC
  qc.coverage.plots = T,
  qc.coverage.histograms = T,
  qc.coverage.violins = F,
  qc.coverage.threshold.plot = c(1, 2, 3, 5, 10),
  #FILTERING
  filtering.snp = "no",
  filtering.missing.value.quantile = 0.4,
  filtering.coverage.threshold = 5,
  filtering.low.coverage.masking = T,
  filtering.high.coverage.outliers = T, 
  filtering.sex.chromosomes.removal = T,
  #filtering.deviation.threshold = 0,
  #imputation.method = "none"
  #INFERENCE
  inference = T,
  inference.targets.sva = c("b_vs_rest", "b_vs_c", "i_vs_rest", "i_vs_c", "e_vs_c"), ## Not 100% if this is ok
  inference.sva.num.method = "leek",
  inference.age.prediction = F,
  #inference.reference.methylome.column = 
  #EXPLORATORY
  exploratory.columns = c("treatment", "prep_batch", "seq_batch", "b_vs_rest", "b_vs_c", "i_vs_rest", "i_vs_c", "e_vs_c"),
  exploratory.clustering.top.sites = 200,
  #exploratory.gene.symbols  = 
  #exploratory.deviation.plots = 
  #exploratory.custom.loci.bed = 
  #DIFFERENTIAL
  differential.comparison.columns = c("b_vs_rest", "b_vs_c", "i_vs_rest", "i_vs_c", "e_vs_c"),
  differential.enrichment.go = F,
  differential.enrichment.lola = T,
  differential.enrichment.lola.dbs = "/tmp/Analysis/LOLACore17", ##This needs intervention - ?? "mm10/" ??
  differential.adjustment.sva = T,
  covariate.adjustment.columns = c("prep_batch", "seq_batch")
  #min.group.size = 
  #differential.comparison.columns.all.pairwise = 
  #differential.variability = 
  #differential.variability.method = 
  #differential.permutations = 
  #OTHER
  #export.to.ewasher
  )
#SETUP RUN OPTIONS


#STEP-BY-STEP ANALYSIS
## DATA_IMPORT -> (QUALITY_CONTROL, PREPROCESSING), PREPROCESSING -> INFERENCE -> (EXPLORATION, DIFFERENTIAL)
result <- rnb.run.import(data.source=c(dataset_dir, sample_annotation),
                         data.type="bs.bed.dir", 
                         dir.reports=report_dir) ## results in a list with two elements: the dataset (rnb.set) and a report

rnb.run.qc(result$rnb.set, report_dir) ## Outputs just report i think. I think it needs to write in /dev
PP_result <- rnb.run.preprocessing(result$rnb.set, report_dir)
I_result <- rnb.run.inference(PP_result$rnb.set, report_dir)
E_result <- rnb.run.exploratory(I_result$rnb.set, report_dir)
D_result <- rnb.run.differential(I_result$rnb.set, report_dir)
rnb.execute.tnt(I_result$rnb.set, report_dir)

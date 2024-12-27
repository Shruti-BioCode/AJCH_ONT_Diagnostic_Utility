source("Rscripts/utility.R")
reference = "hg19"

# Cohort information
raw_count = fread("data/ont_sample.stats", stringsAsFactors=F, sep="\t")
raw_count = raw_count[Methyl=="y"]

# Get methylation data
methylfiles = c(list.files("data/positive_cohort/methylation/", full.names = T),list.files("data/negative_cohort/methylation/", full.names = T))
samplenames = gsub(".*/|.methyl.cpg.bed.gz","",methylfiles)

# read output of modbam bed files for all the reference
methylationData = generate_methylation_data(methylfiles = methylfiles, samplenames = samplenames, reference = reference)
saveRDS(methylationData,"data/methylationData.rds")
methylationData=readRDS("data/methylationData.rds")



# function to run the analysis
run_complete_methylation_analysis<-function(outdir = "results",methylationData = NULL, cohort_samples=NULL,remove_unstable_episignature=FALSE){

  episign_matrix = create_episign_matrix_from_methylationData(methylationData = methylationData,valuetype = "methylMean")
  sma_methylation_data = methylationData$sma_methylationData[methylationData$sma_methylationData$regionID=="smn1",]
  imd_methylation_data = methylationData$ImD_methylationData
 
  
  # create output directories 
  mndd_outdir = file.path(outdir,"MNDD/")
  sma_outdir = file.path(outdir,"SMA/")
  ImD_outdir = file.path(outdir,"ImD/")
  dir.create(outdir, showWarnings = FALSE)
  dir.create(mndd_outdir, showWarnings = FALSE)
  dir.create(sma_outdir, showWarnings = FALSE)
  dir.create(ImD_outdir, showWarnings = FALSE)
  
  # step 1
  if(!is.null(cohort_samples)){
    episignSample_names = gsub("-","_",cohort_samples)
    episign_matrix = episign_matrix[colnames(episign_matrix) %in% episignSample_names]
    sma_methylation_data = sma_methylation_data[sma_methylation_data$sampleID %in% cohort_samples,]
    imd_methylation_data = imd_methylation_data[imd_methylation_data$sampleID %in% cohort_samples,]
  }
 
  ### epimarker
  
  # print(outdir)
  print("Performing profiling samples for MNDD episignature....")
  processedData = prepare_data(episign_matrix = episign_matrix,reference = reference)
  cl <- makeCluster(4)
  registerDoParallel(cl)
  epimarked_result = episignature_analysis(episign_matrix=processedData$episign_matrix, norm_data = processedData$norm_data, reference = reference, outdir = mndd_outdir, remove_unstable_episignature = remove_unstable_episignature)
  stopCluster(cl)
  fwrite(epimarked_result,file.path(mndd_outdir,"complete_epimarker_results.txt"),row.names = F,sep = "\t")
  
  fwrite(epimarked_result[is_Cluster==TRUE,c("Sample","ProbeGroup","pvalue","is_Cluster","norm_correlation","MI","is_NNcluster","is_hcluster","is_normCorr","is_MI",
                                             "Sample_NN_norm_dist_value","Episign_NN_norm_dist_value","Sample_NN_pca_dist","Episign_NN_pca_dist", "hCluster")],file.path(mndd_outdir,"epimarker_results.txt"),row.names = F,sep = "\t")

  ### SMA analysis
  print("Performing analysis of methylation marks for SMA...")
  sma = methylation_profiling_heatmap(methylationData = sma_methylation_data, readdepth =5, scoreth = 98, pdffile=file.path(sma_outdir,"SMA_heatmap.pdf"))
  
  
  ### Angleman analysis
  print("Performing analysis of methylation marks for Angelman/Padre Willi")
  angelman = methylation_profiling_heatmap(methylationData = imd_methylation_data[regionID=="angelman_paderwili",], readdepth =5, scoreth = 98, pdffile=file.path(ImD_outdir,"Angelma_heatmap.pdf"))
  
  
  return(list("sma"=sma, "angelman"=angelman, "epimarked_result"=epimarked_result))
}

# analysis for negative cohort
negative_cohort_results = run_complete_methylation_analysis(outdir = "results/negative_cohort", methylationData = methylationData, cohort_samples=raw_count$SampleName[raw_count$Cohort=="Negative Cohort Samples"],remove_unstable_episignature = TRUE)

# analysis for positive cohort
positive_cohort_results = run_complete_methylation_analysis(outdir = "results/positive_cohort", methylationData = methylationData, cohort_samples=raw_count$SampleName[raw_count$Cohort!="Negative Cohort Samples"])



b = negative_cohort_results$epimarked_result[is_Cluster==TRUE,c("Sample","ProbeGroup","is_NNcluster","is_hcluster","is_normCorr","is_MI","norm_correlation","is_Cluster","pvalue","adjust_pval_fdr")]
b[order(Sample)]

c = positive_cohort_results$epimarked_result[is_Cluster==TRUE,c("Sample","ProbeGroup","is_NNcluster","is_hcluster","is_normCorr","is_MI","norm_correlation","is_Cluster","pvalue","adjust_pval_fdr")]
c[order(Sample)]



source("Rscripts/utility.R")


reference = "hg19"
outdir = "results/Epimarker_evaluation"
datafiles = "data/methylation_data_illumina_nanopore_samples_controls.xlsx"


### Read methylation data value from the excel file and create matrix for using in the Epimarker
read_excel_allsheets <- function(filename, tibble = FALSE) {
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) data.table(read_excel(filename, sheet = X), check.names = F, stringsAsFactors=F))
  names(x) <- sheets
  x
}


# get episignature specific probes
episignature_probes = extract_episignature_probes(reference = reference)
episignature_probes$Chr = gsub("chr","",episignature_probes$Chr)
episignature_probes$Pos = episignature_probes$`Position (hg38)`

# read methylation data from the file
episign = unique(rbindlist(read_excel_allsheets(datafiles),idcol = F,use.names = T,fill = T))
colnames(episign) = gsub("Control\\.\\.\\.","Control_",colnames(episign))
colnames(episign)[3:37] = paste0("Episign_",colnames(episign)[3:37])
colnames(episign) = gsub("Episign_Control.*","Episign_Control",colnames(episign))
episign = unique(merge.data.frame(unique(episignature_probes[,c("Probes","Chr","Pos")]),episign))
episign$Episign_ADCADN[episign$Episign_ADCADN=="NA"] = 0.000 
episign$Episign_ADCADN = as.numeric(episign$Episign_ADCADN)

# remove samples with non-confirmed clinical testing as reported in the paper Geysens et al,
unconfirmed_methylation=c("EPI_13", "EPI_18", "EPI_20")
episign = episign[!colnames(episign) %in% unconfirmed_methylation]

# generate methylation value matrix for the disease set and the samples
episign_matrix = episign
probename = episign_matrix$Probes
episign_matrix[,c("Chr","Pos","Probes")] = NULL
episign_matrix <- as.data.frame(sapply(episign_matrix, as.numeric))
row.names(episign_matrix) = probename
episign_matrix[is.na(episign_matrix)] = 0
episign_norm_matrix = normalize_data(episign_matrix = episign_matrix)


cl <- makeCluster(4)
registerDoParallel(cl)
# calculate statistics based on PCA distance and correlations
epimarked_result = episignature_analysis(episign_matrix=episign_matrix, reference = reference, outdir = outdir)
stopCluster(cl)

a = epimarked_result[is_Cluster==TRUE,c("Sample","ProbeGroup","is_NNcluster","is_hcluster","is_normCorr","is_MI","norm_correlation","is_Cluster","pvalue","proportion")]
a[order(Sample)]

a = epimarked_result[adjust_pval_bonferroni<=0.05,c("Sample","ProbeGroup","Episign","is_NNcluster","is_hcluster","is_Euc_Corr_pca","is_Euc_MI_pca","norm_correlation","is_Cluster","pvalue","adjust_pval_fdr","adjust_pval_BH","adjust_pval_bonferroni")]
a[order(Sample)]

a = epimarked_result[is_Cluster!=FALSE,c("Sample","ProbeGroup","is_Cluster")]

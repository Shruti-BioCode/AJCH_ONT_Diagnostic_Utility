
setwd(gsub("Rscripts","",dirname(rstudioapi::getSourceEditorContext()$path)))

source("Rscripts/utility.R")
reference = "hg19"

# Cohort information
raw_count = fread("data/ont_sample.stats", stringsAsFactors=F, sep="\t")
raw_count$sampleID = gsub("-","_",raw_count$SampleName)
raw_count = raw_count[Methyl=="y"]

# Get methylation data
methylfiles = list.files("data/sma_methylation/", full.names = T)
samplenames = gsub(".*/|.methyl.cpg.bed.gz","",methylfiles)

# read output of modbam bed files for all the reference
methylationData = generate_methylation_data(methylfiles = methylfiles, samplenames = samplenames, reference = reference)
saveRDS(methylationData,"data/sma_methylationData.rds")
methylationData=readRDS("data/sma_methylationData.rds")

# get methylation data for smn1 / smn2 
smn1_methylation_data = methylationData$sma_methylationData[grepl("smn1",methylationData$sma_methylationData$sampleID) & methylationData$sma_methylationData$regionID=="smn1"]
smn1_methylation_data$sampleID = gsub("_smn1","",smn1_methylation_data$sampleID)

smn2_methylation_data = methylationData$sma_methylationData[grepl("smn2",methylationData$sma_methylationData$sampleID) & methylationData$sma_methylationData$regionID=="smn2"]
smn2_methylation_data$sampleID = gsub("_smn2","",smn2_methylation_data$sampleID)



# function to run the analysis
run_sma_deconvulated_methylation_analysis<-function(sma_outdir = "results", sma_outfile = "sma.pdf" ,sma_methylation_data = NULL, cohort_samples=NULL){

   
   
  # create output directories
   dir.create(sma_outdir, showWarnings = FALSE)

  # step 1
  if(!is.null(cohort_samples)){
    sma_methylation_data = sma_methylation_data[sma_methylation_data$sampleID %in% cohort_samples,]
    missed_samples = cohort_samples[! cohort_samples %in% sma_methylation_data$sampleID]
    missed_data = unique(sma_methylation_data[,-c("sampleID")])
    missed_data = missed_data[,c("modification","score","read_depth","methyl","methylNdistinct","methylNbases","methylMin","methylMax","methylMean","methylMedian","methylSum","methylStdev","methylN0count","methylN100count","methylP0count","methylP100count"):=NA]
    missed_data = unique(missed_data)
    if(length(missed_samples)>0){
       for(sn in missed_samples){
         sn_missed_data = missed_data
         sn_missed_data$sampleID = sn
         sma_methylation_data = rbind(sma_methylation_data,sn_missed_data)
       }
    }
  }

  
  ### SMA analysis
  print("Performing analysis of methylation marks for SMA...")
  sma = methylation_profiling_heatmap(methylationData = sma_methylation_data, readdepth = 1, scoreth = 98, pdffile=file.path(sma_outdir,sma_outfile))
  
  return(sma)
}

# analysis for negative cohort
all_results = run_sma_deconvulated_methylation_analysis(sma_outdir = "results/SMA", methylationData = methylationData, cohort_samples=raw_count$SampleName)


# analysis for negative cohort
negative_cohort_smn1_results = run_sma_deconvulated_methylation_analysis(sma_outdir = "results/negative_cohort/SMA", sma_outfile = "negative_cohort_smn1.pdf", sma_methylation_data = smn1_methylation_data, cohort_samples=raw_count$SampleName[raw_count$Cohort=="Negative Cohort Samples"])
negative_cohort_smn2_results = run_sma_deconvulated_methylation_analysis(sma_outdir = "results/negative_cohort/SMA", sma_outfile = "negative_cohort_smn2.pdf", sma_methylation_data = smn2_methylation_data, cohort_samples=raw_count$SampleName[raw_count$Cohort=="Negative Cohort Samples"])

# analysis for positive cohort
positive_cohort_smn1_results = run_sma_deconvulated_methylation_analysis(sma_outdir = "results/positive_cohort/SMA", sma_outfile = "positive_cohort_smn1.pdf", sma_methylation_data = smn1_methylation_data, cohort_samples=raw_count$SampleName[raw_count$Cohort!="Negative Cohort Samples"])
positive_cohort_smn2_results = run_sma_deconvulated_methylation_analysis(sma_outdir = "results/positive_cohort/SMA", sma_outfile = "positive_cohort_smn2.pdf", sma_methylation_data = smn2_methylation_data, cohort_samples=raw_count$SampleName[raw_count$Cohort!="Negative Cohort Samples"])


# figure for ddpcr
sma = raw_count[!is.na(SMN1_ddpcr),c("SampleName","SMN1_ddpcr","SMN2_ddpcr")]
sma = melt.data.table(sma,id.vars = c("SampleName"),variable.name = "gene",value.name = "CN",variable.factor = F,value.factor = F)
sma$SampleName = factor(sma$SampleName,levels = rev(c("OXN-007","OXN-008","OXN-010","OXN-068","OXN-069","OXN-070","OXN-071","OXN-018","OXN-021","OXN-060","OXN-031","OXN-035","OXN-038","OXN-047","OXN-062","OXN-065")))
p1 = ggplot(sma,aes(x=SampleName,y=CN, fill=gene)) + geom_bar(position = "dodge",stat = "identity") + 
  theme_classic() + coord_flip()+  
  theme(text = element_text(size=12,colour = "black"),axis.text = element_text(size=12,colour = "black"))

pdf("results/sma/SMA_ddpcr.pdf")
p1
dev.off()


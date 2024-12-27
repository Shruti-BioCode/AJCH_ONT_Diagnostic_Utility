

options(width=150)

source("Rscripts/utility.R")

sample_metadata = fread("data/ont_sample.stats", stringsAsFactors=F, sep="\t")
sample_metadata = unique(sample_metadata[,c("SampleName","MRN","All_AGC_IDs","Clinical_indication","Gender","Ethinicity")])
colnames(sample_metadata)[1] = "Indiv"

analyse_splicing_variants<-function(files_dir=NULL, file_pattern="_spliceaiSNP.vcf", outdir="results", outfile_prefix="splice",annotation_data=annotationData, splice_th=0.7, ptrans_th=95){

  # create output directories
  dir.create(outdir, showWarnings = FALSE)
  
  vcffiles = list.files(path=files_dir, pattern=file_pattern,full.names=T)
  samplename = gsub(paste0(".*/|",file_pattern),"",vcffiles)
  spliceAI_results = lapply(vcffiles,function(x, annotation_data,splice_th,ptrans_th){
    
    
    # extracting all splice ai variants from VCF
    oxndata = extract_splicevariants(spliceai_out=x)
    
    # annotating splice ai variants with gene information, splicing score and splice distance
    a = annotate_assign_impact(splice_data=oxndata, annotation_data=annotation_data,spliceAi_th = splice_th)
    
    # print(str(a))
    # get high confidence spliceAI variants of interest
    a$VariantOfInterest = ifelse(a$pExons>ptrans_th & (a$ES=="TRUE" | a$alt3p=="TRUE"  | a$alt5p==TRUE ) & a$has_a_CS=="FALSE", TRUE,FALSE)
    
    # get high confidence spliceAI variants of disease causing 
    a$DiseaseVariantOfInterest = ifelse(a$pExons>ptrans_th & (a$ES=="TRUE" | a$alt3p=="TRUE"  | a$alt5p==TRUE ) & a$has_a_CS=="FALSE" & a$is_disease_causing=="TRUE", TRUE,FALSE)
    
    a$DiseaseVoI_other = ifelse(a$pExons>ptrans_th & a$other==TRUE & a$has_a_CS=="FALSE" & a$is_disease_causing=="TRUE", TRUE,FALSE)
    
    return(a)
    
  }, annotation_data=annotation_data,splice_th=splice_th,ptrans_th=ptrans_th)
  
  names(spliceAI_results) = samplename
  
  print(str(spliceAI_results))
  merged_spliceAI_results = rbindlist(spliceAI_results)
  merged_spliceAI_results = merged_spliceAI_results[,occurence:=uniqueN(Indiv),by=c("POS","chrom","genes")]
  
  outfile = file.path(outdir,paste0(outfile_prefix,"_spliceAI_allresults.txt"))
  fwrite(merged_spliceAI_results, outfile,sep="\t",quote=F,row.names=F)
  
  statsfile = file.path(outdir,paste0(outfile_prefix,"_spliceAI_stats.txt"))
  merged_spliceAI_results$sample_variant_ID = paste0(merged_spliceAI_results$Indiv,"_",merged_spliceAI_results$POS,"_",merged_spliceAI_results$chrom)
  merged_spliceAI_results = merged_spliceAI_results[,nsamples:=uniqueN(Indiv),by=c("POS","chrom")]
  
  merged_spliceAI_results = merged_spliceAI_results[,count_passing_qual:=uniqueN(sample_variant_ID[gt_GQ>=10 & gt_DP>=30 & QUAL>=10 & FILTER=="PASS"]),by=c("Indiv")]
  merged_spliceAI_results = merged_spliceAI_results[,count_passing_qual_disease_gene:=uniqueN(sample_variant_ID[gt_GQ>=10 & gt_DP>=30 & QUAL>=10 & FILTER=="PASS" & is_disease_causing=="TRUE"]),by=c("Indiv")]
  merged_spliceAI_results = merged_spliceAI_results[,count_passing_qual_disease_gene_potential:=uniqueN(sample_variant_ID[gt_GQ>=10 & gt_DP>=30 & QUAL>=10 & FILTER=="PASS" & is_disease_causing=="TRUE" & (DiseaseVariantOfInterest=="TRUE" | DiseaseVoI_other=="TRUE")]),by=c("Indiv")]
  merged_spliceAI_results = merged_spliceAI_results[,count_passing_qual_disease_gene_potential_filtered:=uniqueN(sample_variant_ID[gt_GQ>=10 & gt_DP>=30 & QUAL>=10 & FILTER=="PASS" & is_disease_causing=="TRUE" & (DiseaseVariantOfInterest=="TRUE" | DiseaseVoI_other=="TRUE") & nsamples<=2]),by=c("Indiv")]
  stats_sp = unique(merged_spliceAI_results[,c("Indiv","count_passing_qual","count_passing_qual_disease_gene","count_passing_qual_disease_gene_potential","count_passing_qual_disease_gene_potential_filtered")])
  fwrite(stats_sp, statsfile, sep="\t",quote=F,row.names=F)
  
  return(merged_spliceAI_results)
  
}


#-------------------------------- running for positive cohort --------------------------------
pos_cohort_files_dir="data/positive_cohort/spliceai_variants/"
pos_results = analyse_splicing_variants(files_dir=pos_cohort_files_dir, file_pattern="_spliceaiSNP.vcf", outdir="results/positive_cohort/spliceAI", outfile_prefix="positive_cohort")
pos_results = pos_results[,occurence:=uniqueN(Indiv),by=c("POS","chrom","genes")]
pos_variant_of_interest = unique(pos_results[occurence<=2 & gt_GQ>=10 & gt_DP>=30 & QUAL>=10 & FILTER=="PASS" & 
                                               (DiseaseVariantOfInterest=="TRUE" | DiseaseVoI_other=="TRUE") &  nsamples<=2,c("Indiv","POS","chrom","genes","DS_AG","DS_AL","DS_DG","DS_DL","gt_GT_alleles","gt_GT","REF","ALT","inheritance","combined_phenotype","combined_inheritance","nsamples")])
pos_variant_of_interest = unique(pos_results[occurence<=2 & gt_GQ>=10 & gt_DP>=30 & QUAL>=10 & FILTER=="PASS" & 
                                               (DiseaseVariantOfInterest=="TRUE" | DiseaseVoI_other=="TRUE") &  nsamples<=2 & 
                                               (cDP_AG_distEXstart<3 | cDP_AG_distEXend<3 | cDP_AL_distEXstart<3 | cDP_AL_distEXend<3 | cDP_DG_distEXstart<3 | cDP_DG_distEXend<3 | cDP_DL_distEXstart<3 | cDP_DL_distEXend < 3) ,c("Indiv","POS","chrom","genes","DS_AG","DS_AL","DS_DG","DS_DL","gt_GT_alleles","gt_GT","REF","ALT","inheritance","combined_phenotype","combined_inheritance")])
pos_variant_of_interest = merge.data.table(pos_variant_of_interest,sample_metadata,by=c("Indiv"))
fwrite(x = unique(pos_variant_of_interest), file = "results/positive_cohort/spliceAI/positive_cohort_spliceAI_variants_of_interest.txt",row.names = F,quote = F,sep = "\t")

#-------------------------------- running for negative cohort --------------------------------
neg_cohort_files_dir="data/negative_cohort/spliceai_variants/"
neg_results = analyse_splicing_variants(files_dir=neg_cohort_files_dir, file_pattern="_spliceaiSNP.vcf", outdir="results/negative_cohort/spliceAI", outfile_prefix="neagtive_cohort")
neg_results = neg_results[,occurence:=uniqueN(Indiv),by=c("POS","chrom","genes")]

neg_variant_of_interest = unique(neg_results[occurence<=2 & gt_GQ>=10 & gt_DP>=30 & QUAL>=10 & FILTER=="PASS" & 
                                               (DiseaseVariantOfInterest=="TRUE" | DiseaseVoI_other=="TRUE")  &  nsamples<=2 ,c("Indiv","POS","chrom","genes","DS_AG","DS_AL","DS_DG","DS_DL","gt_GT_alleles","gt_GT","REF","ALT","inheritance","combined_phenotype","combined_inheritance")])

neg_variant_of_interest = unique(neg_results[occurence<=2 & gt_GQ>=10 & gt_DP>=30 & QUAL>=10 & FILTER=="PASS" & 
                                               (DiseaseVariantOfInterest=="TRUE" | DiseaseVoI_other=="TRUE") &  nsamples<=2 & 
                                               (cDP_AG_distEXstart<3 | cDP_AG_distEXend<3 | cDP_AL_distEXstart<3 | cDP_AL_distEXend<3 | cDP_DG_distEXstart<3 | cDP_DG_distEXend<3 | cDP_DL_distEXstart<3 | cDP_DL_distEXend < 3) ,c("Indiv","POS","chrom","genes","DS_AG","DS_AL","DS_DG","DS_DL","gt_GT_alleles","gt_GT","REF","ALT","inheritance","combined_phenotype","combined_inheritance")])
neg_variant_of_interest = merge.data.table(neg_variant_of_interest,sample_metadata,by=c("Indiv"))
fwrite(x = unique(neg_variant_of_interest), file = "results/negative_cohort/spliceAI/neagtive_cohort_spliceAI_variants_of_interest.txt",row.names = F,quote = F,sep = "\t")




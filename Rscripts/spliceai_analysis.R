

options(width=150)

source("Rscripts/utility.R")

sample_metadata = fread("data/ont_sample.stats", stringsAsFactors=F, sep="\t")
sample_metadata = unique(sample_metadata[,c("SampleName","MRN","All_AGC_IDs","Clinical_indication","Gender","Ethinicity")])
colnames(sample_metadata)[1] = "Indiv"



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




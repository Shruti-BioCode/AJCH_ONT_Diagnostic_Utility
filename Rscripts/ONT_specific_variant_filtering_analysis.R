source("Rscripts/utility.R")
library(httr)
library(jsonlite)

patient_info = fread("data/ont_sample.stats", stringsAsFactors=F, sep="\t")
patient_info = patient_info[,c("SampleName","MRN","All_AGC_IDs","Clinical_indication")]
colnames(patient_info)[1] = "sampleName"

missed_exome_variants_file = "data/negative_cohort/missed_exons_variants/missed_exome_variants_varseq_annotated.tsv.gz" 
missed_exome_variants = fread(missed_exome_variants_file, header = T,sep="\t",stringsAsFactors = F)
# a = data.frame(colnames(missed_exome_variants))
colIndex = which(colnames(missed_exome_variants)=="Gene Names")
if(length(colIndex)>1){
  colnames(missed_exome_variants)[colIndex[2]] = "Gene Names clinvar"   
}



missed_exome_variants = add_gene_annotation(annotationData = annotationData, data2annotated = missed_exome_variants)


sample_disease_variant = filter_variants(missed_exome_variants = missed_exome_variants)
sample_disease_variant$variant_category = paste0(sample_disease_variant$CLASS,";",sample_disease_variant$`Effect (Combined)`,";",sample_disease_variant$Classification)
sample_disease_variant$variant_category = gsub("^;|;$","",sample_disease_variant$variant_category)
sample_disease_variant = merge.data.table(sample_disease_variant,patient_info,by = c("sampleName"),all.x = T)
fwrite(sample_disease_variant,"results/negative_cohort/SNV/negative_cohort_missed_exonic_variants.txt",row.names = F,sep = "\t",quote = F)

full_stats = data.table(ftable(sample_disease_variant$variant_category~sample_disease_variant$sampleName))
colnames(full_stats) = c("Sample","Variant_category","Frequency")
full_stats = dcast.data.table(full_stats, Sample ~ Variant_category,value.var = c("Frequency"))

variant_count =  data.frame(table(sample_disease_variant$sampleName))
colnames(variant_count) = c("Sample","Total_variants")

uniq_variant_count =  data.frame(table(sample_disease_variant$sampleName[sample_disease_variant$`# Samples`==1]))
colnames(uniq_variant_count) = c("Sample","unique_variants")

statistics_data = Reduce(function(...) merge.data.table(...,by=c("Sample"), all=T), list(variant_count,uniq_variant_count,full_stats))
fwrite(statistics_data,"results/negative_cohort/SNV/negative_cohort_missed_exonic_variants_stats.txt",row.names = F,sep = "\t",quote = F)


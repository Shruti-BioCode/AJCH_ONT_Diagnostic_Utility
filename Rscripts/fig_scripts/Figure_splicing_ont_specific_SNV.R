library(data.table)
library(ggplot2)

pdffile = "results/splicing_ont_only_exonic_variants.pdf"
outfile = "results/splicing_ont_only_exonic_variants_stats.txt"

snvs = fread("data/ont_sample.stats",header=T,sep="\t",quote=F,stringsAsFactors = F,check.names = F)
snvs = snvs[!is.na(`SNVs (raw call)`)]
snvs$ONT_specific_exonic_variants_perc = snvs$ONT_specific_exonic_variants/snvs$`SNVs (raw call)`
snvs$Splice_Variants_perc = snvs$Splice_Variants/snvs$`SNVs (raw call)`


variant_count = unique(snvs[,c("SampleName","SNVs (raw call)","ONT_specific_exonic_variants","Splice_Variants")])
variant_count = melt.data.table(variant_count,id.vars = c("SampleName","SNVs (raw call)"),variable.name = "Category",value.name = "count_variants",variable.factor = F,value.factor = F)
variant_count$Category = ifelse(grepl("splice",variant_count$Category,ignore.case = T,perl = T),"Splice Variants","Exonic Variants(ONT specific)")


variants_perc = unique(snvs[,c("SampleName","ONT_specific_exonic_variants_perc","Splice_Variants_perc")])
variants_perc = melt.data.table(variants_perc,id.vars = c("SampleName"),variable.name = "Category",value.name = "percentage_variants",variable.factor = F,value.factor = F)
variants_perc$Category = ifelse(grepl("splice",variants_perc$Category,ignore.case = T,perl = T),"Splice Variants","Exonic Variants(ONT specific)")

filtered_variants = unique(snvs[,c("SampleName","ONT_specific_exoninc_variants_passing_filter","Splice_Variants_passing_qual_disease_gene_potential_filtered")])
filtered_variants = melt.data.table(filtered_variants,id.vars = c("SampleName"),variable.name = "Category",value.name = "filtered_variants",variable.factor = F,value.factor = F)
filtered_variants$Category = ifelse(grepl("splice",filtered_variants$Category,ignore.case = T,perl = T),"Splice Variants","Exonic Variants(ONT specific)")

snvs_data4plotting = merge.data.table(variants_perc,filtered_variants,by=c("SampleName","Category"))
ratio1 <-  max(snvs_data4plotting$percentage_variants,na.rm = T) / (max(snvs_data4plotting$filtered_variants,na.rm = T)-0.70)
snvs_data4plotting$ValueScaled = snvs_data4plotting$filtered_variants * ratio1
snvs_data4plotting$SampleName = factor(snvs_data4plotting$SampleName, level=snvs$SampleName,ordered=T)

p1 = ggplot() + geom_bar(data=snvs_data4plotting, aes(SampleName, percentage_variants, fill=Category ), stat="identity") + 
# scale_fill_manual(values=c("Splice Variants"="#b18faf", "Exonic Variants(ONT specific)"="#D4BED3")) +  
scale_fill_manual(values=c("Splice Variants"="#D4BED3", "Exonic Variants(ONT specific)"="#89CDC2")) +  
geom_line(data = snvs_data4plotting, aes(x = SampleName, y = ValueScaled, group=Category, color=Category)) + 
geom_point(aes(x=SampleName, y=ValueScaled, color=Category), data = snvs_data4plotting) +   
scale_colour_manual(values=c("Splice Variants"="black", "Exonic Variants(ONT specific)"="white")) +  
scale_y_continuous(name="SNVs (% of all Variants)",limits = c(0,0.03), labels = scales::label_percent(), sec.axis = sec_axis(~./ratio1, name="Filtered Variants(n)"), expand=c(0,0)) + 
theme(panel.background=element_blank(), axis.line=element_line(colour="black"), axis.text.y=element_text(colour="black"),axis.text.x=element_text(colour="black",angle = 35,hjust = 1),legend.position = "bottom",legend.direction = "horizontal")

pdf(pdffile, family="ArialMT",width = 15)
print(p1)
dev.off()

variant_count = merge.data.table(variant_count, snvs_data4plotting,by=c("SampleName","Category") )
fwrite(variant_count, outfile, row.names = F,sep = "\t",quote = F)
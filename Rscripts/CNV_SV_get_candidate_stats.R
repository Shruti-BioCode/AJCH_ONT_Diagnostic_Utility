
setwd(gsub("Rscripts","",dirname(rstudioapi::getSourceEditorContext()$path)))

source("Rscripts/utility.R")

raw_count = fread("data/ont_sample.stats", stringsAsFactors=F, sep="\t")
raw_count_data = melt.data.table(raw_count[,c("SampleName","QDNAseq","sniffles","cuteSV")],id.vars=c("SampleName"),variable.name="Method",variable.factor=F,value.name="raw_count")

plot_summarize_data <- function(alldataset=NULL, pdffile="rplot.pdf", intervals_val = 100){
  
  print("plotting all summarized data...")
#  alldataset=allsummary
  
  # bringing cuteSV & sniffles to QDNAseq range. 
  alldataset$true_raw_count = alldataset$raw_count
  
  # for positive cohort
  alldataset$raw_count[! alldataset$Method %in% c("QDNAseq","SVs")] = alldataset$true_raw_count[! alldataset$Method %in% c("QDNAseq","SVs")]/intervals_val
  
  # for negative cohort
  # alldataset$raw_count[! alldataset$Method %in% c("QDNAseq","SVs")] = alldataset$true_raw_count[! alldataset$Method %in% c("QDNAseq","SVs")]/100
  
  p1 = ggplot(alldataset[!is.na(true_raw_count)], aes(x=raw_count, y=SampleName, colour=factor(Method), shape=Method)) +
    geom_point(size = 4) + geom_segment(aes(xend = 0, yend = SampleName), colour="grey50", linewidth = 0.2, lineend = "butt") + 
    scale_colour_manual(values=c("cuteSV"="#582766","sniffles"="#EDB081","QDNAseq"="#CE4763"))+
    scale_shape_manual(values=c("cuteSV"=16,"sniffles"=15,"QDNAseq"=17))+
    scale_x_continuous(name="# CNVs/SVs calls", breaks=seq(0,max(alldataset$raw_count,na.rm=T),by=intervals_val)) + 
    theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.y=element_blank(), axis.line.x=element_line(colour="black"),  axis.text=element_text(colour="black"),legend.position="bottom" ,legend.direction="horizontal", legend.key.size = unit(1, 'cm'), legend.text=element_text(size=12)) + ggtitle("Raw calls")
  
  barcolor="#5975A4"
  p2 = ggplot(data=alldataset[Method=="SVs"], aes(x=SampleName, y=passing_disrupting_disease_unique_calls)) + geom_bar(stat="identity", fill=barcolor) + 
    geom_line(data=alldataset[Method=="QDNAseq"], aes(x = SampleName, y = passing_disrupting_disease_unique_calls, group=1), color="white") + 
    geom_point(aes(x=SampleName, y=passing_disrupting_disease_unique_calls), data = alldataset[Method=="QDNAseq"], color="white") + 
    scale_y_continuous(name="# CNVs/SVs HC calls ", breaks=seq(0,max(alldataset$passing_disrupting_disease_unique_calls[alldataset$Method=="SVs"],na.rm=T),by=5)) + 
    coord_flip() + 
    theme(panel.background=element_blank(), axis.line=element_line(colour="black"), axis.text=element_text(colour="black")) + 
    ggtitle("passing_disrupting_disease_unique_calls")
  
  
  p2.1 = ggplot(alldataset[!is.na(true_raw_count)] , aes(x=true_raw_count, y=SampleName, colour=factor(Method), shape=Method)) +
    geom_point(size = 4) + geom_segment(aes(xend = 0, yend = SampleName), colour="grey50", linewidth = 0.2, lineend = "butt") + 
    scale_colour_manual(values=c("cuteSV"="#582766","sniffles"="#EDB081","QDNAseq"="#CE4763"))+
    scale_shape_manual(values=c("cuteSV"=16,"sniffles"=15,"QDNAseq"=17))+
    xlab("# CNVs/SVs calls") + 
    theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.y=element_blank(), axis.line.x=element_line(colour="black"),  axis.text=element_text(colour="black"),legend.position="bottom" ,legend.direction="horizontal", legend.key.size = unit(1, 'cm'), legend.text=element_text(size=12)) + ggtitle("Raw calls")
  p2.1 + scale_x_break(c(350, 19999),scale=3.25)
  
  filtering_step = alldataset[!is.na(true_raw_count),c("SampleName", "Method", "true_raw_count" , "passing_calls", "passing_gene_calls", "passing_disrupting_gene_calls", "passing_disrupting_disease_calls", "passing_disrupting_disease_unique_calls")] 
  filtering_step = melt.data.table(filtering_step, id.vars=c("Method","SampleName"), variable.factor=F, variable.name="Filtering_Steps",value.name="count")
  filtering_step$Filtering_Steps = factor(filtering_step$Filtering_Steps, levels=c( "passing_disrupting_disease_unique_calls", "passing_disrupting_disease_calls", "passing_disrupting_gene_calls", "passing_gene_calls", "passing_calls", "true_raw_count" ), ordered=T)
  filtering_step = filtering_step[,average_count:=round(mean(count),0), by=c("Method","Filtering_Steps")]
  filtering_step$Method = factor(filtering_step$Method, levels=c("QDNAseq","sniffles","cuteSV"),ordered=T)
  a = unique(filtering_step[!is.na(Method),c("Method","Filtering_Steps","average_count")])
  
  px = ggplot(a[Method!="SVs"] , aes(x=average_count, y=Filtering_Steps, fill=Method, shape=Method)) +
    geom_bar(stat = "identity", position=position_dodge2(width=0.5),alpha=0.75, width=0.70) + 
    scale_fill_manual(values=c("cuteSV"="#582766","sniffles"="#EDB081","QDNAseq"="#CE4763"))+
    scale_shape_manual(values=c("cuteSV"=16,"sniffles"=15,"QDNAseq"=17))+
    theme_bw() + theme( panel.border = element_blank(), axis.line.y=element_blank(), axis.line.x=element_line(colour="black"),  axis.text=element_text(colour="black", size = 12),legend.position="bottom" ,legend.direction="horizontal", legend.key.size = unit(1, 'cm'), legend.text=element_text(size=12)) 
  
  px + scale_x_break(c(20, 50),scale=2.25) + scale_x_break(c(400, 15000),scale=3.25) 
  
  
  
  pdf(pdffile, family="ArialMT", width=13, height=9)
  print(p1)
  print(p2)
  print(p2.1 + scale_x_break(c(400, 19999),scale=3.25))	
  # print(px + scale_x_break(c(20, 50),scale=2.25) + scale_x_break(c(400, 15000),scale=3.25))
  print(px + scale_x_break(c(20, 50),scale=2.25) + scale_x_break(c(400, 5000),scale=3.25))	
  dev.off()
  
  
}

analyse_data <-function(qdnaseq_dir = NULL, cuteSV_dir=NULL, snifflesSV_dir=NULL, outdir=NULL, outfile_prefix="SV",
                        pattern_qdnaseq="annotated_", pattern_cutesv="annotated_cutesv_",pattern_sniffles="annotated_",raw_count_data=NULL,intervals_val=100){
  ### read the files ####
  qdnaseq = read_SV_CNV_files(datadir = qdnaseq_dir, pattern=pattern_qdnaseq)
  cuteSV = read_SV_CNV_files(datadir = cuteSV_dir, pattern=pattern_cutesv)
  snifflesSV = read_SV_CNV_files(datadir = snifflesSV_dir, pattern=pattern_sniffles)
  disrupting_SVfile = file.path(outdir,paste0(outfile_prefix,"_SV_disrupting_genes.txt"))
  unique_disrupting_SVfile = file.path(outdir,paste0(outfile_prefix,"_SV_unique_disrupting_genes.txt")) 
  all_CNV_file = file.path(outdir,paste0(outfile_prefix,"_all_CNVs_annotated.txt")) 
  outfile = file.path(outdir,paste0(outfile_prefix,"_summary.txt")) 
  pdffile = file.path(outdir,paste0(outfile_prefix,"_summary.pdf"))
  rdsfile = file.path(outdir,paste0(outfile_prefix,"_data.RDS"))
  
  ### add annotation on unique, disease causing , dirupting genes..
  qdnaseq = add_CNV_SV_filtering_annotation(full_data=qdnaseq, method="qdnaseq")
  cuteSV = add_CNV_SV_filtering_annotation(full_data=cuteSV, method="cuteSV")
  snifflesSV = add_CNV_SV_filtering_annotation(full_data=snifflesSV, method="sniffles")
  
  
  ### SVs disrupting genes
  cols2get = colnames(snifflesSV)[colnames(snifflesSV) %in% colnames(cuteSV)]
  cols2get = cols2get[! cols2get %in% c("REF","ALT")]
  a = cuteSV[is_disrupting_gene=="TRUE", cols2get, with=FALSE]
  b = snifflesSV[is_disrupting_gene=="TRUE", cols2get, with=FALSE]
  disrupting_SVs = data.table(rbind.data.frame(a,b))
  unique_disrupting_SVs = get_unique_CNV_SV_within_sample(full_data = disrupting_SVs)
  
  ### get statistics
  print("getting summary statistics...")
  summary_qdnaseq = get_CNV_SV_stats(full_data = qdnaseq, Method = "QDNAseq")
  summary_cuteSV = get_CNV_SV_stats(full_data = cuteSV, Method = "cuteSV")
  summary_snifflesSV = get_CNV_SV_stats(full_data = snifflesSV, Method = "sniffles")
  summary_SVs = get_CNV_SV_stats(full_data = disrupting_SVs, Method = "SVs")
  usummary_SVs = get_CNV_SV_stats(full_data = unique_disrupting_SVs, Method = "uniqueSVs")
  
  allsummary = data.table(rbind.data.frame(summary_cuteSV,summary_snifflesSV,summary_qdnaseq,summary_SVs,usummary_SVs))
  allsummary = merge.data.table(allsummary,raw_count_data, by=c("SampleName","Method"), all.x=T)
  
  print("writing output files...")
  fwrite(allsummary,file=outfile,quote=F,sep="\t",row.names=F)
  fwrite(qdnaseq,file=all_CNV_file,quote=F,sep="\t",row.names=F)
  fwrite(disrupting_SVs,disrupting_SVfile,quote=F,sep="\t",row.names=F)
  fwrite(unique_disrupting_SVs[passing_disrupting_disease_unique_calls=="TRUE",],unique_disrupting_SVfile,quote=F,sep="\t",row.names=F)
  
  plot_summarize_data(alldataset=allsummary, pdffile=pdffile,intervals_val=intervals_val)
  
  saveRDS(list("qdnaseq" = qdnaseq,
               "cuteSV" = cuteSV,
               "snifflesSV" = snifflesSV,
               "disrupting_SVs" = disrupting_SVs,
               "unique_disrupting_SVs" = unique_disrupting_SVs
               ), file = rdsfile)
}

############# Negative Cohort analysis ###############
qdnaseq_dir = "data/negative_cohort/annotated_qdnaseq/"
cuteSV_dir = "data/negative_cohort/annotated_cuteSV/"
snifflesSV_dir = "data/negative_cohort/annotated_sniffles/"
outdir = "results/negative_cohort/CNV_SV"
outfile_prefix = "negative_cohort"
analyse_data(qdnaseq_dir = qdnaseq_dir, cuteSV_dir=cuteSV_dir, snifflesSV_dir=snifflesSV_dir, outdir=outdir, outfile_prefix=outfile_prefix, raw_count_data=raw_count_data, intervals_val=100)
neg_cohort_result = readRDS("results/negative_cohort/CNV_SV/negative_cohort_data.RDS")

cols2select = c("SampleName","Chromosome","Start","End","variant_unique_id","GT","SV_type","disrupting_no_of_disease_genes_unique_in_variant","disrupting_unique_gene_list","all_disease_phenotype","all_inheritance","all_inheritance_combined","all_classification","all_sv_combined_effect","all_clingen_annotation","passing_disrupting_disease_unique_calls")
svs = unique(neg_cohort_result$disrupting_SVs[disrupting_no_of_disease_genes_unique_in_variant>0 & passing_disrupting_disease_unique_calls=="TRUE",cols2select, with=FALSE])
svs = merge.data.table(svs,unique(raw_count[,c("SampleName","MRN","All_AGC_IDs","Clinical_indication")]),by = c("SampleName"),all.x = T)
cnvs = unique(neg_cohort_result$qdnaseq[disrupting_no_of_disease_genes_unique_in_variant>0 & passing_disrupting_disease_unique_calls=="TRUE",cols2select, with=FALSE])
cnvs = merge.data.table(cnvs,unique(raw_count[,c("SampleName","MRN","All_AGC_IDs","Clinical_indication")]),by = c("SampleName"),all.x = T)

cols2select = c("SampleName","Chromosome","outer_start","outer_end","VariantID","GT","SV_type","disrupting_no_of_disease_genes_unique_in_variant","disrupting_unique_gene_list","all_disease_phenotype","all_inheritance","all_inheritance_combined","all_classification","all_sv_combined_effect","all_clingen_annotation","passing_disrupting_disease_unique_calls","all_variant_unique_id","all_variant_info","all_methods")
usvs = unique(neg_cohort_result$unique_disrupting_SVs[disrupting_no_of_disease_genes_unique_in_variant>0 & passing_disrupting_disease_unique_calls=="TRUE",cols2select, with=FALSE])
usvs = merge.data.table(usvs,unique(raw_count[,c("SampleName","MRN","All_AGC_IDs","Clinical_indication")]),by = c("SampleName"),all.x = T)


fwrite(svs,"results/negative_cohort/CNV_SV/SVs_passing_disrupting_disease_unique_calls.txt",sep = "\t",quote = F,row.names = F)
fwrite(usvs,"results/negative_cohort/CNV_SV/SVs_unique_passing_disrupting_disease_unique_calls.txt",sep = "\t",quote = F,row.names = F)
fwrite(cnvs,"results/negative_cohort/CNV_SV/CNVs_passing_disrupting_disease_unique_calls.txt",sep = "\t",quote = F,row.names = F)



############# Negative Cohort analysis ###############
qdnaseq_dir = "data/positive_cohort/annotated_qdnaseq/"
cuteSV_dir = "data/positive_cohort/annotated_cuteSV/"
snifflesSV_dir = "data/positive_cohort/annotated_sniffles/"
outdir = "results/positive_cohort/CNV_SV"
outfile_prefix = "positive_cohort"
analyse_data(qdnaseq_dir = qdnaseq_dir, cuteSV_dir=cuteSV_dir, snifflesSV_dir=snifflesSV_dir, outdir=outdir, outfile_prefix=outfile_prefix, raw_count_data=raw_count_data,intervals_val=100)

allsummary = fread("results/negative_cohort/CNV_SV/negative_cohort_summary.txt")
plot_summarize_data(alldataset=allsummary, pdffile="results/negative_cohort/CNV_SV/negative_cohort_summary.pdf",intervals_val=intervals_val)


allsummary = fread("results/positive_cohort/CNV_SV/positive_cohort_summary.txt")
plot_summarize_data(alldataset=allsummary, pdffile="results/positive_cohort/CNV_SV/positive_cohort_summary.pdf",intervals_val=intervals_val)

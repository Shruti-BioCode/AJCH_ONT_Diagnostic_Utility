library(GenomicRanges)
library(data.table)
library(ggfortify)
library(viridis)
library(pheatmap)
library(gridExtra)
library(grid)
library(caret)
library(factoextra)
library(ggplot2)
library(doParallel)
library(dplyr)
library(tidyr)
library(valr)
library(ggbreak)
library(WGCNA)
library(stringr)
library(vcfR)
library(readxl)

packageList = c("GenomicRanges","data.table","ggfortify","viridis","pheatmap","gridExtra","grid","caret","factoextra","ggplot2","doParallel","dplyr","WGCNA")


options(digits = 9) 
options(scipen=999)

referenceDB = list()

# bed file for the episign methylation, region for SMA and angleman/padrewili for hg19
referenceDB[["hg19"]] = list(
  "probeBed" = "reference/hg19/methylation/regions.bed",
  # episign values and probes used for profile
  "Episign_Marker" =  "reference/episign/methylation_42_disease_markers.txt",
  "Episign_Values" = "reference/episign/methylation_42_disease.txt",
  
  # sma values consists of methylation modification values for the discovered region of SMN1 from 2 positive and 2 carriers and angleman as negative/control.
 # "SMA" = "reference/hg19/sma.values",
  
  # Currently ImD values consists of methylation modification values for Angleman positive and smna +ve cases as negative/control.
  #"ImD_values" = "reference/hg19/imd.values",
  #"ImD_bed" = "reference/hg19/imd.bed"
) 

# bed file for the episign methylation, region for SMA and angleman/padrewili for hg19
referenceDB[["hg38"]] = list(
  "probeBed" = "reference/hg38/methylation/regions.bed",
  # episign values and probes used for profile
  "Episign_Marker" =  "reference/episign/methylation_42_disease_markers.txt",
  "Episign_Values" = "reference/episign/methylation_42_disease.txt",
  
  # sma values consists of methylation modification values for the discovered region of SMN1 from 2 positive and 2 carriers and angleman as negative/control.
  #"SMA" = "reference/hg38/sma.values",
  
  # Currently ImD values consists of methylation modification values for Angleman positive and smna +ve cases as negative/control.
  #"ImD_values" = "reference/hg38/imd.values",
  #"ImD_bed" = "reference/hg38/imd.bed"
) 

annotation_file_list = list()
annotation_file_list[["omim"]] = "reference/hg19/Gene-based/OMIM/20230317_OMIM-2-annotations.tsv.gz"
annotation_file_list[["genecc"]] = "reference/hg19/Gene-based/GenCC/20230317_GenCC.tsv.gz"
annotation_file_list[["hgmd"]] = "reference/hg19/genomic/HGMD_Pro_2023.2_hg19.vcf.gz"
annotation_file_list[["exonCords"]] = "reference/hg19/genomic/hg19.ncbiRefSeq.primaryChr.slop50.annotation.bed.gz"
annotation_file_list[["inheritance"]] = "reference/hg19/Gene-based/gene_inheritance.txt"
annotation_file_list[["DGV"]] = "reference/hg19/Gene-based/DGV/DGV_gold.gz"
annotation_file_list[["Clingen"]] = "reference/hg19/Gene-based/ClinGen/ClinGen_gene_curation_list_GRCh37.tsv"

# ---------------------------------------- Annotations extraction functions ---------------------------------
prepare_gene_annotation <-function(annotation_file_list=NULL){
  
  print("loading gene annotation...")
  annotation_data = list()
  for(fn in names(annotation_file_list)){
    if(grepl("genecc", fn,ignore.case=T)){
      fd = fread(annotation_file_list[[fn]],stringsAsFactors=F,check.names=F,sep="\t")
      fd$GenCC_moi = NULL  
    }
    if(grepl("omim", fn,ignore.case=T)){
      fd = fread(annotation_file_list[[fn]],stringsAsFactors=F,check.names=F,sep="\t")
      fd$OMIM_inheritance = NULL  
    }
    if(grepl("hgmd", fn,ignore.case=T)){
      fd = fread(annotation_file_list[[fn]],stringsAsFactors=F,check.names=F,sep="\t")
      fd$phenotype = gsub(".*PHEN=\"","",fd$INFO)
      fd$phenotype = gsub("\".*","",fd$phenotype)
      fd$genes = gsub(".*;GENE=","",fd$INFO)
      fd$genes = gsub(";.*","",fd$genes)
      fd = fd[, hgmd_phenotype:=paste(unique(phenotype),sep="",collapse=";"),by="genes"]
      fd = unique(fd[, c("genes","hgmd_phenotype")])
    }
    if(grepl("exonCords", fn,ignore.case=T)){
      fd = fread(annotation_file_list[[fn]],stringsAsFactors=F,check.names=F,sep="\t")
      fd = fd[, Exno:=ifelse(strand=="-",(max(EXnumber)+1)-EXnumber,EXnumber),by=c("genes","transcript")]
      fd = fd[, nIsoforms:=uniqueN(transcript),by=c("genes")]
      fd = fd[, nExons:=uniqueN(transcript),by=c("genes","EXstart","EXend")]
      fd$pExons = fd$nExons/fd$nIsoforms * 100
    }

    if(grepl("inheritance", fn,ignore.case=T)){
      fd = fread(annotation_file_list[[fn]],stringsAsFactors=F,check.names=F,sep="\t")
      fd$inheritance = ifelse(! is.na(fd$HGMD_inheritance),fd$HGMD_inheritance,
                              ifelse(!is.na(fd$OMIM_inheritance),fd$OMIM_inheritance,
                                     ifelse(!is.na(fd$GeneCC_inheritance),fd$GeneCC_inheritance,"Unknown")))
      fd$combined_inheritance = paste0(fd$genes,"(GeneCC-", fd$GeneCC_inheritance,";OMIM-", fd$OMIM_inheritance,";HGMD-", fd$HGMD_inheritance,")")
      fd$combined_inheritance = gsub("GeneCC\\-;|HGMD\\-;|OMIM\\-;|GeneCC\\-NA;|HGMD\\-NA;|OMIM\\-NA;","",fd$combined_inheritance)
    }
    
    if(grepl("DGV", fn,ignore.case=T)){
      fd = unique(fread(annotation_file_list[[fn]],sep="\t",stringsAsFactors=F,check.names=F,header = F,select = c(1:4,8,9,13,14)))
      colnames(fd) = c("chrom","start","end","genes","dg_start","dg_end","annot","overlap")
      fd$DGV_SV = gsub(".*variant_sub_type=|;outer_start.*|\\.","",fd$annot)
      fd$DGV_SV = gsub("Gain","DUP",fd$DGV_SV)
      fd$DGV_SV = gsub("Loss","DEL",fd$DGV_SV)
      fd$DGV_SV_frequency = gsub(".*Frequency=|%;PopulationSummary.*|\\.","",fd$annot)
      fd$DGV_SV_frequency = as.numeric(fd$DGV_SV_frequency)
      fd$DGV_ID = gsub("^ID=|;.*|\\.","",fd$annot)
      fd$DGV_SV[fd$DGV_SV == ""] = NA
      fd$DGV_ID[fd$DGV_ID == ""] = NA
      
      total_bases = unique(fd[ ,c("genes","chrom","start","end")])
      total_bases = total_bases[ ,total_bases:=sum((end-start)) ,c("genes","chrom")]
      total_bases = unique(total_bases[,c("genes","chrom","total_bases")])
      
      fd = merge.data.table(fd,total_bases,by=c("genes","chrom"))
      fd = fd[ ,covered_bases:=sum(overlap), by = c("genes","DGV_ID","chrom")]
      fd$dgv_perc_overlap =  fd$covered_bases/fd$total_bases*100
      fd = unique(fd[,c("genes","DGV_SV","DGV_ID","dgv_perc_overlap","DGV_SV_frequency")])
    }
    
    if(grepl("Clingen", fn,ignore.case=T)){
      fd = fread(annotation_file_list[[fn]],stringsAsFactors=F,check.names=F,sep="\t",skip = 5)
      fd = fd[,c("#Gene Symbol","Haploinsufficiency Description","Haploinsufficiency Score","Triplosensitivity Description","Triplosensitivity Score")]
      colnames(fd) = c("genes","hd","hs","td","ts")
      fd$Clingen_annotation = paste0(fd$genes,"(HS:",fd$hd,",",fd$hs,"|TS:",fd$td,",",fd$ts)
      fd = unique(fd[,c("genes","Clingen_annotation")])
    }
    
    annotation_data[[fn]] = fd
  }
  annotation_data = Reduce(function(...) merge.data.table(...,by=c("genes"), all=T, allow.cartesian=TRUE), annotation_data)
  annotation_data = annotation_data[!is.na(EXstart)]
  annotation_data$combined_phenotype = paste0(annotation_data$genes,"(GeneCC-", annotation_data$GenCC_disease,";HGMD-", annotation_data$hgmd_phenotype,";OMIM-", annotation_data$OMIM_phenotype,";)")
  annotation_data$combined_phenotype = gsub("%2C",",",annotation_data$combined_phenotype)
  annotation_data$combined_phenotype = gsub("%2D","-",annotation_data$combined_phenotype)
  annotation_data$combined_phenotype = gsub("%2E",".",annotation_data$combined_phenotype)
  annotation_data$combined_phenotype = gsub("%3A",":",annotation_data$combined_phenotype)
  annotation_data$combined_phenotype = gsub("%3B",";",annotation_data$combined_phenotype)
  annotation_data$combined_phenotype = gsub("GeneCC\\-;|HGMD\\-;|OMIM\\-;|GeneCC\\-NA;|HGMD\\-NA;|OMIM\\-NA;","",annotation_data$combined_phenotype)
  annotation_data$combined_phenotype = gsub("NA;|;NA","",annotation_data$combined_phenotype)
  
  annotation_data$combined_classification = paste0(annotation_data$genes,"(GeneCC-", annotation_data$GenCC_classification,")")
  annotation_data$combined_classification = gsub("GeneCC\\-;|HGMD\\-;|OMIM\\-;|GeneCC\\-NA;|HGMD\\-NA;|OMIM\\-NA;","",annotation_data$combined_classification)
  annotation_data$combined_classification[is.na(annotation_data$GenCC_classification) | annotation_data$GenCC_classification==""]=NA
  
  annotation_data$is_disease_causing = ifelse(!is.na(annotation_data$GenCC_disease) | !is.na(annotation_data$hgmd_phenotype) | !is.na(annotation_data$OMIM_phenotype), TRUE, FALSE)
  annotation_data$combined_phenotype[annotation_data$is_disease_causing==FALSE]=NA
  
  return(annotation_data)
}

annotationData = prepare_gene_annotation(annotation_file_list=annotation_file_list)


# ---------------------------------------------------------- Methylation Analysis functions --------------------------------------------------

# normalize methylation data 
normalize_data <-function(episign_matrix = NULL){
  
  # convert all na vlaues to zero
  episign_matrix = as.matrix(episign_matrix)
  episign_matrix[is.na(episign_matrix)] = 0
  dat.sample.cs = apply(episign_matrix,2,function(x){ 
    mean_diff = (mean(x)-x) / sd(x)
    normvalue = (mean_diff - min(mean_diff)) / (max(mean_diff) - min(mean_diff))
    return(normvalue)
    }
    )
  # epidata = as.matrix(t(episign_matrix))
  # preproc = preProcess(epidata, method = c("center","scale"))
  # dat.sample.cs = predict(preproc, epidata)
  # dat.sample.cs = dat.sample.cs
  # print(dat.sample.cs[row.names(dat.sample.cs) %in% c("cg11807105","cg27314558"),])
  return(data.frame(dat.sample.cs))
}

# compute principal component
pca_analysis_data <-function(episign_matrix = NULL, to_normalize=TRUE){
  
  dat.sample.cs = episign_matrix
  if(to_normalize){
    dat.sample.cs = normalize_data(episign_matrix = episign_matrix)  
  }
  # normalize the values & generate pca plot
  norm_pca_res <- prcomp(t(dat.sample.cs))

  return(norm_pca_res)
  
}

# compute mutual information
compute_mutual_information <-function(mat=NULL, samplenames = NULL){
  

  mu_dist = mutualInfoAdjacency(mat)
  sample_dist = mu_dist$AdjacencySymmetricUncertainty[! row.names(mu_dist$AdjacencySymmetricUncertainty) %in% samplenames,samplenames]
  dist_zscore = (sample_dist - mean(sample_dist)) / sd(sample_dist)
  dist_pvalue = 2*pnorm(dist_zscore,lower.tail = F)
  dist_pvaluel = pnorm(dist_zscore,lower.tail = T)
  dist_pvalueg = pnorm(dist_zscore,lower.tail = F)
  
  nn_sample_dist = sample_dist[which.max(sample_dist)]
  nn_sample = names(sample_dist[which.max(sample_dist)])
  
  dist_data = data.table("Diseases" = names(sample_dist), "MI"= sample_dist, "MI_zscore"=dist_zscore,"MI_pvalue"=dist_pvalue,"MI_pvaluel"=dist_pvaluel,"MI_pvalueg"=dist_pvalueg,"MI_nn"= nn_sample, "MI_nn_value"=nn_sample_dist, stringsAsFactors = F, keep.rownames = T)
  
  return(dist_data)
}

# compute distance
compute_distance <-function(mat=NULL, samplenames = NULL){
  
  
  pca_dist = as.matrix(dist(mat))
  sample_dist = pca_dist[! row.names(pca_dist) %in% samplenames,samplenames]
  dist_zscore = (sample_dist - mean(sample_dist)) / sd(sample_dist)
  dist_pvalue = 2*pnorm(dist_zscore,lower.tail = F)
  dist_pvaluel = pnorm(dist_zscore,lower.tail = T)
  dist_pvalueg = pnorm(dist_zscore,lower.tail = F)
  
  nn_sample_dist = sample_dist[which.min(sample_dist)]
  nn_sample = names(sample_dist[which.min(sample_dist)])
  
  dist_data = data.table("Diseases" = names(sample_dist), "distance"= sample_dist, "dist_zscore"=dist_zscore,"dist_pvalue"=dist_pvalue,"dist_pvaluel"=dist_pvaluel,"dist_pvalueg"=dist_pvalueg,"dist_nn"= nn_sample, "dist_nn_value"=nn_sample_dist, stringsAsFactors = F, keep.rownames = T)
  
  return(dist_data)
}

# compute pca distance
compute_pca_distance <-function(sample_methyl_matrix = NULL, samplenames=NULL, to_normalize=TRUE){
  
  # calculate similarity based on distance from PCA
  norm_pca_res = pca_analysis_data(episign_matrix = sample_methyl_matrix, to_normalize=to_normalize)
  pca_dist = compute_distance(mat = norm_pca_res$x[,c(1:ncol(norm_pca_res$x))], samplenames = samplenames)
  
  summ = data.frame(summary(norm_pca_res)$importance)[,1:2]*100
  pca_dist$PC1 = summ[2,1]
  pca_dist$PC2 = summ[2,2]
  pca_dist$cumulativePC1 = summ[3,1]
  pca_dist$cumulativePC2 = summ[3,2]
  
  return(pca_dist)
}

# compute correlation
compute_correlation <-function(sample_methyl_matrix = NULL, samplenames=NULL){
  
  sample_corr = cor(sample_methyl_matrix, method = "spearman")
  corr_dist = compute_distance(mat=abs(sample_corr), samplenames = samplenames)
  
  sample_corr = sample_corr[! row.names(sample_corr) %in% samplenames, samplenames]
  corr_zscore = (sample_corr - mean(sample_corr)) / sd(sample_corr)
  corr_pvalue = 2*pnorm(corr_zscore,lower.tail = F) 
  corr_pvaluel = pnorm(corr_zscore,lower.tail = T) 
  corr_pvalueg = pnorm(corr_zscore,lower.tail = F) 
  
  nn_sample_corr = sample_corr[which.max(sample_corr)]
  nn_sample = names(sample_corr[which.max(sample_corr)])
  
  colnames(corr_dist) = paste0("corr_",colnames(corr_dist))
  colnames(corr_dist) = gsub("corr_Diseases","Diseases",colnames(corr_dist))
  
  corr_data = data.table("Diseases" = names(sample_corr), "correlation"= sample_corr, "corr_zscore"=corr_zscore,"corr_pvalue"=corr_pvalue,"corr_pvaluel"=corr_pvaluel,"corr_pvalueg"=corr_pvalueg,"corr_nn"= nn_sample, "corr_nn_value"=nn_sample_corr, stringsAsFactors = F, keep.rownames = T)
  corr_data = merge.data.table(corr_data,corr_dist,by=c("Diseases"))
  return(corr_data)
}

# extract probe information
extract_episignature_probes <-function(reference = "hg19"){
  
  # reference datasets
  referencesData =  referenceDB[[ reference ]]
  
  # reading the methylation marker file of the paper
  markerdata = fread(referencesData$Episign_Marker, header=T, sep="\t", stringsAsFactors=F, check.names=F)
  cntrl_markerdata = unique(markerdata[,c("Probes","Chr","Position (hg19)","Position (hg38)")])
  cntrl_markerdata$Disease = "Episign_Control"
  cntrl_markerdata$is_marker = TRUE

  markerdata = melt.data.table(markerdata, id.vars=c("Probes","Chr","Position (hg19)","Position (hg38)"), variable.factor=F, variable.name="Disease", value.name="is_marker")
  markerdata = markerdata[is_marker==TRUE]
  markerdata$Disease = paste0("Episign_",markerdata$Disease)
  markerdata = data.table(rbind(markerdata,cntrl_markerdata))

  # reading probe bed file  
  marker_probeBed = fread(referencesData$probeBed, header=T, sep="\t", stringsAsFactors=F, check.names=F)
  colnames(marker_probeBed) = gsub("probeID","Probes",colnames(marker_probeBed))
  
  # merging episign specific probes with probe coordinates(bed files)
  markerdata = merge.data.table(markerdata,marker_probeBed,by = c("Probes"))
  
  # get all episign related values from the paper
  epivalues = fread(referencesData$Episign_Values,header=T, sep="\t", stringsAsFactors=F)
  colnames(epivalues)[-c(1)] = paste0("Episign_",colnames(epivalues)[-c(1)])
  matrix_rownames = epivalues$Probe
  epivalues = as.matrix(epivalues[,-c("Probe")])
  row.names(epivalues) = matrix_rownames
  
  # standardize & normalize methylation values 
  norm_data = normalize_data(episign_matrix = epivalues)
  
  # get optimal clusters for each Episign
  cluster_info=NULL
  for(d in unique(markerdata$Disease)){
    
    # get normalized methylation matrix of probes specific to Episign
    dnormData = t(norm_data[row.names(norm_data) %in% markerdata$Probes[markerdata$Disease==d & markerdata$is_marker==TRUE],grepl("Episign",colnames(norm_data))])
    
    # Perform hierarchical clustering
    res <- pheatmap(dnormData,silent = T)
    for(nclust in seq(2,35,1)){
      clust <- data.table("Diseases"=names(cutree(res$tree_row, k = nclust)),"Cluster"=as.vector(cutree(res$tree_row, k = nclust)),check.names = F)
      n = nrow(clust[Cluster==Cluster[Diseases==d]])
      # If the episignature forms unique cluster within all episign markers
      if(n==1){
        cluster_info = rbind(cluster_info,data.frame("Disease"=d,"optimal_cluster"= nclust))
        break
      }
    }
  }
  
  # merging episign specific probes with optimal number of clusters
  markerdata = merge.data.table(markerdata,cluster_info,by = c("Disease"),all.x = T)
  
  return(markerdata)
}

# Function to read modbam input files and intersect with the disease specific probes and region 
extract_methylationtags <-function(modbamBed = NULL, samplename=NULL, reference = "hg19"){
  
  referencesData =  referenceDB[[ reference ]]
  
  tryCatch(
    {
      probeBed = fread( referencesData$probeBed, header=T,sep="\t",stringsAsFactors = F)
      probeBed = probeBed[! grepl("_|un",chrom,ignore.case = T)]
      probeBed$chrom = gsub("chr","",probeBed$chrom)
      
      methylationBed = fread(modbamBed,header = F,sep="\t",stringsAsFactors = F)
      methylationBed = methylationBed[! grepl("_|un",V1,ignore.case = T)]
      methylationBed$V1 = gsub("chr","",methylationBed$V1)
      methylationBed = methylationBed[,c(1,2,3,11,6,5,10)]
      colnames(methylationBed) = c("chrom","start","end","mod_perc","strand","score","read_depth")
      
      # methData = nrow(methylationBed[! (is.na(mod_perc) & mod_perc==0)])
      methData = nrow(methylationBed[! (is.na(mod_perc))])
      
      episign_methylationData = NULL
      sma_methylationData = NULL
      ImD_methylationData = NULL
      
      
      if(methData > 0){
        epimarks = GRanges(seqnames = probeBed$chrom,
                           ranges = IRanges(start = probeBed$start, end = probeBed$end),
                           region_chrom = probeBed$chrom,
                           region_start = probeBed$start,
                           region_end = probeBed$end,
                           regionID  = probeBed$probeID,
                           methyl_signature = probeBed$dataset
        )
        
        methylations = GRanges(seqnames = methylationBed$chrom,
                               ranges = IRanges(start = methylationBed$start, end = methylationBed$end),
                               #  position = ifelse(methylationBed$strand == "+", methylationBed$start, methylationBed$end),
                               position =  methylationBed$start,
                               modification = methylationBed$mod_perc,
                               score = methylationBed$score,
                               read_depth = methylationBed$read_depth
        )
        
        #get overlaps		
        overlaps <- findOverlaps(query = epimarks, subject = methylations)
        methylationData = data.table(data.frame(cbind(mcols(epimarks[queryHits(overlaps)]),mcols(methylations[subjectHits(overlaps)]))))
        colnames(methylationData) = c("chrom","start","end","regionID","methyl_signature","position","modification","score","read_depth")
        methylationData$region_coordinate = paste0(methylationData$chrom,"_", methylationData$start,"_", methylationData$end,"_",methylationData$regionID)
        methylationData = methylationData[, methyl := mean(modification,na.rm=T), by=c("position","chrom")]
        methylationData = methylationData[, methylNdistinct:=uniqueN(methyl), by =c("chrom","position")]
        methylationData = methylationData[, methylNbases:=uniqueN(position), by=c("region_coordinate")]
        
        methylationData = unique(methylationData[methylNdistinct!=0])
        
        methylationData = methylationData[, methylMin:=min(methyl,na.rm = T), by=c("region_coordinate")]
        methylationData = methylationData[, methylMax:=max(methyl,na.rm = T), by=c("region_coordinate")]
        methylationData = methylationData[, methylMean:=mean(methyl,na.rm = T), by=c("region_coordinate")]
        methylationData = methylationData[, methylMedian:=median(methyl,na.rm = T), by=c("region_coordinate")]
        methylationData = methylationData[, methylSum:=sum(methyl,na.rm = T), by=c("region_coordinate")]
        methylationData = methylationData[, methylStdev:=sd(methyl,na.rm = T), by=c("region_coordinate")]
        methylationData = methylationData[, methylN0count:=length(!is.na(methyl) & methyl==0), by=c("region_coordinate")]
        methylationData = methylationData[, methylN100count:=length(!is.na(methyl) & methyl==100), by=c("region_coordinate")]
        methylationData = methylationData[, methylP0count:=methylN0count/length(methyl)*100, by=c("region_coordinate")]
        methylationData = methylationData[, methylP100count:=methylN100count/length(methyl)*100, by=c("region_coordinate")]
        methylationData$sampleID = samplename 
        episign_methylationData = unique(methylationData[methyl_signature=="Episign",-c("position","modification","methyl","score","read_depth")])
        sma_methylationData = unique(methylationData[methyl_signature=="SMA",])
        ImD_methylationData = unique(methylationData[methyl_signature=="ImD",])  
        
        gc()
        return(list("episign_methylationData"=episign_methylationData,"sma_methylationData"=sma_methylationData,"ImD_methylationData"=ImD_methylationData))
      }else{
        print(paste0("Sample ",samplename," has no methylation value. Please check the file"))
      }
    },
    error = function(e) {
      message(paste("File is empty: ", modbamBed))
      message(conditionMessage(e))
    }
  )
   
}

# Function to generate methylation extracted information for all the samples
generate_methylation_data <- function(methylfiles = NULL, samplenames = NULL, reference = "hg19"){
 
  ### check if the file & sample input is correct
  nfiles = length(unique(methylfiles))
  nsamples = length(unique(samplenames))
  if(!is.null(samplenames)){
    dups = length(samplenames[duplicated(samplenames)])
    if(dups > 0){
      print("Same sample name cannot be used. 
            Exiting the program.
            Please enter unique name for each file")
    }
  if(nsamples!=nfiles){
    print("Number of samples differ than the number of files. 
            Exiting the program.
            Please ensure the number of samples is same as the number of files")
  }  
  }else{
    samplenames = gusb(".methyl.cpg.bed.gz","",basename(methylfiles))
    dups = length(samplenames[duplicated(samplenames)])
    if(dups > 0){
        samplenames = paste0(samplenames,"-",c(1:length(samplenames)))
    }
  }
  
  # declare variables
  episign_methylationData = list()
  sma_methylationData = list()
  ImD_methylationData = list()
  
  for(i in 1:length(methylfiles)){
    print(paste0("Processing sample ", samplenames[i]," ...."))
    methylationData = extract_methylationtags(modbamBed = methylfiles[i], samplename=samplenames[i], reference = reference)
    episign_methylationData[[samplenames[i]]] = methylationData$episign_methylationData
    sma_methylationData[[samplenames[i]]] = methylationData$sma_methylationData
    ImD_methylationData[[samplenames[i]]] = methylationData$ImD_methylationData
    rm(methylationData)
    gc()
  } 
  episign_methylationData = rbindlist(episign_methylationData,use.names = T,fill = T)
  sma_methylationData = rbindlist(sma_methylationData,use.names = T,fill = T)
  ImD_methylationData = rbindlist(ImD_methylationData,use.names = T,fill = T)
  
  # step 1 find ambiguous probes (with same name with different coordinates)
  dup_probes = unique(episign_methylationData[,c("regionID","region_coordinate")])
  dup_probes = dup_probes$regionID[duplicated(dup_probes$regionID)]
  
  episign_methylationData = unique(episign_methylationData[!regionID %in% dup_probes,])
  
  return(list("episign_methylationData"=episign_methylationData,"sma_methylationData"=sma_methylationData,"ImD_methylationData"=ImD_methylationData))
}

# Function to create methylation matrix for all the episignature diseases and the user samples 
generate_episign_matrix <- function(methylationData=NULL, valuetype = "mean" , reference = "hg19", probes_of_interest=NULL){

  epivalues = NULL
  
  if(!is.null(reference)){
    referencesData =  referenceDB[[ reference ]]
    
    # reading the methylation value file of the paper
    epivalues = fread(referencesData$Episign_Values,header=T, sep="\t", stringsAsFactors=F)
    colnames(epivalues)[-c(1)] = paste0("Episign_",colnames(epivalues)[-c(1)])  
    probes_of_interest = epivalues$Probe
  }
  
  
  
  # methylation modification value from the samples based on user choice of mean/min/max/median/sum of the modifications across the probes
  cols2select = c("regionID","methylMean","sampleID")
  
  if(valuetype=="median"){
    cols2select = c("regionID","methylMedian","sampleID")
  }
  if(valuetype=="min"){
    cols2select = c("regionID","methylMin","sampleID")
  }
  if(valuetype=="max"){
    cols2select = c("regionID","methylMax","sampleID")
  }
  if(valuetype=="sum"){
    cols2select = c("regionID","methylSum","sampleID")
  }
  
  methylvalues = unique(methylationData$episign_methylationData[methylationData$episign_methylationData$regionID %in% probes_of_interest,cols2select, with=FALSE])
  colnames(methylvalues) = c("Probe","Modification","Sample")
  methylvalues$Modification = methylvalues$Modification/100
  
  
  if(length(unique(methylvalues$Sample))>1){
    methylvalues = unique(dcast.data.table(methylvalues, Probe ~ Sample, fill = 0, value.var = "Modification"))
  }else{
    new_colnames = c("Probe",methylvalues$Sample[1])
    methylvalues = unique(methylvalues[,c("Probe","Modification")])
    colnames(methylvalues) = new_colnames
  }
  
  # merge with the episignature data
  if(!is.null(reference)){
    methylvalues = merge.data.table(epivalues,methylvalues, by=c("Probe"),all.x = T)
  }
  matrix_rownames = methylvalues$Probe
  methylvalues = as.matrix(methylvalues[,-c("Probe")])
  row.names(methylvalues) = matrix_rownames
  return(methylvalues)
}

# Function to create the matrix from the methylation data
create_episign_matrix_from_methylationData <- function(methylationData=NULL, valuetype = "methylMean"){
  episign_matrix = unique(methylationData$episign_methylationData[,c("regionID","sampleID",valuetype),with=FALSE])
  episign_matrix = data.frame(dcast.data.table(episign_matrix,regionID~sampleID,value.var = valuetype),check.names = F)
  row.names(episign_matrix) = episign_matrix$regionID
  episign_matrix$regionID = NULL
  colnames(episign_matrix)=gsub("-","_",colnames(episign_matrix))
  return(episign_matrix)
}
  
# Function to create methylation matrix for all the episignature diseases and the user samples 
create_episign_matrix <- function(methylvalues=NULL, reference = "hg19", combine_normalize_data=FALSE){
  
  if(is.null(reference)){
     print("Reference was not provided")
     exit
  }
  
  
  # reading the methylation value file of the paper
  referencesData =  referenceDB[[ reference ]]
  epivalues = fread(referencesData$Episign_Values,header=T, sep="\t", stringsAsFactors=F)
  colnames(epivalues)[-c(1)] = paste0("Episign_",colnames(epivalues)[-c(1)])  
  probes_of_interest = unique(epivalues$Probe)
  
  # converting from percentage to decimal
  if( any(methylvalues>1)){
    methylvalues = methylvalues / 100  
  }
  
  # if merging normalized data
  if(combine_normalize_data==TRUE){
    
    # normalize the episign values
    epivalues = data.frame(epivalues)
    row.names(epivalues) = epivalues$Probe
    epivalues$Probe = NULL
    epivalues = normalize_data(episign_matrix = epivalues)
    epivalues$Probe = row.names(epivalues)
    
    # normalize the ONT methylation probes
    methylvalues = normalize_data(episign_matrix = methylvalues)
    
  }
  
  # merge with the episignature data
  methylvalues$Probe = row.names(methylvalues)
  methylvalues = methylvalues[methylvalues$Probe %in% probes_of_interest,]
  methylvalues = data.frame(merge.data.table(epivalues,methylvalues, by=c("Probe"),all.x = T))
  row.names(methylvalues) = methylvalues$Probe
  methylvalues$Probe = NULL
  methylvalues = methylvalues[!is.na(row.names(methylvalues)),]
  methylvalues[is.na(methylvalues)] = 0
  return(methylvalues)
}

# Function to create pca plot
episignature_plot_pca <- function(episign_matrix=NULL, title="PCA plot", return = FALSE, to_normalize=TRUE, colourcode=c("Control"="#4CCE76","Episign"="grey35","FBXO22"="#DF2F62","OXN"="#722FDF")){
  
  norm_pca_res = pca_analysis_data(episign_matrix = episign_matrix, to_normalize = to_normalize)
  
  
  db = norm_pca_res$x
  if(ncol(norm_pca_res$x)>1){
    db = data.frame(norm_pca_res$x[,1:2])  
  }
  
  db$Dataset = row.names(db)
  
  if(is.null(colourcode)){
    colourcode=c("Control"="#4CCE76","Episign"="grey35","FBXO22"="#DF2F62","OXN"="#722FDF","ONT"="#2F9FDF")
  }
  for(i in 1:length(colourcode)){
    db$Dataset[grepl(names(colourcode)[i],db$Dataset)] = names(colourcode)[i]
  }
  # fill in remaining colourgroups
  namesOthers = db$Dataset[!grepl(paste0(names(colourcode),collapse = "|"),db$Dataset)]
  otherGroup = rep(x="#6ABB23",length(namesOthers))
  names(otherGroup) = namesOthers
  colourcode = c(colourcode,otherGroup)
  
  alphacode = rep(1,length(colourcode))
  names(alphacode) = unique(db$Dataset)
  alphacode[names(alphacode)=="Episign"] = 0.3
  
  pcData = data.frame(norm_pca_res$x)
  pcData$Disorders = row.names(pcData)
  pcData = melt.data.table(data.table(pcData),id.vars = c("Disorders"),variable.name = "PC",value.name = "Methylation",variable.factor = F,value.factor = F)
  pcData$PC = factor(pcData$PC, levels=unique(pcData$PC), ordered=T)
  pcData$Disorders = factor(pcData$Disorders, levels=rev(unique(pcData$Disorders)), ordered=T)
  
  columnorder = colnames(norm_pca_res$x)

   p1 = autoplot(norm_pca_res, x=1, y=2, label = TRUE, label.size=5, label.repel=T, max.overlaps=100, size = 9, data=db, color='Dataset',alpha = 'Dataset') +
    geom_vline(xintercept = 0, colour="grey") + geom_hline(yintercept = 0,colour="grey") + 
   # scale_y_continuous(expand = expansion(mult = c(.05, .05))) + scale_x_continuous(expand = expansion(mult = c(.05, .05))) + 
    scale_colour_manual(values=colourcode) + scale_alpha_manual(values=alphacode) + ggtitle(title) +
    theme_classic() + theme(axis.text = element_text(colour = "black",size = 17),axis.title = element_text(colour = "black",size = 18),aspect.ratio = 0.9, plot.title = element_text(hjust = 0.5, size = 21))
   p2 = ggplot(pcData,aes(y=Disorders,x=PC,fill=Methylation)) + geom_tile(stat = "identity")+
     scale_fill_gradient2() + ggtitle(title) + 
     theme_classic() + theme(axis.text = element_text(colour = "black",size = 12),aspect.ratio = 0.9, legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5, size = 21))
   
  
  if(return == TRUE){
    return(list(p1,p2))   
  }else{
    print(p1)
    print(p2)
  }
  
  
}

# Function to create heatmap
episignature_plot_heatmap<-function(episign_matrix=NULL, title="heatmap", return = FALSE, clustering_method = "ward.D2", scale_by="none"){

 # convert all na vlaues to zero
 episign_matrix[is.na(episign_matrix)] = 0
 if(return == TRUE){
   pmap = pheatmap(t(episign_matrix),clustering_distance_cols = "euclidean",clustering_distance_rows = "euclidean", clustering_method = clustering_method, show_rownames = T, show_colnames = F,fontsize = 10, main = title, treeheight_row = 40, treeheight_col = 40, silent = T, scale = scale_by)
#   grid.arrange(pmap$gtable, vp=viewport(width=0.9, height=0.75)) 
   return(pmap)   
 }else{
   pmap = pheatmap(t(episign_matrix),clustering_distance_cols = "euclidean",clustering_distance_rows = "euclidean", clustering_method = clustering_method, show_rownames = T, show_colnames = F,fontsize = 10, main = title, treeheight_row = 40, treeheight_col = 40, scale = scale_by, silent = T)
   grid.arrange(pmap$gtable, vp=viewport(width=0.9, height=0.75)) 
 }
 
}

# calculate distance & correlation
measure_metric <-function(sample_methyl_norm=NULL, samplenames=NULL,to_normalize = FALSE){
  
  # calculate similarity based on distance from PCA
  pca_dist_data = compute_pca_distance(sample_methyl_matrix = sample_methyl_norm, samplenames=samplenames,to_normalize = to_normalize)
  colnames(pca_dist_data) = paste0("pca_",colnames(pca_dist_data))
  colnames(pca_dist_data) = gsub("pca_Diseases","Diseases",colnames(pca_dist_data))
  
  # calculate similarity based on distance from PCA
  norm_dist_data = compute_distance(mat = t(sample_methyl_norm), samplenames=samplenames)
  colnames(norm_dist_data) = paste0("norm_",colnames(norm_dist_data))
  colnames(norm_dist_data) = gsub("norm_Diseases","Diseases",colnames(norm_dist_data))
  
  # calculate similarity based on correlation for normalized & standardized methylation values
  norm_corr_data = compute_correlation(sample_methyl_matrix = sample_methyl_norm, samplenames=samplenames)
  colnames(norm_corr_data) = paste0("norm_",colnames(norm_corr_data))
  colnames(norm_corr_data) = gsub("norm_Diseases","Diseases",colnames(norm_corr_data))
  
  # mutual information
  mi_data = compute_mutual_information(mat = sample_methyl_norm,samplenames = samplenames)
  
  # merge data 
  sample_dataset = Reduce(function(...) merge(..., all=T), list("norm_corr_data"=norm_corr_data,"norm_dist_data"=norm_dist_data,"pca_dist_data"=pca_dist_data,"mi_data" = mi_data ))
  return(sample_dataset)
  
}

# perform hierarchial_clustering_analysis to get the episign sample is clustering with
hierarchial_clustering_analysis <-function(sample_methyl_norm=NULL, samplename = NULL, episign=NULL){
  
  max_clusters = ncol(sample_methyl_norm)-1
  sclust = NA
  res <- pheatmap(t(sample_methyl_norm), clustering_distance_cols = "euclidean",clustering_distance_rows = "euclidean", clustering_method = "ward.D2",silent = T)
  
  # assess if at any level of the hierarchial tree, the sample and the episign are in the same branch and cluster together
  for(nclust in seq(2,max_clusters,1)){
    clust <- data.table("Groups"=names(cutree(res$tree_row, k = nclust)),"ClusterNo"=as.vector(cutree(res$tree_row, k = nclust)),check.names = F)
    epiCluster = clust$ClusterNo[clust$Groups==episign]
    sampCluster = clust$ClusterNo[clust$Groups==samplename]
    clust = clust[ClusterNo==ClusterNo[Groups==episign]]
    if((epiCluster == sampCluster) & nrow(clust)==2){
      sclust = paste0(clust$Groups,"(",clust$ClusterNo,")",collapse = ";")
      break
    }
  }
  
  return(sclust)
  
}


# Function to compute all measure values
compute_all_measures <-function(sample_methyl_norm=NULL, samplenames=NULL, probeset="All Episign probes", episign = NULL){
  
  # Hierarchial clustering: each disease specific episign should form a cluster separated from the others for its corresponding set of Episigns 
  # hence, 2 clusters would indicate that the Episign is clustered separtely from other Episigns and the sample in question if is similar to the Episign will cluster with it.
   
  sample_methyl_norm[is.na(sample_methyl_norm)] = 0
  
  if(! is.null(episign)){
    sample_dataset = measure_metric(sample_methyl_norm=sample_methyl_norm, samplenames=episign)
  }else{
    sample_dataset = measure_metric(sample_methyl_norm=sample_methyl_norm, samplenames=samplenames)
  }
  sample_dataset$Sample = samplenames
  sample_dataset$Probeset = probeset
  
  return(sample_dataset)
}

# merge data set for sample and episign to get nearest neighbour for
merge_datasets <-function(d_dataset=NULL,s_dataset=NULL, samplenames=NULL, episignname=NULL){
  
  d_dataset = d_dataset[d_dataset$Diseases==samplenames,]
  d_names = colnames(d_dataset)
  d_names[grepl("_nn",d_names)] = paste0("Episign_NN_",d_names[grepl("_nn",d_names)])
  d_names = gsub("_nn","",d_names,ignore.case = F,fixed = T)
  colnames(d_dataset) = d_names
  d_dataset$Diseases = NULL
  
  s_dataset = s_dataset[s_dataset$Diseases== episignname,c("Probeset","pca_dist_nn","pca_dist_nn_value","norm_dist_nn","norm_dist_nn_value","norm_corr_nn","norm_corr_nn_value","norm_corr_dist_nn","norm_corr_dist_nn_value","MI_nn","MI_nn_value")]
  colnames(s_dataset)[2:ncol(s_dataset)] = paste0("Sample_NN_",colnames(s_dataset)[2:ncol(s_dataset)])
  colnames(s_dataset) = gsub("_nn","",colnames(s_dataset),ignore.case = F,fixed = T)
  
  db = merge.data.frame(d_dataset,s_dataset,by=c("Probeset"))
  columnOrder = colnames(db)
  db$Episign = episignname
  columnOrder=c("Probeset","Sample","Episign",columnOrder[grepl("Episign_",columnOrder)],columnOrder[grepl("Sample_",columnOrder)],columnOrder[!grepl("Episign|Sample|Probeset",columnOrder)])
  db = db[,columnOrder]
  
  return(db)
}

# Function to plot top Episigns with high correlation
plot_all_results<-function(sample_methyl_matrix=NULL,norm_data=NULL,markerdata=NULL, outfile=NULL, colourcode=NULL, episign_matrix=NULL){
  
  episign_disease =  unique(markerdata$Disease)
  sample_methyl_matrix = as.matrix(sample_methyl_matrix)
  norm_data = as.matrix(norm_data)
  sample_methyl_matrix[is.na(sample_methyl_matrix)] = 0
  ### generating output ###
  pdf(paste0(outfile,".pdf"),width = 15,height = 12,onefile=T)
  episignature_plot_pca(episign_matrix = norm_data,title = paste0("PCA on standardized & normalized methylation values","(#",nrow(norm_data),")"),to_normalize = FALSE,colourcode = colourcode)
  episignature_plot_heatmap(episign_matrix=sample_methyl_matrix, title=paste0("Heatmap on raw methylation values","(#",nrow(norm_data),")"), return = FALSE)
  episignature_plot_heatmap(episign_matrix=norm_data, title=paste0("Heatmap on standardized & normalized methylation values","(#",nrow(norm_data),")"), return = FALSE)
  for( episign in episign_disease){
    norm_disorder_matrix = norm_data[row.names(norm_data) %in% markerdata$Probes[markerdata$Disease==episign],]
    disorder_matrix = sample_methyl_matrix[row.names(sample_methyl_matrix) %in% markerdata$Probes[markerdata$Disease==episign],]
    
    cols2select = colnames(episign_matrix)
    cols2select = unique(c( cols2select[! grepl("Episign",cols2select)], cols2select[cols2select %in% c(episign,"Episign_Control")]))
    cohort_matrix = episign_matrix[row.names(episign_matrix) %in% markerdata$Probes[markerdata$Disease==episign], colnames(episign_matrix) %in% cols2select]
 
    if(nrow(norm_disorder_matrix)>0){
      title = paste0(episign," Episign specific plot on standardized & normalized methylation values","(#",nrow(norm_disorder_matrix),")")
      episignature_plot_pca(episign_matrix=norm_disorder_matrix, title = title, return = FALSE, to_normalize = FALSE,colourcode = colourcode)  
      episignature_plot_heatmap(episign_matrix=norm_disorder_matrix, title=title, return = FALSE)
      title = paste0(episign," Episign specific plot on raw methylation values","(#",nrow(disorder_matrix),")")
      episignature_plot_heatmap(episign_matrix=disorder_matrix, title=title, return = FALSE)
      episignature_plot_heatmap(episign_matrix=cohort_matrix, title=title, return = FALSE, scale_by = "column")
      
    }else{
      print(paste0("No probes for: ", episign))
    }
  }
  
  dev.off()
}

# Function to analyse the methylation data for the episignature datasets
# outputs pdf file for each user sample and a file with the pca distance, correlation and pvalue

prepare_data <-function(episign_matrix = NULL, reference = "hg19"){
  
  # reading the methylation marker file of the paper
  markerdata = extract_episignature_probes(reference = reference)
  
  # convert all na vlaues to zero
  episign_matrix[is.na(episign_matrix)] = 0
  
  # standardize and normalize the methylation matrix with MNDD methylation values
  print("Processing data for standardization & normalization...")
  norm_data = create_episign_matrix(methylvalues=episign_matrix, reference = reference, combine_normalize_data=TRUE)
  print("Processing completed")
  
  # raw methylation values matrix with MNDD 
  episign_matrix = create_episign_matrix(methylvalues=episign_matrix, reference = reference, combine_normalize_data=FALSE)
  
  return(list("episign_matrix"=episign_matrix,"norm_data"=norm_data))
}


episignature_analysis<-function(episign_matrix=NULL, norm_data = NULL , reference = "hg19", outdir = NULL, remove_unstable_episignature=FALSE){
  
  # create output directory
  dir.create(outdir, showWarnings = FALSE)
  
  if(is.null(norm_data)){
    print("no normalized data has been provided...")
    print("standardization & normalization will be performed on the raw data")
    norm_data = normalize_data(episign_matrix = episign_matrix) # normaliztion on the 3,643 probes
  }
  
  episign_matrix = as.matrix(episign_matrix)
  norm_data = as.matrix(norm_data)
  
  episign_matrix[is.na(episign_matrix)] = 0
  
  # get the MNDD probes - 3,643 
  markerdata = extract_episignature_probes(reference = reference)

  # get all user samples
  cohort_samples = colnames(episign_matrix)[!colnames(episign_matrix) %in% markerdata$Disease]

  if(remove_unstable_episignature==TRUE){
    # removing low reproducibility episignatures
    markerdata = markerdata[!Disease %in% c("Episign_CdLS","Episign_AUTS18","Episign_RSTS")]
    
  }
  
  # overiew analysis with all samples
  ov_outfile = file.path(outdir,"Epimarker_Overview.pdf")
  pmapRaw = episignature_plot_heatmap(episign_matrix = episign_matrix, title="All Episign Specific Probe", return = TRUE, clustering_method = "ward.D")
  pmapnorm = episignature_plot_heatmap(episign_matrix = norm_data, title="All Episign Specific Probe on standardized & normalized", return = TRUE, clustering_method = "ward.D")
  
  pdf(ov_outfile,".pdf",width = 15,height = 12)
  episignature_plot_pca(episign_matrix = episign_matrix)
  grid.arrange(pmapRaw$gtable, vp=viewport(width=0.9, height=0.75))  
  grid.arrange(pmapnorm$gtable, vp=viewport(width=0.9, height=0.75))  
  dev.off()
  

  ## Standardized and Normalized data values will be used for all the downstram analysis
   results =  foreach(index=1:length(cohort_samples), .inorder = FALSE, .combine = rbind, .multicombine = TRUE, .export = ls(envir=globalenv()),.packages = packageList) %dopar% {  

    samplenames = cohort_samples[index]
    sample_methyl_matrix = episign_matrix[,grepl(paste0("Episign|",samplenames),colnames(episign_matrix))]
    sample_methyl_norm = norm_data[,grepl(paste0("Episign|",samplenames),colnames(norm_data))]
  
    episign_disease =  unique(markerdata$Disease)
    
    print(samplenames)
    # computing distance and correlation for each disease specific probes on raw and normalized & standardized methylation values
    sample_dataset = rbindlist(lapply(episign_disease, function(d, sample_methyl_matrix, sample_methyl_norm, samplenames, markerdata, compute_all_measures, merge_datasets) {
       
      dball = NULL
      db = NULL
      
      # probe specific Values for Episign & Samples
      probeset = paste0(d," specific probes") 
      
      # print(probeset)
      norm_disorder_matrix = sample_methyl_norm[row.names(sample_methyl_norm) %in% markerdata$Probes[markerdata$Disease==d],]
      d_dataset = compute_all_measures(sample_methyl_norm=norm_disorder_matrix, samplenames=samplenames, probeset=probeset, episign = d)
      s_dataset = compute_all_measures(sample_methyl_norm=norm_disorder_matrix, samplenames=samplenames, probeset=probeset)
      db = merge_datasets(d_dataset = d_dataset, s_dataset = s_dataset, samplenames=samplenames, episignname=d)
      clust = hierarchial_clustering_analysis(sample_methyl_norm=norm_disorder_matrix, samplename = samplenames, episign=d)
      db$hCluster = clust
      
      # for all episign probes
      # probeset = "All Episign Probes"
      # norm_disorder_matrix = sample_methyl_norm[row.names(sample_methyl_norm),]
      # d_dataset = compute_all_measures(sample_methyl_norm=norm_disorder_matrix, samplenames=samplenames, probeset=probeset, episign = d)
      # s_dataset = compute_all_measures(sample_methyl_norm=norm_disorder_matrix, samplenames=samplenames, probeset=probeset)
      # dball = merge_datasets(d_dataset = d_dataset, s_dataset = s_dataset, samplenames=samplenames, episignname=d)
      # clust = hierarchial_clustering_analysis(sample_methyl_norm=norm_disorder_matrix, samplename = samplenames, episign=d)
      # dball$hCluster = clust
      
      # merge dataset for all episign probes and episign specific probes
      dball = data.table(rbind(dball,db),check.names = F, stringsAsFactors = F)
      
      # Assess if the sample matches with any Episign
      dball$ProbeGroup = gsub(" .*","",dball$Probeset)
    
      return(dball)       
      
    }, sample_methyl_matrix = sample_methyl_matrix, sample_methyl_norm = sample_methyl_norm, samplenames = samplenames, markerdata = markerdata, compute_all_measures = compute_all_measures, merge_datasets = merge_datasets))
  

    #outfile
    outfile = file.path(outdir,samplenames)

    # sample_dataset = sample_dataset[Probeset!="All Episign Probes" &  norm_correlation>=0.7, is_NNcluster:= ifelse(Sample_NN_norm_dist == Episign  & Episign_NN_norm_dist==Sample & ProbeGroup==Episign, TRUE, FALSE) , by=c("Sample")]
    # sample_dataset = sample_dataset[Probeset!="All Episign Probes" &  norm_correlation>=0.7, is_hcluster:= ifelse(!is.na(hCluster), TRUE, FALSE) , by=c("Sample")]
    # sample_dataset = sample_dataset[Probeset!="All Episign Probes" &  norm_correlation>=0.75, is_Euc_Corr_pca:= ifelse(norm_correlation == max(norm_correlation[ Sample_NN_norm_dist==Sample_NN_norm_corr & ProbeGroup==Sample_NN_norm_dist & ProbeGroup==Sample_NN_pca_dist]), TRUE, FALSE) , by=c("Sample")]
    # sample_dataset = sample_dataset[Probeset!="All Episign Probes" &  norm_correlation>=0.75, is_Euc_MI_pca:= ifelse(MI == max(MI[Sample_NN_norm_dist==Sample_NN_MI & ProbeGroup==Sample_NN_norm_dist & ProbeGroup==Sample_NN_pca_dist]), TRUE, FALSE) , by=c("Sample")]
    # sample_dataset = sample_dataset[,is_Cluster := ifelse((is_NNcluster==TRUE & is_hcluster==TRUE) | is_Euc_Corr_pca==TRUE | is_Euc_MI_pca==TRUE , TRUE,FALSE), by= c("Sample")]

    # sample_dataset = sample_dataset[, is_NNcluster:= ifelse(Sample_NN_norm_dist == Episign  & Episign_NN_norm_dist==Sample & ProbeGroup==Episign, TRUE, FALSE) , by=c("Sample")]
    # sample_dataset = sample_dataset[, is_hcluster:= ifelse(!is.na(hCluster), TRUE, FALSE) , by=c("Sample")]
    # sample_dataset = sample_dataset[, is_Euc_Corr_pca:= ifelse(norm_correlation == max(norm_correlation[ Sample_NN_norm_dist==Sample_NN_norm_corr & ProbeGroup==Sample_NN_norm_dist & ProbeGroup==Sample_NN_pca_dist]), TRUE, FALSE) , by=c("Sample")]
    # sample_dataset = sample_dataset[, is_Euc_MI_pca:= ifelse(MI == max(MI[Sample_NN_norm_dist==Sample_NN_MI & ProbeGroup==Sample_NN_norm_dist & ProbeGroup==Sample_NN_pca_dist]), TRUE, FALSE) , by=c("Sample")]
    # sample_dataset$nCriteria = rowSums(sample_dataset[,c("is_NNcluster","is_hcluster","is_Euc_Corr_pca","is_Euc_MI_pca")]==TRUE,na.rm = T)

    # is nearest neighbour based on eucledian distance and PCA
    sample_dataset$is_NNcluster = ifelse(((sample_dataset$Sample_NN_norm_dist == sample_dataset$Episign  & sample_dataset$Episign_NN_norm_dist == sample_dataset$Sample) | 
                                            (sample_dataset$Sample_NN_pca_dist == sample_dataset$Episign  & sample_dataset$Episign_NN_pca_dist == sample_dataset$Sample)) & 
                                            sample_dataset$ProbeGroup==sample_dataset$Episign, TRUE, FALSE) 
    # is nearest neighbour in heirarical clustering
    sample_dataset$is_hcluster = ifelse(!is.na(sample_dataset$hCluster), TRUE, FALSE) 
    
    # has maximum correlation
    max_corr = max(sample_dataset$norm_correlation[sample_dataset$Sample_NN_norm_corr == sample_dataset$Episign  & sample_dataset$ProbeGroup == sample_dataset$Episign],na.rm = T)
    print(max_corr)
    if(is.finite(max_corr)){
      sample_dataset$is_normCorr = ifelse(sample_dataset$norm_correlation == max_corr, TRUE, FALSE)   
    }else{
      sample_dataset$is_normCorr = FALSE
    }
    
    # has maximum mutual information
    max_mi = max(sample_dataset$MI[sample_dataset$Sample_NN_MI==sample_dataset$Episign & sample_dataset$ProbeGroup==sample_dataset$Episign])
    print(max_mi)
    if(is.finite(max_mi)){
      sample_dataset$is_MI = ifelse(sample_dataset$MI == max_mi, TRUE, FALSE) 
    }else{
      sample_dataset$is_MI = FALSE
    }
    
    #calculate number of criteria 
    sample_dataset$nCriteria = rowSums(sample_dataset[,c("is_NNcluster","is_hcluster","is_normCorr","is_MI")]==TRUE,na.rm = T)
    
    # calculate fisher test
    stats_test = rbindlist(lapply(episign_disease,function(x,sample_dataset){
      match_sample_episign = sum(sample_dataset[,c("Sample_NN_pca_dist","Sample_NN_norm_dist","Sample_NN_MI")]==x)
      not_match_sample_episign = sum(sample_dataset[,c("Sample_NN_pca_dist","Sample_NN_norm_dist","Sample_NN_MI")]!=x)
      match_episign_episign = sum(sample_dataset[,c("Episign_NN_pca_dist","Episign_NN_norm_dist","Episign_NN_MI")]==x)
      not_match_episign_episign = sum(sample_dataset[,c("Episign_NN_pca_dist","Episign_NN_norm_dist","Episign_NN_MI")]!=x)
      proportion = match_sample_episign/(match_sample_episign+not_match_sample_episign)*100
      pvalue = fisher.test(matrix(c(match_sample_episign,not_match_sample_episign,match_episign_episign,not_match_episign_episign),byrow = T,nrow = 2),simulate.p.value = T ,alternative = "two.sided")$p.value
      stats_table = data.frame("Episign"=x,"pvalue"=pvalue,"proportion"=proportion,"Episign_matches"=match_sample_episign,"Episign_mismatches"=not_match_sample_episign)
      return(stats_table)
    },sample_dataset))
    
    # merge fisher statistical result with the data
    sample_dataset = merge.data.table(sample_dataset, stats_test, by=c("Episign"))
    
    max_match_samples = paste0(sample_dataset$Episign[max(sample_dataset$proportion)==sample_dataset$proportion],collapse = ";")
    
    sample_dataset$max_match = max_match_samples
    sample_dataset$adjust_pval_fdr = p.adjust(sample_dataset$pvalue, method = "fdr")
    sample_dataset$adjust_pval_BH = p.adjust(sample_dataset$pvalue, method = "BH")
    sample_dataset$adjust_pval_bonferroni = p.adjust(sample_dataset$pvalue, method = "bonferroni")
    
    # sample_dataset$is_Cluster = ifelse(sample_dataset$nCriteria==max(sample_dataset$nCriteria) & sample_dataset$norm_correlation>=0.7 & sample_dataset$nCriteria>0 & sample_dataset$pvalue<0.005, TRUE,FALSE)
    sample_dataset$is_Cluster = ifelse(sample_dataset$nCriteria==1 & sample_dataset$norm_correlation>0.7 & sample_dataset$adjust_pval_fdr<0.001 , TRUE,
                                       ifelse( sample_dataset$nCriteria>1 & sample_dataset$norm_correlation>0.7  & sample_dataset$adjust_pval_fdr<0.5, TRUE,
                                         FALSE))
    
    # write all measure value for the sample
    fwrite(sample_dataset,paste0(outfile,".txt"),row.names = F,sep = "\t",quote = F)
    
    # plot all heatmaps & pca
    print(paste0("Plotting PCA & heatmaps for all Disorders for sample",samplenames))
    plot_all_results(sample_methyl_matrix=sample_methyl_matrix, norm_data=sample_methyl_norm, markerdata=markerdata, outfile=outfile,episign_matrix = episign_matrix)
    print("Plotting completed")
    gc()
    return(sample_dataset)
  } 

     
  return(results)
  
}


episignature_feature_selection_optimization <-function(epidata = NULL, condition="FBXO22", to_normalize=TRUE){

  
  samplenames = condition
  # get pca and extract probes and their contribution to individual PCs
  norm_data = epidata
  if(to_normalize){
    norm_data = normalize_data(episign_matrix = epidata)  
  }
  norm_pca_res = pca_analysis_data(episign_matrix = norm_data, to_normalize = FALSE)

  episign_diff = data.frame("probeName" = row.names(norm_data), "norm_diff" = rowMeans(norm_data[,colnames(norm_data)!= samplenames]) - norm_data[,samplenames],stringsAsFactors = F, check.names = F)
  episign_diff$abs_norm_diff = abs(episign_diff$norm_diff)
 
  episign_diff = episign_diff[order(episign_diff$abs_norm_diff,decreasing = TRUE),]
  totProbes = nrow(episign_diff)
  optimization_probes = list()
  optimization_result = list()
  
  #foreach(nprobes=1:totProbes) %do% {
  for(nprobes in c(5:totProbes)){
    probeset=episign_diff$probeName[c(1:nprobes)]
     probesetName = paste0("probeGroup",nprobes)
    probe_matrix = norm_data[row.names(norm_data) %in% probeset,]

    probe_dataset = compute_all_measures(sample_methyl_norm=probe_matrix, samplenames=samplenames, probeset=probesetName, episign = NULL)
    optimization_probes[[probesetName]] = probeset
    optimization_result[[probesetName]] = probe_dataset
    
  }
  optimization_result = rbindlist(optimization_result,use.names = TRUE,idcol = "iteration")
  gc()
  return(list("optimization_result"=optimization_result,"optimization_probes"=optimization_probes))
  
}


# ----------------------------------- Imprinting diseases & SMA analysis -----------------------------------------------

methylation_profiling_heatmap <-function(methylationData = NULL, readdepth = 5, scoreth = 98, pdffile="result/smn1.pdf"){


  # average by base for +ve and -ve strand
  methylationData$modification[methylationData$read_depth <= readdepth & (methylationData$score/1000) <= scoreth] = 0
  methylationData$modification[!is.finite(methylationData$modification) ] = 0
  methylationData = unique(methylationData[, c("chrom","position","sampleID","modification","read_depth")])
  methylationData = methylationData[,smn1_sampleavg:=mean(read_depth),by=c("sampleID")]
  methylationData = methylationData[,smn1_posavg:=mean(read_depth),by=c("position")]
  methylationData = methylationData[,smn1_nomod:=uniqueN(position[modification==0])/uniqueN(position)*100, by=c("sampleID")]
  methylationData = methylationData[,smn1_mod:=uniqueN(position[modification>0])/uniqueN(position)*100, by=c("sampleID")]
  
  #	p1 = ggplot(methylationData, aes(x=pos, y=read_depth, group = sample)) + geom_smooth(method="auto" , fill="#69b3a2", se=TRUE, level=0.95, alpha=0.2) + 
  #	theme_classic() + theme(axis.text = element_text(size=15, colour="black")) +
  #	scale_y_continuous(limits=c(0,50), breaks=seq(0,50,by=10),expand=c(0,0))+
  #	ggrepel::geom_label_repel(data = a , aes(x=70240513, y=smn1_avg, label=sample), nudge_x=0.75, nudge_y=0.25, box.padding=0, label.size = 0.1, max.overlaps = 15, direction="both", size=2.5)
  
  
  # make unique positions for each cpg 
  
  fdata = data.frame(dcast(methylationData, chrom + position ~ sampleID, value.var = "modification", fill=0),check.names=F,stringsAsFactors=F)
  row.names(fdata) = paste0(fdata$chrom,":",fdata$position)
  fdata$chrom = NULL
  fdata$position = NULL
  fdata=as.matrix(fdata)
  
  color = hcl.colors(20, "YlGnBu")
  phm = pheatmap(fdata,clustering_distance_cols = "euclidean",clustering_distance_rows = "euclidean", cluster_rows = F, cluster_cols = T, show_rownames = F, show_colnames = T, color = color, scale="none")
  dev.off()
  
  color = hcl.colors(25, "YlGnBu")
  phm1 = pheatmap(fdata,clustering_distance_cols = "euclidean",clustering_distance_rows = "euclidean", cluster_rows = F, cluster_cols = T, show_rownames = F, show_colnames = T, color = color, scale="none")
  dev.off()
  
  color = hcl.colors(50, "YlGnBu")
  phm2 = pheatmap(fdata,clustering_distance_cols = "euclidean",clustering_distance_rows = "euclidean", cluster_rows = F, cluster_cols = T, show_rownames = F, show_colnames = T, color = color, scale="none")
  dev.off()
  
  
  color = hcl.colors(10, "YlGnBu")
  print(color)
  phm3 = pheatmap(fdata,clustering_distance_cols = "euclidean",clustering_distance_rows = "euclidean", cluster_rows = F, cluster_cols = T, show_rownames = F, show_colnames = T, color = color, scale="none")
  dev.off()
  
  colorder = phm$tree_col$order
  methylationDataSummary = unique(methylationData[,c("sampleID","smn1_nomod","smn1_sampleavg","smn1_mod")])
  methylationDataSummary$sampleID = factor(methylationDataSummary$sampleID, levels=rev(colnames(fdata)[colorder]),ordered=T)
  p2 = ggplot(methylationDataSummary, aes(x=sampleID, y=smn1_nomod)) + geom_bar(stat="identity") + 
    theme_classic() + theme(axis.text = element_text(size=15, colour="black")) + coord_flip()
  
  p3 = ggplot(methylationDataSummary, aes(x=sampleID, y=smn1_mod)) + geom_bar(stat="identity") + 
    theme_classic() + theme(axis.text = element_text(size=15, colour="black")) + coord_flip()
  
  p4 = ggplot(methylationDataSummary, aes(x=sampleID,y=0, size=smn1_nomod)) + geom_point() + 
    theme_classic() + theme(axis.text = element_text(size=15, colour="black")) + coord_flip()
  
  pdf(pdffile)
  print(phm)
  print(p2)
  plot.new()
  print(phm1)
  print(p3)
  plot.new()
  print(phm2)
  plot.new()
  print(phm3)
  print(p3)
  print(p4)	
  dev.off()
  
  return(list("summary_data"=methylationDataSummary,"methylation_data"=fdata))
  
  
}

# ----------------------------------------------------------- CNV/SV Analysis functions ---------------------------------------------------------
# Functions for CNV/SV/SNV


read_SV_CNV_files <-function(datadir = NULL, pattern=""){
  filelist = data.frame("file_path" = list.files(datadir, pattern = "",full.names = T),
                        "samplename" = gsub(paste0(pattern,"|.txt|.tsv"),"", basename(list.files(datadir, pattern = "",full.names = T)))
  )
  print(paste0("reading files from the directory ",datadir))
  genomic_data = list()
  for(index in 1:nrow(filelist)){
    
    print(paste0("reading sample ",filelist$samplename[index],"..."))
    
    sample_data = fread(file = filelist$file_path[index],check.names = F,sep = "\t",header = T)

      if(any(grepl("Chromosome",colnames(sample_data)))){
        sample_data = sample_data %>% rename(SV_type=Type)  
      }
      if(any(grepl("SV_chrom",colnames(sample_data)))){
        sample_data = sample_data %>% rename(Chromosome=SV_chrom, Start=SV_start, End=SV_end, VariantID=AnnotSV_ID)
        colnames(sample_data)[15] = "Format_value"
        sample_data$GT = gsub(":.*","",sample_data$Format_value)
        sample_data$INFO = gsub("RNAMES=.*;STDEV_LEN","STDEV_LEN",sample_data$INFO)
        
        
        variantids_with_genes = sample_data$VariantID[!sample_data$Location==""]
        sample_data = sample_data[,nGene_count:=unique(Gene_count[!(is.na(Gene_count) | Gene_count=="")]),by=c("VariantID")]
        sample_data = sample_data[!(VariantID %in% variantids_with_genes & Location==""),]
        sample_data$Gene_count = sample_data$nGene_count
        sample_data$nGene_count = NULL
      }
    
    genomic_data[[filelist$samplename[index]]] = sample_data
  }
  genomic_data = data.table(rbindlist(genomic_data, use.names = T,idcol = "SampleName", fill = T),check.names = F)
  if(any(grepl("Samples_ID",colnames(genomic_data)))){
    genomic_data$Samples_ID = NULL  
  }
  genomic_data = unique(genomic_data)
  return(genomic_data)
}


intersect_SV_CNV_calls <- function(fulldataset = NULL){
  
  overlapdata = list()
  
  dataset = unique(fulldataset[fulldataset$SV_type!="BND",c("Chromosome","Start","End","SampleName","GT","SV_type", "Location", "ALT")])
  
  dataset$Chrom = dataset$Chromosome
  
  dataset$Chromosome = paste0("chr",dataset$Chromosome)
  dataset$Chromosome = gsub("chrchr","chr",dataset$Chromosome)
  
  for(samp in unique(dataset$SampleName)){
    
    print(paste0("procesing ",samp))
    sampleSV = dataset[dataset$SampleName==samp]
    otherSV = dataset[dataset$SampleName!=samp]
    # make Granges
    sampleSV = GRanges(seqnames = sampleSV$Chromosome,
                       ranges = IRanges(start = sampleSV$Start, end = sampleSV$End),
                       Chromosome = sampleSV$Chrom,
                       Start = sampleSV$Start,
                       End = sampleSV$End,
                       SampleName = sampleSV$SampleName,
                       SV_type = sampleSV$SV_type,
                       Location = sampleSV$Location,
                       Sequence = sampleSV$ALT,
                       GT = sampleSV$GT
    )
    
    otherSV =GRanges(seqnames = otherSV$Chromosome,
                     ranges = IRanges(start = otherSV$Start, end = otherSV$End),
                     other_samp = otherSV$SampleName,
                     other_SV = otherSV$SV_type,
                     other_Loc = otherSV$Location,
                     other_seq = otherSV$ALT,                     
                     other_GT = otherSV$GT)
    #get overlaps		
    overlaps <- findOverlaps(query = sampleSV, subject = otherSV)
    sample_SV_data = data.table(data.frame(cbind(mcols(sampleSV[queryHits(overlaps)]),mcols(otherSV[subjectHits(overlaps)])),stringsAsFactors=F,check.names = F))

    sample_SV_data = sample_SV_data[,nsamplesGT:=uniqueN(other_samp[other_GT==GT & other_SV==SV_type & other_Loc==Location]),by=c("Chromosome","Start","End","SampleName")]
    sample_SV_data = sample_SV_data[,other_samplesGT:=paste0(unique(other_samp[other_GT==GT & other_SV==SV_type & other_Loc==Location]),collapse = ";"),by=c("Chromosome","Start","End","SampleName")]
    
    sample_SV_data = sample_SV_data[,nsamplesSeq:=uniqueN(other_samp[other_GT==GT & other_SV==SV_type & other_Loc==Location & other_seq==Sequence]),by=c("Chromosome","Start","End","SampleName")]
    sample_SV_data = sample_SV_data[,other_samplesSeq:=paste0(unique(other_samp[other_GT==GT & other_SV==SV_type & other_Loc==Location & other_seq==Sequence]),collapse = ";"),by=c("Chromosome","Start","End","SampleName")]
    
    sample_SV_data = sample_SV_data[,nsamples:=uniqueN(other_samp),by=c("Chromosome","Start","End","SampleName")]
    sample_SV_data = sample_SV_data[,other_samples:=paste0(unique(other_samp),collapse = ";"),by=c("Chromosome","Start","End","SampleName")]
    
    sample_SV_data = unique(sample_SV_data[,c("Chromosome","Start","End","SampleName","nsamplesGT","other_samplesGT","nsamplesSeq","other_samplesSeq","nsamples","other_samples")])
    sample_SV_data$nsamples = sample_SV_data$nsamples 
    sample_SV_data$nsamplesGT = sample_SV_data$nsamplesGT 
    sample_SV_data$nsamplesSeq = sample_SV_data$nsamplesSeq 
    
    overlapdata[[samp]] = sample_SV_data
    
  }
  overlapdata =  data.table(rbindlist(overlapdata, use.names = T, fill = T),check.names = F, stringsAsFactors = F)
  return(overlapdata)
}


intersect_SV_CNV_gene_exon_wise_calls <- function(full_data = NULL, phenotype_annotation = annotationData, col_prefix=NULL, dgv_perc_overlap_th=70, DGV_SV_frequency_th=10){
  
  phenotype_annotation = unique(phenotype_annotation[phenotype_annotation$is_disease_causing==TRUE,c("genes","combined_phenotype","combined_inheritance","inheritance","combined_classification","Clingen_annotation","DGV_SV","dgv_perc_overlap","DGV_SV_frequency")]) 
  colnames(phenotype_annotation) = c("Gene_name","disease_phenotype","combined_inheritance","inheritance","classification","Clingen_annotation","DGV_SV","dgv_perc_overlap","DGV_SV_frequency")
  
  dataset = NULL
  if(any(grepl("All protein coding genes", colnames(full_data)))){
    dataset = unique(full_data[, c("Chromosome","Start","End","variant_unique_id","SampleName","GT","SV_type", "Location", "ALT","All protein coding genes")])  
  }else{
    full_data$is_disrupting_gene = ifelse(grepl("x",full_data$Location),"TRUE",ifelse(full_data$l1!=full_data$l2,"TRUE","FALSE")) 
    dataset = unique(full_data[, c("Chromosome","Start","End","variant_unique_id","SampleName","GT","SV_type", "Location", "ALT","Gene_name")])  
  }
  colnames(dataset)=c("Chromosome","Start","End","variant_unique_id","SampleName","GT","SV_type", "Location", "ALT","Gene_name")
  dataset$Gene_name = gsub(", ",",",dataset$Gene_name)
  
  print("checking for unique genes within SVs")
  
  # get gene related annotations
  split_dataset = unique(data.table(separate_rows(dataset, Gene_name, sep = "[,;]+", convert = TRUE)))
  split_dataset = merge.data.table(split_dataset,phenotype_annotation,by=c("Gene_name"),all.x = T,allow.cartesian=TRUE) 
  split_dataset = split_dataset[,all_disease_phenotype:=paste0(unique(disease_phenotype[! (is.na(disease_phenotype) | disease_phenotype=="")]),collapse = "::"),by=c("variant_unique_id","Chromosome","Start","End","SampleName")]
  split_dataset = split_dataset[,all_inheritance:=paste0(unique(inheritance[! (is.na(inheritance) | inheritance=="")]),collapse = "::"),by=c("variant_unique_id","Chromosome","Start","End","SampleName")]
  split_dataset = split_dataset[,all_inheritance_combined:=paste0(unique(combined_inheritance[! (is.na(combined_inheritance) | combined_inheritance=="")]),collapse = "::"),by=c("variant_unique_id","Chromosome","Start","End","SampleName")]
  split_dataset = split_dataset[,all_classification:=paste0(unique(classification[! (is.na(classification) | classification=="")]),collapse = "::"),by=c("variant_unique_id","Chromosome","Start","End","SampleName")]
  split_dataset = split_dataset[,all_clingen_annotation:=paste0(unique(Clingen_annotation[! (is.na(Clingen_annotation) | Clingen_annotation=="")]),collapse = "::"),by=c("variant_unique_id","Chromosome","Start","End","SampleName")]
  
  split_dataset$all_gene_location = paste0(split_dataset$Gene_name,"(",split_dataset$Location,")")
  split_dataset = split_dataset[,all_genes:=paste0(unique(all_gene_location),collapse = ";"),by=c("variant_unique_id","Chromosome","Start","End","SampleName")]
  split_dataset$all_genes = gsub("\\(genic\\)|\\(\\)","",split_dataset$all_genes)
  
  split_dataset$sv_gt_effect = ifelse(is.na(split_dataset$inheritance),"unknown", ifelse(split_dataset$GT=="0/1" & split_dataset$inheritance=="AR","unaffected","affected") )
  common_variant = split_dataset$Gene_name[!is.na(split_dataset$DGV_SV) & split_dataset$SV_type==split_dataset$DGV_SV & split_dataset$dgv_perc_overlap>=dgv_perc_overlap_th & split_dataset$DGV_SV_frequency>=DGV_SV_frequency_th]
  
  split_dataset$sv_dgv_effect = ifelse(split_dataset$Gene_Name %in% common_variant,"common_DGV","not_DGV")
  split_dataset$sv_combined_effect = paste0(split_dataset$Gene_name,"(GT_effect-",split_dataset$sv_gt_effect,",DGV_effect-",split_dataset$sv_dgv_effect,")")
  split_dataset = split_dataset[,all_sv_combined_effect:=paste0(unique(sv_combined_effect),collapse = ";"),by=c("variant_unique_id","Chromosome","Start","End","SampleName")]
  
  # get sample names with same alterations to the genes
  split_dataset = split_dataset[,samples:=paste0(sort(unique(SampleName)), collapse = ","),by=c("Gene_name","SV_type","Location")]
  split_dataset$gene_sample_svtype = paste0(split_dataset$Gene_name,"(",split_dataset$samples,")")
  split_dataset = split_dataset[,no_of_samples_for_gene:=uniqueN(SampleName),by=c("Gene_name","SV_type","Location")]
  split_dataset = split_dataset[,no_of_samples_for_gene_GT:=uniqueN(SampleName),by=c("Gene_name","SV_type","Location","GT")]
  
  split_dataset$is_unique = ifelse(!is.na(split_dataset$DGV_SV) & split_dataset$SV_type==split_dataset$DGV_SV & split_dataset$dgv_perc_overlap>=70 & split_dataset$DGV_SV_frequency>=10,FALSE,
                            ifelse(split_dataset$no_of_samples_for_gene_GT==1 & split_dataset$GT=="1/1",TRUE,
                            ifelse(split_dataset$no_of_samples_for_gene>1,FALSE,TRUE)))
  
  
  # get gene stats within the variants  
  split_dataset = split_dataset[,no_of_genes_in_variant:=uniqueN(Gene_name),by=c("variant_unique_id","Chromosome","Start","End","SampleName")]
  split_dataset = split_dataset[,no_of_genes_unique_in_variant:=uniqueN(Gene_name[is_unique==TRUE]),by=c("variant_unique_id","Chromosome","Start","End","SampleName")]

  split_dataset = split_dataset[,no_of_disease_genes_in_variant:=uniqueN(Gene_name[!is.na(disease_phenotype)]),by=c("variant_unique_id","Chromosome","Start","End","SampleName")]
  split_dataset = split_dataset[,no_of_disease_genes_unique_in_variant:=uniqueN(Gene_name[!is.na(disease_phenotype) & is_unique==TRUE]),by=c("variant_unique_id","Chromosome","Start","End","SampleName")]

  # get gene list information within variants
  split_dataset = split_dataset[,common_gene_list:=paste0(sort(unique(Gene_name[is_unique!=TRUE])),collapse = ";"), by=c("variant_unique_id","Chromosome","Start","End","SampleName")]
  split_dataset = split_dataset[,unique_gene_list:=paste0(sort(unique(Gene_name[is_unique==TRUE])),collapse = ";"), by=c("variant_unique_id","Chromosome","Start","End","SampleName")]

  split_dataset = split_dataset[,gene_list_with_samples_svtype:=paste0(sort(unique(gene_sample_svtype)),collapse = ";"), by=c("variant_unique_id","Chromosome","Start","End","SampleName")]
  
  # make unique set for each variant
  split_dataset = unique(split_dataset[,c("Chromosome","Start","End","variant_unique_id","SampleName","all_genes",
                                          "all_disease_phenotype","all_inheritance","all_inheritance_combined","all_classification","all_sv_combined_effect","all_clingen_annotation",
                                          "no_of_genes_in_variant","no_of_genes_unique_in_variant",
                                          "no_of_disease_genes_in_variant", "no_of_disease_genes_unique_in_variant",
                                          "unique_gene_list","common_gene_list",
                                          "gene_list_with_samples_svtype")]) 

  if(!is.null(col_prefix)){
   colnames(split_dataset) = c(colnames(split_dataset)[c(1:5)],paste0(col_prefix,"_",colnames(split_dataset)[c(6:ncol(split_dataset))]))
  }
  
  return(split_dataset)
}


annotate_for_unique_genes <-function(full_data=NULL, dgv_perc_overlap_th=70, DGV_SV_frequency_th=10){
  
  print("checking for genes within SVs...")
  all_gene_full_data = intersect_SV_CNV_gene_exon_wise_calls(full_data = full_data[is_genecontaining==TRUE,], dgv_perc_overlap_th=dgv_perc_overlap_th, DGV_SV_frequency_th=DGV_SV_frequency_th)
  
  print("filter 3: Affects CDS")
  print("checking for genes within SVs affected in CDS region...")
  disrupting_full_data = intersect_SV_CNV_gene_exon_wise_calls(full_data = full_data[is_disrupting_gene=="TRUE",], col_prefix="disrupting", dgv_perc_overlap_th=dgv_perc_overlap_th, DGV_SV_frequency_th=DGV_SV_frequency_th )
  disrupting_full_data = disrupting_full_data[,c("SampleName","Chromosome","Start","End","variant_unique_id","disrupting_all_genes",
                                                 "disrupting_no_of_genes_in_variant","disrupting_no_of_disease_genes_in_variant","disrupting_no_of_disease_genes_unique_in_variant",
                                                 "disrupting_unique_gene_list")]
  
  print("merging all gene annotations within SVs")
  
  full_data = Reduce(function(...) merge.data.table(...,by=c("SampleName","Chromosome","Start","End","variant_unique_id"), all.x =T), list(full_data,all_gene_full_data,disrupting_full_data))
 
  # re-annotating full annotation rows of AnnotSV with the split gene data values
  variants_disrupting_genes = unique(full_data$variant_unique_id[full_data$is_disrupting_gene=="TRUE"])
  full_data$is_disrupting_gene[full_data$variant_unique_id %in% variants_disrupting_genes]="TRUE"
  full_data$is_disrupting_gene[!full_data$variant_unique_id %in% variants_disrupting_genes]="FALSE"
  
  full_data$no_of_genes_in_variant[is.na(full_data$no_of_genes_in_variant)] = 0
  full_data$no_of_genes_unique_in_variant[is.na(full_data$no_of_genes_unique_in_variant)] = 0
  full_data$no_of_disease_genes_in_variant[is.na(full_data$no_of_disease_genes_in_variant)] = 0
  full_data$no_of_disease_genes_unique_in_variant[is.na(full_data$no_of_disease_genes_unique_in_variant)] = 0
  
  full_data$disrupting_no_of_genes_in_variant[is.na(full_data$disrupting_no_of_genes_in_variant)] = 0
  full_data$disrupting_no_of_disease_genes_in_variant[is.na(full_data$disrupting_no_of_disease_genes_in_variant)] = 0
  full_data$disrupting_no_of_disease_genes_unique_in_variant[is.na(full_data$disrupting_no_of_disease_genes_unique_in_variant)] = 0

  allcolumnNames = c("SampleName","Chromosome","Start","End","variant_unique_id","GT","SV_type","SV_length","LOG2CNT","SCORE","Classification","QUAL","FILTER","AF","SUPPORT","RE","is_passing","is_genecontaining","is_disrupting_gene",
                     "is_disease_gene","is_unique","is_disease_disrupting_unique","all_genes","all_disease_phenotype","all_inheritance","all_inheritance_combined","all_classification","all_sv_combined_effect","all_clingen_annotation",
                     "no_of_genes_in_variant","no_of_disease_genes_in_variant","no_of_genes_unique_in_variant","no_of_disease_genes_unique_in_variant",
                     "unique_gene_list","common_gene_list",
                     "disrupting_all_genes","disrupting_no_of_genes_in_variant","disrupting_no_of_disease_genes_in_variant","disrupting_no_of_disease_genes_unique_in_variant",
                     "disrupting_unique_gene_list",
                     "VariantID","ID", "ALT","INFO","FORMAT","Format_value")
  cols2select = allcolumnNames[allcolumnNames %in% colnames(full_data)]
  full_data = unique(full_data[,cols2select, with=FALSE])
 
  
  # filter 3: encompassing genes filters
  print("filter 4: Affects disease associated gene list filters")
  full_data$is_disease_gene = ifelse(full_data$no_of_disease_genes_in_variant > 0,"TRUE","FALSE") 
  full_data$is_disrupting_disease_gene = ifelse(full_data$disrupting_no_of_disease_genes_in_variant > 0,"TRUE","FALSE") 
  
  print("filter 5: unique to the sample")
  full_data$is_unique = ifelse(full_data$no_of_genes_unique_in_variant > 0,"TRUE","FALSE")
  
  full_data$is_disease_unique = ifelse(full_data$no_of_disease_genes_unique_in_variant > 0,"TRUE","FALSE")
  full_data$is_disease_disrupting_unique = ifelse(full_data$disrupting_no_of_disease_genes_unique_in_variant > 0,"TRUE","FALSE")
  
  full_data$passing_calls = ifelse(full_data$is_passing == "TRUE", "TRUE","FALSE")
  full_data$passing_gene_calls = ifelse(full_data$is_genecontaining == "TRUE" & full_data$is_passing == "TRUE", "TRUE","FALSE")
  full_data$passing_disrupting_gene_calls = ifelse(full_data$is_disrupting_gene == "TRUE" & full_data$is_passing == "TRUE", "TRUE","FALSE")  
  full_data$passing_disrupting_disease_calls = ifelse(full_data$is_disrupting_disease_gene == "TRUE" & full_data$is_passing == "TRUE", "TRUE","FALSE")  
  full_data$passing_disrupting_disease_unique_calls = ifelse(full_data$is_disease_disrupting_unique == "TRUE" & full_data$is_passing == "TRUE", "TRUE","FALSE") 
  
  allcolumnNames = c("SampleName","Chromosome","Start","End","variant_unique_id","GT","SV_type","SV_length","LOG2CNT","SCORE","Classification","QUAL","FILTER","AF","SUPPORT","RE",
                     "passing_disrupting_disease_unique_calls","passing_disrupting_disease_calls","passing_disrupting_gene_calls","passing_gene_calls", "passing_calls",
                     "no_of_genes_in_variant", "all_genes","all_disease_phenotype","all_inheritance","all_inheritance_combined","all_classification","all_sv_combined_effect","all_clingen_annotation","unique_gene_list",
                     "no_of_genes_unique_in_variant","no_of_disease_genes_in_variant","no_of_disease_genes_unique_in_variant","disrupting_no_of_genes_in_variant","disrupting_no_of_disease_genes_in_variant","disrupting_no_of_disease_genes_unique_in_variant",
                     "disrupting_all_genes","disrupting_unique_gene_list","common_gene_list",
                     "is_genecontaining","is_disrupting_gene","is_disrupting_disease_gene","is_disease_disrupting_unique","is_disease_gene","is_unique",
                     "VariantID","ID", "ALT","INFO","FORMAT","Format_value")
  cols2select = allcolumnNames[allcolumnNames %in% colnames(full_data)]
  full_data = unique(full_data[,cols2select, with=FALSE])

  return(full_data)
}
  

add_CNV_SV_filtering_annotation <-function( full_data = NULL, method = "qdnaseq", AF_thresh = 0.3, Support_read = 5, Log2CNT = 0.5, phenotype_annotation = annotationData, dgv_perc_overlap_th=70, DGV_SV_frequency_th=10){
  
  
  print("Adding gene annotations...")
  full_data$variant_unique_id = paste0(full_data$SampleName,"_",full_data$VariantID)
  
  if(! any(grepl("Location", colnames(full_data)))){
    full_data$Location = ifelse(full_data$`All protein coding genes` == "","genomic","genic")
  }
  if(! any(grepl("ALT", colnames(full_data)))){
    full_data$ALT = ""
  }

  print("Checking on filtering criterias...")
  if(method == "qdnaseq"){
    
    # filter 1: passing quality filters
    print("filter 1: passing quality filters")
    full_data$is_passing = ifelse(abs(full_data$LOG2CNT) >= Log2CNT, "TRUE","FALSE")
    
    # filter 2: encompassing genes filters
    print("filter 2: encompassing genes filters")
    full_data$is_genecontaining = ifelse(! full_data$`All protein coding genes` == "","TRUE","FALSE") 
    full_data$is_disrupting_gene = ifelse(! full_data$`All protein coding genes` == "","TRUE","FALSE")
  } else{
    
      print("filter 1: passing quality filters")
      # filter 1: passing quality filters
      if(method=="sniffles"){
        full_data$is_passing = ifelse(full_data$AF >= AF_thresh  & full_data$SUPPORT >= Support_read, "TRUE","FALSE")
      }
      if(method=="cuteSV"){
        full_data$is_passing = ifelse(full_data$AF >= AF_thresh  & full_data$RE >= Support_read, "TRUE","FALSE")
      }
    
      # filter 2: encompassing genes filters
      print("filter 2: encompassing genes filters")
      full_data$is_genecontaining = ifelse(full_data$Gene_count>0 | ! (is.na(full_data$Gene_name) | full_data$Gene_name=="" ),"TRUE","FALSE") 
      setDT(full_data)[, paste0("l", 1:2) := tstrsplit(Location, "-")]
      full_data$is_disrupting_gene = ifelse(grepl("x",full_data$Location),TRUE,ifelse(full_data$l1!=full_data$l2,TRUE,FALSE)) 
  }

  full_data = annotate_for_unique_genes(full_data = full_data, dgv_perc_overlap_th=70, DGV_SV_frequency_th=10)
  full_data$method = method
  gc()
  return(unique(full_data))
}


get_unique_CNV_SV_within_sample <- function(full_data = NULL,overlap_th = 0.8){
  
  print("getting CNV summary statistics...")
  full_data$chrom = paste0("chr",full_data$Chromosome)
  full_data$start = full_data$Start
  full_data$end = full_data$End
  
  x <- group_by(full_data, SampleName)
  y <- group_by(unique(full_data[,c("chrom","start","end","SampleName","method","ID","SV_type","GT")]), SampleName)
  result = data.table(bed_intersect(x, y))
 # print(head(result))
  
  result = result[SV_type.x==SV_type.y & GT.x==GT.y & (.overlap/(End.x-Start.x)>=overlap_th),]
  result$sv_value = paste0(result$variant_unique_id.x,"(",result$INFO.x,";QUAL=",result$QUAL.x,";FILTER=",result$FILTER.x,";FORMAT=",result$FORMAT.x,"=",result$Format_value.x,")")
  result$method_value = paste0(result$method.x,"(",result$`.overlap`,")")
  
  result = result[,ID:=paste(sort(unique(ID.y,ID.x)),sep="",collapse=","),by=c("ID.x","SampleName.x","SV_type.x")]
  result = result[,all_variant_unique_id:=paste(sort(unique(variant_unique_id.x)),sep="",collapse=","),by=c("ID","SampleName.x","SV_type.x")]
  result = result[,all_variant_info:=paste(sort(unique(sv_value)),sep="",collapse=","),by=c("ID","SampleName.x","SV_type.x")]
  result = result[,all_methods:=paste(sort(unique(method_value)),sep="",collapse=","),by=c("ID","SampleName.x","SV_type.x")]
  result = result[,outer_start:=min(Start.x),by=c("ID","SampleName.x","SV_type.x")]
  result = result[,outer_end:=max(End.x),by=c("ID","SampleName.x","SV_type.x")]
  
  #update all columns
  
  result = result[,outer_end:=max(End.x),by=c("ID","SampleName.x","SV_type.x")]
  
  result = unique(result[,c("Start.x","End.x","chrom","start.y","end.y","SampleName.y","method.y","ID.y","SV_type.y","ID.x","sv_value","method_value",".overlap","variant_unique_id.x","QUAL.x","FILTER.x","FORMAT.x","Format_value.x","SV_length.x","AF.x",
                            "VariantID.x","INFO.x","method.x","start.x","end.x","GT.y"):=NULL])
  colnames(result) = gsub("\\.x","",colnames(result))
  
  
  colorder1 = c("ID","SampleName","Chromosome","GT","outer_start","outer_end")
  colorder2 = colnames(result)[!colnames(result) %in% colorder1]
  result = unique(result[,c(colorder1,colorder2),with=FALSE])
  colnames(result)[1]="VariantID"
  
  return(unique(result))
}


get_CNV_SV_stats <-function(full_data = NULL, Method="cuteSV"){
  
  print(paste0("getting ",Method," SV summary statistics..."))

  # Calculating Freq
  summary_data <- full_data %>% 
    group_by(SampleName) %>% # For each sample 
    mutate(
      calls = n_distinct(VariantID), 
      passing_calls = n_distinct(VariantID[passing_calls == "TRUE"]), 
      passing_gene_calls = n_distinct(VariantID[passing_gene_calls == "TRUE"]), 
      passing_disrupting_gene_calls = n_distinct(VariantID[passing_disrupting_gene_calls == "TRUE"]),
      passing_disrupting_disease_calls = n_distinct(VariantID[passing_disrupting_disease_calls == "TRUE"]),
      passing_disrupting_disease_unique_calls = n_distinct(VariantID[passing_disrupting_disease_unique_calls == "TRUE"]),
    )
  summary_data = unique(summary_data[,c("SampleName","calls","passing_calls","passing_gene_calls","passing_disrupting_gene_calls","passing_disrupting_disease_calls","passing_disrupting_disease_unique_calls")])
  
  summary_data$Method = Method
  return(summary_data)
  
}

# --------------------------------------- Filter for SNVs -----------------------------------
make_long_datatable <-function(Disease_Variants=NULL){
  
  # make table in long format
  b = colnames(Disease_Variants)
  samples = unique(gsub(" .*","",b[grepl("OXN",colnames(Disease_Variants))]))
  sample_disease_variant = list()
  annotationcols = b[!grepl("OXN",colnames(Disease_Variants))]
  for(sn in samples){                                                                                                                                                                                                 
    snCols = b[grepl(sn,colnames(Disease_Variants))]
    newsnCols = gsub(".* |\\(|\\)","",snCols)
    
    sample_disease_variant[[sn]] = Disease_Variants[,c(annotationcols,snCols),with=FALSE]
    colnames(sample_disease_variant[[sn]]) = c(annotationcols,newsnCols)
    
  }
  return(sample_disease_variant)
}

add_gene_annotation<-function(annotationData = NULL, data2annotated = NULL){
  phenotype_annotation = unique(annotationData[annotationData$is_disease_causing==TRUE,c("genes","combined_phenotype")]) 
  colnames(phenotype_annotation) = c("Gene Names","disease_phenotype")
  
  columnorder = c(colnames(data2annotated),"disease_phenotype")
  
  single_no_gene = data2annotated[! grepl(",",data2annotated$`Gene Names`),]
  single_no_gene = merge.data.table(single_no_gene,phenotype_annotation,by=c("Gene Names"),all.x = T) 
  single_no_gene = single_no_gene[,columnorder,with=FALSE]
  
  multigene = data2annotated[grepl(",",data2annotated$`Gene Names`),]
  multigene$disease_phenotype = 
    apply(multigene, 1, function(x, annotation){
      # print(x[names(x)=="Gene Names"])
      a = str_split_1(x[names(x)=="Gene Names"],",")
      paste0(annotation$disease_phenotype[annotation$`Gene Names` %in% a],collapse = "; ")
    }, phenotype_annotation)
  
  data2annotated = data.table(rbind(single_no_gene,multigene)) 
  return(data2annotated)
  
}

get_splicing_variants<-function(missed_exome_variants = NULL){
  
  spliceai = unique(missed_exome_variants[`Effect (Combined)`=="LoF" & `Sequence Ontology (Combined)` %in% c("splice_donor_variant","splice_acceptor_variant","splice_region_variant"), c("unique_identifier","HGVS c. (Clinically Relevant)")])
  spliceai = data.table(separate_rows(spliceai, `HGVS c. (Clinically Relevant)`, sep = "[,]", convert = TRUE))
  spliceai$nm = gsub("\\:.*","",spliceai$`HGVS c. (Clinically Relevant)`)
  spliceai$position = "splice_site"
  spliceai$position[grepl("c\\.\\-|n\\.\\-|\\*|\\-[0-9]{2}+|\\+[0-9]{2}+|\\-[3-9]|\\+[3-9]",spliceai$`HGVS c. (Clinically Relevant)`,perl = T)] = "introninc"
  
  spliceai$ensembl_url = paste0("https://grch37.rest.ensembl.org/vep/human/hgvs/",spliceai$`HGVS c. (Clinically Relevant)`,"?content-type=application/json&vcf_string=1")
  
  
  # get spliceAI scores for splicing variants
  spliceai = rbindlist(apply(spliceai,1,function(x){
    x = data.frame(as.list(x))
    si = data.frame("t_refseq_ids"=NA,"g_name"=NA,"t_id"=NA,"t_type"=NA,"t_priority"=NA,"DS_AG"=NA,"DS_AL"=NA,"DS_DG"=NA,"DS_DL"=NA,"DP_AG"=NA,"DP_AL"=NA,"DP_DG"=NA,"DP_DL"=NA)
    db = data.frame(x,si)
    
    try(
      {
        ensembl_data = jsonlite::read_json(x$ensembl_url[1],simplifyVector = T)
        Sys.sleep(5)
        spliceai_url = paste0("https://spliceai-37-xwkwwwxdwq-uc.a.run.app/spliceai/?hg=37&variant=chr",ensembl_data$vcf_string[1])
        splice_data = jsonlite::read_json(spliceai_url,simplifyVector = T)
        if(any(grepl("scores",names(splice_data)))){
          si = splice_data$scores[,c("t_refseq_ids","g_name","t_id","t_type","t_priority","DS_AG","DS_AL","DS_DG","DS_DL","DP_AG","DP_AL","DP_DG","DP_DL")]
          db = data.frame(x,si)
        }else{
          si = data.frame("t_refseq_ids"=NA,"g_name"=NA,"t_id"=NA,"t_type"=NA,"t_priority"=NA,"DS_AG"=NA,"DS_AL"=NA,"DS_DG"=NA,"DS_DL"=NA,"DP_AG"=NA,"DP_AL"=NA,"DP_DG"=NA,"DP_DL"=NA)
          db = data.frame(x,si)
        }
      },silent = T)
    Sys.sleep(5)
    return(db)
  }))
  
  spliceai$t_priority_value = ifelse(is.na(spliceai$t_priority) | spliceai$t_priority=="N",0,
                                     ifelse(spliceai$t_priority=="MS",3,ifelse(spliceai$t_priority=="MP",2,ifelse(spliceai$t_priority=="C",1,0))))
  spliceai$spliceai_score = paste0(spliceai$g_name,"(",spliceai$t_refseq_ids,", ",spliceai$t_id,", ",spliceai$t_type,":",
                                   "DS_AG:",spliceai$DS_AG,"|DS_AL:",spliceai$DS_AL,"|DS_DG:",spliceai$DS_DG,"|DS_DL:",spliceai$DS_DL,
                                   "|DP_AG:",spliceai$DP_AG,"|DP_AL:",spliceai$DP_AL,"|DP_DG:",spliceai$DP_DG,"|DP_DL:",spliceai$DP_DL,")")
  
  #  MS: MANE Select transcript
  #  MP: MANE Plus Clinical transcript
  #  C: Canonical transcript
  
  
  # select spliceai passing variants
  passing_filter_DS_0.2 = spliceai[(DS_AG>=0.2  | DS_AL>=0.2  | DS_DG>=0.2  | DS_DL>=0.2) & position == "splice_site" ]
  passing_filter_DS_0.2 = passing_filter_DS_0.2[,spliceai_scores:= paste0(unique(spliceai_score),collapse = "; "),by=c("unique_identifier")]
  passing_filter_DS_0.2 = passing_filter_DS_0.2[,splice_transcript_types:= paste0(unique(t_type),collapse = "; "),by=c("unique_identifier")]
  passing_filter_DS_0.2 = passing_filter_DS_0.2[,splice_positions:= paste0(unique(position),collapse = "; "),by=c("unique_identifier")]
  passing_filter_DS_0.2$splice_notes = "passing splice score"
  passing_filter_DS_0.2 = unique(passing_filter_DS_0.2[,c("unique_identifier","spliceai_scores","splice_positions","splice_transcript_types","splice_notes")])
  
  # select variants not passing variants and report for the main transcript
  not_passing_filter_DS_0.2 = spliceai[! spliceai$unique_identifier %in% unique(passing_filter_DS_0.2$unique_identifier) & !is.na(DS_AG)]
  not_passing_filter_DS_0.2 = data.table(not_passing_filter_DS_0.2 %>%
                                           group_by(unique_identifier) %>%
                                           filter(t_priority_value == max(t_priority_value, na.rm=TRUE)))
  not_passing_filter_DS_0.2 = not_passing_filter_DS_0.2[,spliceai_scores:= paste0(unique(spliceai_score),collapse = "; "),by=c("unique_identifier")]
  not_passing_filter_DS_0.2 = not_passing_filter_DS_0.2[,splice_transcript_types:= paste0(unique(t_type),collapse = "; "),by=c("unique_identifier")]
  not_passing_filter_DS_0.2 = not_passing_filter_DS_0.2[,splice_positions:= paste0(unique(position),collapse = "; "),by=c("unique_identifier")]
  not_passing_filter_DS_0.2$splice_notes = "failing splice score"
  not_passing_filter_DS_0.2 = unique(not_passing_filter_DS_0.2[,c("unique_identifier","spliceai_scores","splice_positions","splice_transcript_types","splice_notes")])
  
  # variants with no information
  no_spliceai = unique(spliceai[! spliceai$unique_identifier %in% unique(c(passing_filter_DS_0.2$unique_identifier,not_passing_filter_DS_0.2$unique_identifier)),c("unique_identifier","position","nm")])
  no_spliceai$t_type = ifelse(grepl("NM",no_spliceai$nm),"protein_coding","non_coding")
  no_spliceai = no_spliceai[,splice_positions:= paste0(unique(position),collapse = "; "),by=c("unique_identifier")]
  no_spliceai = no_spliceai[,splice_transcript_types:= paste0(unique(t_type),collapse = "; "),by=c("unique_identifier")]
  no_spliceai$spliceai_scores = NA
  no_spliceai$splice_notes = "no_splice_scores"
  no_spliceai = unique(no_spliceai[,c("unique_identifier","spliceai_scores","splice_positions","splice_transcript_types","splice_notes")])
  
  
  # make splice ai score related annotation
  splice_annotation = unique(rbindlist(list(passing_filter_DS_0.2, not_passing_filter_DS_0.2,no_spliceai)))
  
  # get complete list of selected/filtered variants 
  splice_variants_passing_filters = unique(missed_exome_variants[unique_identifier %in% passing_filter_DS_0.2$unique_identifier] )
  splice_variants_failing_filters = unique(missed_exome_variants[unique_identifier %in% not_passing_filter_DS_0.2$unique_identifier] )
  splice_variants_with_noinfo = unique(missed_exome_variants[unique_identifier %in% no_spliceai$unique_identifier])
  
  return(list("splice_variants_passing_filters" = splice_variants_passing_filters,
              "splice_variants_failing_filters" = splice_variants_failing_filters,
              "splice_variants_with_noinfo" = splice_variants_with_noinfo,
              "splice_annotation" = splice_annotation))
  
}

filter_variants <-function(missed_exome_variants = NULL){
  
  missed_exome_variants$unique_identifier = paste0(missed_exome_variants$`Chr:Pos`,"_",missed_exome_variants$`Ref/Alt`)
  
  # filter for DM & Classification
  Disease_Variants =  missed_exome_variants[CLASS == "DM" | Classification %in% c("Pathogenic","Likely Pathogenic") ] 
  
  # filter for LOFs without splicing variants
  lof_selection = c("disruptive_inframe_insertion","disruptive_inframe_deletion","exon_loss_variant","frameshift_variant","inframe_insertion","inframe_deletion","stop_lost","stop_retained_variant","stop_gained")  
  LOFs = missed_exome_variants[(`Effect (Combined)`=="LoF" &  `Sequence Ontology (Combined)` %in% lof_selection) & ! unique_identifier %in% Disease_Variants$unique_identifier]
  LOFs = LOFs[grepl("NM",LOFs$`HGVS c. (Clinically Relevant)`),]
  
  # filter for LOFs with splicing variants
  splice_variants = get_splicing_variants(missed_exome_variants=missed_exome_variants)
  
  splice_variants_passing_filters = splice_variants$splice_variants_passing_filters
  
  Disease_Variants = rbindlist(list(Disease_Variants,LOFs,splice_variants_passing_filters))
  Disease_Variants = merge.data.table(Disease_Variants,splice_variants$splice_annotation,all.x = T)
  
  sample_disease_variant = make_long_datatable(Disease_Variants=Disease_Variants)
  sample_disease_variant = rbindlist(sample_disease_variant,use.names = F,idcol = "sampleName")
  sample_disease_variant = sample_disease_variant[GT!="./."]
  
  notpassing_filter = merge.data.table(splice_variants$splice_variants_failing_filters,splice_variants$splice_annotation,all.x = T)
  notpassing_filter =rbindlist( make_long_datatable(Disease_Variants = notpassing_filter ),use.names = F,idcol = "sampleName")
  notpassing_filter$variant_category = paste0(notpassing_filter$CLASS,";",notpassing_filter$`Effect (Combined)`,";",notpassing_filter$Classification)
  notpassing_filter$variant_category = gsub("^;|;$","",notpassing_filter$variant_category)
  notpassing_filter = merge.data.table(notpassing_filter,patient_info,by = c("sampleName"),all.x = T)
  fwrite(notpassing_filter,"results/negative_cohort/SNV/negative_cohort_missed_exonic_variants_notpassing_filter.txt",row.names = F,sep = "\t",quote = F)
  
  
  return(sample_disease_variant)
}



# --------------------------------------- Splice variant analysis functions ----------------------------------
extract_splicevariants <- function(spliceai_out=""){
  
  splicedata = read.vcfR(spliceai_out, verbose = FALSE) 
  splicedata = data.table(vcfR2tidy(splicedata, single_frame=T )$dat)
  splicedata = splicedata[!is.na(SpliceAI)]
  
  split_splicedata = rbindlist(apply(splicedata, 1, function(x){
    spfull = unique(unlist(strsplit(x[names(x) == "SpliceAI"],split = ',')))
    sites = data.table(t(x[names(x) != "SpliceAI"]),"spfull"=spfull)
    sites = setDT(sites)[, c("ALLELE","SYMBOL","DS_AG","DS_AL","DS_DG","DS_DL","DP_AG","DP_AL","DP_DG","DP_DL") := tstrsplit(spfull, "\\|")]
  }))
  
  columns <-c("POS", "QUAL", "gt_GQ","gt_DP","gt_AF","DS_AG","DS_AL","DS_DG","DS_DL","DP_AG","DP_AL","DP_DG","DP_DL")
  split_splicedata[, columns] <- lapply(columns, function(x) as.numeric(split_splicedata[[x]]))
  split_splicedata$cDP_AG = split_splicedata$POS + split_splicedata$DP_AG
  split_splicedata$cDP_AL = split_splicedata$POS + split_splicedata$DP_AL
  split_splicedata$cDP_DG = split_splicedata$POS + split_splicedata$DP_DG
  split_splicedata$cDP_DL = split_splicedata$POS + split_splicedata$DP_DL
  
  return(split_splicedata)
}

find_overlap <- function(splice_data = NULL, annotation_data = NULL){
  
  sampleSNV = makeGRangesFromDataFrame(splice_data,
                                       keep.extra.columns=TRUE,
                                       ignore.strand=TRUE,
                                       seqinfo=NULL,
                                       seqnames.field="CHROM",
                                       start.field="POS",
                                       end.field="POS",
                                       starts.in.df.are.0based=FALSE)
  
  
  geneExon = makeGRangesFromDataFrame(annotation_data,
                                      keep.extra.columns=TRUE,
                                      ignore.strand=FALSE,
                                      seqinfo=NULL,
                                      seqnames.field="chrom",
                                      start.field="slEXstart", 
                                      end.field = "slEXend",
                                      starts.in.df.are.0based=FALSE)
  
  #get overlaps		
  overlaps <- findOverlaps(query = sampleSNV, subject = geneExon)
  db = data.frame(geneExon[subjectHits(overlaps)],stringsAsFactors=F,check.names=F)
  samp = data.frame(sampleSNV[queryHits(overlaps)],stringsAsFactors=F,check.names=F)
  colnames(samp)[c(1:2)] = c("chrom","POS")
  samp$end = NULL
  samp$strand = NULL
  
  overlapdata = data.table(cbind.data.frame(samp, db),stringsAsFactors=F,check.names=F)
  overlapdata = overlapdata[SYMBOL == genes]
  
  return(overlapdata)
}

annotate_assign_impact <-function(splice_data = NULL, annotation_data = NULL, spliceAi_th = 0.5){
  
  # intersect with the gene / exon data
  overlapdata = find_overlap(splice_data = splice_data, annotation_data = annotation_data)
  overlapdata = overlapdata[,c("DGV_SV","dgv_perc_overlap","DGV_ID","DGV_SV_frequency"):=NULL]
  overlapdata = unique(overlapdata[!grepl("\\|.\\|.\\|.\\|.\\|.\\|.\\|.\\|.",overlapdata$spfull),])
  overlapdata$strand = as.character(overlapdata$strand)
  
  # assigning Exon start & end based on strand information
  overlapdata$EXstart_ = ifelse(overlapdata$strand=="+", overlapdata$EXstart, overlapdata$EXend)
  overlapdata$EXend_ = ifelse(overlapdata$strand=="+", overlapdata$EXend, overlapdata$EXstart)
  overlapdata$EXstart = overlapdata$EXstart_
  overlapdata$EXend = overlapdata$EXend_
  overlapdata$EXstart_ = NULL
  overlapdata$EXend_ = NULL
  
  #calculate distances w.r.t start & end and assign upstream & downstream
  
  locations = c("cDP_AG","cDP_AL","cDP_DG","cDP_DL")
  db = overlapdata[,c("POS","EXstart","EXend","strand","cDP_AG","cDP_AL","cDP_DG","cDP_DL")]
  
  cals = do.call("cbind",lapply(locations, function(x,db){ 
    cols2names = paste0(x,c("_distEXstart","_distEXend","_EXannotation"))
    
    a = as.vector(unlist(db[,c(x),with = FALSE]))
    diffstart = abs(a - db$EXstart)
    diffend = abs(a - db$EXend)
    EXSt = ifelse(db$strand == "+", ifelse(a==db$EXstart,"EXstart", ifelse(a > db$EXstart,"EXStdn","EXStup")), ifelse(a==db$EXstart,"EXstart", ifelse(a < db$EXstart,"EXStdn","EXStup")) )
    EXEn = ifelse(db$strand == "+", ifelse(a==db$EXend,"EXend", ifelse(a > db$EXend,"EXEndn","EXEnup")), ifelse(a==db$EXend,"EXend", ifelse(a < db$EXend,"EXEndn","EXEnup")) )
    EXann = ifelse(db$strand == "+", ifelse(a >= db$EXend, EXEn, EXSt ), ifelse( a <= db$EXend, EXEn, EXSt))
    ab = data.frame(diffstart,diffend,EXann,stringsAsFactors=F)
    colnames(ab) = cols2names
    return(ab)
  },db = db))
  

  overlapdata = data.table(cbind(overlapdata,cals),stringsAsFactors=F)
  
  # assign constitutive splicing event
  overlapdata$CS = ifelse((overlapdata$DS_AG>0 & overlapdata$cDP_AG_EXannotation == "EXstart") | (overlapdata$DS_DG>0 & overlapdata$cDP_DG_EXannotation == "EXend"), TRUE, FALSE)
  
  
  # assign exon acceptor loss events - exon skipping  
  overlapdata$ES = ifelse((overlapdata$DS_AL>spliceAi_th & overlapdata$DS_AG==0 & overlapdata$cDP_AL_EXannotation == "EXstart") | (overlapdata$DS_DG>spliceAi_th & overlapdata$DS_DL==0 & overlapdata$cDP_DG_EXannotation == "EXstart"), TRUE, FALSE)
  
  # assign alternate 5' events 
  overlapdata$alt5p = ifelse( ( overlapdata$DS_DL > spliceAi_th & overlapdata$cDP_DL_EXannotation == "EXend" ) | ( overlapdata$DS_DG > spliceAi_th & ( overlapdata$cDP_DG_EXannotation == "EXStdn" |  overlapdata$cDP_DG_EXannotation == "EXEndn") ), TRUE, FALSE)
  
  # assign alternate 3' events   
  overlapdata$alt3p = ifelse( ( overlapdata$DS_AL > spliceAi_th & overlapdata$cDP_AL_EXannotation == "EXstart" ) | ( overlapdata$DS_AG > spliceAi_th & ( overlapdata$cDP_AG_EXannotation == "EXStup" |  overlapdata$cDP_AG_EXannotation == "EXEnup")), TRUE, FALSE)
  
  
  overlapdata$other = ifelse( ( overlapdata$DS_AG > spliceAi_th & overlapdata$cDP_AG_EXannotation == "EXEndn" ) | ( overlapdata$DS_DG > spliceAi_th & overlapdata$cDP_DG_EXannotation == "EXStup"), TRUE, FALSE)
  
  
  overlapdata = overlapdata[,has_a_CS:=ifelse(CS=="TRUE",TRUE,FALSE),by=c("genes","POS")]
  overlapdata = overlapdata[,pCS:=uniqueN(transcript[CS=="TRUE"])/nExons*100,by=c("genes","POS")]
  overlapdata = overlapdata[,pES:=uniqueN(transcript[ES=="TRUE" & has_a_CS =="FALSE"])/nExons*100,by=c("genes","POS")]
  overlapdata = overlapdata[,palt5p:=uniqueN(transcript[alt5p=="TRUE" & has_a_CS =="FALSE" ])/nExons*100,by=c("genes","POS")]
  overlapdata = overlapdata[,palt3p:=uniqueN(transcript[alt3p=="TRUE" & has_a_CS =="FALSE" ])/nExons*100,by=c("genes","POS")]
  overlapdata = overlapdata[,pother:=uniqueN(transcript[other=="TRUE" & has_a_CS =="FALSE" ])/nExons*100,by=c("genes","POS")]
  
  return(overlapdata)
}



source("Rscripts/utility.R")

cohort_sample = fread("data/ont_sample.stats",header = T, sep = "\t",stringsAsFactors = F)

pos_cohort_sample = cohort_sample[!grepl("Neg",Cohort)]
neg_cohort_sample = cohort_sample[grepl("Neg",Cohort)]
neg_cohort_sample$episignname = gsub("-","_",neg_cohort_sample$SampleName)

methylationData=readRDS("data/methylationData.rds")


# ----------------------------------- HMA ----------------------------------------
reference = "hg19"

episign_matrix = create_episign_matrix_from_methylationData(methylationData = methylationData,valuetype = "methylMean")
episign_matrix = create_episign_matrix(methylvalues=episign_matrix, reference = reference, combine_normalize_data=FALSE)


cols2select = colnames(episign_matrix)
cols2select = unique(c( cols2select[cols2select %in% neg_cohort_sample$episignname], cols2select[cols2select %in% c("Episign_HMA","Episign_Control")]))

markerdata = extract_episignature_probes(reference = "hg19")
  
cohort_matrix = episign_matrix[row.names(episign_matrix) %in% markerdata$Probes[markerdata$Disease=="Episign_HMA"], colnames(episign_matrix) %in% cols2select]

pmap = pheatmap(t(cohort_matrix),clustering_distance_cols = "euclidean",clustering_distance_rows = "euclidean", clustering_method = "ward.D2", show_rownames = T, show_colnames = F,fontsize = 10, main = "HMA_figure", treeheight_row = 40, treeheight_col = 40, scale = "column", silent = T,color = hcl.colors(5, "BluYl"))
grid.arrange(pmap$gtable, vp=viewport(width=0.9, height=0.75)) 

pdf("results/negative_cohort/MNDD/OXN-062_HMA_positive.pdf")
pmap
dev.off()


# ---------------------------------- Angelman ------------------------------------
angelman = methylation_profiling_heatmap(methylationData = methylationData$ImD_methylationData[methylationData$ImD_methylationData$sampleID %in% pos_cohort_sample$SampleName,], readdepth =5, scoreth = 98, pdffile="results/positive_cohort/Angleman_heatmap.pdf")
  
angelman_distribution = angelman$methylation_data
pos = row.names(angelman_distribution)
angelman_distribution = data.table(angelman_distribution)
angelman_distribution$position = pos

angelman_distribution = melt.data.table(angelman_distribution,id.vars = c("position"),value.name = c("modification"),variable.name = "sampleID")
angelman_distribution$Group = ifelse(angelman_distribution$sampleID == "OXN-018","OXN-018\n(Angelman)",ifelse(angelman_distribution$sampleID == "OXN-068","OXN-068","Negative"))
angelman_distribution$sampleID = factor(angelman_distribution$sampleID,levels = c("OXN-018","OXN-068", "OXN-010", "OXN-069", "OXN-012", "OXN-006","OXN-071","OXN-043", "OXN-008","OXN-070", "OXN-020", "OXN-007","OXN-021"),ordered = T)

pairwise.wilcox.test(angelman_distribution$modification, angelman_distribution$Group,
                     p.adjust.method = "BH")

p2 = ggplot(angelman_distribution,aes(sampleID, modification,fill = Group)) + 
  geom_violin(scale = "width",alpha = 0.5) +
  geom_jitter(size = 0.5,height = 0.0, width = 0.07)+
  theme_classic() + ylab("Methylation modification (%)")+ xlab("")+
  theme(axis.text = element_text(size=18,colour = "black"), axis.title.y = element_text(size=18,colour = "black") , axis.text.x = element_text(size=18, angle=45, hjust = 1,colour = "black"))


p1 = ggplot(angelman_distribution,aes(Group,modification,fill = Group)) + 
  geom_violin(scale = "width",alpha = 0.5) +
  geom_jitter(size = 0.5,height = 0.0, width = 0.07)+
  theme_classic() + ylab("Methylation modification (%)")+ xlab("")+
  theme(axis.text = element_text(size=18,colour = "black"), axis.title.y = element_text(size=18,colour = "black") , axis.text.x = element_text(size=18, angle=45, hjust = 1,colour = "black"))

pdf("results/positive_cohort/angleman_oxn068.pdf")
print(p1)
print(p2)
dev.off()


# ------------------------------------------- SMA --------------------------------------

samples2plot = c(cohort_sample$SampleName[cohort_sample$Cohort=="Negative Cohort Samples"],c("OXN-007","OXN-008","OXN-010","OXN-068","OXN-069","OXN-070","OXN-071"))
oxn060sma = methylation_profiling_heatmap(methylationData = methylationData$sma_methylationData[methylationData$sma_methylationData$sampleID %in% samples2plot & methylationData$sma_methylationData$regionID=="smn1",], readdepth =5, scoreth = 98, pdffile="results/negative_cohort/OXN060_SMA_heatmap.pdf")

sma = methylation_profiling_heatmap(methylationData = methylationData$sma_methylationData[methylationData$sma_methylationData$sampleID %in% pos_cohort_sample$SampleName & methylationData$sma_methylationData$regionID=="smn1",], readdepth =5, scoreth = 98, pdffile="results/positive_cohort/SMA_heatmap.pdf")


sma_distribution = sma$methylation_data
pos = row.names(sma_distribution)
sma_distribution = data.table(sma_distribution)
sma_distribution$position = pos

sma_distribution = melt.data.table(sma_distribution,id.vars = c("position"),value.name = c("modification"),variable.name = "sampleID")
sma_distribution$Group = ifelse(sma_distribution$sampleID %in% c("OXN-068","OXN-069","OXN-010","OXN-007","OXN-008"),"SMA +ve",ifelse(sma_distribution$sampleID %in% c("OXN-070","OXN-071"),"SMA Carrier","Negative"))
sma_distribution$sampleID = factor(sma_distribution$sampleID,levels = c("OXN-007","OXN-008","OXN-010","OXN-068", "OXN-069", "OXN-071","OXN-070", "OXN-006","OXN-012","OXN-018", "OXN-020", "OXN-021","OXN-043"),ordered = T)

pairwise.wilcox.test(sma_distribution$modification, sma_distribution$Group,
                     p.adjust.method = "BH")

p2 = ggplot(sma_distribution,aes(sampleID, modification,fill = Group)) + 
  geom_violin(scale = "width",alpha = 0.5) +
  geom_jitter(size = 0.5,height = 0.0, width = 0.07)+
  theme_classic() + ylab("Methylation modification (%)")+ xlab("")+
  theme(axis.text = element_text(size=18,colour = "black"), axis.title.y = element_text(size=18,colour = "black") , axis.text.x = element_text(size=18, angle=45, hjust = 1,colour = "black"))


p1 = ggplot(sma_distribution,aes(Group,modification,fill = Group)) + 
  geom_violin(scale = "width",alpha = 0.5) +
  geom_jitter(size = 0.5,height = 0.0, width = 0.07)+
  theme_classic() + ylab("Methylation modification (%)")+ xlab("")+
  theme(axis.text = element_text(size=18,colour = "black"), axis.title.y = element_text(size=18,colour = "black") , axis.text.x = element_text(size=18, angle=45, hjust = 1,colour = "black"))

pdf("results/positive_cohort/sma_distribution.pdf")
print(p1)
print(p2)
dev.off()

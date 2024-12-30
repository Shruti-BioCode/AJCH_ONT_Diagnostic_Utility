#### This is a composite repository for analysing long read sequencing data from Oxford Nanopore Technologies in a heterogenous disease cohort for the study ####
**"Long read sequencing enhances pathogenic and novel variation discovery in patients with rare diseases"** 

It contains all the scripts used for analysis and plotting in the paper. 
1. Running "Epimarker" for methylation analysis - comparison with known 34 Mendelian Neurological disorders, Angelman Syndrome and SMA from methylation bedfiles. [rscript] (https://github.com/Shruti-BioCode/AJCH_ONT_Diagnostic_Utility/blob/main/Rscripts/Epimarker.R)
2. CNV & SVs analysis using "Funnel Down" approach from processed & annotated files of ClassifyCNV and AnnotSV. [rscript](https://github.com/Shruti-BioCode/AJCH_ONT_Diagnostic_Utility/blob/main/Rscripts/CNV_SV_get_candidate_stats.R)
3. Splicing variants analysis from spliceAi output files. [rscript](https://github.com/Shruti-BioCode/AJCH_ONT_Diagnostic_Utility/blob/main/Rscripts/spliceai_analysis.R)
4. Long read sepcific SNV analysis, post removing variants detected from whole exome sequencing. [rscript] (https://github.com/Shruti-BioCode/AJCH_ONT_Diagnostic_Utility/blob/main/Rscripts/ONT_specific_variant_filtering_analysis.R)

Below are the softwares used for the analysis
```
Softwares used sequencing and base calling:
- MinKnow distribution (version 22.05.7)
- Guppy (version 6.1.5)

Software for processing:
- Aligner - minimap2 (version 2.22-r1101).
- Variant/Methylation caller - Epi2Me workflow wf-human-variation (v1.2.0), CuteSV(v2.0.3)
- Annotation -  AnnotSV(v3.2.3), ClassifyCNV(1.1.1), SpliceAI (1.3.1)

All codes are compatible with R version 4.3.0 and genome build hg19
```

Follow Wiki for instruction on analysing your own dataset  
1. [Epimarker - methylation analysis](https://github.com/Shruti-BioCode/AJCH_ONT_Diagnostic_Utility/wiki/Epimarker-%E2%80%90-methylation-analysis)
2. [Funnel down - CNVs/SVs analysis](https://github.com/Shruti-BioCode/AJCH_ONT_Diagnostic_Utility/wiki/Funnel-Down-Filtering-%E2%80%90-CNVs---SVs)
3. [Splice Variant Analysis](https://github.com/Shruti-BioCode/AJCH_ONT_Diagnostic_Utility/wiki/Splice-Variant-Analysis)

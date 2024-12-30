#### This is a composite repository for analysing long read sequencing data from Oxford Nanopore Technologies in a heterogenous disease cohort for the study ####
**"Long read sequencing enhances pathogenic and novel variation discovery in patients with rare diseases"** 



It provides instructions for
1. Running "Epimarker" for methylation analysis - comparison with known 34 Mendelian Neurological disorders, Angelman Syndrome and SMA from methylation bedfiles.
2. Filtering, CNV & SVs from processed & annotated files of ClassifyCNV and AnnotSV
3. Filtering Splicing variants from spliceAi output files.
4. Filtering SNVs.

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

Follow Wiki for the different modules  
1. [Epimarker - methylation analysis](https://github.com/Shruti-BioCode/AJCH_ONT_Diagnostic_Utility/wiki/Epimarker-%E2%80%90-methylation-analysis)
2. [Funnel down - CNVs/SVs analysis](https://github.com/Shruti-BioCode/AJCH_ONT_Diagnostic_Utility/wiki/Funnel-Down-Filtering-%E2%80%90-CNVs---SVs)

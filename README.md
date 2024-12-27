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


### Epimarker - Methylation Analysis ###

Step 1: process methylation calls using the command below
``` Rscript
# read the list of files (Epi2Me -methyl output) from input directory. 
methylfiles = c("/path/to/your/bedmethyl/files", full.names = T)
samplenames = gsub(".*/|.methyl.cpg.bed.gz","",methylfiles)

# generate the methylation data for the episign probes, SMA & Angleman genomic region
methylationData = generate_methylation_data(methylfiles = methylfiles, samplenames = samplenames, reference = "hg19")
```

Step 2: Analyse the methylation data for Episigns, Angelman & SMA
```Rscript
methylation_results = run_complete_methylation_analysis(outdir = "/path/to/output/directory", methylationData = methylationData, cohort_samples=samplenames,remove_unstable_episignature = TRUE).
remove_unstable_episignature :  options (TRUE/FALSE), default: FALSE.
                                If selected TRUE, will remove Episigns CdLS, AUTS18 and RSTS that have been reported to have weak/overlapping methylation profile.
```

Output files:
```
|---MNDD
        |---Epimarker_Overview.pdf - pca & heatmap plot across all the samples under analysis
        |---epimarker_results.txt - samples matching episign profile
        |---complete_epimarker_results.txt - all profile comparisons compared along with the distance/correlation matrices values
        |---sample.pdf - pca & heatmap for the sample "sample" across all comparison
        |---sample.pdf - all profile comparisons compared along with the distance/correlation matrices values for the sample "sample"
|---SMA |---SMA_heatmap.pdf -  heatmap along with (%) bases modified across all the samples
|---ImD |---Angelman_heatmap.pdf -  heatmap along with (%) bases modified across all the samples
```


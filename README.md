

#  AAV vectors RNAseq preprocessing and analysis

paper_title: AAV vectors trigger DNA damage response-dependent pro-inflammatory signalling in human iPSC-derived CNS models and mouse brain  
doi: doi.org/10.1038/s41467-025-58778-3  
data_source: "GEO GSE253820"  
BioProject PRJNA1300358
SRA runs SRR27673943 to SRR27673956

This repository provides an end to end workflow for RNA seq analysis of human iPSC derived astrocytes comparing AAV vector transduced cells to control samples. It includes two parts:  

- Preprocessing in shell downloads SRA runs, builds sample FASTQs, runs QC and trimming, aligns to GRCh38, and produces gene level counts

- Analysis in R Markdown performs differential expression with DESeq2, PCA and volcano plots, a heatmap, and selected gene expression comparisons

## 1. Data sources

This task is based on publicly available sequencing data from a study of **\[AAV vectors trigger DNA damage response-dependent pro-inflammatory signalling in human iPSC-derived CNS models and mouse brain]**. The dataset includes multiple samples under different conditions (UT, AAV9) and was originally sequenced using **\[Illumina NovaSeq 2×150]**.
The subsampled and cleaned FASTQs are stored in `data/` and are used as the inputs for the workflow.

Public GEO GSE253820 and BioProject PRJNA1300358 human iPSC derived astrocytes paired end RNA seq

Groups and runs used

UT SRR27673943 SRR27673944 SRR27673945 SRR27673949  
AAV9 SRR27673950 SRR27673951 SRR27673956  

Reference

Genome: GRCh38, p13
Transcript annotation: gencode.v35.annotation.gtf
Alignment: STAR
Alignment index: hg38

---

## 2. How to download

Install the tools with conda and fetch FASTQs from SRA using prefetch and fasterq dump. See workflow scripts for exact commands.  


---

## 3. Pre-processing / subsampling

The preprocessing script creates a structured working directory data_pre_processing/ containing:

raw/ # downloaded and processed SRA FASTQs
fastq/ # concatenated FASTQs (one per sample)
aligned/ # hisat-aligned BAMs
counts/ # gene-level count matrix
qc/ # FastQC output
hisat2_index/ # reference index for alignment

You can skip raw FASTQ processing and go straight to analysis using data_RNAcounts/final_counts_symbols.tsv.


---

## 4. How the workflow works

The shell script executes the following steps. Save your script and run it from a terminal. It creates and uses ./data_pre_processing.

Step 1 Quality control

Purpose check base quality and GC content Tools FastQC and optional MultiQC Inputs FASTQs in fastq Outputs per sample HTML and zip in qc

Step 2 Trimming

Purpose remove adapters and low quality bases Tools Trimmomatic with TruSeq3 SE adapters Inputs FASTQs in fastq Outputs trimmed FASTQs in trimmed and trimming logs in logs

Step 3 Reference setup

Purpose obtain the GRCh38 HISAT2 transcriptome aware index Tools curl and tar Outputs index in hisat2_index/grch38_tran

Step 4 Alignment

Purpose align single end reads to GRCh38 and sort BAM Tools HISAT2 and samtools Inputs trimmed FASTQs in trimmed Outputs BAM and BAI in aligned plus logs in logs

Step 5 Post alignment QC

Purpose review mapping statistics Tools samtools

Step 6 Quantification

Purpose count reads per gene using the GENCODE GTF Tools featureCounts Inputs BAM files and GENCODE v35 GTF Outputs counts in counts and a final matrix data_pre_processing/counts/final_counts_symbols.tsv

Step 7 Analysis

Purpose differential expression and figures Tools R DESeq2 ggplot2 pheatmap Inputs featureCounts output and sample metadata Outputs PCA and volcano plots, QC tables, and selected gene comparisons

In R

run workflow/analysis_clean.Rmd to load data and perform the main analysis

Notes You can skip raw FASTQ processing and start from the counts matrix if you already have one.
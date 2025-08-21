# Bulk RNA-Seq Analysis

This repository contains a Bulk RNA-Seq expression analysis pipeline comparing LNCap and PC3 prostate cancer cells under normoxia and hypoxia conditions.

---

## Dataset Information

- GEO Accession: [GSE106305](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE106305)
- Samples:
    - LNCaP, Empty Vector, Normoxia (rep1, rep2)
        - [GSM3145509](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3145509)
        - [GSM3145510](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3145510)
    - LNCaP, Empty Vector, Hypoxia (rep1, rep2)
        - [GSM3145513](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3145513)
        - [GSM3145514](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3145514)
    - PC3, siCtrl, Normoxia (rep1, rep2)
        - [GSM3145517](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3145517)
        - [GSM3145518](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3145518)
    - PC3, siCtrl, Hypoxia (rep1, rep2)
        - [GSM3145521](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3145521)
        - [GSM3145522](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3145522)

---

## Steps

1. **Dataset Download** (SRA file → FASTQ)  
2. **Quality Control and Preprocessing** (FastQC + MultiQC + Trimmomatic)  
3. **Alignment and Post-Alignment QC** (HISAT2 + Samtools + Qualimap)  
4. **Quantification** (featureCounts → gene counts matrix)  
5. **Differential Expression Analysis** (DESeq2 in R)

Note: Scripts for each step are provided in the `scripts/` folder. 

---

## Requirements

### Primary analysis (Linux tools)

- Sra toolkit
- Fastqc
- Multiqc
- Trimmomatic
- HISAT2
- Samtools
- FeatureCounts
- Qualimap

### Differential Expression Analysis (R packages)

- R (v4.5)
- BiocManager
- DESeq2
- tidyverse
- dplyr
- tibble
- EnhancedVolcano
- RColorBrewer
- pheatmap

Note:
Install all dependencies using:  
```bash
./requirements.sh         # Linux tools  
Rscript install_packages.R # R packages

---

## Results

- Quality control reports (FastQC, MultiQC)

- Alignment stats (HISAT2, Qualimap)

- Gene counts matrix (FeatureCounts)

- Differential expression results (DESeq2)

- Plots: PCA, Volcano, Heatmaps

---

## References

The datasets and analysis workflow were adapted from: 
- [Guo H, Ci X, Ahmed M, et al. ONECUT2 is a driver of neuroendocrine prostate cancer.Nat Commun. 2019;10(1):278](https://pmc.ncbi.nlm.nih.gov/articles/PMC6336817/#Sec11)
- [Bulk RNA-sequencing pipeline and differential gene expression analysis](https://erilu.github.io/bulk-rnaseq-analysis/#Obtaining_raw_data_from_GEO)

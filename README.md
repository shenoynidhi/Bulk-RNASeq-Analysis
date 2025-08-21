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

1. Dataset Download
2. Quality control (QC) and Preprocessing
3. Alignment and Post-alignment QC
4. Quantification
5. Differential expression analysis

Note: Scripts for each step can be run from the scripts/ folder. 

---

## Requirements

### Primary analysis:

- Sra toolkit
- Fastqc
- Multiqc
- Trimmomatic
- HISAT2
- Samtools
- FeatureCounts
- Qualimap

### Differential Expression Analysis

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
All the requirements can be directly installed by running the `./requirements.sh` and `Rscript install_packages.R`

---

## Results
  

---

## References

The datasets and analysis workflow were adapted from: 
[Guo H, Ci X, Ahmed M, et al. ONECUT2 is a driver of neuroendocrine prostate cancer. Nat Commun. 2019;10(1):278. Published 2019 Jan 17. doi:10.1038/s41467-018-08133-6](https://pmc.ncbi.nlm.nih.gov/articles/PMC6336817/#Sec11) 

# Bulk RNA-Seq Analysis

This repository contains a Bulk RNA-Seq expression analysis pipeline comparing LNCap and PC3 prostate cancer cells under normoxia and hypoxia conditions.

---

##Dataset Information

- GEO Accession: [GSE106305](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE106305)
- Samples:
    - LNCaP, Empty Vector, Normoxia (rep1, rep2)
        -[GSM3145509](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3145509)
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

##Steps

1. Dataset Download: script add
2. Quality control (QC) and Preprocessing
3. Alignment and Post-alignment QC
4. Quantification
5. Differential expression analysis


---

##Requirements

###Primary analysis:

- sra toolkit
- Fastqc
- Multiqc
- Trimmomatic
- HISAT2
- Samtools
- featureCounts
- qualimap

###Differential Expression Analysis

- R v4.5
- Rstudio
- 

Note: 
All the requiremnts can be directly installed by running the `requirements.sh`
` ./requirements.sh`

---

##Results

- Count matrix (in `results/`)  
- Plots and DEGs list  

---

## References

The datasets and analysis workflow were adapted from:  

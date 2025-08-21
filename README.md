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

1. **Download requirements** 

```bash
# Install Linux tools
bash requirements.sh
# Install R packages
Rscript install_packages.R

2. **Download dataset**

```bash
# Download SRA files
prefetch SRR7179504
# Convert SRA to FASTQ
fastq-dump --outdir fastq --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip SRR7179504.sra
# Use instead to automate multiple SRA downloads
python3 scripts/fastq_download.py

3. **Quality Control and Preprocessing**

```bash
#FASTQC for quality control
fastqc fastq/*.fastq.gz -o fastqc_results/ --threads 8
#MULTIQC for single QC report
multiqc fastqc_results/ -o multiqc_report/
#Trimming (optional for the above listed samples)
java -jar Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 4 fastq/SRR7179504.fastq.gz fastq/SRR7179504_trimmed.fastq.gz TRAILING:10 -phred33
#Concatenate individual sample runs and rename as per GEO
cat fastq/SRR7179504_pass.fastq.gz fastq/SRR7179505_pass.fastq.gz fastq/SRR7179506_pass.fastq.gz fastq/SRR7179507_pass.fastq.gz > fastq/LNCAP_Normoxia_S1.fastq.gz
cat fastq/SRR7179508_pass.fastq.gz fastq/SRR7179509_pass.fastq.gz fastq/SRR7179510_pass.fastq.gz fastq/SRR7179511_pass.fastq.gz > fastq/LNCAP_Normoxia_S2.fastq.gz
cat fastq/SRR7179520_pass.fastq.gz fastq/SRR7179521_pass.fastq.gz fastq/SRR7179522_pass.fastq.gz fastq/SRR7179523_pass.fastq.gz > fastq/LNCAP_Hypoxia_S1.fastq.gz
cat fastq/SRR7179524_pass.fastq.gz fastq/SRR7179525_pass.fastq.gz fastq/SRR7179526_pass.fastq.gz fastq/SRR7179527_pass.fastq.gz > fastq/LNCAP_Hypoxia_S2.fastq.gz
mv fastq/SRR7179536_pass.fastq.gz fastq/PC3_Normoxia_S1.fastq.gz
mv fastq/SRR7179537_pass.fastq.gz fastq/PC3_Normoxia_S2.fastq.gz
mv fastq/SRR7179540_pass.fastq.gz fastq/PC3_Hypoxia_S1.fastq.gz
mv fastq/SRR7179541_pass.fastq.gz fastq/PC3_Hypoxia_S2.fastq.gz
#Clear downloaded SRA files 
rm -rf SRR*
   
4. **Alignment and Post-Alignment QC**

```bash
#Download and extract prebuilt human genome index for alignment using HISAT2
wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
tar -xvzf grch38_genome.tar.gz
#Alignment using HISAT2 and convert to .bam format
hisat2 -q -x grch38/genome -U ../LNCAP_Hypoxia_S1.fastq.gz | samtools sort -o LNCAP_Hypoxia_S1.bam
#Use instead if RAM less than 8GB, uses 2 threads (-@ 2) and 1GB RAM per thread (-m 1G) and ./tmp helps to store intermediate files and automatically clean up incase system's default tmp directory is fast filling (low space)
mkdir -p tmp
hisat2 -q -x grch38/genome -U ../LNCAP_Hypoxia_S1.fastq.gz | samtools sort -m 1G -@ 2 -T ./tmp -o LNCAP_Hypoxia_S1.bam
#Index .bam file
samtools index LNCAP_Hypoxia_S1.bam
#Use instead to automate multiple sample alignments
./scripts/alignment.sh
#Optional, quality Check for aligned .bam files (use --java-mem-size=6G for RAM less than 8GB) 
qualimap rnaseq -bam <bam file> -gtf gencode.v48.primary_assembly.annotation.gtf -outdir rnaseq_qc_results --java-mem-size=6G

5. **Quantification** 

```bash
#Download and extract GTF file for read count estimation
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.primary_assembly.annotation.gtf.gz
tar -xvzf gencode.v48.primary_assembly.annotation.gtf.gz
#Run featureCounts to get read counts from .bam file
mkdir -p quants
featureCounts -S 2 -a Homo_sapiens.GRCh38.114.gtf -o quants/featurecounts.txt tmp.bam
#Use instead for automating for each of the .bam files
./scripts/featurecounts.sh

6. **Differential Expression Analysis** (DESeq2 in R)
 

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

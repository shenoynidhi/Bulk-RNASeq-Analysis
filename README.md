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
```
2. **Download dataset**: Fetch SRA files under each of the above mentioned samples for analysis

```bash
# For a single SRA file
prefetch SRR7179504
# Convert SRA to FASTQ
fastq-dump --outdir fastq --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip SRR7179504/SRR7179504.sra
# Optional: faster alternative
fasterq-dump SRR7179504/SRR7179504.sra --split-files -O fastq --threads 8
pigz -p 8 fastq/SRR7179504.fastq"
# To download all SRA files in one go, use the provided script
python scripts/fastq_download.py
```
3. **Quality Control (QC) and Preprocessing**: Check the quality of the reads followed by trimming (optional) and renaming for analysis

```bash
# FASTQC for quality control
fastqc fastq/*.fastq.gz -o fastqc_results/ --threads 8
# MULTIQC for single QC report
multiqc fastqc_results/ -o multiqc_report/
# Trimming (optional for the above listed samples)
java -jar Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 4 fastq/SRR7179504.fastq.gz fastq/SRR7179504_trimmed.fastq.gz TRAILING:10 -phred33
# Concatenate individual sample runs and rename as per GEO (_pass is appended by fastq-dump, _1 appended by fasterq-dump on FASTQ file generation)
cat fastq/SRR7179504_pass.fastq.gz fastq/SRR7179505_pass.fastq.gz fastq/SRR7179506_pass.fastq.gz fastq/SRR7179507_pass.fastq.gz > fastq/LNCAP_Normoxia_S1.fastq.gz
cat fastq/SRR7179508_pass.fastq.gz fastq/SRR7179509_pass.fastq.gz fastq/SRR7179510_pass.fastq.gz fastq/SRR7179511_pass.fastq.gz > fastq/LNCAP_Normoxia_S2.fastq.gz
cat fastq/SRR7179520_pass.fastq.gz fastq/SRR7179521_pass.fastq.gz fastq/SRR7179522_pass.fastq.gz fastq/SRR7179523_pass.fastq.gz > fastq/LNCAP_Hypoxia_S1.fastq.gz
cat fastq/SRR7179524_pass.fastq.gz fastq/SRR7179525_pass.fastq.gz fastq/SRR7179526_pass.fastq.gz fastq/SRR7179527_pass.fastq.gz > fastq/LNCAP_Hypoxia_S2.fastq.gz
mv fastq/SRR7179536_pass.fastq.gz fastq/PC3_Normoxia_S1.fastq.gz
mv fastq/SRR7179537_pass.fastq.gz fastq/PC3_Normoxia_S2.fastq.gz
mv fastq/SRR7179540_pass.fastq.gz fastq/PC3_Hypoxia_S1.fastq.gz
mv fastq/SRR7179541_pass.fastq.gz fastq/PC3_Hypoxia_S2.fastq.gz
# Clear downloaded SRA files 
rm -rf SRR*
```   
4. **Alignment and Post-Alignment QC**: Map the reads to human reference genome (GRCh38) followed by indexing and additional quality check (optional)

```bash
# Download and extract prebuilt human genome index for alignment
wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
tar -xvzf grch38_genome.tar.gz
# Alignment using HISAT2 and convert to .bam format
mkdir -p alignment
cd alignment/
hisat2 -q -x grch38/genome -U fastq/LNCAP_Hypoxia_S1.fastq.gz | samtools sort -o LNCAP_Hypoxia_S1.bam
# Optional Alternative: For RAM less than 8GB, uses 2 threads (-@ 2) and 1GB RAM per thread (-m 1G) and ./tmp helps to store intermediate files and automatically clean up incase system's default tmp directory is at low space
mkdir -p tmp
hisat2 -q -x grch38/genome -U fastq/LNCAP_Hypoxia_S1.fastq.gz | samtools sort -m 1G -@ 2 -T ./tmp -o LNCAP_Hypoxia_S1.bam
# Index .bam file
samtools index LNCAP_Hypoxia_S1.bam
# Use for multiple files
./scripts/alignment.sh
# Optional: Quality Check for aligned .bam files (use --java-mem-size=6G for RAM less than 8GB) 
qualimap rnaseq -bam LNCAP_Hypoxia_S1.bam -gtf gencode.v48.primary_assembly.annotation.gtf -outdir rnaseq_qc_results --java-mem-size=6G
```
5. **Quantification**: Estimate read counts for each gene/transcripts

```bash
#Download and extract GTF file for read count estimation
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.primary_assembly.annotation.gtf.gz
tar -xvzf gencode.v48.primary_assembly.annotation.gtf.gz
#Run featureCounts to get read counts from .bam file
mkdir -p quants
featureCounts -S 2 -a gencode.v48.primary_assembly.annotation.gtf -o quants/LNCAP_Hypoxia_S1_featurecounts.txt LNCAP_Hypoxia_S1.bam
#To process all .bam files in one go
./scripts/counts.sh
```
6. **Differential Expression Analysis (DEA)**: Identify and study genes that exhibit significantly different expression levels between the two conditions

    - **WorkFlow Overview** 

        1. **Initial DESeq analysis** (01_initial_DESeq.R)

            - Load count and annotation files – Input your raw RNA-seq metadata and count matrix generated in the previous step and gene annotation to start the analysis.

            - Filter genes by type and zero counts – Remove genes that are unlikely to be informative (non-coding or mostly zero across samples).

            - Run DESeq normalization – Normalize counts to account for differences in sequencing depth and prepare data for downstream analyses.

        2. **Exploratory visualizations** (02_visualizations.R)

            - PCA and sample distance heatmaps – Explore how samples cluster and identify any batch effects or outliers.

            - Density plots for raw and normalized counts – Compare distribution of counts before and after normalization.

            - Heatmap of top variable genes – Identify genes with the most variation across all samples, highlighting biologically meaningful differences.

        3. **Cell line-specific DEA** (03_cell_line_DE.R)

            - Subset data for each cell line (LNCaP or PC3) – Analyze each cell line separately to account for major differences in gene expression.

            - Differential expression analysis – Identify genes significantly up- or down-regulated between conditions (Hypoxia vs Normoxia).

            - Volcano plots and heatmaps for top DE genes – Visualize significant DE genes and their expression patterns.

            - Gene set enrichment analysis (GSEA) – Discover pathways or functional categories enriched in the DE genes.

            - Optional: boxplots for specific genes – Inspect expression of individual genes of interest across conditions.

- Note: `GRCh38annotation.csv` and `h.all.v2025.1.Hs.symbols.gmt` are provided in the repository
---

## Requirements

### Primary analysis (Linux tools)

- sra-toolkit
- pigz
- FASTQC
- MULTIQC
- Trimmomatic
- HISAT2
- samtools
- subread
- qualimap

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
- org.Hs.eg.db
- ggplot2
- ggrepel
- readr
- stats
- fgsea
- stringr
- stats

- Note: Qualimap installation should be done via conda or manually, refer [Documentation](http://qualimap.conesalab.org/), [Conda installation](https://anaconda.org/bioconda/qualimap)
---

## Results

- Raw counts matrix (`GSE106305_counts_matrix`)
- Annotated Count matrix (`biotype_count_matrix`)
- Filtered Count matrix (`filtered_biotype_nozero_count_matrix`)
- Normalized Count matrix (`normalized_counts`)
- Sample Variability Visualizations: PCA, Distance plot, Variable gene heatmap and Density plot for raw vs normalized counts
- LNCAP Results (Hypoxia vs Normoxia): DEGs list, Volcano Plot, Top 20 DE genes heatmap and GSEA plots
- PC3 Results (Hypoxia vs Normoxia): DEGs list, Volcano Plot, Top 20 DE genes heatmap and GSEA plots
- Boxplot for comparison of normalized counts of IGFBP1 gene across samples and conditions  

---

## References

The datasets and analysis workflow were adapted from: 
- [Guo H, Ci X, Ahmed M, et al. ONECUT2 is a driver of neuroendocrine prostate cancer.Nat Commun. 2019;10(1):278](https://pmc.ncbi.nlm.nih.gov/articles/PMC6336817/#Sec11)
- [Bulk RNA-sequencing pipeline and differential gene expression analysis](https://erilu.github.io/bulk-rnaseq-analysis/#Obtaining_raw_data_from_GEO)

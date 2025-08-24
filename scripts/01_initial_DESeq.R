#Load required libraries
library(DESeq2)
library(tidyverse)
library(dplyr)
library(tibble)

#Loading required data files

raw_counts <- read.csv("./fastq/alignment/quants/GSE106305_counts_matrix.csv", header = TRUE, row.names = "Geneid", stringsAsFactors = FALSE)
head(raw_counts)
raw_counts <- raw_counts[,sort(colnames(raw_counts))]
colSums(raw_counts)

condition <- c(rep("LNCAP_Hypoxia", 2), rep("LNCAP_Normoxia", 2), rep("PC3_Hypoxia", 2), rep("PC3_Normoxia", 2))
print(condition)

my_colData <- as.data.frame(condition)
rownames(my_colData) <- colnames(raw_counts)
head(my_colData)

#creating dds object
dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = my_colData,
                              design = ~condition)
dds
head(counts(dds))
dim(counts(dds))

count_matrix <- counts(dds)
dim(count_matrix)

#Annotation file can be downloaded from Ensembl BioMart (http://uswest.ensembl.org/biomart/martview/)
#Loading annotation file and filtering for specific biotypes
annotation_file <- "GRCh38annotation.csv"
annotation <- read.csv(annotation_file, header=TRUE, stringsAsFactors = FALSE)

counts_gse <- read.csv("./fastq/alignment/quants/GSE106305_counts_matrix.csv",
                     header = TRUE,
                     stringsAsFactors = FALSE)
#remove version numbers on gene from both dfs
counts_gse$Geneid <- sub("\\..*$", "", counts_gse$Geneid)
annotation$Geneid <- sub("\\..*$", "", annotation$Geneid)
annotated_counts <- left_join(counts_gse, annotation, by = "Geneid") %>%
  select(Geneid, Genesymbol, Genetype, 
         LNCAP_Hypoxia_S1, LNCAP_Hypoxia_S2, LNCAP_Normoxia_S1, LNCAP_Normoxia_S2, 
         PC3_Hypoxia_S1, PC3_Hypoxia_S2, PC3_Normoxia_S1, PC3_Normoxia_S2)

biotypes_to_keep <- c("protein_coding", "IG_J_gene", "IG_V_gene", "IG_C_gene", "IG_D_gene", "TR_D_gene", "TR_C_gene", "TR_V_gene", "TR_J_gene")

filtered_counts <- annotated_counts %>%
  filter(Genetype %in% biotypes_to_keep)
head(filtered_counts, n = 3)

write.csv(filtered_counts, file = "biotype_count_matrix.csv", sep = ",", row.names = FALSE)
zero_counts1 <- rowSums(filtered_counts[, 4:11] == 0)
zero_summary2 <- table(zero_counts1)
print(zero_summary2)

#Filtering genes based on zero counts
#keep genes that have non zero counts in 2 samples 
keep_genes <- zero_counts1 < 7
filtered_counts_nozero <- filtered_counts[keep_genes, ]
cat("Number of genes after filtering (zeros in <7 samples):", nrow(filtered_counts_nozero), "\n")

#display how many genes fall in 0-6 zero count samples after filtering
new_zero_counts <- rowSums(filtered_counts_nozero[, 4:11] == 0)
cat("New zero counts distribution:\n")
print(table(new_zero_counts))

output_file <- "filtered_biotype_nozero_count_matrix.csv"
write.csv(filtered_counts_nozero, file = output_file, sep = ",", row.names = FALSE)

head(filtered_counts_nozero, n = 3)
#remove version numbers from dds
rownames(dds) <- sub("\\..*$", "", rownames(dds))
#filter dds object to keep only genes that we have saved after filtering zero counts
dds_filtered <- dds[rownames(dds) %in% filtered_counts_nozero$Geneid, ]
cat("Dimensions of filtered DESeqDataSet:", dim(dds_filtered), "\n")

#run DESeq analysis and normalization 
dds <- DESeq(dds_filtered)
dds
normalized_counts <- counts(dds, normalized = T)
normalized_counts_df <- as.data.frame(normalized_counts)
write.csv(normalized_counts_df, file = "normalized_counts.csv", row.names = TRUE)

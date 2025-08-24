#Differential expression analysis for each cell line separately
#based on exploratory analysis, we can see that samples cluster by cell line (PC3 vs LNCAP) and then by condition (Hypoxia vs Normoxia) so we run DEA separately for each cell line
#for lncap
dds_lncap <- dds[, grepl("LNCAP", colnames(dds))]
dds_lncap
dds_lncap$condition <- droplevels(dds_lncap$condition)
dds_lncap$condition <- relevel(dds_lncap$condition, ref = "LNCAP_Normoxia")
dds_lncap <- DESeq(dds_lncap)
#extract results for lncap DEA
res_lncap <- results(dds_lncap, contrast = c("condition", "LNCAP_Hypoxia", "LNCAP_Normoxia"))
res_lncap
summary(res_lncap)

reslncapOrdered <- res_lncap[order(res_lncap$padj), ]

sum(reslncapOrdered$padj < 0.05, na.rm = TRUE)
head(reslncapOrdered)
summary(reslncapOrdered)
write.csv(as.data.frame(reslncapOrdered), file = "DEGs_lncap.csv")

#volcano plot for lncap
res_df <- as.data.frame(reslncapOrdered)
res_df <- na.omit(res_df)
res_df$gene <- rownames(res_df)

res_df$regulation <- "Not Significant"
res_df$regulation[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Upregulated"
res_df$regulation[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Downregulated"

qp <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = regulation)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Upregulated" = "#FEA405", 
                                "Downregulated" = "purple", 
                                "Not Significant" = "gray")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  annotate("text", x = min(res_df$log2FoldChange), y = -log2(0.05) + 0.5,
           label = "padj = 0.05", hjust = 0, size = 3) +
  theme_minimal() +
  labs(title = "Volcano Plot", 
       x = "Log2 Fold Change", 
       y = "-Log10 Adjusted P-Value") +
  theme(plot.title = element_text(hjust = 0.5))

v_plot <- "vp_lncap.png"
ggsave(v_plot, plot = qp, width = 8, height = 6, dpi = 300)

#heatmap for top 20 DE genes in lncap
library(pheatmap)
library(RColorBrewer)
library(DESeq2)
library(dplyr)

DE_gene_heatmap <- function(res_lncap, count_matrix, padj_cutoff = 0.0001, ngenes = 20) {
  # Generate the color palette
  brewer_palette <- "RdBu"
  ramp <- colorRampPalette(RColorBrewer::brewer.pal(11, brewer_palette))
  mr <- ramp(256)[256:1]  # Reversed palette (blue to red)

  # Convert DESeqResults to data frame and get significant genes
  significant_genes <- as.data.frame(res_lncap) %>%
    filter(padj < padj_cutoff) %>%
    arrange(desc(log2FoldChange)) %>%
    head(ngenes)

  # Extract count data for significant genes
  heatmap_values <- count_matrix[rownames(significant_genes), ]  # Use row names (Ensembl IDs)
  
  #map ensembl IDs in res_lncap to gene symbols
  gene_symbols <- annotation$Genesymbol[match(rownames(significant_genes), annotation$Geneid)]
  rownames(heatmap_values) <- gene_symbols
  
  # Scale rows for heatmap (z-score normalization)
  heatmap_values <- t(scale(t(heatmap_values)))

  # Plot the heatmap using pheatmap
  p <- pheatmap::pheatmap(heatmap_values,
                          color = mr,
                          scale = "none",  # Already scaled
                          cluster_rows = TRUE,
                          cluster_cols = TRUE,
                          fontsize_col = 10,
                          fontsize_row = max(6, 200/ngenes),  # Minimum font size of 6
                          border_color = NA,
                          main = paste("Top", ngenes, "DE Genes (padj <", padj_cutoff, ")"))

  # Return the pheatmap object
  return(invisible(p))
}

count_matrix <- assay(dds_lncap)  # Replace with your count matrix if different
png("de_gene_heatmap.png",
    width = 1000, height = 1700, res = 200)  # Adjust dimensions and resolution
d <- DE_gene_heatmap(res, count_matrix, padj_cutoff = 0.001, ngenes = 30)
dev.off()

#Gene set enrichment analysis for Lncap DEGs
#Downloaded hallmark gene sets from MSigDB (https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp)
library(fgsea)
hallmark_pathway <- gmtPathways("h.all.v2025.1.Hs.symbols.gmt")
head(names(hallmark_pathway))
head(hallmark_pathway$HALLMARK_HYPOXIA, 20)
library(org.Hs.eg.db)

res_lncap$symbol <- mapIds(org.Hs.eg.db,
                           keys = rownames(res_lncap),
                           column = "SYMBOL",
                           keytype = "ENSEMBL",
                           multiVals = "first")
res_clean <- res_lncap[!is.na(res_lncap$symbol), ]        # remove NAs
res_clean <- res_clean[!duplicated(res_clean$symbol), ]   # remove duplicates

# named numeric vector: names = symbols, values = log2FC
lncap_ranked_list <- res_clean$log2FoldChange
names(lncap_ranked_list) <- res_clean$symbol

# sort decreasing (most up first)
lncap_ranked_list <- sort(lncap_ranked_list, decreasing = TRUE)

fgsea_results <- fgsea(pathways = hallmark_pathway,
                       stats = lncap_ranked_list,
                       minSize = 15,
                       maxSize = 500,
                       nperm = 1000)
fgsea_results_ordered <- fgsea_results[order(-NES)]
head(fgsea_results_ordered[, .(pathway, padj, NES)])
plotEnrichment(hallmark_pathway[["HALLMARK_HYPOXIA"]], lncap_ranked_list)

#waterfall plot for fgsea results
waterfall_plot <- function (fsgea_results, graph_title) {
  fgsea_results %>% 
    mutate(short_name = str_split_fixed(pathway, "_",2)[,2])%>% # removes 'HALLMARK_' from the pathway title 
    ggplot( aes(reorder(short_name,NES), NES)) +
      geom_bar(stat= "identity", aes(fill = padj<0.05))+
      coord_flip()+
      labs(x = "Hallmark Pathway", y = "Normalized Enrichment Score", title = graph_title)+
      theme(axis.text.y = element_text(size = 7), 
            plot.title = element_text(hjust = 1))
}
library(stringr)
waterfall_plot(fgsea_results, "Hallmark pathways altered by hypoxia in LNCaP cells")
ggsave("GSEA_hallmark_lncap.png", p, width = 8, height = 6, dpi = 300)

#GSEA using Reactome pathways for lncap DEGs
res_lncap <- read.csv("DEGs_lncap.csv", row.names = 1)
head(res_lncap)

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(stats)

ncbi_list <- clusterProfiler::bitr(
  geneID = rownames(res_lncap),        # use Ensembl IDs from row names
  fromType = "ENSEMBL",          
  toType = "ENTREZID", 
  OrgDb = org.Hs.eg.db
)

res_lncap$ENSEMBL <- rownames(res_lncap)

res_mapped <- res_lncap %>%
  left_join(ncbi_list, by = "ENSEMBL") %>%
  filter(!is.na(ENTREZID)) %>%
  distinct(ENTREZID, .keep_all = TRUE)

ngenes <- res_mapped$log2FoldChange
names(ngenes) <- res_mapped$ENTREZID
ngenes <- sort(ngenes, decreasing = TRUE)

library(ReactomePA)
enp_gsea <- gsePathway(
  ngenes,
  organism = "human",
  #pvalueCutoff = 0.05,
  verbose = FALSE
)

pathways <- enp_gsea@result
pathways <- pathways[order(pathways$p.adjust), ]  # Sort by FDR (adjusted p-value)
top_pathways <- pathways[order(abs(pathways$NES), decreasing = TRUE), ]  # Sort by NES

library(dplyr)
library(forcats)

top20 <- top_pathways[1:20, ] %>%
  mutate(Description = fct_reorder(Description, NES))  # Reorder factor for y-axis

library(ggplot2)

r1 <- ggplot(top20, aes(x = NES,
                        y = Description,
                        color = p.adjust,
                        size = setSize)) +
  geom_point(alpha = 0.9) +
  scale_color_gradient(low = "#0072B2", high = "#D55E00", name = "FDR (p.adjust)") +
  scale_size(range = c(3, 10), name = "Gene Set Size") +
  labs(
    title = "Top 20 Enriched Pathways",
    subtitle = "Gene Set Enrichment Analysis (GSEA)",
    x = "Normalized Enrichment Score (NES)",
    y = NULL,
    caption = "Data source: clusterProfiler::gsePathway"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y = element_text(size = 5),
    axis.text.x = element_text(size = 5),
    plot.title = element_text(face = "bold", size = 9),
    plot.subtitle = element_text(size = 7),
    legend.position = "right"
  )
ggsave("GSEA_reactome_lncap.png", r1, width = 8, height = 6, dpi = 300)

#Results can similarly be generated for PC3 cell line by subsetting dds object for PC3 samples and running DEA, GSEA and visualizations as shown above for LNCAP

#Additionally, boxplot for comparison of normalized counts of a specific gene (IGFBP1) in different samples and conditions
annotation <- read.csv("GRCh38annotation.csv", header = T, stringsAsFactors = F)
annotation$Geneid <- sub("\\..*$", "", annotation$Geneid)
normalized_data <- counts(dds, normalized = T) 
condition <- dds@colData$condition
ensembl_id <- annotation$Geneid[which(annotation$Genesymbol == "IGFBP1")]
expression <- normalized_data[ensembl_id,]
gene_name <- annotation$Genesymbol[which(annotation$Geneid == ensembl_id)]
gene_tib <- tibble(condition = condition, expression = expression)
ggplot(gene_tib, aes(x = condition, y = expression))+
geom_boxplot(outlier.size = NULL)+
geom_point()+
labs (title = paste0("Expression of ", gene_name, " - ", ensembl_id), x = "group", y = paste0("Normalized     expression"))+
theme(axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11))
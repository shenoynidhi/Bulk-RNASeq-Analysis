#Visualizations for sample variability

#Variance stabilizing transformation
vsd <- vst(dds_filtered, blind = TRUE)  # blind=TRUE for exploratory PCA
plot_PCA = function (vsd.obj) {
  pcaData <- plotPCA(vsd.obj,  intgroup = c("condition"), returnData = T)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  ggplot(pcaData, aes(PC1, PC2, color=condition)) +
    geom_point(size=3) +
    labs(x = paste0("PC1: ",percentVar[1],"% variance"),
         y = paste0("PC2: ",percentVar[2],"% variance"),
         title = "PCA Plot colored by condition") +
    ggrepel::geom_text_repel(aes(label = name), color = "black")
}

png(filename = "pcab.png", 
    width = 2000, height = 2000, res = 300)  # adjust width/height as needed
plot_PCA(vsd)
dev.off()

plotDists = function (vsd.obj) {
  sampleDists <- dist(t(assay(vsd.obj)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(vsd.obj$condition)
  colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(55)
  pheatmap::pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists,  clustering_distance_cols = sampleDists, col = colors, fontsize_row = 4, fontsize_col = 4, fontsize_legend = 4, fontsize = 4)
}
png(filename = "distance_plot.png", width = 1000, height = 900, res = 300)  # adjust width/height as needed
plotDists(vsd)
dev.off()

#heatmap for top 500 most variable genes across all samples
variable_gene_heatmap <- function (vsd.obj, num_genes = 500, annotation, title = "Heatmap") {
  brewer_palette <- "RdBu"
  ramp <- colorRampPalette( RColorBrewer::brewer.pal(11, brewer_palette))
  mr <- ramp(256)[256:1]
  stabilized_counts <- assay(vsd.obj)
  row_variances <- rowVars(stabilized_counts)
  top_variable_genes <- stabilized_counts[order(row_variances, decreasing=T)[1:num_genes],]
  top_variable_genes <- top_variable_genes - rowMeans(top_variable_genes, na.rm=T)
  gene_names <- annotation$Genesymbol[match(rownames(top_variable_genes), annotation$Geneid)]
  rownames(top_variable_genes) <- gene_names
  coldata <- as.data.frame(vsd.obj@colData)
  coldata$sizeFactor <- NULL
  pheatmap::pheatmap(top_variable_genes, color = mr, annotation_col = coldata, fontsize_col = 8, fontsize_row = 250/num_genes, border_color = NA, main = title)
}

png(filename = "variable_gene_heatmap.png", 
    width = 1000, height = 1000, res = 300)  # adjust width/height as needed
variable_gene_heatmap(vsd, num_genes = 40, annotation = annotation)
dev.off()

#comparison of density plots for raw and vst normalized counts for all samples
raw_counts <- assay(dds)
vst_counts <- assay(vsd)

png("density_plots_raw_vst.png",
    width = 4000, height = 4000, res = 300)  # Adjust width, height (pixels), and resolution (dpi)

par(mfrow = c(4, 4), mar = c(3, 3, 2, 1))  # mar adjusts margins (bottom, left, top, right)

for (i in 1:8) {
  # Raw counts density
  plot(density(raw_counts[, i]),
       main = paste("Raw - Sample", colnames(raw_counts)[i]),
       xlab = "Expression",
       col = "red",
       lwd = 2,
       ylim = c(0, max(sapply(1:8, function(j) max(density(raw_counts[, j])$y, na.rm = TRUE))))  # Uniform y-axis
       )
  
  # VST counts density (next panel)
  plot(density(vst_counts[, i]),
       main = paste("VST - Sample", colnames(vst_counts)[i]),
       xlab = "Expression",
       col = "blue",
       lwd = 2,
       ylim = c(0, max(sapply(1:8, function(j) max(density(vst_counts[, j])$y, na.rm = TRUE))))  # Uniform y-axis
       )
}
dev.off()

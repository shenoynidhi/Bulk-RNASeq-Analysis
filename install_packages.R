# Ensure BiocManager is present
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# CRAN packages
cran_pkgs <- c("tidyverse", "dplyr", "tibble",
               "EnhancedVolcano", "RColorBrewer", "pheatmap")

for (pkg in cran_pkgs) {
    if (!require(pkg, character.only = TRUE)) {
        install.packages(pkg)
    }
}

# Bioconductor packages
bioc_pkgs <- c("DESeq2")

for (pkg in bioc_pkgs) {
    if (!require(pkg, character.only = TRUE)) {
        BiocManager::install(pkg, ask = FALSE, update = TRUE)
    }
}

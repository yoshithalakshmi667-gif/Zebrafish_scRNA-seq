# 05_normalization.R
# Normalize zebrafish scRNA-seq data

library(Seurat)

# 1. Load QC-filtered Seurat object

seurat_obj <- readRDS("zebrafish_seurat_QC_filtered.rds")

cat("Loaded filtered Seurat object:\n")
print(seurat_obj)

# 2. Normalize RNA expression

seurat_obj <- NormalizeData(
    seurat_obj,
    normalization.method = "LogNormalize",
    scale.factor = 10000
)

cat("\nNormalization completed successfully.\n")

# 3. Check RNA assay layers

cat("\nRNA assay information:\n")
print(seurat_obj[["RNA"]])

# 4. Save normalized object

saveRDS(
    seurat_obj,
    file = "zebrafish_seurat_normalized.rds"
)

cat("\nNormalized Seurat object saved successfully.\n")

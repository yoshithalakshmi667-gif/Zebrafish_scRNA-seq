# 07_PCA.R
# Principal Component Analysis

library(Seurat)
# 1. Load variable-gene Seurat object

seurat_obj <- readRDS("zebrafish_seurat_variable_genes.rds")

cat("Loaded Seurat object:\n")
print(seurat_obj)

# 2. Scale the data

seurat_obj <- ScaleData(
    seurat_obj,
    features = VariableFeatures(seurat_obj)
)

cat("\nScaling completed successfully.\n")

# 3. Run PCA

seurat_obj <- RunPCA(
    seurat_obj,
    features = VariableFeatures(seurat_obj),
    npcs = 50
)

cat("\nPCA completed successfully.\n")

# 4. Display PCA information

print(
    seurat_obj[["pca"]]
)

# 5. PCA elbow plot

pdf(
    "PCA_elbow_plot.pdf",
    width = 8,
    height = 6
)

print(
    ElbowPlot(
        seurat_obj,
        ndims = 50
    )
)

dev.off()

cat("\nPCA elbow plot saved as PCA_elbow_plot.pdf\n")

# 6. Save PCA object

saveRDS(
    seurat_obj,
    file = "zebrafish_seurat_PCA.rds"
)

cat("\nPCA Seurat object saved successfully.\n")

library(Seurat)

# Load Seurat object after PCA
seurat_obj <- readRDS("zebrafish_seurat_PCA.rds")

# Use PCs 1-15 based on the PCA elbow plot
seurat_obj <- FindNeighbors(
    seurat_obj,
    dims = 1:15
)

# Perform clustering
seurat_obj <- FindClusters(
    seurat_obj,
    resolution = 0.5
)

# Display clustering information
cat("\nNumber of clusters:\n")
print(length(levels(Idents(seurat_obj))))

cat("\nCluster sizes:\n")
print(table(Idents(seurat_obj)))

# Save clustered object
saveRDS(
    seurat_obj,
    file = "zebrafish_seurat_clustered.rds"
)

cat("\nClustering completed successfully.\n")
cat("Seurat object saved as zebrafish_seurat_clustered.rds\n")

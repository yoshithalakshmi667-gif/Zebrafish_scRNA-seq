library(Seurat)

# Load clustered Seurat object
seurat_obj <- readRDS("zebrafish_seurat_clustered.rds")

# Run UMAP using PCs 1-15
seurat_obj <- RunUMAP(
    seurat_obj,
    dims = 1:15
)

# Save UMAP plot colored by cluster
pdf(
    "UMAP_clusters.pdf",
    width = 9,
    height = 7
)

print(
    DimPlot(
        seurat_obj,
        reduction = "umap",
        group.by = "seurat_clusters",
        label = TRUE,
        repel = TRUE
    )
)

dev.off()

# Save updated object
saveRDS(
    seurat_obj,
    file = "zebrafish_seurat_UMAP.rds"
)

cat("\nUMAP completed successfully.\n")
cat("UMAP plot saved as UMAP_clusters.pdf\n")

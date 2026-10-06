# Verify final cell-type annotations

library(Seurat)

cat("Loading annotated Seurat object...\n")

seurat_obj <- readRDS("zebrafish_seurat_annotated.rds")

cat("Object loaded successfully.\n\n")

# Number of cells
cat("Number of cells:", ncol(seurat_obj), "\n")

# Number of clusters
cat("Number of clusters:", length(unique(Idents(seurat_obj))), "\n")

# Check cell-type annotations
cat("\nCell-type counts:\n")
print(table(seurat_obj$cell_type))

# Check for missing annotations
cat("\nNumber of cells with missing cell type:",
    sum(is.na(seurat_obj$cell_type)), "\n")

# Check unique cell types
cat("\nNumber of unique cell types:",
    length(unique(seurat_obj$cell_type)), "\n")

# Show cluster-to-cell-type mapping
cat("\nCluster-to-cell-type mapping:\n")

annotation_check <- unique(
  seurat_obj@meta.data[, c("seurat_clusters", "cell_type")]
)

annotation_check <- annotation_check[
  order(annotation_check$seurat_clusters),
]

print(annotation_check)

# Save verification table
write.csv(
  annotation_check,
  "verified_cluster_annotations.csv",
  row.names = FALSE
)

cat("\nVerification completed successfully.\n")
cat("Saved: verified_cluster_annotations.csv\n")


library(Seurat)

cat("Loading annotated Seurat object...\n")

seurat_obj <- readRDS("zebrafish_seurat_annotated.rds")

cat("Object loaded successfully.\n\n")

# Show metadata columns
cat("Metadata columns:\n")
print(colnames(seurat_obj@meta.data))

# Check sample_name
cat("\nUnique sample_name values:\n")
print(unique(seurat_obj$sample_name))

cat("\nNumber of cells in each sample:\n")
print(table(seurat_obj$sample_name))

# Check ClusterNames
cat("\nOriginal atlas ClusterNames:\n")
print(length(unique(seurat_obj$ClusterNames)))

# Check original biological annotation
cat("\nOriginal atlas Cluster labels:\n")
print(length(unique(seurat_obj$Cluster)))

cat("\nDevelopmental-stage check completed.\n")

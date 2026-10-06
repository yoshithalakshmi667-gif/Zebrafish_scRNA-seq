# 04_QC_filtering.R
# Filter low-quality cells from zebrafish scRNA-seq data

library(Seurat)
library(Matrix)
# 1. Load QC object

seurat_obj <- readRDS("zebrafish_seurat_QC_before_filtering.rds")

cat("Cells before filtering:", ncol(seurat_obj), "\n")

# 2. Define QC thresholds

min_genes <- 500
max_genes <- 6000
max_mito <- 10

cat("\nQC thresholds:\n")
cat("Minimum genes:", min_genes, "\n")
cat("Maximum genes:", max_genes, "\n")
cat("Maximum mitochondrial %:", max_mito, "\n")

# 3. Count cells that will be removed

low_genes <- seurat_obj$nFeature_RNA < min_genes
high_genes <- seurat_obj$nFeature_RNA > max_genes
high_mito <- seurat_obj$percent.mt > max_mito

cat("\nCells with too few genes:", sum(low_genes), "\n")
cat("Cells with too many genes:", sum(high_genes), "\n")
cat("Cells with high mitochondrial %:", sum(high_mito), "\n")

# 4. Apply filtering

seurat_obj_filtered <- subset(
    seurat_obj,
    subset =
        nFeature_RNA >= min_genes &
        nFeature_RNA <= max_genes &
        percent.mt <= max_mito
)

# 5. Report filtering results

cat("\nCells after filtering:", ncol(seurat_obj_filtered), "\n")

cells_removed <- ncol(seurat_obj) - ncol(seurat_obj_filtered)

cat("Cells removed:", cells_removed, "\n")

cat(
    "Percentage of cells retained:",
    round(100 * ncol(seurat_obj_filtered) / ncol(seurat_obj), 2),
    "%\n"
)

# 6. Save filtered object

saveRDS(
    seurat_obj_filtered,
    file = "zebrafish_seurat_QC_filtered.rds"
)

cat("\nFiltered Seurat object saved successfully.\n")

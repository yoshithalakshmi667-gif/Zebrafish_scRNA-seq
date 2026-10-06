library(Seurat)
library(ggplot2)

# Load final Seurat object
seurat_obj <- readRDS("zebrafish_seurat_developmental.rds")

# Check available metadata
cat("Metadata columns:\n")
print(colnames(seurat_obj@meta.data))

# Check cell-type annotations
cat("\nCell-type counts:\n")
print(table(seurat_obj$cell_type))

# Create UMAP colored by biological cell type
p <- DimPlot(
    seurat_obj,
    reduction = "umap",
    group.by = "cell_type",
    label = TRUE,
    repel = TRUE,
    raster = FALSE
) +
    ggtitle("Zebrafish Developmental scRNA-seq: Cell-Type Annotations") +
    theme_classic() +
    theme(
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10)
    )

# Display plot
print(p)

# Save PDF
ggsave(
    "UMAP_cell_type_annotations.pdf",
    plot = p,
    width = 14,
    height = 10
)

cat("\nSaved: UMAP_cell_type_annotations.pdf\n")

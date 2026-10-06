library(Seurat)
library(ggplot2)

# Load final developmental Seurat object
seurat_obj <- readRDS("zebrafish_seurat_developmental.rds")

# Check developmental-stage metadata
cat("Developmental stage counts:\n")
print(table(seurat_obj$developmental_stage))

# Create UMAP colored by developmental stage
p <- DimPlot(
    seurat_obj,
    reduction = "umap",
    group.by = "developmental_stage",
    raster = FALSE
) +
    ggtitle("Zebrafish Developmental scRNA-seq by Developmental Stage") +
    theme_classic() +
    theme(
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10)
    )

# Display plot
print(p)

# Save PDF
ggsave(
    "UMAP_developmental_stage.pdf",
    plot = p,
    width = 10,
    height = 8
)

cat("\nSaved: UMAP_developmental_stage.pdf\n")

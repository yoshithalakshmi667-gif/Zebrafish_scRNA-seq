# 06_variable_genes.R
# Identify highly variable genes

library(Seurat)

# 1. Load normalized Seurat object

seurat_obj <- readRDS("zebrafish_seurat_normalized.rds")

cat("Loaded normalized Seurat object:\n")
print(seurat_obj)

# 2. Find highly variable genes

seurat_obj <- FindVariableFeatures(
    seurat_obj,
    selection.method = "vst",
    nfeatures = 2000
)

cat("\nNumber of variable genes:\n")
print(length(VariableFeatures(seurat_obj)))

# 3. Display top variable genes

cat("\nTop 20 variable genes:\n")

print(
    head(
        VariableFeatures(seurat_obj),
        20
    )
)

# 4. Plot highly variable genes

pdf(
    "variable_genes.pdf",
    width = 10,
    height = 7
)

print(
    VariableFeaturePlot(seurat_obj)
)

dev.off()

cat("\nVariable gene plot saved as variable_genes.pdf\n")

# 5. Save variable-gene object

saveRDS(
    seurat_obj,
    file = "zebrafish_seurat_variable_genes.rds"
)

cat("\nVariable-gene Seurat object saved successfully.\n")

# Quality control for zebrafish scRNA-seq data

library(Seurat)
library(Matrix)

# 1. Load raw Seurat object

seurat_obj <- readRDS("zebrafish_seurat_raw.rds")

cat("Loaded Seurat object:\n")
print(seurat_obj)

# 2. Identify mitochondrial genes

mito_genes <- grep(
    "^ENSDARG.*-mt-",
    rownames(seurat_obj),
    ignore.case = TRUE,
    value = TRUE
)

cat("\nNumber of mitochondrial genes found:", length(mito_genes), "\n")

cat("\nMitochondrial genes:\n")
print(mito_genes)

# 3. Calculate mitochondrial percentage

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
    seurat_obj,
    features = mito_genes
)

# 4. Display QC statistics

cat("\nQC statistics - genes detected per cell:\n")
print(summary(seurat_obj$nFeature_RNA))

cat("\nQC statistics - UMIs per cell:\n")
print(summary(seurat_obj$nCount_RNA))

cat("\nQC statistics - mitochondrial percentage:\n")
print(summary(seurat_obj$percent.mt))

# 5. QC violin plots

pdf(
    "QC_violin_plots_before_filtering.pdf",
    width = 12,
    height = 5
)

print(
    VlnPlot(
        seurat_obj,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3,
        pt.size = 0
    )
)

dev.off()

cat("\nQC violin plot saved as QC_violin_plots_before_filtering.pdf\n")

# 6. QC scatter plots

pdf(
    "QC_scatter_plots_before_filtering.pdf",
    width = 10,
    height = 5
)

print(
    FeatureScatter(
        seurat_obj,
        feature1 = "nCount_RNA",
        feature2 = "nFeature_RNA"
    ))

print(
    FeatureScatter(
        seurat_obj,
        feature1 = "nCount_RNA",
        feature2 = "percent.mt"
    ))

dev.off()

cat("\nQC scatter plots saved as QC_scatter_plots_before_filtering.pdf\n")

# 5. Save QC object before filtering

saveRDS(
    seurat_obj,
    file = "zebrafish_seurat_QC_before_filtering.rds"
)

cat("\nQC object saved successfully.\n")

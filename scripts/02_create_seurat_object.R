# Create Seurat object from expression matrix and metadata
library(Seurat)
library(Matrix)

# 1. Load expression matrix

expr_matrix <- readRDS(
    "zebrafish_expression_sparse.rds"
)

# 2. Read metadata

metadata <- read.delim(
    "meta.tsv",
    header = TRUE,
    stringsAsFactors = FALSE
)

# 3. Verify dimensions

cat("Expression matrix dimensions:\n")
print(dim(expr_matrix))

cat("\nMetadata dimensions:\n")
print(dim(metadata))

# 4. Verify cell IDs

if (!all(colnames(expr_matrix) == metadata$Cell)) {
    stop("ERROR: Expression matrix and metadata cell IDs do not match.")
}

cat("\nCell IDs match successfully.\n")

# 5. Create Seurat object

seurat_obj <- CreateSeuratObject(
    counts = expr_matrix,
    meta.data = metadata,
    project = "Zebrafish_scRNA"
)

# 6. Check object

print(seurat_obj)

# 7. Save Seurat object

saveRDS(
    seurat_obj,
    file = "zebrafish_seurat_raw.rds"
)

cat("\nSeurat object saved successfully.\n")

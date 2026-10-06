# 12_developmental_metadata.R
# Create developmental-stage metadata

library(Seurat)

cat("Loading annotated Seurat object...\n")

seurat_obj <- readRDS("zebrafish_seurat_annotated.rds")

cat("Object loaded successfully.\n\n")

# Create developmental stage

seurat_obj$developmental_stage <- sub(
  ".*_",
  "",
  seurat_obj$sample_name
)

# Convert to ordered factor
seurat_obj$developmental_stage <- factor(
  seurat_obj$developmental_stage,
  levels = c("1dpf", "2dpf", "5dpf")
)

# Check developmental stages

cat("Developmental stages:\n")
print(table(seurat_obj$developmental_stage))

cat("\nSample by developmental stage:\n")
print(table(
  seurat_obj$sample_name,
  seurat_obj$developmental_stage
))

# Check cell types by developmental stage

cat("\nCell types by developmental stage:\n")

cell_type_stage_table <- table(
  seurat_obj$cell_type,
  seurat_obj$developmental_stage
)

print(cell_type_stage_table)

# Save summary table

write.csv(
  as.data.frame(cell_type_stage_table),
  "cell_type_by_developmental_stage_counts.csv",
  row.names = FALSE
)

# Save updated Seurat object

saveRDS(
  seurat_obj,
  "zebrafish_seurat_developmental.rds"
)

cat("\nDevelopmental metadata created successfully.\n")
cat("Saved:\n")
cat("  cell_type_by_developmental_stage_counts.csv\n")
cat("  zebrafish_seurat_developmental.rds\n")

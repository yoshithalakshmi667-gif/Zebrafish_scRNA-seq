# 13_cell_type_composition.R
# Cell-type composition across developmental stages

library(Seurat)
library(ggplot2)

cat("Loading developmental Seurat object...\n")

seurat_obj <- readRDS("zebrafish_seurat_developmental.rds")

cat("Object loaded successfully.\n\n")

# Create count table

count_table <- as.data.frame(
  table(
    seurat_obj$cell_type,
    seurat_obj$developmental_stage
  )
)

colnames(count_table) <- c(
  "CellType",
  "DevelopmentalStage",
  "CellCount"
)

# Calculate percentage within each developmental stage

count_table$Percentage <- ave(
  count_table$CellCount,
  count_table$DevelopmentalStage,
  FUN = function(x) {
    100 * x / sum(x)
  }
)

# Print results

cat("Cell-type composition:\n\n")

print(count_table)

# Save complete table

write.csv(
  count_table,
  "cell_type_composition_by_developmental_stage.csv",
  row.names = FALSE
)

# Create percentage plot

p <- ggplot(
  count_table,
  aes(
    x = DevelopmentalStage,
    y = Percentage,
    fill = CellType
  )
) +
  geom_bar(
    stat = "identity",
    position = "fill"
  ) +
  labs(
    title = "Cell-Type Composition Across Zebrafish Development",
    x = "Developmental Stage",
    y = "Proportion of Cells"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 8)
  )

# Save plot

ggsave(
  "cell_type_composition_developmental_stages.pdf",
  plot = p,
  width = 12,
  height = 8
)

cat("\nAnalysis completed successfully.\n")

cat("Saved:\n")
cat("  cell_type_composition_by_developmental_stage.csv\n")
cat("  cell_type_composition_developmental_stages.pdf\n")

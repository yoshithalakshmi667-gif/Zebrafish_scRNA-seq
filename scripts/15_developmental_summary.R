# Summarize cell-type proportions across development

library(Seurat)
library(dplyr)
library(ggplot2)

cat("Loading developmental Seurat object...\n")

seurat_obj <- readRDS("zebrafish_seurat_developmental.rds")

# Calculate cell-type counts

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

# Calculate percentages

count_table <- count_table %>%
  group_by(DevelopmentalStage) %>%
  mutate(
    Percentage = 100 * CellCount / sum(CellCount)
  ) %>%
  ungroup()

# Save complete table

write.csv(
  count_table,
  "developmental_cell_type_proportions.csv",
  row.names = FALSE
)

# Print major populations

major_cell_types <- c(
  "Neural progenitor cells",
  "Neurogenic progenitor cells",
  "Neurons",
  "Neural / glial cells",
  "Neural crest / mesenchymal cells",
  "Perivascular / cranial mesenchymal cells",
  "Perivascular cells",
  "Fibroblast / ECM cells",
  "Epidermal epithelial cells",
  "Gill / ion-transporting epithelial cells",
  "Chondrocytes / cartilage cells",
  "Hepatic cells",
  "Macrophages / immune cells",
  "Photoreceptors"
)

major_table <- count_table %>%
  filter(CellType %in% major_cell_types)

cat("\nMajor developmental populations:\n\n")
print(major_table)

p <- ggplot(
  major_table,
  aes(
    x = DevelopmentalStage,
    y = Percentage,
    group = CellType,
    linetype = CellType
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(
    title = "Major Cell-Type Proportions Across Development",
    x = "Developmental Stage",
    y = "Percentage of Cells",
    linetype = "Cell Type"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold"
    ),
    legend.text = element_text(size = 8)
  )

ggsave(
  "major_cell_type_developmental_changes.pdf",
  plot = p,
  width = 11,
  height = 8
)

cat("\nAnalysis completed successfully.\n")

cat("Saved:\n")
cat("  developmental_cell_type_proportions.csv\n")
cat("  major_cell_type_developmental_changes.pdf\n")

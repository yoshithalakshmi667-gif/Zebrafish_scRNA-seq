library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)

cat("========================================\n")
cat("Creating developmental DEG heatmap\n")
cat("========================================\n")

# --------------------------------------------------
# Load Seurat object
# --------------------------------------------------

seurat_obj <- readRDS(
    "zebrafish_seurat_developmental.rds"
)

cat(
    "Cells:",
    ncol(seurat_obj),
    "\n"
)

# --------------------------------------------------
# Read DEG results
# --------------------------------------------------

deg_1_2 <- read.csv(
    "DEG_1dpf_vs_2dpf.csv",
    stringsAsFactors = FALSE
)

deg_2_5 <- read.csv(
    "DEG_2dpf_vs_5dpf.csv",
    stringsAsFactors = FALSE
)

deg_1_5 <- read.csv(
    "DEG_1dpf_vs_5dpf.csv",
    stringsAsFactors = FALSE
)

# --------------------------------------------------
# Select top genes from each comparison
# --------------------------------------------------

top_1_2 <- deg_1_2 %>%
    arrange(desc(abs(avg_log2FC))) %>%
    pull(gene) %>%
    head(15)

top_2_5 <- deg_2_5 %>%
    arrange(desc(abs(avg_log2FC))) %>%
    pull(gene) %>%
    head(15)

top_1_5 <- deg_1_5 %>%
    arrange(desc(abs(avg_log2FC))) %>%
    pull(gene) %>%
    head(15)

top_genes <- unique(
    c(
        top_1_2,
        top_2_5,
        top_1_5
    )
)

cat(
    "Number of unique genes selected:",
    length(top_genes),
    "\n"
)

# --------------------------------------------------
# Keep genes present in Seurat object
# --------------------------------------------------

available_genes <- top_genes[
    top_genes %in% rownames(seurat_obj)
]

cat(
    "Genes found in Seurat object:",
    length(available_genes),
    "\n"
)

# --------------------------------------------------
# Set developmental stage
# --------------------------------------------------

Idents(seurat_obj) <- "developmental_stage"

# --------------------------------------------------
# Average expression
# --------------------------------------------------

avg_expression <- AverageExpression(
    seurat_obj,
    assays = "RNA",
    features = available_genes,
    group.by = "developmental_stage",
    slot = "data"
)$RNA

cat("\nAverage expression matrix dimensions:\n")
print(dim(avg_expression))

cat("\nAverage expression column names:\n")
print(colnames(avg_expression))

# --------------------------------------------------
# Automatically identify stage columns
# --------------------------------------------------

stage_columns <- colnames(avg_expression)

cat("\nStage columns detected:\n")
print(stage_columns)

# We expect exactly three developmental groups.
if (length(stage_columns) != 3) {

    stop(
        paste(
            "Expected 3 developmental-stage columns, but found",
            length(stage_columns)
        )
    )
}

# --------------------------------------------------
# Rename columns to clean stage names
# --------------------------------------------------

colnames(avg_expression) <- c(
    "1dpf",
    "2dpf",
    "5dpf"
)

# --------------------------------------------------
# Z-score each gene
# --------------------------------------------------

heatmap_matrix <- t(
    scale(
        t(avg_expression)
    )
)

# Replace NA values
heatmap_matrix[
    is.na(heatmap_matrix)
] <- 0

# --------------------------------------------------
# Save heatmap matrix
# --------------------------------------------------

heatmap_df <- as.data.frame(
    heatmap_matrix
)

heatmap_df$gene <- rownames(
    heatmap_matrix
)

heatmap_df <- heatmap_df %>%
    select(
        gene,
        everything()
    )

write.csv(
    heatmap_df,
    "developmental_DEG_heatmap_matrix.csv",
    row.names = FALSE
)

# --------------------------------------------------
# Convert to long format
# --------------------------------------------------

heatmap_long <- heatmap_df %>%
    pivot_longer(
        cols = c(
            "1dpf",
            "2dpf",
            "5dpf"
        ),
        names_to = "developmental_stage",
        values_to = "z_score"
    )

# --------------------------------------------------
# Order genes
# --------------------------------------------------

heatmap_long$gene <- factor(
    heatmap_long$gene,
    levels = rev(
        unique(
            heatmap_df$gene
        )
    )
)

# --------------------------------------------------
# Create heatmap
# --------------------------------------------------

p <- ggplot(
    heatmap_long,
    aes(
        x = developmental_stage,
        y = gene,
        fill = z_score
    )
) +

    geom_tile() +

    labs(
        title = "Top Developmental Differentially Expressed Genes",
        x = "Developmental Stage",
        y = "Gene",
        fill = "Z-score"
    ) +

    theme_minimal() +

    theme(
        axis.text.x = element_text(
            angle = 45,
            hjust = 1
        ),

        axis.text.y = element_text(
            size = 6
        ),

        plot.title = element_text(
            hjust = 0.5,
            face = "bold"
        )
    )

# --------------------------------------------------
# Save plot
# --------------------------------------------------

ggsave(
    "developmental_DEG_heatmap.pdf",
    plot = p,
    width = 9,
    height = 12
)

cat("\n========================================\n")
cat("DEG HEATMAP COMPLETED SUCCESSFULLY\n")
cat("========================================\n")

cat(
    "Saved: developmental_DEG_heatmap.pdf\n"
)

cat(
    "Saved: developmental_DEG_heatmap_matrix.csv\n"
)

library(Seurat)
library(dplyr)

# Load clustered Seurat object
seurat_obj <- readRDS("zebrafish_seurat_UMAP.rds")

# Create cluster vs original annotation table
comparison <- as.data.frame(
    table(
        SeuratCluster = seurat_obj$seurat_clusters,
        OriginalCluster = seurat_obj$ClusterNames
    )
)

# Remove combinations with zero cells
comparison <- comparison %>%
    filter(Freq > 0)

# Find the most common original annotation for each of our clusters
dominant_annotation <- comparison %>%
    group_by(SeuratCluster) %>%
    slice_max(order_by = Freq, n = 1, with_ties = FALSE) %>%
    ungroup()

# Calculate percentage of each cluster represented by its dominant annotation
cluster_sizes <- comparison %>%
    group_by(SeuratCluster) %>%
    summarise(
        TotalCells = sum(Freq),
        .groups = "drop"
    )

dominant_annotation <- dominant_annotation %>%
    left_join(cluster_sizes, by = "SeuratCluster") %>%
    mutate(
        Percent = round(100 * Freq / TotalCells, 2)
    )

# Rename columns for clarity
dominant_annotation <- dominant_annotation %>%
    select(
        OurCluster = SeuratCluster,
        OriginalAnnotation = OriginalCluster,
        Cells = Freq,
        TotalClusterCells = TotalCells,
        Percent
    )

# Save concise summary
write.csv(
    dominant_annotation,
    "dominant_original_annotation_by_cluster.csv",
    row.names = FALSE
)

# Print only the concise summary
cat("\nDominant original annotation for each of our clusters:\n\n")
print(dominant_annotation, row.names = FALSE)

cat("\nSummary saved as:\n")
cat("dominant_original_annotation_by_cluster.csv\n")

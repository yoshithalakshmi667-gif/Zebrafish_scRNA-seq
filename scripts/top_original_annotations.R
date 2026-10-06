library(Seurat)
library(dplyr)

# Load Seurat object
seurat_obj <- readRDS("zebrafish_seurat_UMAP.rds")

# Create cluster vs original annotation table
comparison <- as.data.frame(
    table(
        OurCluster = seurat_obj$seurat_clusters,
        OriginalAnnotation = seurat_obj$ClusterNames
    )
)

# Remove zero-count combinations
comparison <- comparison %>%
    filter(Freq > 0)

# Calculate percentage within each of our clusters
comparison <- comparison %>%
    group_by(OurCluster) %>%
    mutate(
        TotalCells = sum(Freq),
        Percent = round(100 * Freq / TotalCells, 2)
    ) %>%
    ungroup()

# Keep top 5 original annotations for every cluster
top5 <- comparison %>%
    group_by(OurCluster) %>%
    slice_max(
        order_by = Freq,
        n = 5,
        with_ties = FALSE
    ) %>%
    arrange(
        as.numeric(as.character(OurCluster)),
        desc(Freq)
    ) %>%
    ungroup()

# Save results
write.csv(
    top5,
    "top5_original_annotations_by_cluster.csv",
    row.names = FALSE
)

cat("\nTop 5 original annotations for each cluster:\n\n")
print(top5, row.names = FALSE)

cat("\nSaved as:\n")
cat("top5_original_annotations_by_cluster.csv\n")

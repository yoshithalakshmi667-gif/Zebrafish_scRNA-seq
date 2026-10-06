library(Seurat)

# Load UMAP/clustering object
seurat_obj <- readRDS("zebrafish_seurat_UMAP.rds")

# Find marker genes for every cluster
markers <- FindAllMarkers(
    seurat_obj,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25
)

# Save complete marker table
write.csv(
    markers,
    file = "cluster_markers.csv",
    row.names = FALSE
)

# Get top 10 markers from each cluster
top10_markers <- markers |>
    dplyr::group_by(cluster) |>
    dplyr::slice_max(
        order_by = avg_log2FC,
        n = 10
    )

# Save top markers
write.csv(
    top10_markers,
    file = "top10_cluster_markers.csv",
    row.names = FALSE
)

# Heatmap of top 10 markers from each cluster
pdf(
    "top10_cluster_marker_heatmap.pdf",
    width = 14,
    height = 20
)

print(
    DoHeatmap(
        seurat_obj,
        features = unique(top10_markers$gene),
        group.by = "seurat_clusters"
    )
)

dev.off()

cat("\nMarker heatmap saved as top10_cluster_marker_heatmap.pdf\n")

cat("\nMarker analysis completed successfully.\n")

cat("\nNumber of marker genes identified:\n")
print(nrow(markers))

cat("\nTop markers by cluster:\n")
print(top10_markers)

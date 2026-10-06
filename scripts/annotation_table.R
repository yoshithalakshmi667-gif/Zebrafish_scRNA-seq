library(Seurat)
library(dplyr)

# Load Seurat object
seurat_obj <- readRDS("zebrafish_seurat_UMAP.rds")

# Load marker results from Step 10
markers <- read.csv(
    "cluster_markers.csv",
    stringsAsFactors = FALSE
)

# Make sure cluster is treated as a number for sorting
markers$cluster <- as.character(markers$cluster)

# Select top 10 positive markers for each cluster
top10 <- markers %>%
    group_by(cluster) %>%
    slice_max(
        order_by = avg_log2FC,
        n = 10,
        with_ties = FALSE
    ) %>%
    arrange(
        as.numeric(cluster),
        desc(avg_log2FC)
    ) %>%
    ungroup()

# Create a compact table with marker genes
annotation_table <- top10 %>%
    group_by(cluster) %>%
    summarise(
        Top_Markers = paste(gene, collapse = ", "),
        .groups = "drop"
    )

# Rename cluster column
colnames(annotation_table)[1] <- "OurCluster"

# Add dominant original atlas annotation
original_labels <- seurat_obj@meta.data %>%
    mutate(
        OurCluster = as.character(seurat_clusters)
    ) %>%
    count(OurCluster, ClusterNames, name = "Cells") %>%
    group_by(OurCluster) %>%
    slice_max(
        order_by = Cells,
        n = 1,
        with_ties = FALSE
    ) %>%
    ungroup() %>%
    select(
        OurCluster,
        OriginalAnnotation = ClusterNames,
        OriginalCells = Cells
    )

# Combine marker information and original annotation
annotation_table <- annotation_table %>%
    left_join(
        original_labels,
        by = "OurCluster"
    ) %>%
    arrange(as.numeric(OurCluster))

# Save results
write.csv(
    annotation_table,
    "cluster_annotation_table.csv",
    row.names = FALSE
)

# Print results
cat("\nCluster annotation table:\n\n")
print(annotation_table, row.names = FALSE)

cat("\nSaved as:\n")
cat("cluster_annotation_table.csv\n")

library(dplyr)

# Load marker results
markers <- read.csv(
    "cluster_markers.csv",
    stringsAsFactors = FALSE
)

# Select top 10 markers for each cluster
top10 <- markers %>%
    group_by(cluster) %>%
    slice_max(
        order_by = avg_log2FC,
        n = 10,
        with_ties = FALSE
    ) %>%
    arrange(
        as.numeric(as.character(cluster)),
        desc(avg_log2FC)
    ) %>%
    ungroup()

# Keep useful columns
top10 <- top10 %>%
    select(
        cluster,
        gene,
        avg_log2FC,
        pct.1,
        pct.2
    )

# Rename columns
colnames(top10) <- c(
    "OurCluster",
    "MarkerGene",
    "AvgLog2FC",
    "PctCluster",
    "PctOtherClusters"
)

# Save complete table
write.csv(
    top10,
    "top10_markers_for_annotation.csv",
    row.names = FALSE
)

cat("\nTop 10 markers for all clusters:\n\n")
print(top10, n = Inf, width = Inf)

cat("\nSaved as:\n")
cat("top10_markers_for_annotation.csv\n")

library(Seurat)

# Load final Seurat object

seurat_obj <- readRDS("zebrafish_seurat_developmental.rds")

# Check developmental stages

cat("Developmental stage counts:\n")
print(table(seurat_obj$developmental_stage))

# Set developmental stage as the active identity

Idents(seurat_obj) <- "developmental_stage"

cat("\nActive identities:\n")
print(table(Idents(seurat_obj)))

# Function for DEG analysis

run_stage_DE <- function(object, stage1, stage2, output_file) {

    cat("Comparing:", stage1, "vs", stage2, "\n")

    # Check that identities exist
    if (!stage1 %in% levels(Idents(object))) {
        stop(paste("Stage not found:", stage1))
    }

    if (!stage2 %in% levels(Idents(object))) {
        stop(paste("Stage not found:", stage2))
    }

    cat(stage1, "cells:",
        sum(Idents(object) == stage1), "\n")

    cat(stage2, "cells:",
        sum(Idents(object) == stage2), "\n")

    # Differential expression
    markers <- FindMarkers(
        object,
        ident.1 = stage1,
        ident.2 = stage2,
        logfc.threshold = 0.25,
        min.pct = 0.10,
        test.use = "wilcox"
    )

    # Add gene names
    markers$gene <- rownames(markers)

    # Put gene column first
    markers <- markers[, c(
        "gene",
        setdiff(colnames(markers), "gene")
    )]

    # Sort by adjusted p-value
    markers <- markers[
        order(markers$p_val_adj),
    ]

    # Save results
    write.csv(
        markers,
        output_file,
        row.names = FALSE
    )

    cat("Number of DE genes:", nrow(markers), "\n")

    cat("Saved:", output_file, "\n")

    cat("\nTop 20 DE genes:\n")
    print(head(markers, 20))

    return(markers)
}

# 1 dpf vs 2 dpf

DE_1_vs_2 <- run_stage_DE(
    seurat_obj,
    "1dpf",
    "2dpf",
    "DEG_1dpf_vs_2dpf.csv"
)

# 2 dpf vs 5 dpf

DE_2_vs_5 <- run_stage_DE(
    seurat_obj,
    "2dpf",
    "5dpf",
    "DEG_2dpf_vs_5dpf.csv"
)

# 1 dpf vs 5 dpf

DE_1_vs_5 <- run_stage_DE(
    seurat_obj,
    "1dpf",
    "5dpf",
    "DEG_1dpf_vs_5dpf.csv"
)

cat("ALL DEVELOPMENTAL DE ANALYSES COMPLETED\n")

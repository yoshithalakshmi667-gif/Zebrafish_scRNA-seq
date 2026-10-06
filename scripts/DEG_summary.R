library(dplyr)

# Function to summarize DEG results
summarize_DEG <- function(
    input_file,
    output_prefix,
    comparison
) {

    cat("Processing:", comparison, "\n")
    deg <- read.csv(
        input_file,
        stringsAsFactors = FALSE
    )

    # Remove genes without adjusted p-values
    deg <- deg %>%
        filter(!is.na(p_val_adj))

    # Significant genes
    significant <- deg %>%
        filter(
            p_val_adj < 0.05,
            abs(avg_log2FC) >= 0.25
        )

    # Upregulated in ident.1
    upregulated <- significant %>%
        filter(avg_log2FC > 0) %>%
        arrange(desc(avg_log2FC))

    # Downregulated in ident.1
    downregulated <- significant %>%
        filter(avg_log2FC < 0) %>%
        arrange(avg_log2FC)

    # Top 20 upregulated
    top_up <- head(upregulated, 20)

    # Top 20 downregulated
    top_down <- head(downregulated, 20)

    # Save
    write.csv(
        top_up,
        paste0(output_prefix, "_top20_upregulated.csv"),
        row.names = FALSE
    )

    write.csv(
        top_down,
        paste0(output_prefix, "_top20_downregulated.csv"),
        row.names = FALSE
    )

    write.csv(
        significant,
        paste0(output_prefix, "_significant.csv"),
        row.names = FALSE
    )

    # Print summary
    cat("Total genes tested:", nrow(deg), "\n")
    cat("Significant DE genes:", nrow(significant), "\n")
    cat("Upregulated in ident.1:", nrow(upregulated), "\n")
    cat("Downregulated in ident.1:", nrow(downregulated), "\n")

    cat("\nTop 20 UPREGULATED genes:\n")
    print(
        top_up[, c(
            "gene",
            "avg_log2FC",
            "pct.1",
            "pct.2",
            "p_val_adj"
        )]
    )

    cat("\nTop 20 DOWNREGULATED genes:\n")
    print(
        top_down[, c(
            "gene",
            "avg_log2FC",
            "pct.1",
            "pct.2",
            "p_val_adj"
        )]
    )
}

# 1 dpf vs 2 dpf

summarize_DEG(
    "DEG_1dpf_vs_2dpf.csv",
    "DEG_1dpf_vs_2dpf",
    "1 dpf vs 2 dpf"
)

# 2 dpf vs 5 dpf

summarize_DEG(
    "DEG_2dpf_vs_5dpf.csv",
    "DEG_2dpf_vs_5dpf",
    "2 dpf vs 5 dpf"
)

# 1 dpf vs 5 dpf

summarize_DEG(
    "DEG_1dpf_vs_5dpf.csv",
    "DEG_1dpf_vs_5dpf",
    "1 dpf vs 5 dpf"
)

cat("DEG SUMMARY COMPLETED\n")

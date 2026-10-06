library(ggplot2)
library(dplyr)

create_volcano <- function(
    input_file,
    output_file,
    comparison
) {

    cat("Creating volcano plot:", comparison, "\n")

    deg <- read.csv(
        input_file,
        stringsAsFactors = FALSE
    )

    # Remove missing values
    deg <- deg %>%
        filter(
            !is.na(avg_log2FC),
            !is.na(p_val_adj)
        )

    # Avoid -log10(0) = Inf
    # Replace zero adjusted p-values with the smallest
    # positive adjusted p-value in this dataset
    min_positive_p <- min(
        deg$p_val_adj[deg$p_val_adj > 0],
        na.rm = TRUE
    )

    deg <- deg %>%
        mutate(
            plot_p_adj = ifelse(
                p_val_adj == 0,
                min_positive_p,
                p_val_adj
            )
        )

    # Classify genes
    deg <- deg %>%
        mutate(
            significance = case_when(
                p_val_adj < 0.05 &
                    avg_log2FC >= 0.25 ~ "Upregulated",

                p_val_adj < 0.05 &
                    avg_log2FC <= -0.25 ~ "Downregulated",

                TRUE ~ "Not significant"
            )
        )

    # Count categories
    cat(
        "Upregulated:",
        sum(deg$significance == "Upregulated"),
        "\n"
    )

    cat(
        "Downregulated:",
        sum(deg$significance == "Downregulated"),
        "\n"
    )

    cat(
        "Not significant:",
        sum(deg$significance == "Not significant"),
        "\n"
    )

    # Volcano plot
    p <- ggplot(
        deg,
        aes(
            x = avg_log2FC,
            y = -log10(plot_p_adj)
        )
    ) +

        geom_point(
            aes(
                shape = significance
            ),
            size = 1.5,
            alpha = 0.6
        ) +

        geom_vline(
            xintercept = c(-0.25, 0.25),
            linetype = "dashed"
        ) +

        geom_hline(
            yintercept = -log10(0.05),
            linetype = "dashed"
        ) +

        labs(
            title = comparison,
            x = "Average log2 Fold Change",
            y = "-log10 Adjusted P-value",
            shape = "Significance"
        ) +

        theme_classic() +

        theme(
            plot.title = element_text(
                hjust = 0.5,
                face = "bold"
            ),
            legend.position = "right"
        )

    ggsave(
        output_file,
        plot = p,
        width = 9,
        height = 7
    )

    cat(
        "Saved:",
        output_file,
        "\n"
    )
}



create_volcano(
    "DEG_1dpf_vs_2dpf.csv",
    "volcano_1dpf_vs_2dpf.pdf",
    "1 dpf vs 2 dpf"
)


create_volcano(
    "DEG_2dpf_vs_5dpf.csv",
    "volcano_2dpf_vs_5dpf.pdf",
    "2 dpf vs 5 dpf"
)


# 1 dpf vs 5 dpf

create_volcano(
    "DEG_1dpf_vs_5dpf.csv",
    "volcano_1dpf_vs_5dpf.pdf",
    "1 dpf vs 5 dpf"
)


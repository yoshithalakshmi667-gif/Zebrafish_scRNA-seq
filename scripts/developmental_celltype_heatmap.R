library(ggplot2)

# Read developmental cell-type proportions
data <- read.csv(
    "developmental_cell_type_proportions.csv",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

# Check the data
cat("Columns in input file:\n")
print(colnames(data))

cat("\nFirst rows:\n")
print(head(data))

# Identify the percentage column

percentage_col <- "Percentage"

if (!percentage_col %in% colnames(data)) {
    stop("Could not find the Percentage column in the input file.")
}

# Make sure developmental stage is a factor

data$DevelopmentalStage <- factor(
    data$DevelopmentalStage,
    levels = c("1 dpf", "2 dpf", "5 dpf")
)

# Convert percentage to numeric

data[[percentage_col]] <- as.numeric(data[[percentage_col]])

# Create heatmap

p <- ggplot(
    data,
    aes(
        x = DevelopmentalStage,
        y = CellType,
        fill = Percentage
    )
) +
    geom_tile(color = "white") +
    geom_text(
        aes(label = sprintf("%.1f%%", Percentage)),
        size = 3
    ) +
    scale_y_discrete(limits = rev(unique(data$CellType))) +
    labs(
        title = "Cell-Type Composition Across Developmental Stages",
        x = "Developmental Stage",
        y = "Cell Type",
        fill = "Percentage"
    ) +
    theme_classic() +
    theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 8),
        legend.title = element_text(size = 10)
    )

# Display
print(p)

# Save
ggsave(
    "developmental_celltype_heatmap.pdf",
    plot = p,
    width = 10,
    height = 14
)

cat("\nSaved: developmental_celltype_heatmap.pdf\n")

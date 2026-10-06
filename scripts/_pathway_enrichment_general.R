library(Seurat)
library(dplyr)
library(clusterProfiler)
library(org.Dr.eg.db)
library(AnnotationDbi)
library(ggplot2)

# ============================================================
# USER SETTINGS
# ============================================================

deg_file <- "DEG_2dpf_vs_5dpf.csv"
output_dir <- "pathway_enrichment_2dpf_vs_5dpf"

# Positive log2FC = higher in 2 dpf
# Negative log2FC = higher in 5 dpf

padj_cutoff <- 0.05
logfc_cutoff <- 0.25

# ============================================================
# CREATE OUTPUT DIRECTORY
# ============================================================

dir.create(output_dir, showWarnings = FALSE)

# ============================================================
# READ DEG RESULTS
# ============================================================

deg <- read.csv(deg_file, stringsAsFactors = FALSE)

cat("Total genes:", nrow(deg), "\n")

# ============================================================
# SIGNIFICANT DEGs
# ============================================================

sig <- deg %>%
    filter(
        p_val_adj < padj_cutoff,
        abs(avg_log2FC) >= logfc_cutoff
    )

up <- sig %>%
    filter(avg_log2FC > 0)

down <- sig %>%
    filter(avg_log2FC < 0)

cat("Significant DEGs:", nrow(sig), "\n")
cat("Up in 2dpf:", nrow(up), "\n")
cat("Up in 5dpf:", nrow(down), "\n")

# ============================================================
# EXTRACT ENSEMBL IDs
# ============================================================

deg$ENSEMBL <- sub("-.*$", "", deg$gene)
sig$ENSEMBL <- sub("-.*$", "", sig$gene)
up$ENSEMBL <- sub("-.*$", "", up$gene)
down$ENSEMBL <- sub("-.*$", "", down$gene)

# ============================================================
# FUNCTION TO CONVERT ENSEMBL → ENTREZ
# ============================================================

convert_ids <- function(x) {

    x <- unique(x)

    mapping <- AnnotationDbi::select(
        org.Dr.eg.db,
        keys = x,
        columns = c("ENTREZID", "SYMBOL"),
        keytype = "ENSEMBL"
    )

    mapping <- mapping %>%
        filter(!is.na(ENTREZID)) %>%
        distinct(ENSEMBL, ENTREZID, .keep_all = TRUE)

    return(mapping)
}

# ============================================================
# BACKGROUND / UNIVERSE
# ============================================================

background_mapping <- convert_ids(deg$ENSEMBL)

universe <- unique(background_mapping$ENTREZID)

cat("Background mapped genes:", length(universe), "\n")

# ============================================================
# MAP SIGNIFICANT GENES
# ============================================================

sig_mapping <- convert_ids(sig$ENSEMBL)
up_mapping <- convert_ids(up$ENSEMBL)
down_mapping <- convert_ids(down$ENSEMBL)

write.csv(
    sig_mapping,
    file.path(output_dir, "significant_gene_mapping.csv"),
    row.names = FALSE
)

write.csv(
    up_mapping,
    file.path(output_dir, "up_2dpf_gene_mapping.csv"),
    row.names = FALSE
)

write.csv(
    down_mapping,
    file.path(output_dir, "down_2dpf_gene_mapping.csv"),
    row.names = FALSE
)

# ============================================================
# ENRICHMENT FUNCTION
# ============================================================

run_enrichment <- function(gene_ids, direction) {

    gene_ids <- unique(gene_ids)

    cat("\nRunning enrichment for:", direction, "\n")
    cat("Genes:", length(gene_ids), "\n")

    # -------------------------
    # GO Biological Process
    # -------------------------

    ego_bp <- enrichGO(
        gene = gene_ids,
        universe = universe,
        OrgDb = org.Dr.eg.db,
        keyType = "ENTREZID",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        readable = TRUE
    )

    # -------------------------
    # GO Cellular Component
    # -------------------------

    ego_cc <- enrichGO(
        gene = gene_ids,
        universe = universe,
        OrgDb = org.Dr.eg.db,
        keyType = "ENTREZID",
        ont = "CC",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        readable = TRUE
    )

    # -------------------------
    # GO Molecular Function
    # -------------------------

    ego_mf <- enrichGO(
        gene = gene_ids,
        universe = universe,
        OrgDb = org.Dr.eg.db,
        keyType = "ENTREZID",
        ont = "MF",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        readable = TRUE
    )

    # -------------------------
    # KEGG
    # -------------------------

    ekegg <- enrichKEGG(
        gene = gene_ids,
        universe = universe,
        organism = "dre",
        keyType = "ncbi-geneid",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05
    )

    # -------------------------
    # Save results
    # -------------------------

    write.csv(
        as.data.frame(ego_bp),
        file.path(output_dir, paste0(direction, "_GO_BP.csv")),
        row.names = FALSE
    )

    write.csv(
        as.data.frame(ego_cc),
        file.path(output_dir, paste0(direction, "_GO_CC.csv")),
        row.names = FALSE
    )

    write.csv(
        as.data.frame(ego_mf),
        file.path(output_dir, paste0(direction, "_GO_MF.csv")),
        row.names = FALSE
    )

    write.csv(
        as.data.frame(ekegg),
        file.path(output_dir, paste0(direction, "_KEGG.csv")),
        row.names = FALSE
    )

    # -------------------------
    # GO BP dotplot
    # -------------------------

    if (nrow(as.data.frame(ego_bp)) > 0) {

        pdf(
            file.path(output_dir, paste0(direction, "_GO_BP_dotplot.pdf")),
            width = 10,
            height = 8
        )

        print(dotplot(ego_bp, showCategory = 15))

        dev.off()
    }

    # -------------------------
    # KEGG dotplot
    # -------------------------

    if (nrow(as.data.frame(ekegg)) > 0) {

        pdf(
            file.path(output_dir, paste0(direction, "_KEGG_dotplot.pdf")),
            width = 10,
            height = 8
        )

        print(dotplot(ekegg, showCategory = 15))

        dev.off()
    }
}

# ============================================================
# RUN ENRICHMENT
# ============================================================

run_enrichment(
    up_mapping$ENTREZID,
    "up_2dpf"
)

run_enrichment(
    down_mapping$ENTREZID,
    "down_2dpf"
)

cat("\n========================================\n")
cat("Pathway enrichment completed.\n")
cat("Output directory:", output_dir, "\n")
cat("========================================\n")

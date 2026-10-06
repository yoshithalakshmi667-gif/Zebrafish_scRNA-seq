# ============================================================
# 22_pathway_enrichment.R
# Zebrafish developmental scRNA-seq
#
# GO and KEGG enrichment analysis
# Comparison: 1 dpf vs 2 dpf
# ============================================================

library(clusterProfiler)
library(org.Dr.eg.db)
library(AnnotationDbi)
library(ggplot2)


input_file <- "DEG_1dpf_vs_2dpf.csv"
output_dir <- "pathway_enrichment_1dpf_vs_2dpf"

dir.create(output_dir, showWarnings = FALSE)

# 2. Read DEG results
deg <- read.csv(
    input_file,
    stringsAsFactors = FALSE
)

cat("Total genes in DEG file:", nrow(deg), "\n")

# 3. Define significant DEGs

sig <- deg[
    !is.na(deg$p_val_adj) &
    deg$p_val_adj < 0.05 &
    abs(deg$avg_log2FC) >= 0.25,
]

cat("Significant DEGs:", nrow(sig), "\n")

# 4. Separate UP and DOWN genes
#
# Positive FC = higher in 1 dpf
# Negative FC = higher in 2 dpf

up_1dpf <- sig[sig$avg_log2FC > 0, ]
down_1dpf <- sig[sig$avg_log2FC < 0, ]

cat("Up in 1 dpf:", nrow(up_1dpf), "\n")
cat("Down in 1 dpf / Up in 2 dpf:", nrow(down_1dpf), "\n")


extract_ensembl <- function(x) {
    sub("-.*$", "", x)
}

sig$ENSEMBL <- extract_ensembl(sig$gene)
up_1dpf$ENSEMBL <- extract_ensembl(up_1dpf$gene)
down_1dpf$ENSEMBL <- extract_ensembl(down_1dpf$gene)

# 6. Convert Ensembl IDs to Entrez IDs
convert_ids <- function(ensembl_ids) {

    ensembl_ids <- unique(ensembl_ids)

    mapping <- AnnotationDbi::select(
        org.Dr.eg.db,
        keys = ensembl_ids,
        keytype = "ENSEMBL",
        columns = c("ENSEMBL", "ENTREZID", "SYMBOL", "GENENAME")
    )

    mapping <- mapping[
        !is.na(mapping$ENTREZID),
    ]

    mapping <- mapping[
        !duplicated(mapping$ENSEMBL),
    ]

    return(mapping)
}

sig_map <- convert_ids(sig$ENSEMBL)
up_map <- convert_ids(up_1dpf$ENSEMBL)
down_map <- convert_ids(down_1dpf$ENSEMBL)

cat("Mapped significant genes:", nrow(sig_map), "\n")
cat("Mapped UP genes:", nrow(up_map), "\n")
cat("Mapped DOWN genes:", nrow(down_map), "\n")

# 7. Save gene mappings

write.csv(
    sig_map,
    file.path(output_dir, "significant_gene_mapping.csv"),
    row.names = FALSE
)

write.csv(
    up_map,
    file.path(output_dir, "up_1dpf_gene_mapping.csv"),
    row.names = FALSE
)

write.csv(
    down_map,
    file.path(output_dir, "down_1dpf_gene_mapping.csv"),
    row.names = FALSE
)

# 8. Create background/universe
# All genes tested in the DEG analysis
background_ensembl <- unique(
    extract_ensembl(deg$gene)
)

background_map <- convert_ids(background_ensembl)

background_entrez <- unique(
    background_map$ENTREZID
)

cat("Background genes mapped:", length(background_entrez), "\n")

# 9. Function for GO enrichment
run_GO <- function(entrez_ids, label) {

    if (length(entrez_ids) < 5) {
        cat("Too few genes for GO:", label, "\n")
        return(NULL)
    }

    # Biological Process
    ego_BP <- enrichGO(
        gene = entrez_ids,
        universe = background_entrez,
        OrgDb = org.Dr.eg.db,
        keyType = "ENTREZID",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.20,
        readable = TRUE
    )

    # Molecular Function
    ego_MF <- enrichGO(
        gene = entrez_ids,
        universe = background_entrez,
        OrgDb = org.Dr.eg.db,
        keyType = "ENTREZID",
        ont = "MF",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.20,
        readable = TRUE
    )

    # Cellular Component
    ego_CC <- enrichGO(
        gene = entrez_ids,
        universe = background_entrez,
        OrgDb = org.Dr.eg.db,
        keyType = "ENTREZID",
        ont = "CC",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.20,
        readable = TRUE
    )

    # Save tables
    if (!is.null(ego_BP) && nrow(as.data.frame(ego_BP)) > 0) {
        write.csv(
            as.data.frame(ego_BP),
            file.path(output_dir, paste0(label, "_GO_BP.csv")),
            row.names = FALSE
        )

        pdf(
            file.path(output_dir, paste0(label, "_GO_BP_dotplot.pdf")),
            width = 10,
            height = 7
        )
        print(dotplot(ego_BP, showCategory = 20))
        dev.off()
    }

    if (!is.null(ego_MF) && nrow(as.data.frame(ego_MF)) > 0) {
        write.csv(
            as.data.frame(ego_MF),
            file.path(output_dir, paste0(label, "_GO_MF.csv")),
            row.names = FALSE
        )
    }

    if (!is.null(ego_CC) && nrow(as.data.frame(ego_CC)) > 0) {
        write.csv(
            as.data.frame(ego_CC),
            file.path(output_dir, paste0(label, "_GO_CC.csv")),
            row.names = FALSE
        )
    }

    return(list(
        BP = ego_BP,
        MF = ego_MF,
        CC = ego_CC
    ))
}

# 10. GO enrichment
cat("\nRunning GO enrichment for genes UP in 1 dpf...\n")

GO_up <- run_GO(
    unique(up_map$ENTREZID),
    "up_1dpf"
)

cat("\nRunning GO enrichment for genes DOWN in 1 dpf...\n")

GO_down <- run_GO(
    unique(down_map$ENTREZID),
    "down_1dpf"
)

# 11. KEGG enrichment
run_KEGG <- function(entrez_ids, label) {

    if (length(entrez_ids) < 5) {
        cat("Too few genes for KEGG:", label, "\n")
        return(NULL)
    }

    ekegg <- enrichKEGG(
        gene = unique(entrez_ids),
        universe = background_entrez,
        organism = "dre",
        keyType = "ncbi-geneid",
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        qvalueCutoff = 0.20
    )

    if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {

        write.csv(
            as.data.frame(ekegg),
            file.path(output_dir, paste0(label, "_KEGG.csv")),
            row.names = FALSE
        )

        pdf(
            file.path(output_dir, paste0(label, "_KEGG_dotplot.pdf")),
            width = 10,
            height = 7
        )

        print(dotplot(ekegg, showCategory = 20))

        dev.off()
    }

    return(ekegg)
}

cat("\nRunning KEGG enrichment for genes UP in 1 dpf...\n")

KEGG_up <- run_KEGG(
    unique(up_map$ENTREZID),
    "up_1dpf"
)

cat("\nRunning KEGG enrichment for genes DOWN in 1 dpf...\n")

KEGG_down <- run_KEGG(
    unique(down_map$ENTREZID),
    "down_1dpf"
)

# 12. Final summary


cat("Input:", input_file, "\n")
cat("Total tested genes:", nrow(deg), "\n")
cat("Significant genes:", nrow(sig), "\n")
cat("Up in 1 dpf:", nrow(up_1dpf), "\n")
cat("Down in 1 dpf:", nrow(down_1dpf), "\n")
cat("Background mapped genes:", length(background_entrez), "\n")
cat("Output directory:", output_dir, "\n")

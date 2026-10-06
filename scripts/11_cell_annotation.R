# Final cell-type annotation for zebrafish scRNA-seq

library(Seurat)

cat("Loading clustered Seurat object...\n")

# Load Seurat object after clustering
seurat_obj <- readRDS("zebrafish_seurat_UMAP.rds")

cat("Object loaded successfully.\n")
cat("Number of cells:", ncol(seurat_obj), "\n")
cat("Number of clusters:", length(unique(Idents(seurat_obj))), "\n")


# Final cluster annotations

cluster_annotations <- c(
  "0"  = "Neural progenitor cells",
  "1"  = "Neurogenic progenitor cells",
  "2"  = "Neural crest / mesenchymal cells",
  "3"  = "Neurons",
  "4"  = "Perivascular / cranial mesenchymal cells",
  "5"  = "Fibroblast / ECM cells",
  "6"  = "Neural / glial cells",
  "7"  = "Perivascular cells",
  "8"  = "Gill / ion-transporting epithelial cells",
  "9"  = "Epidermal epithelial cells",
  "10" = "Basal epithelial cells",
  "11" = "Erythrocytes",
  "12" = "Ion-transporting / pharyngeal epithelial cells",
  "13" = "Photoreceptors",
  "14" = "Melanocytes / pigment cells",
  "15" = "Macrophages / immune cells",
  "16" = "Intestinal epithelial cells",
  "17" = "Mucosal / olfactory epithelial cells",
  "18" = "Skeletal muscle cells",
  "19" = "Endothelial cells",
  "20" = "Erythroid cells",
  "21" = "Fast skeletal muscle cells",
  "22" = "Myogenic progenitor cells",
  "23" = "Hepatic cells",
  "24" = "Cardiac / slow skeletal muscle cells",
  "25" = "Chondrocytes / cartilage cells",
  "26" = "Ciliated / otic epithelial cells",
  "27" = "Erythroid progenitor cells",
  "28" = "Lens cells",
  "29" = "Notochord cells",
  "30" = "Cephalic skeletal muscle cells",
  "31" = "Pancreatic acinar cells"
)


# Add annotation to Seurat metadata

seurat_obj$cell_type <- unname(
  cluster_annotations[as.character(Idents(seurat_obj))]
)


# Check annotation results

cat("\nCell-type counts:\n")
print(table(seurat_obj$cell_type))


cat("\nCluster-to-cell-type mapping:\n")

annotation_table <- data.frame(
  Cluster = names(cluster_annotations),
  CellType = unname(cluster_annotations),
  stringsAsFactors = FALSE
)

print(annotation_table)


# Save annotation table

write.csv(
  annotation_table,
  "final_cluster_annotations.csv",
  row.names = FALSE
)


# Save annotated Seurat object

saveRDS(
  seurat_obj,
  "zebrafish_seurat_annotated.rds"
)


cat("\nAnnotation completed successfully.\n")
cat("Saved:\n")
cat("  final_cluster_annotations.csv\n")
cat("  zebrafish_seurat_annotated.rds\n")

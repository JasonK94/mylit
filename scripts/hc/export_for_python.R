#!/usr/bin/env Rscript
# Export HC Seurat embeddings + metadata for Python analyses (MELD, scCODA)

.libPaths(c(
  "/home/user3/GJC_KDW_250721/renv/library/R-4.3/x86_64-pc-linux-gnu",
  "/home/user3/R/x86_64-pc-linux-gnu-library/4.3",
  .libPaths()
))

library(Seurat)
library(SeuratObject)
library(qs)
library(data.table)

cat("=== Export HC data for Python ===\n")

sobj <- qread("/data/user3/sobj/hc_only_v1/2_hc_annotated.qs")
cat(sprintf("  Loaded: %d cells, %d patients\n", ncol(sobj), length(unique(sobj$name))))

out_meld <- "/data/user3/sobj/hc_only_v1/meld"
out_sccoda <- "/data/user3/sobj/hc_only_v1/sccoda"
dir.create(out_meld, recursive = TRUE, showWarnings = FALSE)
dir.create(out_sccoda, recursive = TRUE, showWarnings = FALSE)

# Ensure age_group
if (!"age_group" %in% colnames(sobj@meta.data)) {
  sobj$age_group <- ifelse(sobj$age < 35, "Young",
                           ifelse(sobj$age <= 50, "Middle", "Old"))
}

# 1. scVI embeddings
cat("Exporting scVI embeddings...\n")
scvi_emb <- Embeddings(sobj, reduction = "integrated.scvi")
fwrite(as.data.frame(scvi_emb), file.path(out_meld, "scvi_embeddings.csv"), row.names = TRUE)
cat(sprintf("  scVI: %d × %d\n", nrow(scvi_emb), ncol(scvi_emb)))

# 2. UMAP (try multiple names)
cat("Exporting UMAP...\n")
cat(sprintf("  Available reductions: %s\n", paste(names(sobj@reductions), collapse = ", ")))
umap_name <- NULL
for (nm in c("umap.scvi", "umap", "UMAP")) {
  if (nm %in% names(sobj@reductions)) { umap_name <- nm; break }
}
if (!is.null(umap_name)) {
  umap_emb <- Embeddings(sobj, reduction = umap_name)
  fwrite(as.data.frame(umap_emb), file.path(out_meld, "umap_embeddings.csv"), row.names = TRUE)
  cat(sprintf("  UMAP (%s): %d × %d\n", umap_name, nrow(umap_emb), ncol(umap_emb)))
} else {
  cat("  WARNING: No UMAP reduction found, skipping\n")
}

# 3. Metadata
cat("Exporting metadata...\n")
key_cols <- c("name", "GEM", "age", "sex", "age_group", "anno1", "anno2")
meta <- sobj@meta.data[, key_cols, drop = FALSE]
fwrite(meta, file.path(out_meld, "metadata.csv"), row.names = TRUE)
cat(sprintf("  Metadata: %d cells, %d cols\n", nrow(meta), ncol(meta)))

# 4. scCODA composition table (patient × cell type)
cat("Exporting composition table for scCODA...\n")
# anno1
ct_counts <- as.data.frame.matrix(table(sobj$name, sobj$anno1))
patient_meta <- meta[!duplicated(meta$name), c("name", "GEM", "age", "sex", "age_group")]
rownames(patient_meta) <- patient_meta$name
patient_meta <- patient_meta[rownames(ct_counts), ]
comp_df <- cbind(patient_meta, ct_counts)
fwrite(comp_df, file.path(out_sccoda, "composition_anno1.csv"), row.names = FALSE)
cat(sprintf("  anno1: %d patients × %d cell types\n", nrow(ct_counts), ncol(ct_counts)))

# anno2
ct_counts2 <- as.data.frame.matrix(table(sobj$name, sobj$anno2))
comp_df2 <- cbind(patient_meta, ct_counts2)
fwrite(comp_df2, file.path(out_sccoda, "composition_anno2.csv"), row.names = FALSE)
cat(sprintf("  anno2: %d patients × %d compartments\n", nrow(ct_counts2), ncol(ct_counts2)))

cat("Done!\n")

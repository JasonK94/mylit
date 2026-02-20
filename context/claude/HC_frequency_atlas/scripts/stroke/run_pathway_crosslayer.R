#!/usr/bin/env Rscript
# ==============================================================================
# Pathway Analysis + Cross-layer Concordance
# Depends on: DEG consensus results + FGS results
# ==============================================================================

.libPaths(c(
  "/home/user3/R/x86_64-pc-linux-gnu-library/4.3",
  "/home/user3/GJC_KDW_250721/renv/library/R-4.3/x86_64-pc-linux-gnu",
  .libPaths()
))

suppressPackageStartupMessages({
  library(qs)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(msigdbr)
  library(UpSetR)
})

source("/home/user3/data_user3/git_repo/mylit/myR/R/plots.R")

cat("=== Pathway Analysis + Cross-layer Concordance ===\n")
cat("Started:", format(Sys.time()), "\n\n")

BASE <- "/data/user3/sobj/stroke_hc_v8_2"
OUT <- file.path(BASE, "pathway")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---- Load DEG results ----
cat("[1/5] Loading DEG results...\n")

load_deg <- function(dir_path) {
  files <- list.files(dir_path, pattern = "consensus.*\\.csv$|results.*\\.csv$",
                      full.names = TRUE, recursive = TRUE)
  if (length(files) == 0) {
    cat("  No DEG results found in:", dir_path, "\n")
    return(NULL)
  }
  results <- list()
  for (f in files) {
    df <- read.csv(f)
    ct <- gsub("(consensus_|results_|\\.csv)", "", basename(f))
    df$celltype <- ct
    results[[ct]] <- df
  }
  do.call(rbind, results)
}

deg_l1 <- load_deg(file.path(BASE, "deg_consensus/L1_cohort_anno1"))
deg_l2 <- load_deg(file.path(BASE, "deg_consensus/L2_g3_anno1"))

if (is.null(deg_l1) && is.null(deg_l2)) {
  cat("  No DEG results available yet. Skipping DEG-based analyses.\n")
  cat("  Run DEG consensus first: scripts/run_deg_and_milo.sh\n")
} else {
  cat("  L1 DEG:", if(!is.null(deg_l1)) nrow(deg_l1) else 0, "results\n")
  cat("  L2 DEG:", if(!is.null(deg_l2)) nrow(deg_l2) else 0, "results\n")
}

# ---- Load FGS results ----
cat("\n[2/5] Loading FGS results...\n")

fgs_n50 <- read.csv(file.path(BASE, "fgs/whole_IS_50_cmgi_genes.csv"))
cat("  FGS n=50:", nrow(fgs_n50), "genes\n")

# Check n=100 if available
fgs_n100_file <- file.path(BASE, "fgs/whole_IS_100_cmgi_genes.csv")
if (file.exists(fgs_n100_file)) {
  fgs_n100 <- read.csv(fgs_n100_file)
  cat("  FGS n=100:", nrow(fgs_n100), "genes\n")
} else {
  fgs_n100 <- NULL
  cat("  FGS n=100: not yet available\n")
}

# ---- Pathway: FGS genes ----
cat("\n[3/5] Pathway enrichment of FGS signature genes...\n")

# HALLMARK gene sets
h_gene_sets <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)

# KEGG
kegg_sets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>%
  dplyr::select(gs_name, gene_symbol)

# GO BP
gobp_sets <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>%
  dplyr::select(gs_name, gene_symbol)

# FGS enrichment (top 50 genes)
fgs_top <- head(fgs_n50$gene, 50)

# Positive contribution genes (higher in g3==2 Bad)
fgs_bad <- fgs_n50 %>% filter(contribution < 0) %>% head(30) %>% pull(gene)
fgs_good <- fgs_n50 %>% filter(contribution > 0) %>% head(30) %>% pull(gene)

for (gs_name in c("HALLMARK", "KEGG", "GOBP")) {
  gs_data <- switch(gs_name,
    "HALLMARK" = h_gene_sets,
    "KEGG" = kegg_sets,
    "GOBP" = gobp_sets
  )

  # Enrichment: all FGS genes
  enr_all <- enricher(fgs_top, TERM2GENE = gs_data, pvalueCutoff = 0.1)
  if (!is.null(enr_all) && nrow(as.data.frame(enr_all)) > 0) {
    p <- dotplot(enr_all, showCategory = 20) +
      ggtitle(paste0("FGS top 50 genes: ", gs_name))
    ggsave(file.path(OUT, paste0("01_fgs_enrichment_", gs_name, ".pdf")), p,
           width = 10, height = 8)
    cat("  Saved:", gs_name, "enrichment\n")
  }

  # Enrichment: Bad-associated genes (negative contribution)
  enr_bad <- enricher(fgs_bad, TERM2GENE = gs_data, pvalueCutoff = 0.1)
  if (!is.null(enr_bad) && nrow(as.data.frame(enr_bad)) > 0) {
    p <- dotplot(enr_bad, showCategory = 15) +
      ggtitle(paste0("Bad-outcome genes (g3==2): ", gs_name))
    ggsave(file.path(OUT, paste0("02_fgs_bad_enrichment_", gs_name, ".pdf")), p,
           width = 10, height = 7)
  }
}

# ---- Cross-layer: FGS vs FAM (FindAllMarkers) overlap ----
cat("\n[4/5] Cross-layer concordance...\n")

# FindAllMarkers results
fam <- read.csv(file.path(BASE, "results/FindMarkers/1_fam.csv"))
fam_sig <- fam %>% filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.5)

# Compare FGS genes vs top cluster markers
fgs_gene_set <- fgs_n50$gene

# Check which anno1 clusters are enriched for FGS genes
cluster_enrichment <- fam_sig %>%
  mutate(is_fgs = gene %in% fgs_gene_set) %>%
  group_by(cluster) %>%
  summarise(
    n_markers = n(),
    n_fgs_overlap = sum(is_fgs),
    fgs_pct = round(n_fgs_overlap / n_markers * 100, 1),
    fgs_genes = paste(gene[is_fgs], collapse = ","),
    .groups = "drop"
  ) %>%
  arrange(desc(fgs_pct))

write.csv(cluster_enrichment,
          file.path(OUT, "03_fgs_cluster_enrichment.csv"), row.names = FALSE)
cat("  Saved: 03_fgs_cluster_enrichment.csv\n")
cat("  Top clusters by FGS enrichment:\n")
print(head(cluster_enrichment, 10))

# If DEG results exist, do cross-layer analysis
if (!is.null(deg_l1) && !is.null(deg_l2)) {
  cat("\n  Cross-layer DEG overlap...\n")

  # Identify column names for gene, log2FC, p-value
  gene_col <- intersect(c("gene", "Gene", "feature"), colnames(deg_l1))[1]
  lfc_col <- intersect(c("avg_log2FC", "logFC", "log2FoldChange", "estimate"), colnames(deg_l1))[1]
  p_col <- intersect(c("p_val_adj", "padj", "FDR", "fdr"), colnames(deg_l1))[1]

  if (!is.na(gene_col) && !is.na(lfc_col) && !is.na(p_col)) {
    # Significant genes per layer
    l1_sig <- deg_l1 %>% filter(.data[[p_col]] < 0.05) %>% pull(gene_col) %>% unique()
    l2_sig <- deg_l2 %>% filter(.data[[p_col]] < 0.05) %>% pull(gene_col) %>% unique()

    # UpSet plot
    gene_lists <- list(
      "HC_vs_IS" = l1_sig,
      "Good_vs_Bad" = l2_sig,
      "FGS_top50" = fgs_gene_set
    )

    png(file.path(OUT, "04_upset_crosslayer.png"), width = 1200, height = 800)
    upset(fromList(gene_lists), order.by = "freq", nsets = 3,
          mainbar.y.label = "Intersection Size",
          sets.x.label = "Set Size")
    dev.off()
    cat("  Saved: 04_upset_crosslayer.png\n")

    # Scatter plot: L1 logFC vs L2 logFC (shared genes)
    shared_genes <- intersect(
      deg_l1 %>% filter(!is.na(.data[[lfc_col]])) %>% pull(gene_col),
      deg_l2 %>% filter(!is.na(.data[[lfc_col]])) %>% pull(gene_col)
    )
    if (length(shared_genes) > 10) {
      l1_lfc <- deg_l1 %>%
        filter(.data[[gene_col]] %in% shared_genes) %>%
        group_by(across(all_of(gene_col))) %>%
        summarise(lfc_l1 = mean(.data[[lfc_col]], na.rm = TRUE), .groups = "drop")

      l2_lfc <- deg_l2 %>%
        filter(.data[[gene_col]] %in% shared_genes) %>%
        group_by(across(all_of(gene_col))) %>%
        summarise(lfc_l2 = mean(.data[[lfc_col]], na.rm = TRUE), .groups = "drop")

      merged <- inner_join(l1_lfc, l2_lfc, by = gene_col)
      r <- cor(merged$lfc_l1, merged$lfc_l2, use = "complete.obs")

      p_scatter <- ggplot(merged, aes(x = lfc_l1, y = lfc_l2)) +
        geom_point(alpha = 0.5, size = 1) +
        geom_smooth(method = "lm", se = TRUE, color = "red") +
        geom_hline(yintercept = 0, linetype = "dashed") +
        geom_vline(xintercept = 0, linetype = "dashed") +
        labs(x = "log2FC (HC vs Stroke)", y = "log2FC (Good vs Bad)",
             title = paste0("Cross-layer concordance (r = ", round(r, 3), ")")) +
        theme_bw()

      ggsave(file.path(OUT, "05_crosslayer_lfc_scatter.pdf"), p_scatter,
             width = 8, height = 8)
      cat("  Saved: 05_crosslayer_lfc_scatter.pdf (r =", round(r, 3), ")\n")
    }
  }
}

# ---- Summary ----
cat("\n[5/5] Summary...\n")

summary_text <- c(
  paste("FGS n=50 top genes:", paste(head(fgs_n50$gene, 10), collapse = ", ")),
  paste("FGS bad-outcome genes:", paste(head(fgs_bad, 10), collapse = ", ")),
  paste("FGS good-outcome genes:", paste(head(fgs_good, 10), collapse = ", ")),
  paste("FAM total significant:", nrow(fam_sig)),
  paste("Cluster most enriched for FGS:", cluster_enrichment$cluster[1],
        "(", cluster_enrichment$fgs_pct[1], "%)"),
  ""
)
writeLines(summary_text, file.path(OUT, "analysis_summary.txt"))

cat("\n=== PATHWAY + CROSS-LAYER COMPLETE ===\n")
cat("Finished:", format(Sys.time()), "\n")
cat("Output:", OUT, "\n")

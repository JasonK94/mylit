#!/usr/bin/env Rscript
# ==============================================================================
# Fix interaction p-values for gene dynamics v3
# Re-fits ANOVA comparison with corrected column name extraction
# Reads existing summary CSVs, re-computes p-values, updates CSVs
# ==============================================================================

.libPaths(c(
  "/home/user3/R/x86_64-pc-linux-gnu-library/4.3",
  "/home/user3/GJC_KDW_250721/renv/library/R-4.3/x86_64-pc-linux-gnu",
  .libPaths()
))

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(monocle3)
  library(S4Vectors)
  library(SummarizedExperiment)
  library(Seurat)
  library(SeuratObject)
  library(qs)
  library(mgcv)
  library(dplyr)
})

cat("=== Fix Gene Dynamics v3 P-values ===\n")
cat("Started:", format(Sys.time()), "\n\n")

BASE <- "/data/user3/sobj/stroke_hc_v8_2"
TRAJ_V2 <- file.path(BASE, "trajectory_v2")
TRAJ_V3 <- file.path(BASE, "trajectory_v3")

compartments <- c("mono", "cd4", "cd8")
conditions <- list(
  list(name = "g3", col = "g3", filter_col = "index_injury_hand", filter_val = "IS"),
  list(name = "cohort", col = "cohort", filter_col = NULL, filter_val = NULL)
)

for (comp in compartments) {
  cat(sprintf("\n=== Compartment: %s ===\n", comp))

  sobj_path <- file.path(TRAJ_V2, comp, "subset.qs")
  cds_path <- file.path(TRAJ_V2, comp, "monocle3_cds.rds")

  if (!file.exists(sobj_path) || !file.exists(cds_path)) {
    cat("  SKIP: Missing files\n")
    next
  }

  sobj <- qread(sobj_path)
  cds <- readRDS(cds_path)

  pt <- pseudotime(cds)
  pt[is.infinite(pt)] <- NA

  cells_with_pt <- intersect(names(pt), colnames(sobj))
  sobj$pseudotime <- NA_real_
  sobj$pseudotime[match(cells_with_pt, colnames(sobj))] <- pt[cells_with_pt]

  for (cond_info in conditions) {
    cond_name <- cond_info$name
    cond_col <- cond_info$col
    summary_dir <- file.path(TRAJ_V3, comp, paste0("gene_dynamics_", cond_name))
    summary_csv <- file.path(summary_dir, "gene_dynamics_v3_summary.csv")

    if (!file.exists(summary_csv)) {
      cat(sprintf("  SKIP %s/%s: No summary CSV\n", comp, cond_name))
      next
    }

    cat(sprintf("  --- %s/%s ---\n", comp, cond_name))
    summary_df <- read.csv(summary_csv)

    # Prepare data
    meta <- sobj@meta.data
    meta$pseudotime <- sobj$pseudotime
    meta <- meta[!is.na(meta$pseudotime), ]

    if (!is.null(cond_info$filter_col)) {
      meta <- meta[meta[[cond_info$filter_col]] == cond_info$filter_val, ]
    }
    meta <- meta[!is.na(meta[[cond_col]]), ]

    dat <- data.frame(
      pseudotime = meta$pseudotime,
      cond = factor(meta[[cond_col]]),
      log_lib_size = log(meta$nCount_RNA),
      percent.mt = meta$percent.mt,
      batch = factor(meta$GEM),
      stringsAsFactors = FALSE
    )

    # Re-extract p-values for each gene
    counts <- GetAssayData(sobj, slot = "counts")
    new_pvals <- rep(NA_real_, nrow(summary_df))
    k_val <- 8

    for (j in seq_len(nrow(summary_df))) {
      gene <- summary_df$gene[j]
      if (!(gene %in% rownames(counts))) next

      gene_expr <- counts[gene, rownames(meta)]
      gene_dat <- dat
      gene_dat$expr <- as.numeric(gene_expr)

      tryCatch({
        fit <- bam(
          expr ~ s(pseudotime, k = k_val, bs = "cr") + cond +
            s(pseudotime, by = cond, k = k_val, bs = "cr") +
            offset(log_lib_size) + percent.mt + s(batch, bs = "re"),
          data = gene_dat, family = nb(), method = "fREML",
          discrete = TRUE, nthreads = 2
        )

        fit0 <- bam(
          expr ~ s(pseudotime, k = k_val, bs = "cr") + cond +
            offset(log_lib_size) + percent.mt + s(batch, bs = "re"),
          data = gene_dat, family = nb(), method = "fREML",
          discrete = TRUE, nthreads = 2
        )

        anova_res <- anova(fit0, fit, test = "Chisq")
        p_col <- if ("Pr(>Chi)" %in% names(anova_res)) "Pr(>Chi)" else "P(>|Chi|)"
        if (!is.null(anova_res[[p_col]]) && length(anova_res[[p_col]]) >= 2) {
          new_pvals[j] <- anova_res[[p_col]][2]
        }
        cat(sprintf("    %s: p=%s\n", gene, format(new_pvals[j], digits = 3)))
      }, error = function(e) {
        cat(sprintf("    %s: ERROR - %s\n", gene, e$message))
      })
    }

    # Update summary
    summary_df$interaction_p <- new_pvals
    summary_df$interaction_padj <- p.adjust(new_pvals, method = "BH")

    write.csv(summary_df, summary_csv, row.names = FALSE)
    cat(sprintf("  Updated: %s (%d genes, %d significant at padj<0.05)\n",
                summary_csv, nrow(summary_df),
                sum(summary_df$interaction_padj < 0.05, na.rm = TRUE)))
  }

  rm(sobj, cds); gc(verbose = FALSE)
}

cat("\n=== P-VALUE FIX COMPLETE ===\n")
cat("Finished:", format(Sys.time()), "\n")

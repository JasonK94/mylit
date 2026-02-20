#!/usr/bin/env Rscript
# HC-only MILO: Age (continuous) neighbourhood differential abundance
# Reuses cached nhood object from sex analysis to avoid rebuilding nhoods
#
# Usage: Rscript scripts/hc/run_milo_age.R
# Output: /data/user3/sobj/hc_only_v1/milo_age/

.libPaths(c(
  "/home/user3/GJC_KDW_250721/renv/library/R-4.3/x86_64-pc-linux-gnu",
  "/home/user3/R/x86_64-pc-linux-gnu-library/4.3",
  .libPaths()
))

suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(miloR)
  library(SingleCellExperiment)
  library(ggplot2)
  library(edgeR)
  library(dplyr)
})

cat("=== HC-only MILO: Age Analysis ===\n")
cat(sprintf("Time: %s\n\n", Sys.time()))

out_dir <- "/data/user3/sobj/hc_only_v1/milo_age"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Load cached MILO object (with nhoods + distances) ----
milo_cache <- "/data/user3/sobj/hc_only_v1/milo/milo_sex_02_distances_calculated.qs"
if (!file.exists(milo_cache)) {
  milo_cache <- "/data/user3/sobj/hc_only_v1/milo/milo_sex_02_distances_calculated_1.qs"
}
cat(sprintf("Loading cached Milo: %s\n", milo_cache))
milo <- qread(milo_cache)
cat(sprintf("  %d cells, %d nhoods\n", ncol(milo), ncol(nhoods(milo))))
cat(sprintf("  Reductions: %s\n", paste(reducedDimNames(milo), collapse = ", ")))

# ---- Ensure metadata ----
meta <- colData(milo)
cat(sprintf("  Patients: %d\n", length(unique(meta$name))))

# Age as continuous variable
if (!"age" %in% colnames(meta)) {
  stop("'age' column not in Milo colData")
}

# Ensure sex is factor
if (is.character(meta$sex)) meta$sex <- factor(meta$sex)

# Scale age for numerical stability
meta$age_scaled <- scale(as.numeric(meta$age))[,1]
colData(milo) <- meta

# ---- Count cells per nhood per patient ----
cat("\n--- Counting cells per neighbourhood ---\n")
milo <- countCells(milo, meta.data = as.data.frame(colData(milo)), sample = "name")
cat(sprintf("  Nhood count matrix: %d nhoods x %d patients\n",
            nrow(nhoodCounts(milo)), ncol(nhoodCounts(milo))))

# ---- DA Testing: ~ age_scaled + sex ----
cat("\n--- Testing: ~ age_scaled + sex ---\n")
tryCatch({
  design_df <- data.frame(
    name = meta$name,
    age_scaled = meta$age_scaled,
    sex = meta$sex
  )
  # Unique per patient
  design_df <- design_df[!duplicated(design_df$name), ]
  rownames(design_df) <- design_df$name

  # Build design matrix
  design <- model.matrix(~ age_scaled + sex, data = design_df)

  da_results <- testNhoods(
    milo,
    design = design,
    design.df = design_df,
    model.contrasts = "age_scaled",
    reduced.dim = "GRAPH",
    fdr.weighting = "graph-overlap"
  )

  # Add anno1 labels
  nhood_sce <- nhoodExpression(milo)
  # Get the major cell type per nhood
  nhood_anno <- data.frame(Nhood = seq_len(ncol(nhoods(milo))))
  nhood_members <- nhoods(milo)
  anno1_vec <- as.character(meta$anno1)

  major_ct <- sapply(seq_len(ncol(nhood_members)), function(i) {
    members <- which(nhood_members[, i] == 1)
    if (length(members) == 0) return(NA)
    ct_tab <- table(anno1_vec[members])
    names(ct_tab)[which.max(ct_tab)]
  })
  da_results$anno1 <- major_ct

  cat(sprintf("  Total nhoods: %d\n", nrow(da_results)))
  n_sig_01 <- sum(da_results$SpatialFDR < 0.1, na.rm = TRUE)
  n_sig_05 <- sum(da_results$SpatialFDR < 0.05, na.rm = TRUE)
  cat(sprintf("  Significant (SpatialFDR<0.1): %d\n", n_sig_01))
  cat(sprintf("  Significant (SpatialFDR<0.05): %d\n", n_sig_05))

  # Direction
  sig <- da_results[da_results$SpatialFDR < 0.1 & !is.na(da_results$SpatialFDR), ]
  if (nrow(sig) > 0) {
    cat(sprintf("  logFC>0 (increase with age): %d\n", sum(sig$logFC > 0)))
    cat(sprintf("  logFC<0 (decrease with age): %d\n", sum(sig$logFC < 0)))

    cat("\n  Per cell type (SpatialFDR<0.1):\n")
    ct_summary <- sig %>%
      group_by(anno1) %>%
      summarise(
        n_nhoods = n(),
        n_up = sum(logFC > 0),
        n_down = sum(logFC < 0),
        median_logFC = median(logFC),
        .groups = "drop"
      ) %>%
      arrange(desc(n_nhoods))
    for (i in seq_len(nrow(ct_summary))) {
      r <- ct_summary[i, ]
      cat(sprintf("    %-30s n=%d up=%d down=%d med_logFC=%.2f\n",
                  r$anno1, r$n_nhoods, r$n_up, r$n_down, r$median_logFC))
    }
    write.csv(ct_summary, file.path(out_dir, "milo_age_summary_by_celltype.csv"), row.names = FALSE)
  }

  # Save
  write.csv(da_results, file.path(out_dir, "milo_age_da_results.csv"), row.names = FALSE)
  qsave(da_results, file.path(out_dir, "milo_age_da_results.qs"))

  # ---- Plots ----
  # Volcano
  tryCatch({
    da_results$neg_log10_fdr <- -log10(da_results$SpatialFDR + 1e-300)
    p_vol <- ggplot(da_results, aes(x = logFC, y = neg_log10_fdr)) +
      geom_point(aes(color = anno1), alpha = 0.5, size = 1) +
      geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "red") +
      geom_vline(xintercept = 0, linetype = "dashed") +
      theme_minimal(base_size = 11) +
      theme(panel.background = element_rect(fill = "white", color = NA),
            plot.background = element_rect(fill = "white", color = NA),
            legend.text = element_text(size = 7)) +
      labs(title = "MILO: Age effect on neighbourhood abundance",
           x = "logFC (+ = increase with age)", y = "-log10(SpatialFDR)", color = "Cell Type")
    ggsave(file.path(out_dir, "milo_age_volcano.png"), p_vol,
           width = 14, height = 8, dpi = 200, bg = "white")
    cat("  Saved: milo_age_volcano.png\n")
  }, error = function(e) cat(sprintf("  Volcano failed: %s\n", e$message)))

  # Boxplot by cell type
  tryCatch({
    ct_order <- da_results %>%
      group_by(anno1) %>%
      summarise(med = median(logFC, na.rm = TRUE), .groups = "drop") %>%
      arrange(med) %>%
      pull(anno1)
    da_results$anno1 <- factor(da_results$anno1, levels = ct_order)

    p_box <- ggplot(da_results, aes(x = anno1, y = logFC)) +
      geom_boxplot(aes(fill = anno1), outlier.size = 0.5) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      coord_flip() +
      theme_minimal(base_size = 11) +
      theme(panel.background = element_rect(fill = "white", color = NA),
            plot.background = element_rect(fill = "white", color = NA),
            legend.position = "none") +
      labs(title = "MILO: logFC distribution by cell type (age effect)",
           x = "", y = "logFC (+ = increase with age)")
    ggsave(file.path(out_dir, "milo_age_boxplot.png"), p_box,
           width = 10, height = 10, dpi = 200, bg = "white")
    cat("  Saved: milo_age_boxplot.png\n")
  }, error = function(e) cat(sprintf("  Boxplot failed: %s\n", e$message)))

  # Significant nhoods barplot
  tryCatch({
    if (nrow(sig) > 0) {
      sig_ct <- sig %>%
        mutate(direction = ifelse(logFC > 0, "Increase with age", "Decrease with age")) %>%
        group_by(anno1, direction) %>%
        summarise(n = n(), .groups = "drop")

      p_bar <- ggplot(sig_ct, aes(x = reorder(anno1, n), y = n, fill = direction)) +
        geom_col(position = "dodge") +
        coord_flip() +
        scale_fill_manual(values = c("Increase with age" = "firebrick", "Decrease with age" = "steelblue")) +
        theme_minimal(base_size = 11) +
        theme(panel.background = element_rect(fill = "white", color = NA),
              plot.background = element_rect(fill = "white", color = NA)) +
        labs(title = "MILO: Significant nhoods by cell type (age)",
             x = "", y = "# significant nhoods (SpatialFDR<0.1)", fill = "")
      ggsave(file.path(out_dir, "milo_age_sig_barplot.png"), p_bar,
             width = 10, height = 8, dpi = 200, bg = "white")
      cat("  Saved: milo_age_sig_barplot.png\n")
    }
  }, error = function(e) cat(sprintf("  Barplot failed: %s\n", e$message)))

}, error = function(e) {
  cat(sprintf("\n  MILO Age DA FAILED: %s\n", e$message))
  cat(sprintf("  Call: %s\n", deparse(e$call)))
})

cat(sprintf("\n=== Done: %s ===\n", Sys.time()))

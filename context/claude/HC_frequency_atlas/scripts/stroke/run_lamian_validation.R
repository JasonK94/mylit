#!/usr/bin/env Rscript
# ==============================================================================
# Phase 3: Lamian Validation of Gene Dynamics v3 GAMM Results
# ==============================================================================
# Runs Lamian XDE test (multi-sample differential expression along pseudotime)
# to validate GAMM v3 findings. Lamian properly handles pseudoreplication via
# sample-level random effects in a mixed-effects B-spline framework.
#
# Reference: Hou et al., Nature Communications, 2023
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
  library(Matrix)
  library(Lamian)
  library(dplyr)
})

cat("=== Phase 3: Lamian Validation ===\n")
cat("Started:", format(Sys.time()), "\n\n")

BASE <- "/data/user3/sobj/stroke_hc_v8_2"
TRAJ_V2 <- file.path(BASE, "trajectory_v2")
OUT_BASE <- file.path(BASE, "trajectory_v3", "lamian")
ANALYSIS_DIR <- file.path(BASE, "trajectory_v3", "analysis")
dir.create(OUT_BASE, showWarnings = FALSE, recursive = TRUE)

# Target genes (same as v3)
TARGET_GENES <- c(
  "S100A8", "S100A9", "S100A12", "VCAN", "FCN1", "CD14", "FCGR3A",
  "IL1B", "CXCL8", "CCL2", "CCL3", "TNF", "IL6",
  "ISG15", "IFI6", "IFIT1", "MX1", "OAS1",
  "GZMB", "GZMK", "PRF1", "NKG7", "GNLY", "IFNG",
  "FOXP3", "IL2RA", "CTLA4", "TIGIT", "LAG3", "PDCD1",
  "DDIT4", "TXNIP", "HLA-B", "HLA-C", "CCL4", "HIF1A"
)

compartments <- c("mono", "cd4", "cd8")

# ==============================================================================
# Helper: Prepare Lamian inputs from Monocle3 CDS
# ==============================================================================
prepare_lamian_input <- function(cds, condition_col, subset_expr = NULL,
                                  target_genes, min_pct = 0.01) {
  cell_meta <- as.data.frame(colData(cds))

  # Subset if needed (e.g., g3 requires filtering NAs)
  if (!is.null(subset_expr)) {
    keep_mask <- eval(subset_expr, envir = list(cell_meta = cell_meta))
    keep_cells <- rownames(cell_meta)[keep_mask]
    cds <- cds[, keep_cells]
    cell_meta <- cell_meta[keep_cells, ]
  }

  # Get pseudotime
  pt <- pseudotime(cds)
  pt <- pt[!is.na(pt) & is.finite(pt)]
  cell_meta <- cell_meta[names(pt), ]
  cds <- cds[, names(pt)]

  # Get condition and sample columns
  # Use patient_name as sample (biological replicate), NOT GEM (batch)
  condition <- cell_meta[[condition_col]]
  names(condition) <- rownames(cell_meta)
  sample_id <- cell_meta$patient_name

  cat("    Cells:", length(pt), "\n")
  cat("    Patients:", length(unique(sample_id)), "\n")

  # Filter genes: must be in target list AND expressed in >= min_pct of cells
  available_genes <- intersect(target_genes, rownames(cds))

  # Get expression matrix and log-normalize
  # CDS contains raw counts; Lamian expects library-size-normalized log-transformed
  raw_counts <- as.matrix(counts(cds)[available_genes, , drop = FALSE])
  lib_sizes <- colSums(counts(cds))
  # Log-normalize: log1p(count / lib_size * 10000)
  expr_mat <- log1p(sweep(raw_counts, 2, lib_sizes, "/") * 10000)

  # Filter genes by expression percentage (per condition)
  conds <- sort(unique(condition))
  keep_genes <- sapply(available_genes, function(g) {
    pcts <- sapply(conds, function(cv) {
      cells_c <- names(condition)[condition == cv]
      mean(expr_mat[g, cells_c] > 0)
    })
    all(pcts >= min_pct)
  })
  available_genes <- available_genes[keep_genes]
  expr_mat <- expr_mat[available_genes, , drop = FALSE]

  cat("    Genes after sparsity filter:", length(available_genes), "/",
      length(target_genes), "\n")

  # Filter out patients with very few cells (< 10)
  patient_cell_counts <- table(sample_id)
  good_patients <- names(patient_cell_counts)[patient_cell_counts >= 10]
  keep_cells2 <- names(pt)[sample_id %in% good_patients]
  if (length(keep_cells2) < length(pt)) {
    n_dropped <- length(pt) - length(keep_cells2)
    cat("    Dropped", n_dropped, "cells from",
        sum(patient_cell_counts < 10), "patients with <10 cells\n")
    pt <- pt[keep_cells2]
    cell_meta <- cell_meta[keep_cells2, ]
    expr_mat <- expr_mat[, keep_cells2, drop = FALSE]
    condition <- condition[keep_cells2]
    sample_id <- cell_meta$patient_name
  }

  # cellanno: data.frame(cell = cell_names, sample = sample_names)
  cellanno <- data.frame(
    cell = names(pt),
    sample = as.character(sample_id),
    stringsAsFactors = FALSE
  )

  # Design matrix: intercept + condition (binary 0/1)
  # Each patient has one condition value (uniform within patient)
  samples <- unique(cellanno$sample)
  sample_condition <- sapply(samples, function(s) {
    cells_s <- cellanno$cell[cellanno$sample == s]
    tab <- table(condition[cells_s])
    names(tab)[which.max(tab)]
  })

  design <- data.frame(
    intercept = 1,
    condition = as.numeric(sample_condition == conds[2])  # second level = 1
  )
  rownames(design) <- samples
  design <- as.matrix(design)

  n_cond1 <- sum(design[, "condition"] == 0)
  n_cond2 <- sum(design[, "condition"] == 1)
  cat("    Samples:", length(samples),
      "(", conds[1], "=", n_cond1, ",", conds[2], "=", n_cond2, ")\n")

  list(
    expr = expr_mat,
    cellanno = cellanno,
    pseudotime = pt,
    design = design,
    genes = available_genes,
    cond_levels = conds,
    n_cells = length(pt)
  )
}

# ==============================================================================
# Main loop: each compartment × condition
# ==============================================================================
all_results <- list()

for (comp_id in compartments) {
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("  Compartment:", comp_id, "\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")

  cds_file <- file.path(TRAJ_V2, comp_id, "monocle3_cds.rds")
  if (!file.exists(cds_file)) {
    cat("  CDS not found:", cds_file, "\n  Skipping.\n")
    next
  }

  cat("  Loading CDS...\n")
  cds <- readRDS(cds_file)
  cell_meta <- as.data.frame(colData(cds))
  cat("  Cells:", ncol(cds), "\n")

  # Inject missing columns from Seurat subset if needed
  required <- c("nCount_RNA", "percent.mt", "GEM")
  missing <- setdiff(required, colnames(cell_meta))
  if (length(missing) > 0) {
    sub_file <- file.path(TRAJ_V2, comp_id, "subset.qs")
    if (file.exists(sub_file)) {
      sub <- qread(sub_file)
      sub_meta <- sub@meta.data
      common_cells <- intersect(rownames(cell_meta), rownames(sub_meta))
      for (col in missing) {
        if (col %in% colnames(sub_meta)) {
          vals <- sub_meta[common_cells, col]
          names(vals) <- common_cells
          colData(cds)[[col]] <- vals[rownames(cell_meta)]
          cat("    Injected:", col, "\n")
        }
      }
      rm(sub); gc(verbose = FALSE)
    }
    cell_meta <- as.data.frame(colData(cds))
  }

  # Check conditions
  has_g3 <- "g3" %in% colnames(cell_meta) && sum(!is.na(cell_meta$g3)) > 100
  has_cohort <- "cohort" %in% colnames(cell_meta) && length(unique(cell_meta$cohort)) >= 2

  analyses <- list()
  if (has_g3) {
    analyses[["g3"]] <- list(
      condition_col = "g3",
      label = "g3 (IS: Bad vs Good)",
      subset_expr = quote(!is.na(cell_meta$g3))
    )
  }
  if (has_cohort) {
    analyses[["cohort"]] <- list(
      condition_col = "cohort",
      label = "Cohort (HC vs Stroke)",
      subset_expr = NULL
    )
  }

  for (analysis_name in names(analyses)) {
    analysis <- analyses[[analysis_name]]
    cat("\n  --- Lamian XDE:", analysis$label, "---\n")

    out_dir <- file.path(OUT_BASE, comp_id, analysis_name)
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

    # Check if already completed
    result_file <- file.path(out_dir, "lamian_result.rds")
    if (file.exists(result_file)) {
      cat("    Result already exists, loading...\n")
      lamian_res <- readRDS(result_file)
      all_results[[paste0(comp_id, "_", analysis_name)]] <- lamian_res
      next
    }

    # Prepare input
    cat("  Preparing Lamian input...\n")
    input <- tryCatch(
      prepare_lamian_input(
        cds = cds,
        condition_col = analysis$condition_col,
        subset_expr = analysis$subset_expr,
        target_genes = TARGET_GENES,
        min_pct = 0.01
      ),
      error = function(e) {
        cat("    ERROR preparing input:", e$message, "\n")
        NULL
      }
    )

    if (is.null(input)) next

    cat("    Expression matrix:", nrow(input$expr), "genes ×", ncol(input$expr), "cells\n")
    cat("    Pseudotime range:", round(range(input$pseudotime), 3), "\n")
    cat("    Condition levels:", paste(input$cond_levels, collapse = " vs "), "\n")

    # Run Lamian XDE test
    cat("\n    Running Lamian XDE test (test.type='Variable')...\n")
    cat("    Chi-squared test with sample-level random effects.\n")
    cat("    Start:", format(Sys.time()), "\n")

    lamian_res <- tryCatch({
      lamian_test(
        expr = input$expr,
        cellanno = input$cellanno,
        pseudotime = input$pseudotime,
        design = input$design,
        testvar = 2,
        test.type = "Variable",
        test.method = "chisq",    # chi-squared: ~100x faster than permutation
        ncores = 8,
        ncores.fit = 4,
        return.all.data = TRUE,
        overall.only = FALSE,     # Get trend + mean FDR too
        verbose.output = FALSE
      )
    }, error = function(e) {
      cat("    LAMIAN ERROR:", e$message, "\n")
      NULL
    })

    if (is.null(lamian_res)) {
      cat("    Lamian failed for", comp_id, "/", analysis_name, "\n")
      next
    }

    cat("    Lamian completed:", format(Sys.time()), "\n")

    # Save result
    saveRDS(lamian_res, result_file)
    cat("    Saved:", result_file, "\n")

    # Extract XDE classification
    xde_type <- tryCatch(
      getXDEType(lamian_res, cutoff = 0.05),
      error = function(e) {
        cat("    getXDEType error:", e$message, "\n")
        NULL
      }
    )

    # Extract p-values and FDR
    if (!is.null(lamian_res$statistics)) {
      stats_df <- as.data.frame(lamian_res$statistics)
      stats_df$gene <- rownames(stats_df)

      if (!is.null(xde_type)) {
        stats_df$xde_type <- xde_type[stats_df$gene]
      }

      stats_df$compartment <- comp_id
      stats_df$analysis <- analysis_name

      write.csv(stats_df, file.path(out_dir, "lamian_statistics.csv"), row.names = FALSE)
      cat("    Statistics saved:", nrow(stats_df), "genes\n")

      # Print summary — handle both permutation and chisq column names
      fdr_col <- intersect(
        c("fdr.overall", "fdr.chisq.overall"),
        colnames(stats_df)
      )[1]
      if (!is.na(fdr_col)) {
        n_sig <- sum(stats_df[[fdr_col]] < 0.05, na.rm = TRUE)
        cat("    Significant (FDR<0.05):", n_sig, "/",
            nrow(stats_df), "\n")
      }
      if (!is.null(xde_type)) {
        cat("    XDE types:", table(xde_type), "\n")
      }
    }

    all_results[[paste0(comp_id, "_", analysis_name)]] <- lamian_res
  }

  rm(cds); gc(verbose = FALSE)
}

# ==============================================================================
# Concordance analysis: Lamian vs GAMM v3
# ==============================================================================
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("  Concordance Analysis: Lamian vs GAMM v3\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

# Load GAMM v3 effect sizes
gamm_file <- file.path(ANALYSIS_DIR, "effect_sizes_all.csv")
if (file.exists(gamm_file)) {
  gamm_df <- read.csv(gamm_file, stringsAsFactors = FALSE)
} else {
  cat("GAMM effect sizes not found:", gamm_file, "\n")
  gamm_df <- NULL
}

concordance_results <- list()

for (key in names(all_results)) {
  lamian_res <- all_results[[key]]
  parts <- strsplit(key, "_")[[1]]
  comp_id <- parts[1]
  analysis_name <- parts[2]

  cat("\n  ---", key, "---\n")

  # Extract Lamian statistics
  lamian_stats <- NULL
  if (!is.null(lamian_res$statistics)) {
    lamian_stats <- as.data.frame(lamian_res$statistics)
    lamian_stats$gene <- rownames(lamian_stats)
  }

  if (is.null(lamian_stats) || is.null(gamm_df)) {
    cat("    Skipping concordance (missing data)\n")
    next
  }

  # Get GAMM results for this compartment/analysis
  gamm_sub <- gamm_df %>%
    filter(compartment == comp_id, analysis == analysis_name, status == "success")

  # Merge
  common_genes <- intersect(lamian_stats$gene, gamm_sub$gene)
  cat("    Common genes:", length(common_genes), "\n")

  if (length(common_genes) < 3) {
    cat("    Too few common genes, skipping.\n")
    next
  }

  merged <- merge(
    lamian_stats[lamian_stats$gene %in% common_genes, ],
    gamm_sub[gamm_sub$gene %in% common_genes, c("gene", "interaction_padj", "rABC", "ABC", "signed_area")],
    by = "gene"
  )

  # Concordance metrics — handle both permutation and chisq columns
  fdr_col <- intersect(
    c("fdr.overall", "fdr.chisq.overall"),
    colnames(merged)
  )[1]
  if (!is.na(fdr_col)) {
    lamian_sig <- merged[[fdr_col]] < 0.05
    gamm_sig <- merged$interaction_padj < 0.05

    both_sig <- sum(lamian_sig & gamm_sig, na.rm = TRUE)
    lamian_only <- sum(lamian_sig & !gamm_sig, na.rm = TRUE)
    gamm_only <- sum(!lamian_sig & gamm_sig, na.rm = TRUE)
    neither <- sum(!lamian_sig & !gamm_sig, na.rm = TRUE)

    # Jaccard
    jaccard <- both_sig / (both_sig + lamian_only + gamm_only)
    if (is.nan(jaccard)) jaccard <- NA

    cat("    Concordance table:\n")
    cat("      Both significant:   ", both_sig, "\n")
    cat("      Lamian only:        ", lamian_only, "\n")
    cat("      GAMM v3 only:       ", gamm_only, "\n")
    cat("      Neither:            ", neither, "\n")
    cat("      Jaccard index:      ", round(jaccard, 3), "\n")

    concordance_results[[key]] <- data.frame(
      compartment = comp_id,
      analysis = analysis_name,
      n_genes = length(common_genes),
      both_sig = both_sig,
      lamian_only = lamian_only,
      gamm_only = gamm_only,
      neither = neither,
      jaccard = round(jaccard, 3),
      stringsAsFactors = FALSE
    )
  } else {
    cat("    No FDR column found in merged data\n")
  }

  # Save merged table
  out_dir <- file.path(OUT_BASE, comp_id, analysis_name)
  write.csv(merged, file.path(out_dir, "concordance_merged.csv"), row.names = FALSE)
}

# Save concordance summary
if (length(concordance_results) > 0) {
  conc_df <- do.call(rbind, concordance_results)
  write.csv(conc_df, file.path(OUT_BASE, "concordance_summary.csv"), row.names = FALSE)
  cat("\n\nConcordance summary:\n")
  print(conc_df)
}

# ==============================================================================
# Concordance figure
# ==============================================================================
library(ggplot2)

if (length(concordance_results) > 0) {
  conc_df <- do.call(rbind, concordance_results)

  # Bar plot of Jaccard indices
  p <- ggplot(conc_df, aes(x = paste0(compartment, "/", analysis), y = jaccard)) +
    geom_col(fill = "steelblue", width = 0.6) +
    geom_text(aes(label = sprintf("%.2f", jaccard)), vjust = -0.3, size = 4) +
    ylim(0, 1) +
    labs(
      title = "Lamian vs GAMM v3: Concordance (Jaccard Index)",
      subtitle = "Proportion of genes significant in both methods",
      x = "Analysis", y = "Jaccard Index"
    ) +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(OUT_BASE, "concordance_jaccard.png"), p, width = 8, height = 5, dpi = 300)
  cat("Concordance figure saved.\n")
}

# ==============================================================================
# Scatter: Lamian -log10(FDR) vs GAMM -log10(padj)
# ==============================================================================
scatter_plots <- list()
for (key in names(all_results)) {
  lamian_res <- all_results[[key]]
  parts <- strsplit(key, "_")[[1]]
  comp_id <- parts[1]
  analysis_name <- parts[2]

  out_dir <- file.path(OUT_BASE, comp_id, analysis_name)
  merged_file <- file.path(out_dir, "concordance_merged.csv")
  if (!file.exists(merged_file)) next

  merged <- read.csv(merged_file, stringsAsFactors = FALSE)
  fdr_col <- intersect(
    c("fdr.overall", "fdr.chisq.overall"),
    colnames(merged)
  )[1]
  if (is.na(fdr_col)) next

  merged$lamian_logp <- -log10(pmax(merged[[fdr_col]], 1e-300))
  merged$gamm_logp <- -log10(pmax(merged$interaction_padj, 1e-300))

  # Cap at reasonable values for visualization
  merged$lamian_logp <- pmin(merged$lamian_logp, 20)
  merged$gamm_logp <- pmin(merged$gamm_logp, 60)

  # Spearman correlation
  cor_test <- cor.test(merged$lamian_logp, merged$gamm_logp, method = "spearman")

  p <- ggplot(merged, aes(x = gamm_logp, y = lamian_logp)) +
    geom_point(aes(size = abs(rABC)), alpha = 0.7) +
    geom_text(aes(label = gene), size = 2.5, vjust = -0.8, check_overlap = TRUE) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "red", alpha = 0.5) +
    geom_vline(xintercept = -log10(0.05), linetype = "dotted", color = "red", alpha = 0.5) +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("rho = %.2f\np = %.2g", cor_test$estimate, cor_test$p.value),
             hjust = 1.1, vjust = 1.3, size = 4) +
    labs(
      title = paste0("Lamian vs GAMM v3: ", comp_id, "/", analysis_name),
      x = "GAMM v3 -log10(padj)",
      y = "Lamian -log10(FDR)",
      size = "|rABC|"
    ) +
    theme_minimal(base_size = 12)

  ggsave(file.path(out_dir, "pvalue_scatter.png"), p, width = 7, height = 6, dpi = 300)
  scatter_plots[[key]] <- p
}

cat("\n\n=== Lamian Validation COMPLETE ===\n")
cat("Finished:", format(Sys.time()), "\n")
cat("Output:", OUT_BASE, "\n")

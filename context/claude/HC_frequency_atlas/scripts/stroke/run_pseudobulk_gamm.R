#!/usr/bin/env Rscript
# ==============================================================================
# Phase 4: Pseudobulk GAMM Sensitivity Analysis
# ==============================================================================
# Addresses pseudoreplication by aggregating cells to patient × pseudotime bin.
# Instead of fitting cell-level GAMM (inflate p-values from thousands of cells
# per patient), we aggregate into pseudobulk profiles per patient.
#
# Approach:
#   1. Bin pseudotime into 10 quantile-based bins per compartment
#   2. For each gene × patient × bin: compute mean log-normalized expression
#   3. Fit per-gene: lmer(mean_expr ~ s(pt_bin, by=condition) + (1|patient))
#   4. Test interaction term (condition × pseudotime)
#   5. Compute ABC from pseudobulk-level fitted curves
#   6. Compare rankings with cell-level GAMM v3
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
  library(qs)
  library(Matrix)
  library(mgcv)
  library(ggplot2)
  library(dplyr)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

cat("=== Phase 4: Pseudobulk GAMM ===\n")
cat("Started:", format(Sys.time()), "\n\n")

BASE <- "/data/user3/sobj/stroke_hc_v8_2"
TRAJ_V2 <- file.path(BASE, "trajectory_v2")
OUT_BASE <- file.path(BASE, "trajectory_v3", "pseudobulk")
ANALYSIS_DIR <- file.path(BASE, "trajectory_v3", "analysis")
dir.create(OUT_BASE, showWarnings = FALSE, recursive = TRUE)

TARGET_GENES <- c(
  "S100A8", "S100A9", "S100A12", "VCAN", "FCN1", "CD14", "FCGR3A",
  "IL1B", "CXCL8", "CCL2", "CCL3", "TNF", "IL6",
  "ISG15", "IFI6", "IFIT1", "MX1", "OAS1",
  "GZMB", "GZMK", "PRF1", "NKG7", "GNLY", "IFNG",
  "FOXP3", "IL2RA", "CTLA4", "TIGIT", "LAG3", "PDCD1",
  "DDIT4", "TXNIP", "HLA-B", "HLA-C", "CCL4", "HIF1A"
)

N_BINS <- 10
MIN_CELLS_PER_BIN <- 5   # min cells for a patient×bin to be included

compartments <- c("mono", "cd4", "cd8")

# ==============================================================================
# Create pseudobulk profiles
# ==============================================================================
create_pseudobulk <- function(cds, condition_col, subset_expr = NULL,
                               target_genes, n_bins = 10,
                               min_cells = 5) {
  cell_meta <- as.data.frame(colData(cds))

  if (!is.null(subset_expr)) {
    keep_mask <- eval(subset_expr, envir = list(cell_meta = cell_meta))
    keep_cells <- rownames(cell_meta)[keep_mask]
    cds <- cds[, keep_cells]
    cell_meta <- cell_meta[keep_cells, ]
  }

  pt <- pseudotime(cds)
  pt <- pt[!is.na(pt) & is.finite(pt)]
  cell_meta <- cell_meta[names(pt), ]

  # Bin pseudotime into quantile-based bins
  pt_breaks <- quantile(pt, probs = seq(0, 1, length.out = n_bins + 1))
  pt_breaks[1] <- pt_breaks[1] - 0.001  # ensure all included
  pt_bin <- cut(pt, breaks = pt_breaks, labels = FALSE)
  names(pt_bin) <- names(pt)

  # Bin midpoints (for downstream fitting)
  bin_midpoints <- (pt_breaks[-1] + pt_breaks[-(n_bins + 1)]) / 2

  # Get available genes
  available_genes <- intersect(target_genes, rownames(cds))

  # Log-normalize expression
  raw_counts <- counts(cds)[available_genes, names(pt), drop = FALSE]
  lib_sizes <- colSums(counts(cds[, names(pt)]))
  log_norm <- log1p(sweep(as.matrix(raw_counts), 2, lib_sizes, "/") * 10000)

  # Aggregate per patient × bin
  patient_id <- cell_meta$patient_name
  condition <- cell_meta[[condition_col]]
  names(condition) <- rownames(cell_meta)
  names(patient_id) <- rownames(cell_meta)

  # Create grouping key
  group_key <- paste(patient_id, pt_bin, sep = "_bin")

  pb_list <- list()
  for (gene in available_genes) {
    expr_vec <- log_norm[gene, ]

    # Aggregate
    agg <- data.frame(
      cell = names(pt),
      patient = patient_id,
      bin = pt_bin,
      cond = condition,
      expr = expr_vec,
      stringsAsFactors = FALSE
    )

    # Summarize per patient × bin
    pb <- agg %>%
      group_by(patient, bin, cond) %>%
      summarise(
        mean_expr = mean(expr),
        n_cells = n(),
        pct_expr = mean(expr > 0),
        .groups = "drop"
      ) %>%
      filter(n_cells >= min_cells) %>%
      mutate(
        pt_mid = bin_midpoints[bin],
        gene = gene
      )

    pb_list[[gene]] <- pb
  }

  pb_all <- do.call(rbind, pb_list)
  cat("    Pseudobulk rows:", nrow(pb_all), "\n")
  cat("    Patients:", length(unique(pb_all$patient)), "\n")
  cat("    Genes:", length(unique(pb_all$gene)), "\n")

  pb_all
}

# ==============================================================================
# Fit pseudobulk GAMM per gene
# ==============================================================================
fit_pb_gamm <- function(pb_gene, n_pred = 200) {
  gene_name <- unique(pb_gene$gene)
  conds <- sort(unique(pb_gene$cond))

  # Check: enough data?
  n_patients <- length(unique(pb_gene$patient))
  if (n_patients < 6) {
    return(list(gene = gene_name, status = "too_few_patients",
                p_value = NA, abc = NA))
  }

  # Check sparsity: if > 90% of pseudobulk values are zero, skip
  if (mean(pb_gene$mean_expr == 0) > 0.9) {
    return(list(gene = gene_name, status = "too_sparse",
                p_value = NA, abc = NA))
  }

  pb_gene$cond <- factor(pb_gene$cond)
  pb_gene$patient <- factor(pb_gene$patient)

  # Fit GAMM: mean_expr ~ s(pt_mid) + cond + s(pt_mid, by=cond) + s(patient, bs="re")
  # Using GAM (not mixed-effects) with patient RE via smooth
  full_formula <- mean_expr ~ s(pt_mid, k = 5) + cond +
    s(pt_mid, by = cond, k = 5) + s(patient, bs = "re")
  null_formula <- mean_expr ~ s(pt_mid, k = 5) + cond +
    s(patient, bs = "re")

  fit_full <- tryCatch(
    gam(full_formula, data = pb_gene, method = "REML"),
    error = function(e) NULL
  )
  fit_null <- tryCatch(
    gam(null_formula, data = pb_gene, method = "REML"),
    error = function(e) NULL
  )

  if (is.null(fit_full) || is.null(fit_null)) {
    return(list(gene = gene_name, status = "model_failed",
                p_value = NA, abc = NA))
  }

  # Interaction p-value from anova
  anova_res <- tryCatch({
    a <- anova(fit_null, fit_full, test = "Chisq")
    a[["Pr(>Chi)"]][2]
  }, error = function(e) NA)

  # Deviance explained
  dev_expl <- summary(fit_full)$dev.expl

  # Generate predictions for ABC
  pt_range <- range(pb_gene$pt_mid)
  pt_seq <- seq(pt_range[1], pt_range[2], length.out = n_pred)

  # Predict for each condition (exclude RE)
  pred_list <- list()
  for (cv in conds) {
    nd <- data.frame(
      pt_mid = pt_seq,
      cond = factor(cv, levels = levels(pb_gene$cond)),
      patient = pb_gene$patient[1]  # dummy, excluded
    )
    pred <- predict(fit_full, newdata = nd, type = "response",
                    exclude = "s(patient)")
    pred_list[[as.character(cv)]] <- pred
  }

  pred_df <- data.frame(
    pt = pt_seq,
    cond1 = pred_list[[as.character(conds[1])]],
    cond2 = pred_list[[as.character(conds[2])]]
  )

  # Compute ABC (trapezoidal)
  diff_abs <- abs(pred_df$cond1 - pred_df$cond2)
  diff_sq <- (pred_df$cond1 - pred_df$cond2)^2
  signed_diff <- pred_df$cond1 - pred_df$cond2

  dt <- diff(pt_seq)
  abc <- sum((diff_abs[-1] + diff_abs[-n_pred]) / 2 * dt)
  rmise <- sqrt(sum((diff_sq[-1] + diff_sq[-n_pred]) / 2 * dt))
  mean_expr <- mean(pb_gene$mean_expr)
  rabc <- if (mean_expr > 0) abc / mean_expr else NA
  signed_area <- sum((signed_diff[-1] + signed_diff[-n_pred]) / 2 * dt)
  mpd <- max(diff_abs)
  mpd_loc <- pt_seq[which.max(diff_abs)]

  list(
    gene = gene_name,
    status = "success",
    p_value = anova_res,
    dev_explained = dev_expl,
    abc = abc,
    rabc = rabc,
    rmise = rmise,
    mpd = mpd,
    mpd_loc = mpd_loc,
    mean_expr = mean_expr,
    signed_area = signed_area,
    n_patients = n_patients,
    n_obs = nrow(pb_gene),
    pred_df = pred_df
  )
}

# ==============================================================================
# Main loop
# ==============================================================================
all_pb_results <- list()

for (comp_id in compartments) {
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("  Compartment:", comp_id, "\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")

  cds_file <- file.path(TRAJ_V2, comp_id, "monocle3_cds.rds")
  if (!file.exists(cds_file)) {
    cat("  CDS not found. Skipping.\n")
    next
  }

  cat("  Loading CDS...\n")
  cds <- readRDS(cds_file)
  cell_meta <- as.data.frame(colData(cds))

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
    cat("\n  --- Pseudobulk GAMM:", analysis$label, "---\n")

    out_dir <- file.path(OUT_BASE, comp_id, analysis_name)
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

    # Check for existing results
    result_file <- file.path(out_dir, "pb_gamm_results.csv")
    if (file.exists(result_file)) {
      cat("    Already completed, skipping.\n")
      next
    }

    # Create pseudobulk
    cat("  Creating pseudobulk profiles...\n")
    pb <- tryCatch(
      create_pseudobulk(
        cds = cds,
        condition_col = analysis$condition_col,
        subset_expr = analysis$subset_expr,
        target_genes = TARGET_GENES,
        n_bins = N_BINS,
        min_cells = MIN_CELLS_PER_BIN
      ),
      error = function(e) {
        cat("    ERROR:", e$message, "\n")
        NULL
      }
    )

    if (is.null(pb)) next

    # Fit per gene
    genes <- unique(pb$gene)
    cat("  Fitting pseudobulk GAMM for", length(genes), "genes...\n")

    gene_results <- list()
    for (gene in genes) {
      pb_gene <- pb[pb$gene == gene, ]
      cat("    ", gene, "...")

      res <- tryCatch(
        fit_pb_gamm(pb_gene),
        error = function(e) {
          cat(" ERROR:", e$message, "\n")
          list(gene = gene, status = paste0("error: ", e$message),
               p_value = NA, abc = NA)
        }
      )

      if (res$status == "success") {
        cat(" OK (p=", format.pval(res$p_value, digits = 2),
            ", rABC=", round(res$rabc, 2), ")\n")
      } else {
        cat(" ", res$status, "\n")
      }

      gene_results[[gene]] <- res
    }

    # Compile results
    summary_rows <- lapply(gene_results, function(r) {
      data.frame(
        compartment = comp_id,
        analysis = analysis_name,
        gene = r$gene,
        status = r$status,
        p_value = r$p_value %||% NA,
        abc = r$abc %||% NA,
        rabc = r$rabc %||% NA,
        rmise = r$rmise %||% NA,
        mpd = r$mpd %||% NA,
        mpd_loc = r$mpd_loc %||% NA,
        mean_expr = r$mean_expr %||% NA,
        signed_area = r$signed_area %||% NA,
        n_patients = r$n_patients %||% NA,
        n_obs = r$n_obs %||% NA,
        dev_explained = r$dev_explained %||% NA,
        stringsAsFactors = FALSE
      )
    })
    summary_df <- do.call(rbind, summary_rows)

    # FDR correction
    sig_mask <- summary_df$status == "success"
    summary_df$padj <- NA
    if (sum(sig_mask) > 0) {
      summary_df$padj[sig_mask] <- p.adjust(
        summary_df$p_value[sig_mask], method = "BH"
      )
    }

    write.csv(summary_df, result_file, row.names = FALSE)
    cat("  Results saved:", result_file, "\n")

    n_sig <- sum(summary_df$padj < 0.05, na.rm = TRUE)
    cat("  Significant (padj<0.05):", n_sig, "/",
        sum(sig_mask), "\n")

    key <- paste0(comp_id, "_", analysis_name)
    all_pb_results[[key]] <- summary_df
  }

  rm(cds); gc(verbose = FALSE)
}

# ==============================================================================
# Concordance: Pseudobulk vs Cell-level GAMM v3
# ==============================================================================
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("  Concordance: Pseudobulk GAMM vs Cell-level GAMM v3\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

gamm_file <- file.path(ANALYSIS_DIR, "effect_sizes_all.csv")
gamm_df <- if (file.exists(gamm_file)) read.csv(gamm_file) else NULL

concordance_pb <- list()

for (key in names(all_pb_results)) {
  pb_df <- all_pb_results[[key]]
  parts <- strsplit(key, "_")[[1]]
  comp_id <- parts[1]
  analysis_name <- parts[2]

  cat("\n  ---", key, "---\n")

  if (is.null(gamm_df)) {
    cat("    No GAMM results to compare.\n")
    next
  }

  gamm_sub <- gamm_df %>%
    filter(compartment == comp_id, analysis == analysis_name,
           status == "success")

  pb_sub <- pb_df %>%
    filter(status == "success")

  common_genes <- intersect(pb_sub$gene, gamm_sub$gene)
  if (length(common_genes) < 3) {
    cat("    Too few common genes:", length(common_genes), "\n")
    next
  }

  merged <- merge(
    pb_sub[pb_sub$gene %in% common_genes,
           c("gene", "padj", "rabc", "abc", "signed_area")],
    gamm_sub[gamm_sub$gene %in% common_genes,
             c("gene", "interaction_padj", "rABC", "ABC",
               "signed_area")],
    by = "gene",
    suffixes = c(".pb", ".cell")
  )

  # Concordance table
  pb_sig <- merged$padj < 0.05
  cell_sig <- merged$interaction_padj < 0.05

  both <- sum(pb_sig & cell_sig, na.rm = TRUE)
  pb_only <- sum(pb_sig & !cell_sig, na.rm = TRUE)
  cell_only <- sum(!pb_sig & cell_sig, na.rm = TRUE)
  neither <- sum(!pb_sig & !cell_sig, na.rm = TRUE)
  jaccard <- both / (both + pb_only + cell_only)
  if (is.nan(jaccard)) jaccard <- NA

  cat("    Common genes:", length(common_genes), "\n")
  cat("    Both sig:", both, " PB only:", pb_only,
      " Cell only:", cell_only, " Neither:", neither, "\n")
  cat("    Jaccard:", round(jaccard, 3), "\n")

  # rABC correlation
  cor_rabc <- cor.test(merged$rabc, merged$rABC,
                        method = "spearman", exact = FALSE)
  cat("    rABC correlation (Spearman):", round(cor_rabc$estimate, 3),
      " p=", format.pval(cor_rabc$p.value, digits = 2), "\n")

  concordance_pb[[key]] <- data.frame(
    compartment = comp_id,
    analysis = analysis_name,
    n_genes = length(common_genes),
    both_sig = both,
    pb_only = pb_only,
    cell_only = cell_only,
    neither = neither,
    jaccard = round(jaccard, 3),
    rabc_rho = round(cor_rabc$estimate, 3),
    rabc_p = cor_rabc$p.value,
    stringsAsFactors = FALSE
  )

  # Save merged table
  out_dir <- file.path(OUT_BASE, comp_id, analysis_name)
  write.csv(merged, file.path(out_dir, "concordance_merged.csv"),
            row.names = FALSE)

  # Scatter: rABC (pseudobulk vs cell-level)
  p <- ggplot(merged, aes(x = rABC, y = rabc)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_text(aes(label = gene), size = 2.5, vjust = -0.8,
              check_overlap = TRUE) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "grey50") +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("rho = %.2f\np = %.2g",
                             cor_rabc$estimate, cor_rabc$p.value),
             hjust = 1.1, vjust = 1.3, size = 4) +
    labs(
      title = paste0("rABC: Pseudobulk vs Cell-level — ",
                     comp_id, "/", analysis_name),
      x = "Cell-level GAMM v3 rABC",
      y = "Pseudobulk GAMM rABC"
    ) +
    theme_minimal(base_size = 12)

  ggsave(file.path(out_dir, "rabc_scatter.png"), p,
         width = 7, height = 6, dpi = 300)
}

# Save concordance
if (length(concordance_pb) > 0) {
  conc_df <- do.call(rbind, concordance_pb)
  write.csv(conc_df, file.path(OUT_BASE, "concordance_summary.csv"),
            row.names = FALSE)
  cat("\n\nPseudobulk concordance summary:\n")
  print(conc_df)
}

cat("\n\n=== Pseudobulk GAMM COMPLETE ===\n")
cat("Finished:", format(Sys.time()), "\n")
cat("Output:", OUT_BASE, "\n")

#!/usr/bin/env Rscript
# ==============================================================================
# Trajectory Effect Size Computation
# Re-fits GAMM v3 models and computes ABC/rABC/RMISE/MPD effect sizes
# Also generates cross-layer and cross-compartment comparison figures
#
# Phase 1: Per-gene effect sizes from GAMM predictions
# Phase 2: Cross-layer scatter, compartment summary, gene ranking heatmap
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
  library(ggplot2)
  library(tibble)
  library(mgcv)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

source("/data/user3/git_repo/mylit/myR/R/analysis/pseudotime.R")

cat("=== Trajectory Effect Size Computation ===\n")
cat("Started:", format(Sys.time()), "\n\n")

# ==============================================================================
# Configuration
# ==============================================================================
BASE <- "/data/user3/sobj/stroke_hc_v8_2"
TRAJ_V2 <- file.path(BASE, "trajectory_v2")
OUT_BASE <- file.path(BASE, "trajectory_v3")
ANALYSIS_DIR <- file.path(OUT_BASE, "analysis")
dir.create(ANALYSIS_DIR, showWarnings = FALSE, recursive = TRUE)

N_PRED_POINTS <- 200  # denser than original 100 for better integration

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
# Phase 1: Fit GAMM + compute effect sizes
# ==============================================================================

#' Compute effect size metrics from two prediction vectors on a common grid
#'
#' @param pred_df data.frame with columns: pseudotime, fit_response, cond
#' @return named list of effect size metrics
compute_effect_sizes <- function(pred_df) {
  conds <- unique(pred_df$cond)
  if (length(conds) != 2) {
    return(list(ABC = NA, rABC = NA, RMISE = NA, MPD = NA, MPD_location = NA,
                mean_expr = NA, direction_prop = NA, signed_area = NA))
  }

  p1 <- pred_df %>% filter(cond == conds[1]) %>% arrange(pseudotime)
  p2 <- pred_df %>% filter(cond == conds[2]) %>% arrange(pseudotime)

  # Common pseudotime grid (intersection of ranges)
  pt_min <- max(min(p1$pseudotime), min(p2$pseudotime))
  pt_max <- min(max(p1$pseudotime), max(p2$pseudotime))
  if (pt_max <= pt_min) {
    return(list(ABC = NA, rABC = NA, RMISE = NA, MPD = NA, MPD_location = NA,
                mean_expr = NA, direction_prop = NA, signed_area = NA))
  }

  pt_grid <- seq(pt_min, pt_max, length.out = N_PRED_POINTS)

  # Interpolate both curves to common grid
  f1 <- approx(p1$pseudotime, p1$fit_response, xout = pt_grid)$y
  f2 <- approx(p2$pseudotime, p2$fit_response, xout = pt_grid)$y

  if (any(is.na(f1)) || any(is.na(f2))) {
    return(list(ABC = NA, rABC = NA, RMISE = NA, MPD = NA, MPD_location = NA,
                mean_expr = NA, direction_prop = NA, signed_area = NA))
  }

  diff_vec <- f1 - f2
  abs_diff <- abs(diff_vec)
  dt <- diff(pt_grid)[1]  # uniform spacing

  # ABC (L1 norm): trapezoidal rule
  ABC <- sum((abs_diff[-1] + abs_diff[-N_PRED_POINTS]) / 2 * dt)

  # Signed area (for direction)
  signed_area <- sum((diff_vec[-1] + diff_vec[-N_PRED_POINTS]) / 2 * dt)

  # RMISE (L2 norm)
  sq_diff <- diff_vec^2
  RMISE <- sqrt(sum((sq_diff[-1] + sq_diff[-N_PRED_POINTS]) / 2 * dt))

  # MPD (max pointwise difference)
  MPD <- max(abs_diff)
  MPD_location <- pt_grid[which.max(abs_diff)]

  # Mean expression (for normalization)
  mean_expr <- mean(c(f1, f2))

  # Relative ABC
  rABC <- if (mean_expr > 1e-6) ABC / mean_expr else NA

  # Direction: proportion of pseudotime where cond1 > cond2
  direction_prop <- mean(diff_vec > 0)

  list(
    ABC = ABC,
    rABC = rABC,
    RMISE = RMISE,
    MPD = MPD,
    MPD_location = MPD_location,
    mean_expr = mean_expr,
    direction_prop = direction_prop,
    signed_area = signed_area,
    cond1 = as.character(conds[1]),
    cond2 = as.character(conds[2]),
    pt_range = pt_max - pt_min
  )
}

#' Fit GAMM and extract predictions + effect sizes for one gene
#' (Slimmed version of analyze_gene_dynamics_v2 — predictions only, no plots)
#' Includes timeout for sparse genes and sparsity filter
fit_and_extract <- function(gene_id, cds_obj, condition_col_name,
                            batch_col = "GEM", covariate_cols = c("percent.mt"),
                            offset_col = "nCount_RNA", k_val = 6,
                            timeout_secs = 300) {

  if (!gene_id %in% rownames(counts(cds_obj))) {
    return(list(gene = gene_id, status = "gene_not_found", predictions = NULL, effect_sizes = NULL))
  }

  cell_metadata <- as.data.frame(colData(cds_obj))

  # Extract pseudotime
  pt_values <- tryCatch(monocle3::pseudotime(cds_obj), error = function(e) NULL)
  if (is.null(pt_values)) {
    return(list(gene = gene_id, status = "pseudotime_failed", predictions = NULL, effect_sizes = NULL))
  }

  # Build data frame
  gene_expr <- as.numeric(counts(cds_obj)[gene_id, ])
  names(gene_expr) <- colnames(cds_obj)
  common_cells <- intersect(names(pt_values), names(gene_expr))

  dat <- data.frame(
    expr = gene_expr[common_cells],
    pseudotime = pt_values[common_cells],
    cond = cell_metadata[common_cells, condition_col_name],
    stringsAsFactors = FALSE
  )
  if (!is.null(offset_col)) {
    dat$lib_size <- cell_metadata[common_cells, offset_col]
    dat$log_lib_size <- log(dat$lib_size + 1)
  }
  if (!is.null(batch_col)) {
    dat$batch <- factor(cell_metadata[common_cells, batch_col])
  }
  for (cv in covariate_cols) {
    dat[[cv]] <- cell_metadata[common_cells, cv]
  }

  dat <- dat[is.finite(dat$expr) & is.finite(dat$pseudotime), ]
  if (!is.null(offset_col)) dat <- dat[is.finite(dat$log_lib_size) & dat$lib_size > 0, ]
  dat$cond <- factor(dat$cond)

  if (nrow(dat) < 30 || nlevels(dat$cond) < 2) {
    return(list(gene = gene_id, status = "insufficient_data", predictions = NULL, effect_sizes = NULL))
  }

  # Sparsity filter: skip genes with <1% expressing cells in any condition
  pct_expr <- tapply(dat$expr > 0, dat$cond, mean)
  if (any(pct_expr < 0.01)) {
    return(list(gene = gene_id, status = "too_sparse", predictions = NULL, effect_sizes = NULL))
  }

  # Build formula
  formula_str <- "expr ~ s(pseudotime, k = k_val, bs = \"cr\")"
  formula_str <- paste0(formula_str, " + cond + s(pseudotime, by = cond, k = k_val, bs = \"cr\")")
  if (!is.null(offset_col)) formula_str <- paste0(formula_str, " + offset(log_lib_size)")
  for (cv in covariate_cols) {
    if (cv %in% colnames(dat)) formula_str <- paste0(formula_str, " + ", cv)
  }
  if (!is.null(batch_col) && nlevels(dat$batch) >= 2) {
    formula_str <- paste0(formula_str, " + s(batch, bs = \"re\")")
  }

  # Fit GAMM with timeout
  fit <- tryCatch({
    setTimeLimit(elapsed = timeout_secs, transient = TRUE)
    on.exit(setTimeLimit(elapsed = Inf, transient = FALSE), add = TRUE)
    gam(as.formula(formula_str), family = nb(link = "log"), data = dat, method = "REML")
  }, error = function(e) {
    if (grepl("time limit|timeout|elapsed", tolower(e$message))) {
      message("  TIMEOUT: ", gene_id, " exceeded ", timeout_secs, "s")
    }
    NULL
  })
  setTimeLimit(elapsed = Inf, transient = FALSE)  # ensure reset
  if (is.null(fit)) {
    return(list(gene = gene_id, status = "fit_failed_or_timeout", predictions = NULL, effect_sizes = NULL))
  }

  # Interaction test
  interaction_p <- NA_real_
  null_parts <- c("expr ~ s(pseudotime, k = k_val, bs = \"cr\")", "cond")
  if (!is.null(offset_col)) null_parts <- c(null_parts, "offset(log_lib_size)")
  for (cv in covariate_cols) {
    if (cv %in% colnames(dat)) null_parts <- c(null_parts, cv)
  }
  if (!is.null(batch_col) && nlevels(dat$batch) >= 2) {
    null_parts <- c(null_parts, "s(batch, bs = \"re\")")
  }
  null_str <- paste(null_parts[1], paste(null_parts[-1], collapse = " + "), sep = " + ")

  fit0 <- tryCatch({
    setTimeLimit(elapsed = timeout_secs, transient = TRUE)
    on.exit(setTimeLimit(elapsed = Inf, transient = FALSE), add = TRUE)
    gam(as.formula(null_str), family = nb(link = "log"), data = dat, method = "REML")
  }, error = function(e) NULL)
  setTimeLimit(elapsed = Inf, transient = FALSE)
  if (!is.null(fit0)) {
    anova_res <- tryCatch(anova(fit0, fit, test = "Chisq"), error = function(e) NULL)
    p_col <- if (!is.null(anova_res) && "Pr(>Chi)" %in% names(anova_res)) "Pr(>Chi)" else "P(>|Chi|)"
    if (!is.null(anova_res) && !is.null(anova_res[[p_col]]) && length(anova_res[[p_col]]) >= 2) {
      interaction_p <- anova_res[[p_col]][2]
    }
  }

  # Predictions (N_PRED_POINTS per condition)
  active_conditions <- levels(dat$cond)
  newd_list <- lapply(active_conditions, function(lvl) {
    pt_sub <- dat$pseudotime[dat$cond == lvl]
    if (length(pt_sub) < 2) return(NULL)
    nd <- data.frame(
      pseudotime = seq(min(pt_sub, na.rm = TRUE), max(pt_sub, na.rm = TRUE),
                       length.out = N_PRED_POINTS),
      cond = factor(lvl, levels = active_conditions)
    )
    if (!is.null(offset_col)) nd$log_lib_size <- median(dat$log_lib_size, na.rm = TRUE)
    for (cv in covariate_cols) {
      if (cv %in% colnames(dat)) nd[[cv]] <- mean(dat[[cv]], na.rm = TRUE)
    }
    if (!is.null(batch_col) && "batch" %in% colnames(dat)) {
      nd$batch <- levels(dat$batch)[1]
    }
    nd
  })
  newd <- do.call(rbind, newd_list)
  if (is.null(newd) || nrow(newd) == 0) {
    return(list(gene = gene_id, status = "prediction_failed", predictions = NULL, effect_sizes = NULL))
  }

  exclude_terms <- if (!is.null(batch_col) && nlevels(dat$batch) >= 2) "s(batch)" else NULL
  newd$fit_response <- predict(fit, newdata = newd, type = "response", exclude = exclude_terms)

  # Compute effect sizes
  es <- compute_effect_sizes(newd)

  # Cell density along pseudotime (for density-weighted ABC)
  density_list <- lapply(active_conditions, function(lvl) {
    pt_sub <- dat$pseudotime[dat$cond == lvl]
    data.frame(cond = lvl, n_cells = length(pt_sub))
  })

  gam_summary <- summary(fit)

  list(
    gene = gene_id,
    status = "success",
    predictions = newd[, c("pseudotime", "cond", "fit_response")],
    effect_sizes = c(es, list(
      interaction_p = interaction_p,
      dev_explained = gam_summary$dev.expl,
      n_cells = nrow(dat)
    ))
  )
}


# ==============================================================================
# Main loop: fit all genes × compartments × conditions
# ==============================================================================

# Check for partial results from previous run
SKIP_COMPLETED <- FALSE  # set to TRUE to skip completed analyses

all_results <- list()

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

  # Inject missing columns if needed
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
      cell_meta <- as.data.frame(colData(cds))
    }
  }

  has_g3 <- "g3" %in% colnames(cell_meta) && sum(!is.na(cell_meta$g3)) > 100
  has_cohort <- "cohort" %in% colnames(cell_meta) && length(unique(cell_meta$cohort)) >= 2

  analyses <- list()
  if (has_g3) {
    analyses[["g3"]] <- list(condition_col = "g3", label = "g3")
  }
  if (has_cohort) {
    analyses[["cohort"]] <- list(condition_col = "cohort", label = "cohort")
  }

  for (analysis_name in names(analyses)) {
    analysis <- analyses[[analysis_name]]
    cat("\n  --- ", comp_id, "/", analysis_name, " ---\n")

    # Skip if already completed (check for existing per-analysis CSV)
    existing_csv <- file.path(ANALYSIS_DIR, paste0("effect_sizes_", comp_id, "_", analysis_name, ".csv"))
    if (SKIP_COMPLETED && file.exists(existing_csv)) {
      cat("    SKIP: results already exist at", existing_csv, "\n")
      existing_df <- read.csv(existing_csv, stringsAsFactors = FALSE)
      for (i in seq_len(nrow(existing_df))) {
        key <- paste(existing_df$compartment[i], existing_df$analysis[i], existing_df$gene[i], sep = "/")
        all_results[[key]] <- list(
          compartment = existing_df$compartment[i], analysis = existing_df$analysis[i],
          gene = existing_df$gene[i], status = existing_df$status[i],
          predictions = NULL,
          effect_sizes = as.list(existing_df[i, setdiff(names(existing_df), c("compartment","analysis","gene","status"))])
        )
      }
      next
    }

    # Subset CDS for g3 (filter NAs)
    if (analysis_name == "g3") {
      keep <- !is.na(cell_meta[[analysis$condition_col]])
      cds_sub <- cds[, which(keep)]
    } else {
      cds_sub <- cds
    }
    cat("    Cells:", ncol(cds_sub), "\n")

    available_genes <- intersect(TARGET_GENES, rownames(cds_sub))
    cat("    Genes:", length(available_genes), "\n")

    for (gene in available_genes) {
      cat("    ", gene, "...")
      res <- tryCatch(
        fit_and_extract(gene, cds_sub, analysis$condition_col),
        error = function(e) {
          cat(" ERROR:", e$message, "\n")
          list(gene = gene, status = paste0("error: ", e$message),
               predictions = NULL, effect_sizes = NULL)
        }
      )

      key <- paste(comp_id, analysis_name, gene, sep = "/")
      all_results[[key]] <- list(
        compartment = comp_id,
        analysis = analysis_name,
        gene = gene,
        status = res$status,
        predictions = res$predictions,
        effect_sizes = res$effect_sizes
      )

      if (res$status == "success") {
        cat(" ABC=", round(res$effect_sizes$ABC, 2),
            " rABC=", round(res$effect_sizes$rABC, 3),
            " p=", format.pval(res$effect_sizes$interaction_p, digits = 2), "\n")
      } else {
        cat(" ", res$status, "\n")
      }
    }
  }

  rm(cds); gc(verbose = FALSE)
}

# ==============================================================================
# Save predictions for future use
# ==============================================================================
cat("\n\nSaving all prediction curves...\n")
pred_list <- lapply(all_results, function(r) {
  if (!is.null(r$predictions)) {
    r$predictions$compartment <- r$compartment
    r$predictions$analysis <- r$analysis
    r$predictions$gene <- r$gene
    r$predictions
  }
})
all_predictions <- do.call(rbind, pred_list)
saveRDS(all_predictions, file.path(ANALYSIS_DIR, "all_prediction_curves.rds"))
cat("  Saved:", nrow(all_predictions), "prediction rows\n")

# ==============================================================================
# Build effect size summary table
# ==============================================================================
cat("\nBuilding effect size summary...\n")

es_rows <- lapply(all_results, function(r) {
  if (is.null(r$effect_sizes)) {
    return(data.frame(
      compartment = r$compartment, analysis = r$analysis, gene = r$gene,
      status = r$status,
      ABC = NA, rABC = NA, RMISE = NA, MPD = NA, MPD_location = NA,
      mean_expr = NA, direction_prop = NA, signed_area = NA,
      interaction_p = NA, dev_explained = NA, n_cells = NA,
      cond1 = NA, cond2 = NA, pt_range = NA,
      stringsAsFactors = FALSE
    ))
  }
  es <- r$effect_sizes
  data.frame(
    compartment = r$compartment, analysis = r$analysis, gene = r$gene,
    status = r$status,
    ABC = es$ABC, rABC = es$rABC, RMISE = es$RMISE,
    MPD = es$MPD, MPD_location = es$MPD_location,
    mean_expr = es$mean_expr, direction_prop = es$direction_prop,
    signed_area = es$signed_area,
    interaction_p = es$interaction_p, dev_explained = es$dev_explained,
    n_cells = es$n_cells,
    cond1 = es$cond1, cond2 = es$cond2, pt_range = es$pt_range,
    stringsAsFactors = FALSE
  )
})
es_df <- do.call(rbind, es_rows)

# Add adjusted p-values per analysis group
es_df <- es_df %>%
  group_by(compartment, analysis) %>%
  mutate(interaction_padj = p.adjust(interaction_p, method = "BH")) %>%
  ungroup()

write.csv(es_df, file.path(ANALYSIS_DIR, "effect_sizes_all.csv"), row.names = FALSE)
cat("  Effect sizes saved:", nrow(es_df), "rows\n")

# Per-analysis CSV
for (comp in unique(es_df$compartment)) {
  for (ana in unique(es_df$analysis[es_df$compartment == comp])) {
    sub <- es_df %>% filter(compartment == comp, analysis == ana) %>% arrange(desc(rABC))
    fname <- paste0("effect_sizes_", comp, "_", ana, ".csv")
    write.csv(sub, file.path(ANALYSIS_DIR, fname), row.names = FALSE)
  }
}

# ==============================================================================
# Phase 2: Comparison Figures
# ==============================================================================
cat("\n\n=== Phase 2: Generating comparison figures ===\n")

# Color schemes
comp_colors <- c(mono = "#E64B35", cd4 = "#4DBBD5", cd8 = "#00A087")
ana_colors <- c(g3 = "#3C5488", cohort = "#F39B7F")

# ------------------------------------------------------------------------------
# Figure 1: Gene ranking heatmap (rABC, 36 genes × 6 analyses)
# ------------------------------------------------------------------------------
cat("  Figure 1: Gene ranking heatmap...\n")

heatmap_df <- es_df %>%
  filter(status == "success") %>%
  mutate(
    label = paste0(compartment, "/", analysis),
    sig = ifelse(!is.na(interaction_padj) & interaction_padj < 0.05, "*", "")
  )

# Order genes by mean rABC across all analyses
gene_order <- heatmap_df %>%
  group_by(gene) %>%
  summarize(mean_rABC = mean(rABC, na.rm = TRUE)) %>%
  arrange(desc(mean_rABC)) %>%
  pull(gene)

heatmap_df$gene <- factor(heatmap_df$gene, levels = rev(gene_order))
heatmap_df$label <- factor(heatmap_df$label,
                           levels = c("mono/cohort", "mono/g3",
                                      "cd4/cohort", "cd4/g3",
                                      "cd8/cohort", "cd8/g3"))

p_heatmap <- ggplot(heatmap_df, aes(x = label, y = gene, fill = rABC)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sig), size = 4, vjust = 0.8) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, na.value = "grey90",
                       name = "rABC\n(relative area\nbetween curves)") +
  labs(title = "Gene Dynamics Effect Sizes (rABC)",
       subtitle = "* = interaction padj < 0.05",
       x = "Compartment / Comparison", y = "Gene") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())

ggsave(file.path(ANALYSIS_DIR, "heatmap_rABC.png"), p_heatmap,
       width = 8, height = 12, dpi = 300)
cat("    Saved: heatmap_rABC.png\n")

# ------------------------------------------------------------------------------
# Figure 2: Cross-layer scatter (rABC cohort vs rABC g3, per compartment)
# ------------------------------------------------------------------------------
cat("  Figure 2: Cross-layer scatter...\n")

cross_df <- es_df %>%
  filter(status == "success") %>%
  select(compartment, analysis, gene, rABC, interaction_padj) %>%
  pivot_wider(names_from = analysis, values_from = c(rABC, interaction_padj))

scatter_plots <- list()
for (comp in compartments) {
  sub <- cross_df %>% filter(compartment == comp)
  if (nrow(sub) == 0 || all(is.na(sub$rABC_cohort)) || all(is.na(sub$rABC_g3))) next

  # Significance labels
  sub <- sub %>% mutate(
    sig_cat = case_when(
      !is.na(interaction_padj_cohort) & interaction_padj_cohort < 0.05 &
        !is.na(interaction_padj_g3) & interaction_padj_g3 < 0.05 ~ "Both sig",
      !is.na(interaction_padj_cohort) & interaction_padj_cohort < 0.05 ~ "Cohort only",
      !is.na(interaction_padj_g3) & interaction_padj_g3 < 0.05 ~ "g3 only",
      TRUE ~ "Neither"
    )
  )

  # Correlation
  cor_test <- tryCatch(
    cor.test(sub$rABC_cohort, sub$rABC_g3, method = "spearman", use = "complete.obs"),
    error = function(e) NULL
  )
  cor_label <- if (!is.null(cor_test)) {
    paste0("rho=", round(cor_test$estimate, 2), ", p=", format.pval(cor_test$p.value, digits = 2))
  } else "NA"

  p <- ggplot(sub, aes(x = rABC_cohort, y = rABC_g3, color = sig_cat)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_text(aes(label = gene), size = 2.5, vjust = -0.8, show.legend = FALSE) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c(
      "Both sig" = "#E64B35", "Cohort only" = "#F39B7F",
      "g3 only" = "#3C5488", "Neither" = "grey70"
    )) +
    labs(title = paste0(toupper(comp), " — Cross-layer concordance"),
         subtitle = cor_label,
         x = "rABC (HC vs Stroke)", y = "rABC (g3: Bad vs Good)",
         color = "Significance") +
    theme_classic(base_size = 11) +
    coord_equal()

  scatter_plots[[comp]] <- p
}

if (length(scatter_plots) > 0) {
  p_combined <- wrap_plots(scatter_plots, ncol = length(scatter_plots)) +
    plot_annotation(title = "Cross-Layer Concordance of Trajectory Effect Sizes")
  ggsave(file.path(ANALYSIS_DIR, "crosslayer_scatter.png"), p_combined,
         width = 6 * length(scatter_plots), height = 6, dpi = 300)
  cat("    Saved: crosslayer_scatter.png\n")
}

# ------------------------------------------------------------------------------
# Figure 3: Compartment summary barplot (median rABC)
# ------------------------------------------------------------------------------
cat("  Figure 3: Compartment summary...\n")

comp_summary <- es_df %>%
  filter(status == "success") %>%
  group_by(compartment, analysis) %>%
  summarize(
    median_rABC = median(rABC, na.rm = TRUE),
    mean_rABC = mean(rABC, na.rm = TRUE),
    n_sig = sum(interaction_padj < 0.05, na.rm = TRUE),
    n_total = n(),
    median_ABC = median(ABC, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(label = paste0(compartment, "/", analysis))

write.csv(comp_summary, file.path(ANALYSIS_DIR, "compartment_summary.csv"), row.names = FALSE)

p_comp <- ggplot(comp_summary, aes(x = reorder(label, -median_rABC),
                                    y = median_rABC, fill = compartment)) +
  geom_col(width = 0.7, alpha = 0.85) +
  geom_text(aes(label = paste0(n_sig, "/", n_total, " sig")),
            vjust = -0.5, size = 3) +
  scale_fill_manual(values = comp_colors) +
  labs(title = "Compartment-level Trajectory Perturbation",
       subtitle = "Median rABC (relative area between curves) across 36 target genes",
       x = "", y = "Median rABC") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(ANALYSIS_DIR, "compartment_summary.png"), p_comp,
       width = 8, height = 5, dpi = 300)
cat("    Saved: compartment_summary.png\n")

# ------------------------------------------------------------------------------
# Figure 4: Top genes by rABC (per analysis)
# ------------------------------------------------------------------------------
cat("  Figure 4: Top gene barplots...\n")

for (comp in compartments) {
  for (ana in c("cohort", "g3")) {
    sub <- es_df %>%
      filter(compartment == comp, analysis == ana, status == "success") %>%
      arrange(desc(rABC)) %>%
      head(20)

    if (nrow(sub) == 0) next

    sub$gene <- factor(sub$gene, levels = rev(sub$gene))
    sub$sig_label <- ifelse(!is.na(sub$interaction_padj) & sub$interaction_padj < 0.05, "*", "")

    p <- ggplot(sub, aes(x = gene, y = rABC)) +
      geom_col(fill = ifelse(ana == "cohort", "#F39B7F", "#3C5488"), alpha = 0.85) +
      geom_text(aes(label = sig_label), hjust = -0.5, size = 5) +
      coord_flip() +
      labs(title = paste0(toupper(comp), " / ", ana, " — Top genes by rABC"),
           subtitle = "* = interaction padj < 0.05",
           x = "", y = "rABC (relative area between curves)") +
      theme_classic(base_size = 11)

    ggsave(file.path(ANALYSIS_DIR, paste0("top_genes_", comp, "_", ana, ".png")), p,
           width = 7, height = 6, dpi = 300)
  }
}
cat("    Saved: top_genes_*.png\n")

# ------------------------------------------------------------------------------
# Figure 5: Pointwise difference curves for top genes (monocyte)
# ------------------------------------------------------------------------------
cat("  Figure 5: Difference curves (mono)...\n")

pred_all <- readRDS(file.path(ANALYSIS_DIR, "all_prediction_curves.rds"))

# For monocyte/cohort - top 8 by rABC
top_mono <- es_df %>%
  filter(compartment == "mono", analysis == "cohort", status == "success") %>%
  arrange(desc(rABC)) %>%
  head(8)

if (nrow(top_mono) > 0) {
  diff_plots <- list()
  for (i in seq_len(nrow(top_mono))) {
    g <- top_mono$gene[i]
    pred_sub <- pred_all %>% filter(compartment == "mono", analysis == "cohort", gene == g)
    conds <- unique(pred_sub$cond)
    if (length(conds) != 2) next

    p1 <- pred_sub %>% filter(cond == conds[1]) %>% arrange(pseudotime)
    p2 <- pred_sub %>% filter(cond == conds[2]) %>% arrange(pseudotime)
    pt_min <- max(min(p1$pseudotime), min(p2$pseudotime))
    pt_max <- min(max(p1$pseudotime), max(p2$pseudotime))
    pt_grid <- seq(pt_min, pt_max, length.out = 200)
    f1 <- approx(p1$pseudotime, p1$fit_response, xout = pt_grid)$y
    f2 <- approx(p2$pseudotime, p2$fit_response, xout = pt_grid)$y

    diff_df <- data.frame(pseudotime = pt_grid, diff = f1 - f2,
                          f1 = f1, f2 = f2)

    p_diff <- ggplot(diff_df, aes(x = pseudotime)) +
      geom_ribbon(aes(ymin = pmin(0, diff), ymax = pmax(0, diff)),
                  fill = ifelse(diff_df$diff > 0, "#E64B35", "#4DBBD5"),
                  alpha = 0.4) +
      geom_line(aes(y = diff), color = "black", linewidth = 0.8) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      labs(title = g, subtitle = paste0("rABC=", round(top_mono$rABC[i], 3)),
           x = "Pseudotime", y = paste0(conds[1], " - ", conds[2])) +
      theme_classic(base_size = 9)

    diff_plots[[i]] <- p_diff
  }

  if (length(diff_plots) > 0) {
    p_diff_combined <- wrap_plots(diff_plots, ncol = 4) +
      plot_annotation(title = "Monocyte / Cohort: Pointwise Difference Curves (Top 8)",
                      subtitle = "Red = HC higher, Blue = Stroke higher")
    ggsave(file.path(ANALYSIS_DIR, "diff_curves_mono_cohort.png"), p_diff_combined,
           width = 16, height = 8, dpi = 300)
    cat("    Saved: diff_curves_mono_cohort.png\n")
  }
}

# Same for mono/g3
top_mono_g3 <- es_df %>%
  filter(compartment == "mono", analysis == "g3", status == "success") %>%
  arrange(desc(rABC)) %>%
  head(8)

if (nrow(top_mono_g3) > 0) {
  diff_plots_g3 <- list()
  for (i in seq_len(nrow(top_mono_g3))) {
    g <- top_mono_g3$gene[i]
    pred_sub <- pred_all %>% filter(compartment == "mono", analysis == "g3", gene == g)
    conds <- unique(pred_sub$cond)
    if (length(conds) != 2) next

    p1 <- pred_sub %>% filter(cond == conds[1]) %>% arrange(pseudotime)
    p2 <- pred_sub %>% filter(cond == conds[2]) %>% arrange(pseudotime)
    pt_min <- max(min(p1$pseudotime), min(p2$pseudotime))
    pt_max <- min(max(p1$pseudotime), max(p2$pseudotime))
    pt_grid <- seq(pt_min, pt_max, length.out = 200)
    f1 <- approx(p1$pseudotime, p1$fit_response, xout = pt_grid)$y
    f2 <- approx(p2$pseudotime, p2$fit_response, xout = pt_grid)$y

    diff_df <- data.frame(pseudotime = pt_grid, diff = f1 - f2)

    p_diff <- ggplot(diff_df, aes(x = pseudotime)) +
      geom_ribbon(aes(ymin = pmin(0, diff), ymax = pmax(0, diff)),
                  fill = ifelse(diff_df$diff > 0, "#3C5488", "#E64B35"),
                  alpha = 0.4) +
      geom_line(aes(y = diff), color = "black", linewidth = 0.8) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      labs(title = g, subtitle = paste0("rABC=", round(top_mono_g3$rABC[i], 3)),
           x = "Pseudotime", y = paste0(conds[1], " - ", conds[2])) +
      theme_classic(base_size = 9)

    diff_plots_g3[[i]] <- p_diff
  }

  if (length(diff_plots_g3) > 0) {
    p_diff_g3 <- wrap_plots(diff_plots_g3, ncol = 4) +
      plot_annotation(title = "Monocyte / g3: Pointwise Difference Curves (Top 8)",
                      subtitle = "Blue = g3=1 (Good) higher, Red = g3=2 (Bad) higher")
    ggsave(file.path(ANALYSIS_DIR, "diff_curves_mono_g3.png"), p_diff_g3,
           width = 16, height = 8, dpi = 300)
    cat("    Saved: diff_curves_mono_g3.png\n")
  }
}

# ==============================================================================
# Final summary
# ==============================================================================
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("TRAJECTORY EFFECT SIZE COMPUTATION COMPLETE\n")
cat("Finished:", format(Sys.time()), "\n\n")

cat("Output directory:", ANALYSIS_DIR, "\n\n")
cat("Files generated:\n")
cat("  - effect_sizes_all.csv (master table)\n")
cat("  - effect_sizes_*_*.csv (per analysis)\n")
cat("  - all_prediction_curves.rds (prediction curves for all genes)\n")
cat("  - compartment_summary.csv\n")
cat("  - heatmap_rABC.png\n")
cat("  - crosslayer_scatter.png\n")
cat("  - compartment_summary.png\n")
cat("  - top_genes_*.png\n")
cat("  - diff_curves_*.png\n")

cat("\n--- Top 5 genes by rABC per analysis ---\n")
es_df %>%
  filter(status == "success") %>%
  group_by(compartment, analysis) %>%
  slice_max(rABC, n = 5) %>%
  select(compartment, analysis, gene, rABC, ABC, interaction_padj) %>%
  print(n = 60)

cat("\n--- Compartment summary ---\n")
print(comp_summary)

#!/usr/bin/env Rscript
# ==============================================================================
# Gene Dynamics v3: Batch-corrected GAMM
# Re-analyses gene dynamics from existing Monocle3 CDS objects (trajectory_v2)
# using analyze_gene_dynamics_v2 with:
#   - offset(log(nCount_RNA)) for library size
#   - percent.mt as fixed covariate
#   - s(GEM, bs="re") for batch random effect
#   - y-axis clipping at 99th percentile
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
})

source("/data/user3/git_repo/mylit/myR/R/analysis/pseudotime.R")

cat("=== Gene Dynamics v3: Batch-corrected GAMM ===\n")
cat("Started:", format(Sys.time()), "\n\n")

BASE <- "/data/user3/sobj/stroke_hc_v8_2"
TRAJ_V2 <- file.path(BASE, "trajectory_v2")
OUT_BASE <- file.path(BASE, "trajectory_v3")
dir.create(OUT_BASE, showWarnings = FALSE, recursive = TRUE)

TARGET_GENES <- c(
  "S100A8", "S100A9", "S100A12", "VCAN", "FCN1", "CD14", "FCGR3A",
  "IL1B", "CXCL8", "CCL2", "CCL3", "TNF", "IL6",
  "ISG15", "IFI6", "IFIT1", "MX1", "OAS1",
  "GZMB", "GZMK", "PRF1", "NKG7", "GNLY", "IFNG",
  "FOXP3", "IL2RA", "CTLA4", "TIGIT", "LAG3", "PDCD1",
  "DDIT4", "TXNIP", "HLA-B", "HLA-C", "CCL4", "HIF1A"
)

compartments <- c("mono", "cd4", "cd8")

# ---- Process each compartment ----
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
  cat("  Available columns:", paste(head(colnames(cell_meta), 20), collapse = ", "), "\n")

  # Verify required columns
  required <- c("nCount_RNA", "percent.mt", "GEM")
  missing <- setdiff(required, colnames(cell_meta))
  if (length(missing) > 0) {
    cat("  WARNING: Missing columns:", paste(missing, collapse = ", "), "\n")
    cat("  Will attempt to load from Seurat subset and inject...\n")

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
    } else {
      cat("  Seurat subset not found either. Skipping compartment.\n")
      next
    }
    cell_meta <- as.data.frame(colData(cds))
  }

  # Check for g3 and cohort
  has_g3 <- "g3" %in% colnames(cell_meta) && sum(!is.na(cell_meta$g3)) > 100
  has_cohort <- "cohort" %in% colnames(cell_meta) && length(unique(cell_meta$cohort)) >= 2

  cat("  g3 available:", has_g3,
      if (has_g3) paste0("(", sum(!is.na(cell_meta$g3)), " cells)") else "", "\n")
  cat("  Cohort available:", has_cohort, "\n")
  cat("  GEM levels:", nlevels(factor(cell_meta$GEM)), "\n")
  cat("  nCount_RNA range:", round(range(cell_meta$nCount_RNA, na.rm = TRUE)), "\n")
  cat("  percent.mt range:", round(range(cell_meta$percent.mt, na.rm = TRUE), 2), "\n")

  # ---- Run for each condition comparison ----
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
    cat("\n  --- Analysis:", analysis$label, "---\n")

    # Create output directory
    out_dir <- file.path(OUT_BASE, comp_id, paste0("gene_dynamics_", analysis_name))
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

    # Subset CDS if needed (e.g., g3 requires filtering NAs)
    if (!is.null(analysis$subset_expr)) {
      keep_mask <- eval(analysis$subset_expr)
      keep_cells <- rownames(cell_meta)[keep_mask]
      cds_sub <- cds[, keep_cells]
      cat("    Subset:", length(keep_cells), "cells\n")
    } else {
      cds_sub <- cds
      cat("    Using all", ncol(cds), "cells\n")
    }

    # Check GEM levels in subset
    sub_meta <- as.data.frame(colData(cds_sub))
    n_gem <- nlevels(factor(sub_meta$GEM))
    cat("    GEM levels in subset:", n_gem, "\n")

    available_genes <- intersect(TARGET_GENES, rownames(cds_sub))
    cat("    Genes:", length(available_genes), "/", length(TARGET_GENES), "\n\n")

    results_list <- list()
    for (gene in available_genes) {
      cat("    ", gene, "...")
      res <- tryCatch({
        analyze_gene_dynamics_v2(
          gene_id = gene,
          cds_obj = cds_sub,
          condition_col_name = analysis$condition_col,
          batch_col = "GEM",
          covariate_cols = c("percent.mt"),
          offset_col = "nCount_RNA",
          output_dir = out_dir,
          k_val = 6,
          min_cells_for_fit = 30,
          ylim_quantile = 0.99,
          plot_width = 7,
          plot_height = 5,
          plot_dpi = 300,
          scale_DR = TRUE,
          plot_split = FALSE
        )
      }, error = function(e) {
        cat(" ERROR:", e$message, "\n")
        list(gene = gene, status = paste0("error: ", e$message),
             metrics = NULL, model_formula = NULL)
      })
      if (!is.null(res$metrics)) {
        cat(" OK (dev.expl=", round(res$metrics$deviance_explained, 3),
            ", P(int)=", format.pval(res$metrics$interaction_p_value, digits = 2), ")\n")
      } else if (res$status != "success") {
        cat(" ", res$status, "\n")
      }
      results_list[[gene]] <- res
    }

    # ---- Save summary ----
    summary_rows <- lapply(results_list, function(r) {
      m <- r$metrics
      if (is.null(m)) {
        return(data.frame(
          gene = r$gene, status = r$status,
          n_cells = NA, dev_explained = NA,
          interaction_p = NA, formula = r$model_formula %||% NA,
          stringsAsFactors = FALSE
        ))
      }
      # Extract TV and DR per condition
      tv_str <- if (is.data.frame(m$TV_per_condition))
        paste(apply(m$TV_per_condition, 1, function(r) paste0(r["cond"], ":", round(as.numeric(r["TV"]), 2))),
              collapse = ";") else NA
      dr_str <- if (is.data.frame(m$DR_per_condition))
        paste(apply(m$DR_per_condition, 1, function(r) paste0(r["cond"], ":", round(as.numeric(r["DR"]), 2))),
              collapse = ";") else NA

      data.frame(
        gene = r$gene, status = r$status,
        n_cells = m$n_cells_fit,
        dev_explained = round(m$deviance_explained, 4),
        interaction_p = m$interaction_p_value,
        TV = tv_str, DR = dr_str,
        formula = r$model_formula %||% NA,
        stringsAsFactors = FALSE
      )
    })
    summary_df <- do.call(rbind, summary_rows)
    write.csv(summary_df, file.path(out_dir, "gene_dynamics_v3_summary.csv"), row.names = FALSE)
    cat("\n    Summary saved:", nrow(summary_df), "genes\n")
    cat("    Output:", out_dir, "\n")
  }

  rm(cds); gc(verbose = FALSE)
}

# ---- Final summary ----
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("GENE DYNAMICS v3 COMPLETE\n")
cat("Finished:", format(Sys.time()), "\n")
cat("Output:", OUT_BASE, "\n")
for (comp_id in compartments) {
  comp_dir <- file.path(OUT_BASE, comp_id)
  if (dir.exists(comp_dir)) {
    cat("\n", comp_id, ":\n")
    dirs <- list.dirs(comp_dir, recursive = FALSE, full.names = FALSE)
    for (d in dirs) {
      n_plots <- length(list.files(file.path(comp_dir, d), pattern = "\\.png$"))
      cat("  ", d, ":", n_plots, "plots\n")
    }
  }
}
cat("\nModel: NB-GAMM with offset(log(nCount_RNA)) + percent.mt + s(GEM, bs='re')\n")
cat("See: docs/gene_dynamics_approaches.md for methodology details\n")

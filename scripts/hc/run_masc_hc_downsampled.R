#!/usr/bin/env Rscript
# HC-only MASC with downsampling: anno1 × sex, anno1 × age_group
# MASC glmer fails with 549K cells. Downsample to ~100K for convergence.
#
# Usage: Rscript scripts/hc/run_masc_hc_downsampled.R
# Output: /data/user3/sobj/hc_only_v1/masc/

.libPaths(c(
  "/home/user3/GJC_KDW_250721/renv/library/R-4.3/x86_64-pc-linux-gnu",
  "/home/user3/R/x86_64-pc-linux-gnu-library/4.3",
  .libPaths()
))

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

cat("=== HC-only MASC (Downsampled) ===\n")
cat(sprintf("Time: %s\n\n", Sys.time()))

MASC_R <- "/home/user3/data_user3/git_repo/mylit/Git_Repo/_wt/masc/myR/R/masc.R"
source(MASC_R)

sobj_path <- "/data/user3/sobj/hc_only_v1/2_hc_annotated.qs"
out_base <- "/data/user3/sobj/hc_only_v1/masc"
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

cat("Loading Seurat object...\n")
sobj <- qread(sobj_path)
cat(sprintf("  Full: %d cells, %d patients\n", ncol(sobj), length(unique(sobj$name))))

# Ensure age_group exists
if (!"age_group" %in% colnames(sobj@meta.data)) {
  sobj$age_group <- ifelse(sobj$age < 35, "Young",
                           ifelse(sobj$age <= 50, "Middle", "Old"))
}

# ---- Downsample: max 1200 cells per patient ----
set.seed(42)
meta <- sobj@meta.data
meta$barcode <- colnames(sobj)

# slice_sample(n=) can't use n() in argument; use split-apply approach
sampled <- do.call(rbind, lapply(split(meta, meta$name), function(df) {
  if (nrow(df) <= 1200) return(df)
  df[sample.int(nrow(df), 1200), , drop = FALSE]
}))

cat(sprintf("  Downsampled: %d cells (from %d)\n", nrow(sampled), ncol(sobj)))
sobj_ds <- sobj[, sampled$barcode]

rm(sobj); gc(verbose = FALSE)

# ---- Run MASC ----
comparisons <- list(
  list(name = "anno1_sex_ds", cluster_var = "anno1",
       contrast_var = "sex", random_effects = "name", fixed_effects = "age"),
  list(name = "anno1_age_group_ds", cluster_var = "anno1",
       contrast_var = "age_group", random_effects = "name", fixed_effects = "sex"),
  list(name = "anno2_sex_ds", cluster_var = "anno2",
       contrast_var = "sex", random_effects = "name", fixed_effects = "age"),
  list(name = "anno2_age_group_ds", cluster_var = "anno2",
       contrast_var = "age_group", random_effects = "name", fixed_effects = "sex")
)

for (comp in comparisons) {
  cat(sprintf("\n--- Running MASC: %s ---\n", comp$name))
  comp_dir <- file.path(out_base, comp$name)
  dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)

  tryCatch({
    res <- run_masc_pipeline(
      seurat_obj = sobj_ds,
      cluster_var = comp$cluster_var,
      contrast_var = comp$contrast_var,
      random_effects = comp$random_effects,
      fixed_effects = comp$fixed_effects,
      save = TRUE,
      output_dir = comp_dir,
      prefix = comp$name,
      force_run = TRUE,
      plotting = TRUE,
      verbose = TRUE
    )

    result <- res$masc_results
    result_sorted <- result[order(result$model.pvalue), ]
    write.csv(result_sorted, file.path(comp_dir, sprintf("%s_results.csv", comp$name)),
              row.names = FALSE)

    n_sig <- sum(result$model.pvalue < 0.05, na.rm = TRUE)
    n_fdr <- sum(result$model.pvalue.fdr < 0.05, na.rm = TRUE)
    cat(sprintf("  %d cell types, %d nom.sig, %d FDR sig\n", nrow(result), n_sig, n_fdr))

    # Print top hits
    top <- head(result_sorted[!is.na(result_sorted$model.pvalue) & result_sorted$model.pvalue < 0.1, ], 10)
    if (nrow(top) > 0) {
      cat("  Top hits:\n")
      or_col <- grep("OR$", colnames(top), value = TRUE)[1]
      for (i in seq_len(nrow(top))) {
        or_val <- if (!is.na(or_col)) top[[or_col]][i] else NA
        cat(sprintf("    %-30s OR=%.3f p=%.2e fdr=%.2e\n",
                    top$cluster[i], or_val, top$model.pvalue[i], top$model.pvalue.fdr[i]))
      }
    }

  }, error = function(e) {
    cat(sprintf("  MASC FAILED for %s: %s\n", comp$name, e$message))
  })
}

cat(sprintf("\n=== MASC (Downsampled) Complete: %s ===\n", Sys.time()))

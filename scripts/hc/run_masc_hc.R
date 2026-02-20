#!/usr/bin/env Rscript
# HC-only MASC: anno1 × sex, anno1 × age_group, anno2 × sex, anno2 × age_group
# Plus age×sex interaction analysis
#
# Uses run_masc_pipeline() from _wt/masc/myR/R/masc.R
#
# Usage: Rscript scripts/hc/run_masc_hc.R
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

cat("=== HC-only MASC Analysis ===\n")
cat(sprintf("Time: %s\n\n", Sys.time()))

# ---- Setup ----
MASC_R <- "/home/user3/data_user3/git_repo/mylit/Git_Repo/_wt/masc/myR/R/masc.R"
source(MASC_R)

sobj_path <- "/data/user3/sobj/hc_only_v1/2_hc_annotated.qs"
out_base <- "/data/user3/sobj/hc_only_v1/masc"
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

cat("Loading Seurat object...\n")
sobj <- qread(sobj_path)
cat(sprintf("  %d cells, %d patients\n", ncol(sobj), length(unique(sobj$name))))

# Ensure age_group exists
if (!"age_group" %in% colnames(sobj@meta.data)) {
  sobj$age_group <- ifelse(sobj$age < 35, "Young",
                           ifelse(sobj$age <= 50, "Middle", "Old"))
  cat(sprintf("  Created age_group: %s\n", paste(table(sobj$age_group), collapse = ", ")))
}

# ---- Define comparisons ----
comparisons <- list(
  list(
    name = "anno1_sex",
    cluster_var = "anno1",
    contrast_var = "sex",
    random_effects = "name",
    fixed_effects = "age"
  ),
  list(
    name = "anno1_age_group",
    cluster_var = "anno1",
    contrast_var = "age_group",
    random_effects = "name",
    fixed_effects = "sex"
  ),
  list(
    name = "anno2_sex",
    cluster_var = "anno2",
    contrast_var = "sex",
    random_effects = "name",
    fixed_effects = "age"
  ),
  list(
    name = "anno2_age_group",
    cluster_var = "anno2",
    contrast_var = "age_group",
    random_effects = "name",
    fixed_effects = "sex"
  )
)

# ---- Run MASC via run_masc_pipeline ----
all_results <- list()

for (comp in comparisons) {
  cat(sprintf("\n--- Running MASC: %s ---\n", comp$name))
  comp_dir <- file.path(out_base, comp$name)
  dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)

  tryCatch({
    res <- run_masc_pipeline(
      seurat_obj = sobj,
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
    cat(sprintf("  MASC completed: %d cell types tested\n", nrow(result)))

    # Summary
    if ("model.pvalue" %in% colnames(result)) {
      n_sig_nom <- sum(result$model.pvalue < 0.05, na.rm = TRUE)
      n_sig_fdr <- sum(result$model.pvalue.adj < 0.05, na.rm = TRUE)
      cat(sprintf("  Nominally significant (p<0.05): %d\n", n_sig_nom))
      cat(sprintf("  FDR significant (padj<0.05): %d\n", n_sig_fdr))

      if (n_sig_nom > 0) {
        cat("  Top hits:\n")
        result_sorted <- result[order(result$model.pvalue), ]
        top <- head(result_sorted[result_sorted$model.pvalue < 0.05, ], 10)
        for (i in seq_len(nrow(top))) {
          or_val <- if ("model.OR" %in% colnames(top)) top$model.OR[i] else NA
          cat(sprintf("    %-30s OR=%.3f p=%.2e\n",
                      rownames(top)[i], or_val, top$model.pvalue[i]))
        }
      }
    }

    all_results[[comp$name]] <- result

  }, error = function(e) {
    cat(sprintf("  MASC FAILED for %s: %s\n", comp$name, e$message))
  })
}

# ---- Age×Sex interaction (special) ----
cat("\n--- Running Age×Sex Interaction Analysis ---\n")
interaction_dir <- file.path(out_base, "anno1_age_sex_interaction")
dir.create(interaction_dir, recursive = TRUE, showWarnings = FALSE)

tryCatch({
  meta <- sobj@meta.data
  meta$cluster <- meta$anno1

  # Per-cell-type proportion analysis with interaction
  patient_meta <- meta %>%
    group_by(name) %>%
    summarise(
      age = first(age),
      sex = first(sex),
      age_group = first(age_group),
      n_cells = n(),
      .groups = "drop"
    )

  cell_types <- sort(unique(meta$cluster))
  interaction_results <- data.frame()

  for (ct in cell_types) {
    ct_counts <- meta %>%
      filter(cluster == ct) %>%
      group_by(name) %>%
      summarise(ct_cells = n(), .groups = "drop")

    pm <- left_join(patient_meta, ct_counts, by = "name")
    pm$ct_cells[is.na(pm$ct_cells)] <- 0
    pm$proportion <- pm$ct_cells / pm$n_cells

    if (sum(pm$ct_cells > 0) < 10) next

    m_full <- tryCatch(lm(proportion ~ age * sex, data = pm), error = function(e) NULL)
    m_reduced <- tryCatch(lm(proportion ~ age + sex, data = pm), error = function(e) NULL)

    if (!is.null(m_full) && !is.null(m_reduced)) {
      anova_res <- anova(m_reduced, m_full)
      interaction_p <- anova_res$`Pr(>F)`[2]

      coefs <- summary(m_full)$coefficients
      age_coef <- coefs["age", "Estimate"]
      sex_coef <- if ("sexM" %in% rownames(coefs)) coefs["sexM", "Estimate"] else NA
      interaction_coef <- if ("age:sexM" %in% rownames(coefs)) coefs["age:sexM", "Estimate"] else NA

      interaction_results <- rbind(interaction_results, data.frame(
        cell_type = ct,
        age_effect = age_coef,
        sex_effect = sex_coef,
        interaction_effect = interaction_coef,
        interaction_p = interaction_p,
        n_patients_nonzero = sum(pm$ct_cells > 0),
        mean_proportion = mean(pm$proportion),
        stringsAsFactors = FALSE
      ))
    }
  }

  interaction_results$interaction_padj <- p.adjust(interaction_results$interaction_p, method = "BH")
  interaction_results <- interaction_results[order(interaction_results$interaction_p), ]

  write.csv(interaction_results, file.path(interaction_dir, "age_sex_interaction.csv"), row.names = FALSE)

  n_sig <- sum(interaction_results$interaction_p < 0.05, na.rm = TRUE)
  n_fdr <- sum(interaction_results$interaction_padj < 0.1, na.rm = TRUE)
  cat(sprintf("  Age×Sex interaction: %d nominally sig, %d FDR<0.1\n", n_sig, n_fdr))

  if (n_sig > 0) {
    cat("  Top interactions:\n")
    top <- head(interaction_results[interaction_results$interaction_p < 0.05, ], 10)
    for (i in seq_len(nrow(top))) {
      cat(sprintf("    %-30s int_coef=%.4f p=%.3e padj=%.3e\n",
                  top$cell_type[i], top$interaction_effect[i],
                  top$interaction_p[i], top$interaction_padj[i]))
    }
  }

  # Interaction heatmap
  tryCatch({
    meta$age_bin <- cut(meta$age, breaks = c(18, 30, 40, 50, 66), include.lowest = TRUE)

    prop_matrix <- meta %>%
      group_by(name, anno1, age_bin, sex) %>%
      summarise(n = n(), .groups = "drop") %>%
      group_by(name) %>%
      mutate(prop = n / sum(n)) %>%
      ungroup() %>%
      group_by(anno1, age_bin, sex) %>%
      summarise(mean_prop = mean(prop), .groups = "drop")

    focus_cts <- interaction_results$cell_type[1:min(15, nrow(interaction_results))]
    prop_sub <- prop_matrix[prop_matrix$anno1 %in% focus_cts, ]
    prop_sub$anno1 <- factor(prop_sub$anno1, levels = rev(focus_cts))

    p_int <- ggplot(prop_sub, aes(x = interaction(age_bin, sex), y = anno1, fill = mean_prop)) +
      geom_tile(color = "white") +
      scale_fill_gradient(low = "white", high = "darkblue") +
      theme_minimal(base_size = 10) +
      theme(panel.background = element_rect(fill = "white", color = NA),
            plot.background = element_rect(fill = "white", color = NA),
            axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Cell type proportion by Age bin x Sex",
           x = "Age bin x Sex", y = "", fill = "Mean proportion")
    ggsave(file.path(interaction_dir, "age_sex_interaction_heatmap.png"),
           p_int, width = 12, height = 8, dpi = 200, bg = "white")
    cat("  Saved: age_sex_interaction_heatmap.png\n")
  }, error = function(e) cat(sprintf("  Interaction heatmap failed: %s\n", e$message)))

}, error = function(e) {
  cat(sprintf("  Interaction analysis FAILED: %s\n", e$message))
})

cat(sprintf("\n=== MASC Complete: %s ===\n", Sys.time()))
cat(sprintf("Results saved to: %s\n", out_base))

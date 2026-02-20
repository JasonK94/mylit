#!/usr/bin/env Rscript
# HC-only MNN (MultiNicheNet) Interpretation: Age and Sex
# Loads MNN results, extracts top DE genes, ligand activities, LR pairs
#
# Usage: Rscript scripts/hc/interpret_mnn.R
# Output: /data/user3/sobj/hc_only_v1/mnn_interpretation/

.libPaths(c(
  "/home/user3/GJC_KDW_250721/renv/library/R-4.3/x86_64-pc-linux-gnu",
  "/home/user3/R/x86_64-pc-linux-gnu-library/4.3",
  .libPaths()
))

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

cat("=== HC MNN Interpretation ===\n")
cat(sprintf("Time: %s\n\n", Sys.time()))

out_dir <- "/data/user3/sobj/hc_only_v1/mnn_interpretation"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Load MNN results ----
interpret_mnn <- function(mnn_path, label, out_path) {
  cat(sprintf("\n--- %s ---\n", label))

  mnn <- tryCatch(qread(mnn_path), error = function(e) {
    cat(sprintf("  FAILED to load: %s\n", e$message))
    return(NULL)
  })
  if (is.null(mnn)) return(NULL)

  cat(sprintf("  Names: %s\n", paste(names(mnn), collapse = ", ")))

  # 1. DE genes summary — celltype_de is a flat data.frame with columns:
  #    gene, cluster_id, logFC, logCPM, F, p_adj.loc, contrast, p_val, p_adj
  de_slot <- mnn$celltype_de
  if (is.data.frame(de_slot) || ("gene" %in% names(de_slot) && "cluster_id" %in% names(de_slot))) {
    # Convert list-of-vectors to data.frame
    de <- as.data.frame(de_slot, stringsAsFactors = FALSE)
    cat(sprintf("  DE genes: %d rows\n", nrow(de)))

    sig_de <- de %>% filter(p_adj < 0.05)
    cat(sprintf("  Significant DE genes (p_adj<0.05): %d\n", nrow(sig_de)))

    if (nrow(sig_de) > 0) {
      ct_summary <- sig_de %>%
        group_by(cluster_id) %>%
        summarise(
          n_de = n(),
          n_up = sum(logFC > 0, na.rm = TRUE),
          n_down = sum(logFC < 0, na.rm = TRUE),
          top_up = paste(head(gene[order(-logFC)], 3), collapse = ","),
          top_down = paste(head(gene[order(logFC)], 3), collapse = ","),
          .groups = "drop"
        ) %>%
        arrange(desc(n_de))

      write.csv(ct_summary, file.path(out_path, sprintf("%s_de_by_celltype.csv", label)),
                row.names = FALSE)
      cat("  Top cell types by DE genes:\n")
      for (i in seq_len(min(10, nrow(ct_summary)))) {
        r <- ct_summary[i, ]
        cat(sprintf("    %-30s %d DE (up=%d, down=%d) | top up: %s\n",
                    r$cluster_id, r$n_de, r$n_up, r$n_down, r$top_up))
      }
    }
  }

  # 2. Ligand activities
  if ("ligand_activities_targets_DEgenes" %in% names(mnn)) {
    la <- mnn$ligand_activities_targets_DEgenes$ligand_activities
    if (!is.null(la)) {
      cat(sprintf("  Ligand activities: %d rows\n", nrow(la)))

      # Top ligands by activity
      activity_col <- intersect(c("activity", "aupr_corrected", "pearson"), colnames(la))[1]
      if (!is.na(activity_col)) {
        top_lig <- la %>%
          arrange(desc(!!sym(activity_col))) %>%
          head(20)
        cat("  Top ligands by activity:\n")
        for (i in seq_len(min(10, nrow(top_lig)))) {
          r <- top_lig[i, ]
          sender <- if ("sender" %in% colnames(r)) r$sender else "?"
          receiver <- if ("receiver" %in% colnames(r)) r$receiver else "?"
          cat(sprintf("    %-12s %s → %s (activity=%.3f)\n",
                      r$test_ligand, sender, receiver, r[[activity_col]]))
        }
        write.csv(top_lig, file.path(out_path, sprintf("%s_top_ligands.csv", label)),
                  row.names = FALSE)
      }
    }
  }

  # 3. Group prioritization
  if ("prioritization_tables" %in% names(mnn)) {
    gp <- mnn$prioritization_tables$group_prioritization_tbl
    if (!is.null(gp) && nrow(gp) > 0) {
      cat(sprintf("  Group prioritization: %d rows\n", nrow(gp)))
      write.csv(gp, file.path(out_path, sprintf("%s_group_prioritization.csv", label)),
                row.names = FALSE)
    } else {
      cat("  Group prioritization: EMPTY (known issue with strict filters)\n")
    }
  }

  # 4. LR pairs
  if ("prioritization_tables" %in% names(mnn)) {
    lr <- mnn$prioritization_tables$prioritization_tbl_ligand_receptor
    if (!is.null(lr) && nrow(lr) > 0) {
      cat(sprintf("  LR pairs: %d rows\n", nrow(lr)))
      # Top LR pairs
      if ("prioritization_score" %in% colnames(lr)) {
        top_lr <- lr %>% arrange(desc(prioritization_score)) %>% head(20)
      } else {
        top_lr <- head(lr, 20)
      }
      write.csv(top_lr, file.path(out_path, sprintf("%s_top_lr_pairs.csv", label)),
                row.names = FALSE)
    }
  }

  return(mnn)
}

# ---- Age MNN ----
age_mnn <- interpret_mnn(
  "/data/user3/sobj/hc_only_v1/mnn/age_v1/multinichenet_results.qs",
  "age", out_dir
)

# ---- Sex MNN ----
sex_mnn <- interpret_mnn(
  "/data/user3/sobj/hc_only_v1/mnn/sex_v1/multinichenet_results.qs",
  "sex", out_dir
)

cat(sprintf("\n=== MNN Interpretation Complete: %s ===\n", Sys.time()))
cat(sprintf("Output: %s\n", out_dir))

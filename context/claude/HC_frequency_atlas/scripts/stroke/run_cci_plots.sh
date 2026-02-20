#!/bin/bash
# ==============================================================================
# CCI Visualization: Run comparison plots for CellChat + MNN
# ==============================================================================

set -e
BASE="/data/user3/sobj/stroke_hc_v8_2"
SCRIPTS="/data/user3/git_repo/mylit/Git_Repo/_wt/cellchat/scripts/cellchat"
MNN_SCRIPTS="/data/user3/git_repo/mylit/Git_Repo/_wt/cci/scripts/cci/mnn"

echo "=== CCI Visualization Started: $(date) ==="

# ---- 1. CellChat Comparison Plots ----

# Layer 2: g3 comparison (X1=good vs X2=bad)
echo ""
echo ">>> [1/4] CellChat L2: g3 comparison plots"
Rscript "$SCRIPTS/run_cellchat_comparison_2a.R" \
  --control "$BASE/cci/cellchat/L2_g3_anno2/merged/X1/cellchat_merged.qs" \
  --disease "$BASE/cci/cellchat/L2_g3_anno2/merged/X2/cellchat_merged.qs" \
  --name1 "g3_Good" --name2 "g3_Bad" \
  --output_dir "$BASE/cci/plots/cellchat_L2_g3" \
  --top_n 20
echo "    >>> CellChat L2 plots done: $(date)"

# Layer 1: HC vs Stroke
echo ""
echo ">>> [2/4] CellChat L1: HC vs Stroke comparison plots"
Rscript "$SCRIPTS/run_cellchat_comparison_2a.R" \
  --control "$BASE/cci/cellchat/L1_cohort_anno2/merged/HC/cellchat_merged.qs" \
  --disease "$BASE/cci/cellchat/L1_cohort_anno2/merged/Stroke/cellchat_merged.qs" \
  --name1 "HC" --name2 "Stroke" \
  --output_dir "$BASE/cci/plots/cellchat_L1_cohort" \
  --top_n 20
echo "    >>> CellChat L1 plots done: $(date)"

# ---- 2. CellChat Individual Plots ----

echo ""
echo ">>> [3/4] CellChat individual plots (per condition)"
for cond_dir in "$BASE/cci/cellchat/L2_g3_anno2/merged/X1" "$BASE/cci/cellchat/L2_g3_anno2/merged/X2" \
                "$BASE/cci/cellchat/L1_cohort_anno2/merged/HC" "$BASE/cci/cellchat/L1_cohort_anno2/merged/Stroke"; do
  cond_name=$(basename "$cond_dir")
  layer_name=$(basename "$(dirname "$(dirname "$cond_dir")")")
  out_name="${layer_name}_${cond_name}"
  echo "  Plotting: $out_name"
  Rscript "$SCRIPTS/plot_cellchat_individual.R" \
    --input "$cond_dir/cellchat_merged.qs" \
    --output_dir "$BASE/cci/plots/individual_${out_name}" \
    --top_links 100 2>&1 | tail -3
done
echo "    >>> Individual plots done: $(date)"

# ---- 3. MNN Results Visualization ----

echo ""
echo ">>> [4/4] MultiNicheNet results visualization"
Rscript -e '
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
})

BASE <- "/data/user3/sobj/stroke_hc_v8_2"

for (layer in c("L2_g3_anno2", "L1_cohort_anno2")) {
  cat("\n=== MNN:", layer, "===\n")
  mnn_file <- file.path(BASE, "cci/mnn", layer, "multinichenet_results.qs")
  if (!file.exists(mnn_file)) { cat("  Not found. Skipping.\n"); next }

  mnn <- qread(mnn_file)
  out_dir <- file.path(BASE, "cci/plots", paste0("mnn_", layer))
  dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

  cat("  Result class:", class(mnn), "\n")
  cat("  Names:", paste(names(mnn), collapse=", "), "\n")

  # Extract prioritized LR interactions
  if ("prioritization_tbl_ligand_receptor" %in% names(mnn)) {
    lr_tbl <- mnn$prioritization_tbl_ligand_receptor
    cat("  LR pairs:", nrow(lr_tbl), "\n")
    write.csv(lr_tbl, file.path(out_dir, "lr_prioritization.csv"), row.names=FALSE)

    # Top LR pairs by priority score
    if ("prioritization_score" %in% colnames(lr_tbl)) {
      top_lr <- lr_tbl %>%
        arrange(desc(prioritization_score)) %>%
        head(30)

      p <- ggplot(top_lr, aes(x=reorder(paste0(ligand, " -> ", receptor), prioritization_score),
                               y=prioritization_score)) +
        geom_col(aes(fill=sender), width=0.7) +
        coord_flip() + theme_minimal() +
        ggtitle(paste("Top 30 LR Pairs -", layer)) +
        xlab("") + ylab("Prioritization Score") +
        theme(axis.text.y=element_text(size=7))
      ggsave(file.path(out_dir, "top_lr_pairs.png"), p, width=10, height=8, dpi=150)
      cat("  Saved top_lr_pairs.png\n")
    }
  }

  # Extract prioritized LR + celltype interactions
  if ("prioritization_tbl_ligand_receptor_target" %in% names(mnn)) {
    lrt_tbl <- mnn$prioritization_tbl_ligand_receptor_target
    cat("  LR-Target triplets:", nrow(lrt_tbl), "\n")
    write.csv(head(lrt_tbl, 500), file.path(out_dir, "lrt_prioritization_top500.csv"),
              row.names=FALSE)
  }

  # DE info
  if ("celltype_de" %in% names(mnn)) {
    de_tbl <- mnn$celltype_de
    cat("  DE results:", nrow(de_tbl), "\n")
    # Top DE genes per celltype
    if ("logFC" %in% colnames(de_tbl) & "p_adj" %in% colnames(de_tbl)) {
      sig_de <- de_tbl %>% filter(p_adj < 0.05) %>%
        group_by(celltype) %>%
        summarise(n_up=sum(logFC > 0.25), n_down=sum(logFC < -0.25), .groups="drop")
      write.csv(sig_de, file.path(out_dir, "de_summary_by_celltype.csv"), row.names=FALSE)
      cat("  Saved DE summary\n")

      p_de <- ggplot(sig_de %>% pivot_longer(c(n_up, n_down)),
                     aes(x=celltype, y=value, fill=name)) +
        geom_col(position="dodge") +
        scale_fill_manual(values=c("n_up"="#B2182B","n_down"="#2166AC"),
                          labels=c("Up","Down"), name="Direction") +
        theme_minimal() + coord_flip() +
        ggtitle(paste("DE genes by celltype -", layer)) +
        ylab("N significant genes (|logFC|>0.25, padj<0.05)")
      ggsave(file.path(out_dir, "de_by_celltype.png"), p_de, width=8, height=5, dpi=150)
      cat("  Saved de_by_celltype.png\n")
    }
  }

  cat("  Done:", layer, "\n")
}
cat("\nMNN visualization complete.\n")
'

echo ""
echo "=== CCI Visualization Finished: $(date) ==="

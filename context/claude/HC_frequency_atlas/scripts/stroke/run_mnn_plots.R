#!/usr/bin/env Rscript
# ==============================================================================
# MultiNicheNet Results Visualization
# ==============================================================================

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

for (layer_info in list(
  list(id="L2_g3_anno2", label="g3 (IS: Bad vs Good)", contrast_label="g3==2 vs g3==1"),
  list(id="L1_cohort_anno2", label="Cohort (Stroke vs HC)", contrast_label="Stroke vs HC")
)) {
  cat("\n", paste(rep("=", 60), collapse=""), "\n")
  cat("MNN:", layer_info$label, "\n")
  cat(paste(rep("=", 60), collapse=""), "\n")

  mnn_file <- file.path(BASE, "cci/mnn", layer_info$id, "multinichenet_results.qs")
  if (!file.exists(mnn_file)) { cat("  Not found.\n"); next }

  mnn <- qread(mnn_file)
  out_dir <- file.path(BASE, "cci/plots", paste0("mnn_", layer_info$id))
  dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

  # ---- 1. DE summary by cell type ----
  de <- mnn$celltype_de
  cat("  DE genes:", nrow(de), "\n")

  sig_de <- de %>%
    filter(p_adj < 0.05) %>%
    group_by(cluster_id) %>%
    summarise(n_up=sum(logFC > 0.25), n_down=sum(logFC < -0.25), .groups="drop")
  write.csv(sig_de, file.path(out_dir, "de_summary_by_celltype.csv"), row.names=FALSE)
  cat("  DE summary:\n")
  print(sig_de)

  sig_long <- sig_de %>% pivot_longer(c(n_up, n_down), names_to="direction", values_to="count")
  p_de <- ggplot(sig_long, aes(x=cluster_id, y=count, fill=direction)) +
    geom_col(position="dodge") +
    scale_fill_manual(values=c("n_up"="#B2182B","n_down"="#2166AC"),
                      labels=c("Down","Up"), name="Direction") +
    theme_minimal() + coord_flip() +
    ggtitle(paste("DE genes -", layer_info$contrast_label)) +
    ylab("N significant (|logFC|>0.25, padj<0.05)") + xlab("")
  ggsave(file.path(out_dir, "de_by_celltype.png"), p_de, width=8, height=5, dpi=150)

  # Volcano per celltype
  for (ct in unique(de$cluster_id)) {
    de_ct <- de %>% filter(cluster_id == ct)
    de_ct$sig <- ifelse(de_ct$p_adj < 0.05 & abs(de_ct$logFC) > 0.25,
                        ifelse(de_ct$logFC > 0, "Up", "Down"), "NS")
    p_vol <- ggplot(de_ct, aes(x=logFC, y=-log10(p_adj), color=sig)) +
      geom_point(size=0.5, alpha=0.5) +
      scale_color_manual(values=c("Up"="#B2182B","Down"="#2166AC","NS"="grey70")) +
      theme_minimal() + ggtitle(paste(ct, "-", layer_info$contrast_label)) +
      geom_hline(yintercept=-log10(0.05), linetype="dashed", alpha=0.3) +
      geom_vline(xintercept=c(-0.25, 0.25), linetype="dashed", alpha=0.3)
    ggsave(file.path(out_dir, paste0("volcano_", gsub("/","_",ct), ".png")),
           p_vol, width=6, height=5, dpi=150)
  }
  cat("  Saved volcano plots\n")

  # ---- 2. Prioritized LR interactions ----
  prio <- mnn$prioritization_tables$group_prioritization_tbl
  cat("  Prioritized LR interactions:", nrow(prio), "\n")
  write.csv(prio, file.path(out_dir, "lr_prioritization.csv"), row.names=FALSE)

  # Create prioritization score (composite)
  if (all(c("scaled_lfc_ligand", "scaled_lfc_receptor", "max_scaled_activity") %in% colnames(prio))) {
    prio <- prio %>%
      mutate(priority_score = scaled_lfc_ligand + scaled_lfc_receptor + max_scaled_activity +
               scaled_pb_ligand + scaled_pb_receptor)

    top30 <- prio %>%
      arrange(desc(priority_score)) %>%
      head(30) %>%
      mutate(lr_label = paste0(sender, ": ", ligand, " -> ", receptor, " (", receiver, ")"))

    p_lr <- ggplot(top30, aes(x=reorder(lr_label, priority_score),
                               y=priority_score, fill=sender)) +
      geom_col(width=0.7) + coord_flip() + theme_minimal() +
      ggtitle(paste("Top 30 LR -", layer_info$contrast_label)) +
      xlab("") + ylab("Priority Score (composite)") +
      theme(axis.text.y=element_text(size=7))
    ggsave(file.path(out_dir, "top_lr_pairs.png"), p_lr, width=10, height=8, dpi=150)
    cat("  Saved top_lr_pairs.png\n")

    # Sender-Receiver heatmap of top interactions
    sr_counts <- prio %>%
      filter(priority_score > quantile(priority_score, 0.9, na.rm=TRUE)) %>%
      count(sender, receiver)
    if (nrow(sr_counts) > 0) {
      p_sr <- ggplot(sr_counts, aes(x=receiver, y=sender, fill=n)) +
        geom_tile() + scale_fill_viridis_c() +
        theme_minimal() +
        ggtitle(paste("Top 10% LR: Sender-Receiver -", layer_info$contrast_label)) +
        theme(axis.text.x=element_text(angle=45, hjust=1))
      ggsave(file.path(out_dir, "sender_receiver_heatmap.png"), p_sr, width=7, height=6, dpi=150)
      cat("  Saved sender_receiver_heatmap.png\n")
    }
  }

  # ---- 3. Ligand activity per receiver ----
  act_tbl <- mnn$prioritization_tables$ligand_activities_target_de_tbl
  if (!is.null(act_tbl) && nrow(act_tbl) > 0) {
    top_act <- act_tbl %>%
      distinct(receiver, ligand, activity_scaled) %>%
      group_by(receiver) %>%
      arrange(desc(activity_scaled)) %>%
      slice_head(n=10) %>%
      ungroup()

    p_act <- ggplot(top_act, aes(x=reorder(ligand, activity_scaled),
                                  y=activity_scaled, fill=receiver)) +
      geom_col() + coord_flip() + facet_wrap(~receiver, scales="free_y") +
      theme_minimal() +
      ggtitle(paste("Top Ligand Activities -", layer_info$contrast_label)) +
      xlab("") + ylab("Scaled Activity") +
      theme(axis.text.y=element_text(size=6), legend.position="none")
    ggsave(file.path(out_dir, "ligand_activity_by_receiver.png"),
           p_act, width=12, height=8, dpi=150)
    cat("  Saved ligand_activity_by_receiver.png\n")
  }

  # ---- 4. Sample-level LR variation ----
  samp_tbl <- mnn$prioritization_tables$sample_prioritization_tbl
  if (!is.null(samp_tbl) && nrow(samp_tbl) > 0) {
    # Top 10 LR pairs by mean production
    top_lr_samp <- samp_tbl %>%
      group_by(sender, receiver, ligand, receptor) %>%
      summarise(mean_prod=mean(ligand_receptor_pb_prod, na.rm=TRUE), .groups="drop") %>%
      arrange(desc(mean_prod)) %>%
      head(10)

    top_lr_ids <- paste0(top_lr_samp$ligand, "_", top_lr_samp$receptor)

    samp_filtered <- samp_tbl %>%
      mutate(lr_id=paste0(ligand, "_", receptor)) %>%
      filter(lr_id %in% top_lr_ids)

    if (nrow(samp_filtered) > 0) {
      p_samp <- ggplot(samp_filtered,
                       aes(x=paste0(ligand, "->", receptor),
                           y=ligand_receptor_pb_prod, fill=group)) +
        geom_boxplot(outlier.size=0.5) + coord_flip() + theme_minimal() +
        ggtitle(paste("Top 10 LR by sample -", layer_info$contrast_label)) +
        xlab("") + ylab("LR pseudobulk product")
      ggsave(file.path(out_dir, "top_lr_sample_boxplot.png"),
             p_samp, width=10, height=6, dpi=150)
      cat("  Saved top_lr_sample_boxplot.png\n")
    }
  }

  cat("  Done:", layer_info$id, "\n")
}

cat("\nMNN visualization complete.\n")

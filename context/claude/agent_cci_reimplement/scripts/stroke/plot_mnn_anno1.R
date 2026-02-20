#!/usr/bin/env Rscript
# ==============================================================================
# MNN anno1 Visualization — DE volcanos, ligand activities, circos from LA table
# ==============================================================================

.libPaths(c(
  "/home/user3/GJC_KDW_250721/renv/library/R-4.3/x86_64-pc-linux-gnu",
  "/home/user3/R/x86_64-pc-linux-gnu-library/4.3",
  .libPaths()
))

suppressPackageStartupMessages({
  library(qs)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(circlize)
  library(RColorBrewer)
  library(patchwork)
  library(ggrepel)
})

cat("=== MNN anno1 Visualization ===\n")
cat("Started:", format(Sys.time()), "\n\n")

BASE <- "/data/user3/sobj/stroke_hc_v8_2"

for (layer in c("L1_cohort", "L2_g3")) {
  mnn_file <- file.path(BASE, "cci/mnn", paste0(layer, "_anno1"), "multinichenet_results.qs")
  if (!file.exists(mnn_file)) {
    cat("  MNN file not found:", mnn_file, "\n")
    next
  }

  cat("\n", strrep("=", 50), "\n")
  cat("  MNN", layer, "anno1\n")
  cat(strrep("=", 50), "\n")

  mnn <- qread(mnn_file)
  out_dir <- file.path(BASE, "cci/plots", paste0("mnn_", layer, "_anno1"))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # ---- 1. DE volcanos ----
  de_genes <- mnn$ligand_activities_targets_DEgenes$de_genes_df
  if (!is.null(de_genes) && nrow(de_genes) > 0) {
    cat("  DE genes:", nrow(de_genes), "\n")

    # Summary
    de_summary <- de_genes %>%
      group_by(receiver) %>%
      summarise(
        n_genes = n(),
        n_sig = sum(p_adj < 0.05, na.rm = TRUE),
        n_up = sum(p_adj < 0.05 & logFC > 0.25, na.rm = TRUE),
        n_down = sum(p_adj < 0.05 & logFC < -0.25, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(desc(n_sig))

    write.csv(de_summary, file.path(out_dir, "de_summary_by_celltype.csv"), row.names = FALSE)
    cat("  Cell types with sig DE:\n")
    print(as.data.frame(head(de_summary, 10)))

    # DE bar summary
    p_bar <- ggplot(de_summary %>% filter(n_sig > 0),
                    aes(x = reorder(receiver, n_sig), y = n_sig)) +
      geom_col(fill = "#d73027", alpha = 0.8) +
      geom_text(aes(label = n_sig), hjust = -0.1, size = 3) +
      coord_flip() +
      theme_minimal(base_size = 11) +
      labs(title = paste("DE Genes per Cell Type -", layer, "anno1"),
           subtitle = paste("Total:", sum(de_summary$n_sig), "significant genes"),
           x = NULL, y = "N significant (p_adj < 0.05)")
    ggsave(file.path(out_dir, "01_de_summary_bar.png"), p_bar, width = 8, height = 6, dpi = 300)

    # Volcano per cell type (top 9)
    top_cts <- head(de_summary$receiver[de_summary$n_sig > 0], 9)
    if (length(top_cts) > 0) {
      plot_list <- list()
      for (ct in top_cts) {
        ct_data <- de_genes %>% filter(receiver == ct)
        ct_data$sig <- ct_data$p_adj < 0.05 & abs(ct_data$logFC) > 0.25
        # Top labels
        top_lab <- ct_data %>% filter(sig) %>%
          arrange(p_adj) %>% head(5)

        p <- ggplot(ct_data, aes(x = logFC, y = -log10(p_adj))) +
          geom_point(aes(color = sig), size = 0.3, alpha = 0.4) +
          scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey70"), guide = "none") +
          geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", linewidth = 0.3) +
          geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "grey50", linewidth = 0.3) +
          theme_minimal(base_size = 8) +
          labs(title = ct, x = "logFC", y = "-log10(p_adj)")

        if (nrow(top_lab) > 0) {
          p <- p + geom_text_repel(
            data = top_lab, aes(label = gene),
            size = 2, max.overlaps = 10, segment.size = 0.2
          )
        }
        plot_list[[ct]] <- p
      }

      n_col <- min(3, length(plot_list))
      n_row <- ceiling(length(plot_list) / n_col)
      combined <- wrap_plots(plot_list, ncol = n_col)
      ggsave(file.path(out_dir, "02_volcano_grid.png"), combined,
             width = 4 * n_col, height = 3.5 * n_row, dpi = 300)
      cat("  Volcano grid saved (", length(top_cts), "cell types)\n")
    }
  }

  # ---- 2. Ligand activities ----
  la <- mnn$ligand_activities_targets_DEgenes$ligand_activities
  if (!is.null(la) && nrow(la) > 0) {
    cat("  Ligand activities:", nrow(la), "rows\n")

    # Top ligands per receiver (by mean activity_scaled)
    la_top <- la %>%
      group_by(receiver, ligand) %>%
      summarise(mean_activity = mean(activity_scaled, na.rm = TRUE), .groups = "drop") %>%
      group_by(receiver) %>%
      slice_max(mean_activity, n = 10) %>%
      ungroup()

    # Pick top 6 receivers by number of high-activity ligands
    top_receivers <- la_top %>%
      group_by(receiver) %>%
      summarise(max_act = max(mean_activity), .groups = "drop") %>%
      arrange(desc(max_act)) %>%
      head(6) %>%
      pull(receiver)

    la_plot <- la_top %>% filter(receiver %in% top_receivers)

    p_la <- ggplot(la_plot, aes(x = reorder(ligand, mean_activity), y = mean_activity)) +
      geom_col(fill = "#4575b4", alpha = 0.8) +
      coord_flip() +
      facet_wrap(~receiver, scales = "free_y", ncol = 3) +
      theme_minimal(base_size = 9) +
      labs(title = paste("Top Ligand Activities -", layer, "anno1"),
           x = NULL, y = "Mean Activity (scaled)")
    ggsave(file.path(out_dir, "03_ligand_activity_grid.png"), p_la, width = 14, height = 8, dpi = 300)
    cat("  Ligand activity grid saved\n")
  }

  # ---- 3. Ligand-activity-target-DE table plots ----
  # NOTE: group_prioritization_tbl is empty (0 rows) because DE was too sparse
  # for the prioritization filter. Use ligand_activities_target_de_tbl instead.
  la_de <- mnn$prioritization_tables$ligand_activities_target_de_tbl
  if (!is.null(la_de) && nrow(la_de) > 0) {
    cat("  LA-DE table:", nrow(la_de), "rows\n")
    cat("  Columns:", paste(colnames(la_de), collapse = ", "), "\n")

    act_col <- intersect(c("activity_scaled", "activity"), colnames(la_de))[1]
    if (!is.na(act_col)) {
      cat("  Using activity column:", act_col, "\n")

      # Top ligand-receiver pairs by activity
      lr_top <- la_de %>%
        group_by(receiver, ligand) %>%
        summarise(mean_act = mean(!!rlang::sym(act_col), na.rm = TRUE),
                  n_targets = n_distinct(target),
                  .groups = "drop") %>%
        arrange(desc(mean_act)) %>%
        head(40)

      if (nrow(lr_top) > 0) {
        p_lr <- ggplot(lr_top, aes(
          x = reorder(paste0(ligand, " -> ", receiver), mean_act),
          y = mean_act, fill = receiver
        )) +
          geom_col() + coord_flip() +
          scale_fill_brewer(palette = "Set3") +
          theme_minimal(base_size = 10) +
          labs(title = paste("Top 40 Ligand -> Receiver -", layer, "anno1"),
               x = NULL, y = "Mean Activity (scaled)")
        ggsave(file.path(out_dir, "04_top_ligand_receiver.png"), p_lr, width = 10, height = 10, dpi = 300)
        cat("  Top ligand-receiver bar saved\n")
      }

      # Ligand -> Receiver circos
      lr_recv_agg <- la_de %>%
        group_by(ligand, receiver) %>%
        summarise(total_act = sum(!!rlang::sym(act_col), na.rm = TRUE), .groups = "drop") %>%
        arrange(desc(total_act)) %>%
        head(60) %>%
        select(from = ligand, to = receiver, weight = total_act)

      if (nrow(lr_recv_agg) > 0) {
        receivers_u <- unique(lr_recv_agg$to)
        n_recv <- length(receivers_u)
        recv_pal <- colorRampPalette(brewer.pal(min(n_recv, 12), "Set3"))(n_recv)
        names(recv_pal) <- receivers_u

        grid_cols <- c()
        for (lig in unique(lr_recv_agg$from)) grid_cols[lig] <- "grey60"
        for (recv in receivers_u) grid_cols[recv] <- recv_pal[recv]

        for (ext in c("pdf", "png")) {
          fname <- file.path(out_dir, paste0("05_circos_ligand_receiver.", ext))
          if (ext == "pdf") pdf(fname, width = 11, height = 11)
          else png(fname, width = 11, height = 11, units = "in", res = 300)

          circos.clear()
          circos.par(gap.after = 2)
          chordDiagram(lr_recv_agg, grid.col = grid_cols, transparency = 0.3,
                       directional = 1, direction.type = c("diffHeight", "arrows"),
                       link.arr.type = "big.arrow",
                       annotationTrack = "grid",
                       preAllocateTracks = list(track.height = 0.08))
          circos.track(track.index = 1, panel.fun = function(x, y) {
            circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                        facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.55)
          }, bg.border = NA)
          title(paste("Top Ligand -> Receiver CCI -", layer, "anno1"))
          legend("topright", legend = receivers_u, fill = recv_pal,
                 title = "Receiver", bty = "n", cex = 0.7)
          dev.off()
        }
        cat("  Circos (ligand->receiver) saved\n")
      }

      # Receiver-level activity heatmap: receiver × ligand
      la_wide <- la_de %>%
        group_by(receiver, ligand) %>%
        summarise(act = mean(!!rlang::sym(act_col), na.rm = TRUE), .groups = "drop")

      # Top ligands overall
      top_ligs <- la_wide %>% group_by(ligand) %>%
        summarise(max_act = max(act), .groups = "drop") %>%
        arrange(desc(max_act)) %>% head(25) %>% pull(ligand)

      # Top receivers
      top_recvs <- la_wide %>% group_by(receiver) %>%
        summarise(max_act = max(act), .groups = "drop") %>%
        arrange(desc(max_act)) %>% head(15) %>% pull(receiver)

      heat_data <- la_wide %>%
        filter(ligand %in% top_ligs, receiver %in% top_recvs)

      if (nrow(heat_data) > 0) {
        p_heat <- ggplot(heat_data, aes(x = receiver, y = ligand, fill = act)) +
          geom_tile() +
          scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027",
                               midpoint = median(heat_data$act), name = "Activity") +
          theme_minimal(base_size = 9) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          labs(title = paste("Ligand Activity Heatmap -", layer, "anno1"),
               x = "Receiver", y = "Ligand")
        ggsave(file.path(out_dir, "06_ligand_receiver_heatmap.png"),
               p_heat, width = 10, height = 8, dpi = 300)
        cat("  Heatmap saved\n")
      }

      # Save top results
      write.csv(head(la_de %>% arrange(desc(!!rlang::sym(act_col))), 500),
                file.path(out_dir, "la_de_top500.csv"), row.names = FALSE)
    }
  }

  rm(mnn); gc(verbose = FALSE)
  cat("  Done for", layer, "\n")
}

# ---- Summary ----
cat("\n", strrep("=", 50), "\n")
cat("MNN anno1 Visualization Complete\n")
cat("Finished:", format(Sys.time()), "\n\n")
for (d in c("mnn_L1_cohort_anno1", "mnn_L2_g3_anno1")) {
  full_path <- file.path(BASE, "cci/plots", d)
  if (dir.exists(full_path)) {
    n_png <- length(list.files(full_path, pattern = "\\.png$", recursive = TRUE))
    n_pdf <- length(list.files(full_path, pattern = "\\.pdf$", recursive = TRUE))
    cat("  ", d, ":", n_png, "PNG +", n_pdf, "PDF\n")
  }
}

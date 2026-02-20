#!/usr/bin/env Rscript
# ==============================================================================
# Comprehensive CCI Visualization
# Part 1: CellChat native plots (bubble, rankNet/river, signaling role, chord)
#          - Build condition-level CellChat from individual sample objects
# Part 2: MNN anno1 circos + DE + ligand activity plots
# ==============================================================================

.libPaths(c(
  "/home/user3/GJC_KDW_250721/renv/library/R-4.3/x86_64-pc-linux-gnu",
  "/home/user3/R/x86_64-pc-linux-gnu-library/4.3",
  .libPaths()
))

suppressPackageStartupMessages({
  library(CellChat)
  library(qs)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(circlize)
  library(RColorBrewer)
})

cat("=== Comprehensive CCI Visualization ===\n")
cat("Started:", format(Sys.time()), "\n\n")

BASE <- "/data/user3/sobj/stroke_hc_v8_2"

# ==============================================================================
# PART 1: CellChat Native Plots
# ==============================================================================

run_cellchat_native_plots <- function(samples_dir, condition_map,
                                       control_name, disease_name,
                                       layer_label, out_dir) {
  cat("\n", strrep("=", 60), "\n")
  cat("  CellChat:", layer_label, "\n")
  cat(strrep("=", 60), "\n")

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  pw_dir <- file.path(out_dir, "pathways")
  dir.create(pw_dir, showWarnings = FALSE)

  # Load all sample CellChat objects
  all_samples <- list.dirs(samples_dir, recursive = FALSE, full.names = FALSE)
  cat("  Total samples:", length(all_samples), "\n")

  cc_all <- list()
  for (sn in all_samples) {
    f <- file.path(samples_dir, sn, "cellchat.qs")
    if (file.exists(f)) {
      cc_all[[sn]] <- qread(f)
    }
  }
  cat("  Loaded:", length(cc_all), "CellChat objects\n")

  # Split by condition
  ctrl_samples <- names(condition_map)[condition_map == control_name]
  dis_samples <- names(condition_map)[condition_map == disease_name]
  cat("  ", control_name, ":", length(ctrl_samples), "samples\n")
  cat("  ", disease_name, ":", length(dis_samples), "samples\n")

  # Pick representative (largest) sample per condition
  pick_rep <- function(sample_list) {
    avail <- intersect(sample_list, names(cc_all))
    if (length(avail) == 0) return(NULL)
    sizes <- sapply(cc_all[avail], function(x) length(x@idents))
    nm <- names(which.max(sizes))
    cat("    Representative:", nm, "(", max(sizes), "cells)\n")
    return(cc_all[[nm]])
  }

  cat("  Picking representative samples...\n")
  cc_ctrl <- pick_rep(ctrl_samples)
  cc_dis <- pick_rep(dis_samples)

  if (is.null(cc_ctrl) || is.null(cc_dis)) {
    cat("  ERROR: Could not find representative samples\n")
    return(invisible(NULL))
  }

  # Compute centrality
  cc_ctrl <- netAnalysis_computeCentrality(cc_ctrl)
  cc_dis <- netAnalysis_computeCentrality(cc_dis)

  # --- Merge ---
  cat("  Merging representative CellChat objects...\n")
  object.list <- list()
  object.list[[control_name]] <- cc_ctrl
  object.list[[disease_name]] <- cc_dis

  cellchat <- tryCatch({
    mergeCellChat(object.list, add.names = names(object.list))
  }, error = function(e) {
    cat("  WARNING: mergeCellChat failed:", e$message, "\n")
    NULL
  })

  # --- Individual condition plots (no merge needed) ---

  # 1. River / signaling role per condition
  cat("  [1] River plots (per condition)...\n")
  for (i in seq_along(object.list)) {
    nm <- names(object.list)[i]
    obj <- object.list[[i]]
    for (pat in c("outgoing", "incoming")) {
      tryCatch({
        fname <- paste0("15_river_", nm, "_", pat, ".png")
        png(file.path(out_dir, fname), width = 10, height = 8, units = "in", res = 300)
        netAnalysis_river(obj, pattern = pat)
        title(paste(nm, "-", pat, "signaling"))
        dev.off()
        cat("    Saved:", fname, "\n")
      }, error = function(e) cat("    river", nm, pat, "failed:", e$message, "\n"))
    }
  }

  # 2. Signaling role heatmap per condition
  cat("  [2] Signaling role heatmap...\n")
  for (i in seq_along(object.list)) {
    nm <- names(object.list)[i]
    obj <- object.list[[i]]
    for (pat in c("outgoing", "incoming")) {
      tryCatch({
        fname <- paste0("16_sigRole_", nm, "_", pat, ".png")
        png(file.path(out_dir, fname), width = 10, height = 6, units = "in", res = 300)
        netAnalysis_signalingRole_heatmap(obj, pattern = pat, title = paste(nm, pat),
                                           width = 8, height = 5)
        dev.off()
        cat("    Saved:", fname, "\n")
      }, error = function(e) cat("    sigRole", nm, pat, "failed:", e$message, "\n"))
    }
  }

  # 3. Pathway circle plots per condition
  cat("  [3] Pathway circle plots...\n")
  pathways_ctrl <- cc_ctrl@netP$pathways
  pathways_dis <- cc_dis@netP$pathways
  common_pathways <- intersect(pathways_ctrl, pathways_dis)
  cat("    Common pathways:", length(common_pathways), "\n")

  for (pw in head(common_pathways, 12)) {
    for (i in seq_along(object.list)) {
      nm <- names(object.list)[i]
      tryCatch({
        fname <- paste0("pathway_", pw, "_", nm, ".png")
        png(file.path(pw_dir, fname), width = 7, height = 7, units = "in", res = 200)
        netVisual_aggregate(object.list[[i]], signaling = pw,
                            layout = "circle", vertex.label.cex = 0.7)
        title(paste(pw, "-", nm), cex.main = 1.2)
        dev.off()
      }, error = function(e) NULL)
    }
    cat("    ", pw, "\n")
  }

  # 4. Contribution plots for key pathways
  cat("  [4] Pathway contribution...\n")
  for (pw in head(common_pathways, 8)) {
    for (i in seq_along(object.list)) {
      nm <- names(object.list)[i]
      tryCatch({
        fname <- paste0("contrib_", pw, "_", nm, ".png")
        p <- netAnalysis_contribution(object.list[[i]], signaling = pw)
        ggsave(file.path(pw_dir, fname), p, width = 6, height = 4, dpi = 200)
      }, error = function(e) NULL)
    }
  }
  cat("    Contribution plots done\n")

  # --- Merged plots (if merge succeeded) ---
  if (!is.null(cellchat)) {
    cat("  [5] Merged comparison plots...\n")

    # rankNet
    tryCatch({
      png(file.path(out_dir, "09_rankNet_comparison.png"), width = 12, height = 8, units = "in", res = 300)
      rankNet(cellchat, mode = "comparison", stacked = TRUE, do.stat = TRUE,
              color.use = c("#4575b4", "#d73027"))
      dev.off()
      cat("    Saved: 09_rankNet_comparison.png\n")
    }, error = function(e) cat("    rankNet failed:", e$message, "\n"))

    tryCatch({
      png(file.path(out_dir, "10_rankNet_unstacked.png"), width = 14, height = 8, units = "in", res = 300)
      rankNet(cellchat, mode = "comparison", stacked = FALSE, do.stat = TRUE,
              color.use = c("#4575b4", "#d73027"))
      dev.off()
      cat("    Saved: 10_rankNet_unstacked.png\n")
    }, error = function(e) cat("    rankNet unstacked failed:", e$message, "\n"))

    # Bubble
    tryCatch({
      png(file.path(out_dir, "11_bubble_comparison.png"), width = 14, height = 10, units = "in", res = 300)
      netVisual_bubble(cellchat, comparison = c(1, 2), angle.x = 45, remove.isolate = TRUE)
      dev.off()
      cat("    Saved: 11_bubble_comparison.png\n")
    }, error = function(e) cat("    bubble failed:", e$message, "\n"))

    # Bubble per sender
    for (src in c("Mono", "DC", "NKc", "Tc", "Bc")) {
      tryCatch({
        fname <- paste0("12_bubble_from_", src, ".png")
        png(file.path(out_dir, fname), width = 12, height = 8, units = "in", res = 300)
        netVisual_bubble(cellchat, sources.use = src, comparison = c(1, 2),
                          angle.x = 45, remove.isolate = TRUE)
        dev.off()
        cat("    Saved:", fname, "\n")
      }, error = function(e) NULL)
    }

    # Differential signaling scatter
    tryCatch({
      png(file.path(out_dir, "14_diff_signalingRole_scatter.png"), width = 10, height = 8, units = "in", res = 300)
      netAnalysis_diff_signalingRole_scatter(cellchat)
      dev.off()
      cat("    Saved: 14_diff_signalingRole_scatter.png\n")
    }, error = function(e) cat("    diff scatter failed:", e$message, "\n"))
  } else {
    cat("  Merged plots skipped (merge failed)\n")
  }

  rm(cc_all, cc_ctrl, cc_dis, object.list); gc(verbose = FALSE)
  cat("  CellChat plots done for", layer_label, "\n")
}

# ---- Build condition map from sample metadata ----
# NOTE: CellChat sample dirs use make.names(patient_name)
#   HC: h58091733 etc. (already valid R names)
#   Stroke: numeric -> X10672798 etc. (make.names adds X prefix)
cat("Building condition maps...\n")

sobj <- qread(file.path(BASE, "5_1_hc_is.qs"))
meta <- sobj@meta.data

# L1: patient_name -> cohort, using make.names() for dir matching
sample_cohort <- meta %>% select(patient_name, cohort) %>% distinct()
sample_cohort$dir_name <- make.names(sample_cohort$patient_name)
cond_map_L1 <- setNames(sample_cohort$cohort, sample_cohort$dir_name)
cat("  L1 map:", sum(cond_map_L1 == "HC"), "HC,", sum(cond_map_L1 == "Stroke"), "Stroke\n")
rm(sobj); gc(verbose = FALSE)

# L2: g3 mapping
sobj2 <- qread(file.path(BASE, "5_2_is_g3.qs"))
meta2 <- sobj2@meta.data
sample_g3 <- meta2 %>% filter(!is.na(g3)) %>% select(patient_name, g3) %>% distinct()
sample_g3$dir_name <- make.names(sample_g3$patient_name)
cond_map_L2 <- setNames(ifelse(sample_g3$g3 == 1, "g3_Good", "g3_Bad"), sample_g3$dir_name)
cat("  L2 map:", sum(cond_map_L2 == "g3_Good"), "Good,", sum(cond_map_L2 == "g3_Bad"), "Bad\n")
rm(sobj2); gc(verbose = FALSE)

# ---- Run L1 ----
run_cellchat_native_plots(
  samples_dir = file.path(BASE, "cci/cellchat/L1_cohort_anno2/samples"),
  condition_map = cond_map_L1,
  control_name = "HC", disease_name = "Stroke",
  layer_label = "L1 (HC vs Stroke)",
  out_dir = file.path(BASE, "cci/plots/cellchat_L1_cohort")
)

# ---- Run L2 ----
run_cellchat_native_plots(
  samples_dir = file.path(BASE, "cci/cellchat/L2_g3_anno2/samples"),
  condition_map = cond_map_L2,
  control_name = "g3_Good", disease_name = "g3_Bad",
  layer_label = "L2 (g3: Good vs Bad)",
  out_dir = file.path(BASE, "cci/plots/cellchat_L2_g3")
)

# ==============================================================================
# PART 2: MNN anno1 Plots
# ==============================================================================

cat("\n", strrep("=", 60), "\n")
cat("  PART 2: MNN anno1 Plots\n")
cat(strrep("=", 60), "\n")

for (layer in c("L1_cohort", "L2_g3")) {
  mnn_file <- file.path(BASE, "cci/mnn", paste0(layer, "_anno1"), "multinichenet_results.qs")
  if (!file.exists(mnn_file)) {
    cat("  MNN file not found:", mnn_file, "\n")
    next
  }

  cat("\n  Loading MNN", layer, "anno1...\n")
  mnn <- qread(mnn_file)

  out_dir <- file.path(BASE, "cci/plots", paste0("mnn_", layer, "_anno1"))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # Check prioritization table
  prio <- mnn$prioritization_tables$group_prioritization_tbl
  cat("  Prioritization table:", if(is.null(prio)) "NULL" else nrow(prio), "rows\n")

  # Circos plots
  if (!is.null(prio) && nrow(prio) > 0) {
    groups <- unique(prio$group)
    cat("  Groups:", paste(groups, collapse = ", "), "\n")

    score_col <- if ("prioritization_score" %in% colnames(prio)) "prioritization_score" else {
      nc <- names(which(sapply(prio, is.numeric)))
      nc[length(nc)]
    }
    cat("  Score column:", score_col, "\n")

    for (grp in groups) {
      grp_data <- prio %>% filter(group == grp)
      cat("  Group", grp, ":", nrow(grp_data), "LR pairs\n")
      if (nrow(grp_data) == 0) next

      top_data <- grp_data %>% arrange(desc(!!rlang::sym(score_col))) %>% head(50)

      # Circos
      tryCatch({
        grp_safe <- gsub("[^a-zA-Z0-9_]", "_", grp)
        links_df <- top_data %>%
          select(ligand, receptor, !!rlang::sym(score_col)) %>%
          rename(from = ligand, to = receptor, weight = !!rlang::sym(score_col))

        senders <- unique(top_data$sender)
        n_col <- max(3, min(length(senders), 12))
        sender_pal <- brewer.pal(n_col, "Set3")[1:length(senders)]
        names(sender_pal) <- senders

        ligand_sender_map <- top_data %>% select(ligand, sender) %>% distinct() %>%
          group_by(ligand) %>% slice(1) %>% ungroup()

        grid_cols <- c()
        for (i in 1:nrow(ligand_sender_map)) {
          grid_cols[ligand_sender_map$ligand[i]] <- sender_pal[ligand_sender_map$sender[i]]
        }
        for (rec in unique(links_df$to)) grid_cols[rec] <- "grey70"

        pdf(file.path(out_dir, paste0("circos_", grp_safe, ".pdf")), width = 10, height = 10)
        circos.clear()
        chordDiagram(links_df, grid.col = grid_cols, transparency = 0.25,
                     directional = 1, direction.type = c("diffHeight", "arrows"),
                     link.arr.type = "big.arrow",
                     annotationTrack = "grid", preAllocateTracks = list(track.height = 0.1))
        circos.track(track.index = 1, panel.fun = function(x, y) {
          circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                      facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.7)
        }, bg.border = NA)
        title(paste("MNN Top 50 LR -", layer, grp))
        legend("topright", legend = names(sender_pal), fill = sender_pal, title = "Sender", bty = "n", cex = 0.8)
        dev.off()
        cat("    Circos saved\n")
      }, error = function(e) cat("    Circos failed:", e$message, "\n"))

      # LR bar plot
      tryCatch({
        top30 <- head(top_data, 30)
        p_lr <- ggplot(top30, aes(x = reorder(paste0(ligand, " -> ", receptor), !!rlang::sym(score_col)),
                                   y = !!rlang::sym(score_col), fill = sender)) +
          geom_col() + coord_flip() + scale_fill_brewer(palette = "Set3") +
          theme_minimal(base_size = 10) +
          labs(title = paste("Top 30 LR -", layer, grp), x = NULL, y = "Priority Score")
        ggsave(file.path(out_dir, paste0("top_lr_", grp_safe, ".png")), p_lr, width = 10, height = 8, dpi = 200)
        cat("    LR bar plot saved\n")
      }, error = function(e) cat("    LR bar failed:", e$message, "\n"))
    }

    write.csv(prio, file.path(out_dir, "lr_prioritization_full.csv"), row.names = FALSE)
  }

  # DE summary and volcanos
  de_tbl <- mnn$celltype_de
  if (!is.null(de_tbl) && is.data.frame(de_tbl) && nrow(de_tbl) > 0) {
    cat("  DE table:", nrow(de_tbl), "rows\n")

    p_adj_col <- intersect(c("p_adj", "p_adj_muscat", "p_val_adj"), colnames(de_tbl))[1]
    lfc_col <- intersect(c("logFC", "avg_log2FC", "log2FoldChange"), colnames(de_tbl))[1]

    if (!is.na(p_adj_col) && "celltype" %in% colnames(de_tbl)) {
      de_summary <- de_tbl %>%
        group_by(celltype) %>%
        summarise(n_genes = n(), n_sig = sum(!!rlang::sym(p_adj_col) < 0.05, na.rm = TRUE), .groups = "drop") %>%
        arrange(desc(n_sig))
      write.csv(de_summary, file.path(out_dir, "de_summary_by_celltype.csv"), row.names = FALSE)

      if (!is.na(lfc_col)) {
        top_cts <- head(de_summary$celltype[de_summary$n_sig > 0], 9)
        for (ct in top_cts) {
          ct_data <- de_tbl %>% filter(celltype == ct)
          p <- ggplot(ct_data, aes(x = !!rlang::sym(lfc_col), y = -log10(!!rlang::sym(p_adj_col)))) +
            geom_point(aes(color = !!rlang::sym(p_adj_col) < 0.05 & abs(!!rlang::sym(lfc_col)) > 0.25),
                       size = 0.5, alpha = 0.5) +
            scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey70"), guide = "none") +
            geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
            theme_minimal(base_size = 10) + labs(title = ct, subtitle = layer)
          ct_safe <- gsub("[/+ .]", "_", ct)
          ggsave(file.path(out_dir, paste0("volcano_", ct_safe, ".png")), p, width = 6, height = 5, dpi = 200)
        }
        cat("  Volcanos for", length(top_cts), "cell types\n")
      }
    }
  }

  # Ligand activities
  la_tbl <- mnn$ligand_activities_targets_DEgenes
  if (!is.null(la_tbl) && is.data.frame(la_tbl) && nrow(la_tbl) > 0 &&
      all(c("receiver", "activity") %in% colnames(la_tbl))) {
    cat("  Ligand activities:", nrow(la_tbl), "rows\n")
    receivers <- unique(la_tbl$receiver)
    for (recv in head(receivers, 6)) {
      recv_data <- la_tbl %>% filter(receiver == recv) %>% slice_max(activity, n = 15)
      if (nrow(recv_data) == 0) next
      p <- ggplot(recv_data, aes(x = reorder(ligand, activity), y = activity)) +
        geom_col(fill = "#d73027", alpha = 0.8) + coord_flip() +
        theme_minimal(base_size = 10) +
        labs(title = paste("Ligand Activity ->", recv), subtitle = layer, x = NULL)
      recv_safe <- gsub("[/+ .]", "_", recv)
      ggsave(file.path(out_dir, paste0("ligand_act_", recv_safe, ".png")), p, width = 7, height = 5, dpi = 200)
    }
    cat("  Ligand activity plots saved\n")
  }

  rm(mnn); gc(verbose = FALSE)
}

# ---- Summary ----
cat("\n", strrep("=", 60), "\n")
cat("CCI Visualization Complete\n")
cat("Finished:", format(Sys.time()), "\n\n")

for (d in c("cellchat_L1_cohort", "cellchat_L2_g3", "mnn_L1_cohort_anno1", "mnn_L2_g3_anno1")) {
  full_path <- file.path(BASE, "cci/plots", d)
  if (dir.exists(full_path)) {
    n_png <- length(list.files(full_path, pattern = "\\.png$", recursive = TRUE))
    n_pdf <- length(list.files(full_path, pattern = "\\.pdf$", recursive = TRUE))
    cat("  ", d, ":", n_png, "PNG +", n_pdf, "PDF\n")
  }
}

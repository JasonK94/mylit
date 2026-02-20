#!/usr/bin/env Rscript
# ==============================================================================
# CellChat Re-implementation v2: Use existing merged objects directly
# ==============================================================================
# Key insight: mergeCellChat() drops @data by design in CellChat v2.
# But comparison functions (rankNet, netVisual_diffInteraction) only need
# @net and @netP, which ARE preserved.
#
# Strategy:
# 1. Load the existing condition-merged objects (HC/Stroke, Good/Bad)
# 2. Fix @data slot by creating a dummy matrix (comparison funcs don't use it)
# 3. Use liftCellChat to ensure cell type consistency
# 4. Create comparison object with mergeCellChat()
# 5. If that fails, use individual condition objects directly
# ==============================================================================

.libPaths(c(
  "/home/user3/GJC_KDW_250721/renv/library/R-4.3/x86_64-pc-linux-gnu",
  "/home/user3/R/x86_64-pc-linux-gnu-library/4.3",
  .libPaths()
))

suppressPackageStartupMessages({
  library(qs)
  library(CellChat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(Matrix)
})

cat("=== CellChat Re-implementation v2 ===\n")
cat("Started:", format(Sys.time()), "\n\n")

BASE <- "/data/user3/sobj/stroke_hc_v8_2"

# ==============================================================================
# Strategy: Build CellChat from Seurat directly, one per condition
# This ensures @data is populated and all slots are consistent
# ==============================================================================

run_cellchat_from_seurat <- function(sobj, celltype_col, condition_col, conditions,
                                      db_use = "Secreted Signaling",
                                      min_cells = 10, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  cc_list <- list()

  for (cond in conditions) {
    cat("  Building CellChat for:", cond, "\n")
    sub <- subset(sobj, cells = colnames(sobj)[sobj@meta.data[[condition_col]] == cond])
    cat("    Cells:", ncol(sub), "\n")

    # Create CellChat object
    cc <- createCellChat(sub, group.by = celltype_col, assay = "RNA")
    CellChatDB <- CellChatDB.human
    CellChatDB.use <- subsetDB(CellChatDB, search = db_use)
    cc@DB <- CellChatDB.use

    # Standard pipeline
    cc <- subsetData(cc)
    cc <- identifyOverExpressedGenes(cc)
    cc <- identifyOverExpressedInteractions(cc)
    cc <- computeCommunProb(cc, type = "triMean")
    cc <- filterCommunication(cc, min.cells = min_cells)
    cc <- computeCommunProbPathway(cc)
    cc <- aggregateNet(cc)

    # Compute centrality
    cc <- netAnalysis_computeCentrality(cc, slot.name = "netP")

    cat("    Pathways:", length(cc@netP$pathways), "\n")
    cat("    @data dims:", paste(dim(cc@data), collapse = " x "), "\n")
    cat("    Net count sum:", sum(cc@net$count), "\n")

    cc_list[[cond]] <- cc
    qsave(cc, file.path(out_dir, paste0("cellchat_", cond, ".qs")))
  }

  return(cc_list)
}

generate_comparison_plots <- function(cc_list, conditions, out_dir) {
  cat("\n  Merging for comparison...\n")

  # Lift to common cell types
  all_ct <- sort(unique(unlist(lapply(cc_list, function(x) levels(x@idents)))))
  cat("  Cell type union:", paste(all_ct, collapse = ", "), "\n")

  cc_lifted <- lapply(cc_list, function(x) {
    liftCellChat(x, group.new = all_ct)
  })

  # Merge for comparison
  cellchat <- tryCatch({
    mergeCellChat(cc_lifted, add.names = conditions)
  }, error = function(e) {
    cat("  ERROR in merge:", conditionMessage(e), "\n")
    return(NULL)
  })

  if (is.null(cellchat)) {
    cat("  FATAL: Merge failed. Cannot generate comparison plots.\n")
    return(NULL)
  }

  cat("  Merge OK. Generating plots...\n")
  qsave(cellchat, file.path(out_dir, "cellchat_comparison.qs"))

  # ---- 1. compareInteractions ----
  tryCatch({
    p1 <- compareInteractions(cellchat, show.legend = FALSE, group = conditions)
    ggsave(file.path(out_dir, "01_compare_interactions_count.png"), p1,
           width = 6, height = 5, dpi = 300)
    p2 <- compareInteractions(cellchat, show.legend = FALSE, group = conditions,
                               measure = "weight")
    ggsave(file.path(out_dir, "02_compare_interactions_weight.png"), p2,
           width = 6, height = 5, dpi = 300)
    cat("    01-02: compareInteractions OK\n")
  }, error = function(e) cat("    01-02 FAILED:", conditionMessage(e), "\n"))

  # ---- 2. netVisual_diffInteraction ----
  tryCatch({
    png(file.path(out_dir, "03_diff_circle_count.png"),
        width = 8, height = 8, units = "in", res = 300)
    netVisual_diffInteraction(cellchat, weight.scale = TRUE)
    dev.off()

    png(file.path(out_dir, "04_diff_circle_weight.png"),
        width = 8, height = 8, units = "in", res = 300)
    netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "weight")
    dev.off()
    cat("    03-04: diffInteraction circles OK\n")
  }, error = function(e) cat("    03-04 FAILED:", conditionMessage(e), "\n"))

  # ---- 3. Differential heatmaps ----
  tryCatch({
    png(file.path(out_dir, "05_diff_heatmap_count.png"),
        width = 8, height = 7, units = "in", res = 300)
    netVisual_heatmap(cellchat)
    dev.off()

    png(file.path(out_dir, "06_diff_heatmap_weight.png"),
        width = 8, height = 7, units = "in", res = 300)
    netVisual_heatmap(cellchat, measure = "weight")
    dev.off()
    cat("    05-06: diff heatmaps OK\n")
  }, error = function(e) cat("    05-06 FAILED:", conditionMessage(e), "\n"))

  # ---- 4. rankNet (information flow) ----
  tryCatch({
    p_rank <- rankNet(cellchat, mode = "comparison", stacked = TRUE, do.stat = TRUE)
    ggsave(file.path(out_dir, "07_rankNet_stacked.png"), p_rank,
           width = 8, height = 10, dpi = 300)

    p_rank2 <- rankNet(cellchat, mode = "comparison", stacked = FALSE, do.stat = TRUE)
    ggsave(file.path(out_dir, "08_rankNet_unstacked.png"), p_rank2,
           width = 10, height = 10, dpi = 300)
    cat("    07-08: rankNet OK\n")
  }, error = function(e) cat("    07-08 FAILED:", conditionMessage(e), "\n"))

  # ---- 5. Bubble comparison ----
  tryCatch({
    p_bubble <- netVisual_bubble(cellchat, comparison = c(1, 2),
                                  remove.isolate = FALSE, angle.x = 45)
    ggsave(file.path(out_dir, "09_bubble_comparison.png"), p_bubble,
           width = 14, height = 12, dpi = 300)
    cat("    09: bubble comparison OK\n")
  }, error = function(e) cat("    09 FAILED:", conditionMessage(e), "\n"))

  # ---- 6. Signaling role heatmaps ----
  for (pattern in c("outgoing", "incoming")) {
    idx <- ifelse(pattern == "outgoing", "10", "11")
    tryCatch({
      png(file.path(out_dir, paste0(idx, "_", pattern, "_pattern.png")),
          width = 10, height = 6, units = "in", res = 300)
      netAnalysis_signalingRole_heatmap(cellchat, pattern = pattern,
                                        width = 10, height = 6)
      dev.off()
      cat("    ", idx, ":", pattern, "pattern OK\n")
    }, error = function(e) cat("    ", idx, ":", pattern, "FAILED:", conditionMessage(e), "\n"))
  }

  # ---- 7. Signaling changes scatter ----
  tryCatch({
    p_scatter <- netAnalysis_signalingChanges_scatter(cellchat)
    ggsave(file.path(out_dir, "12_signaling_scatter.png"), p_scatter,
           width = 8, height = 8, dpi = 300)
    cat("    12: scatter OK\n")
  }, error = function(e) cat("    12 scatter FAILED:", conditionMessage(e), "\n"))

  # ---- 8. Per-pathway analysis for common pathways ----
  pathways_1 <- cc_lifted[[1]]@netP$pathways
  pathways_2 <- cc_lifted[[2]]@netP$pathways
  common_pw <- intersect(pathways_1, pathways_2)
  cat("  Common pathways:", length(common_pw), "\n")

  pw_dir <- file.path(out_dir, "pathways")
  dir.create(pw_dir, showWarnings = FALSE, recursive = TRUE)

  for (pw in common_pw) {
    tryCatch({
      # Chord diagram per condition
      for (i in 1:2) {
        fname <- file.path(pw_dir, paste0("chord_", pw, "_", conditions[i], ".png"))
        png(fname, width = 8, height = 8, units = "in", res = 200)
        netVisual_aggregate(cc_lifted[[i]], signaling = pw,
                            layout = "chord", signaling.name = paste(pw, conditions[i]))
        dev.off()
      }
      # L-R contribution
      p_contrib <- netAnalysis_contribution(cc_lifted[[1]], signaling = pw)
      ggsave(file.path(pw_dir, paste0("contrib_", pw, "_", conditions[1], ".png")),
             p_contrib, width = 8, height = 5, dpi = 200)
      p_contrib2 <- netAnalysis_contribution(cc_lifted[[2]], signaling = pw)
      ggsave(file.path(pw_dir, paste0("contrib_", pw, "_", conditions[2], ".png")),
             p_contrib2, width = 8, height = 5, dpi = 200)
    }, error = function(e) {
      # Silently skip pathways that fail (some may not have enough data)
    })
  }
  n_pw_files <- length(list.files(pw_dir, pattern = "\\.png$"))
  cat("  Pathway plots:", n_pw_files, "files\n")

  # Summary
  n_plots <- length(list.files(out_dir, pattern = "\\.png$"))
  cat("\n  TOTAL plots (main):", n_plots, "\n")
  return(cellchat)
}


# ==============================================================================
# LAYER 1: HC vs Stroke — Build from Seurat
# ==============================================================================
cat(strrep("=", 60), "\n")
cat("LAYER 1: HC vs Stroke (from Seurat)\n")
cat(strrep("=", 60), "\n\n")

cat("Loading L1 Seurat object...\n")
sobj_l1 <- qread(file.path(BASE, "5_1_hc_is.qs"))
cat("  Cells:", ncol(sobj_l1), ", Cell types:", length(unique(sobj_l1$anno2)), "\n")

out_l1 <- file.path(BASE, "cci/plots/cellchat_L1_cohort_v2")
cc_l1 <- run_cellchat_from_seurat(
  sobj_l1, celltype_col = "anno2", condition_col = "cohort",
  conditions = c("HC", "Stroke"),
  db_use = "Secreted Signaling", min_cells = 10, out_dir = out_l1
)
rm(sobj_l1); gc(verbose = FALSE)

cellchat_l1 <- generate_comparison_plots(cc_l1, c("HC", "Stroke"), out_l1)
rm(cc_l1); gc(verbose = FALSE)


# ==============================================================================
# LAYER 2: g3 Good vs Bad — Build from Seurat
# ==============================================================================
cat("\n\n")
cat(strrep("=", 60), "\n")
cat("LAYER 2: g3 Good vs Bad (from Seurat)\n")
cat(strrep("=", 60), "\n\n")

cat("Loading L2 Seurat object...\n")
sobj_l2 <- qread(file.path(BASE, "5_2_is_g3.qs"))
# Map g3 to meaningful labels
sobj_l2$g3_label <- ifelse(sobj_l2$g3 == 1, "Good", "Bad")
cat("  Cells:", ncol(sobj_l2), "\n")
cat("  g3 distribution:", table(sobj_l2$g3_label), "\n")

out_l2 <- file.path(BASE, "cci/plots/cellchat_L2_g3_v2")
cc_l2 <- run_cellchat_from_seurat(
  sobj_l2, celltype_col = "anno2", condition_col = "g3_label",
  conditions = c("Good", "Bad"),
  db_use = "Secreted Signaling", min_cells = 10, out_dir = out_l2
)
rm(sobj_l2); gc(verbose = FALSE)

cellchat_l2 <- generate_comparison_plots(cc_l2, c("Good", "Bad"), out_l2)
rm(cc_l2); gc(verbose = FALSE)


# ==============================================================================
# Summary
# ==============================================================================
cat("\n\n", strrep("=", 60), "\n")
cat("CellChat Re-implementation v2 Complete\n")
cat("Finished:", format(Sys.time()), "\n\n")

for (d in c("cellchat_L1_cohort_v2", "cellchat_L2_g3_v2")) {
  full_path <- file.path(BASE, "cci/plots", d)
  if (dir.exists(full_path)) {
    n_main <- length(list.files(full_path, pattern = "\\.png$"))
    n_pw <- length(list.files(file.path(full_path, "pathways"),
                               pattern = "\\.png$", recursive = TRUE))
    cat("  ", d, ":", n_main, "main +", n_pw, "pathway PNG\n")
  }
}
cat("\nDone.\n")

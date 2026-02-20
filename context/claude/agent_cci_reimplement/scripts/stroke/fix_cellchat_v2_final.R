#!/usr/bin/env Rscript
# ==============================================================================
# Fix remaining CellChat v2 plots: 05-06 heatmaps, 10-11 centrality, 12 scatter
# ==============================================================================
# Strategy for centrality plots:
#   - Use INDIVIDUAL condition objects (which have centrality computed)
#   - Generate side-by-side plots using patchwork
#   - For scatter: use individual objects, not merged
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
  library(ggalluvial)
  library(ComplexHeatmap)
})

BASE <- "/data/user3/sobj/stroke_hc_v8_2"

configs <- list(
  list(
    name = "cellchat_L1_cohort_v2",
    conditions = c("HC", "Stroke"),
    labels = c("HC", "Stroke")
  ),
  list(
    name = "cellchat_L2_g3_v2",
    conditions = c("Good", "Bad"),
    labels = c("Good (g3=1)", "Bad (g3=2)")
  )
)

for (cfg in configs) {
  cat("\n", strrep("=", 60), "\n")
  cat("Fixing:", cfg$name, "\n")
  cat(strrep("=", 60), "\n")

  out_dir <- file.path(BASE, "cci/plots", cfg$name)

  # Load individual objects (with centrality)
  cc1 <- qread(file.path(out_dir, paste0("cellchat_", cfg$conditions[1], ".qs")))
  cc2 <- qread(file.path(out_dir, paste0("cellchat_", cfg$conditions[2], ".qs")))

  cat("  ", cfg$conditions[1], ": @data", paste(dim(cc1@data), collapse="x"),
      ", pathways:", length(cc1@netP$pathways), "\n")
  cat("  ", cfg$conditions[2], ": @data", paste(dim(cc2@data), collapse="x"),
      ", pathways:", length(cc2@netP$pathways), "\n")

  # Also load merged object for plots that work on it
  cc_merged <- qread(file.path(out_dir, "cellchat_comparison.qs"))

  # ---- 05-06: Differential heatmaps ----
  # These should work on merged objects. May have failed due to ComplexHeatmap not loaded.
  for (measure_type in c("count", "weight")) {
    idx <- ifelse(measure_type == "count", "05", "06")
    fname <- file.path(out_dir, paste0(idx, "_diff_heatmap_", measure_type, ".png"))
    if (!file.exists(fname)) {
      tryCatch({
        png(fname, width = 8, height = 7, units = "in", res = 300)
        if (measure_type == "count") {
          netVisual_heatmap(cc_merged)
        } else {
          netVisual_heatmap(cc_merged, measure = "weight")
        }
        dev.off()
        cat("    ", idx, ": diff heatmap", measure_type, "OK\n")
      }, error = function(e) {
        try(dev.off(), silent = TRUE)
        cat("    ", idx, ": diff heatmap", measure_type, "FAILED:", conditionMessage(e), "\n")
      })
    } else {
      cat("    ", idx, ": already exists\n")
    }
  }

  # ---- 10-11: Signaling role heatmap (per condition, side by side) ----
  # Must use INDIVIDUAL objects because merged object loses centrality
  for (pattern in c("outgoing", "incoming")) {
    idx <- ifelse(pattern == "outgoing", "10", "11")
    fname <- file.path(out_dir, paste0(idx, "_", pattern, "_signaling_role.png"))

    tryCatch({
      png(fname, width = 16, height = 7, units = "in", res = 300)
      # Draw side by side using CellChat's built-in comparison mode
      # First approach: pass a list of two objects
      ht1 <- netAnalysis_signalingRole_heatmap(cc1, pattern = pattern,
                                                 title = paste(cfg$labels[1], "-", pattern),
                                                 width = 7, height = 5)
      ht2 <- netAnalysis_signalingRole_heatmap(cc2, pattern = pattern,
                                                 title = paste(cfg$labels[2], "-", pattern),
                                                 width = 7, height = 5)
      draw(ht1 + ht2, ht_gap = unit(1, "cm"))
      dev.off()
      cat("    ", idx, ":", pattern, "signaling role OK (side-by-side)\n")
    }, error = function(e) {
      try(dev.off(), silent = TRUE)
      # Fallback: draw individually
      tryCatch({
        for (i in 1:2) {
          cc_i <- if (i == 1) cc1 else cc2
          fname_i <- file.path(out_dir, paste0(idx, "_", pattern, "_",
                                                 cfg$conditions[i], ".png"))
          png(fname_i, width = 8, height = 6, units = "in", res = 300)
          netAnalysis_signalingRole_heatmap(cc_i, pattern = pattern,
                                             title = paste(cfg$labels[i], "-", pattern),
                                             width = 8, height = 5)
          dev.off()
        }
        cat("    ", idx, ":", pattern, "signaling role OK (individual files)\n")
      }, error = function(e2) {
        try(dev.off(), silent = TRUE)
        cat("    ", idx, ":", pattern, "FAILED completely:", conditionMessage(e2), "\n")
      })
    })
  }

  # ---- 12: Signaling changes scatter (functional similarity) ----
  # netAnalysis_signalingChanges_scatter needs the merged object with centrality
  # Since merged object loses centrality, we'll use a different approach:
  # Compute functional similarity between conditions using individual objects
  tryCatch({
    # Approach: use rankSimilarity or compute outgoing/incoming pattern overlap
    # CellChat's scatter actually shows cell-type-level signaling changes
    # Let's try the merged object but manually add centrality back

    # Re-build merged with centrality: add centrality to netP manually
    # The centrality is stored in cc@netP$centr
    if (!is.null(cc1@netP$centr) && !is.null(cc2@netP$centr)) {
      # Lift individual objects to common cell types
      all_ct <- sort(unique(c(levels(cc1@idents), levels(cc2@idents))))
      cc1_l <- liftCellChat(cc1, group.new = all_ct)
      cc2_l <- liftCellChat(cc2, group.new = all_ct)

      # Re-compute centrality after lift
      cc1_l <- netAnalysis_computeCentrality(cc1_l, slot.name = "netP")
      cc2_l <- netAnalysis_computeCentrality(cc2_l, slot.name = "netP")

      cc_m2 <- mergeCellChat(list(cc1_l, cc2_l), add.names = cfg$conditions)

      p_scatter <- netAnalysis_signalingChanges_scatter(cc_m2, idents.use = all_ct[1])
      ggsave(file.path(out_dir, "12_signaling_scatter.png"), p_scatter,
             width = 8, height = 8, dpi = 300)
      cat("    12: scatter OK (first cell type:", all_ct[1], ")\n")

      # Also generate scatter for each cell type
      scatter_dir <- file.path(out_dir, "scatter")
      dir.create(scatter_dir, showWarnings = FALSE)
      for (ct in all_ct) {
        tryCatch({
          p <- netAnalysis_signalingChanges_scatter(cc_m2, idents.use = ct)
          ct_safe <- gsub("/", "_", ct)
          ggsave(file.path(scatter_dir, paste0("scatter_", ct_safe, ".png")),
                 p, width = 8, height = 8, dpi = 200)
        }, error = function(e) {})
      }
      n_scatter <- length(list.files(scatter_dir, pattern = "\\.png$"))
      cat("    12: generated", n_scatter, "per-celltype scatter plots\n")

      rm(cc_m2, cc1_l, cc2_l)
    } else {
      cat("    12: centrality not available on individual objects\n")
    }
  }, error = function(e) {
    cat("    12 FAILED:", conditionMessage(e), "\n")
  })

  # ---- 13: Information flow comparison (additional) ----
  # This uses the merged object and should already be covered by rankNet (07-08)
  # But let's also add: signaling role scatter (outgoing vs incoming per cell type)
  tryCatch({
    png(file.path(out_dir, "13_signaling_role_scatter.png"),
        width = 12, height = 6, units = "in", res = 300)
    par(mfrow = c(1, 2))
    netAnalysis_signalingRole_scatter(cc1, title = cfg$labels[1])
    netAnalysis_signalingRole_scatter(cc2, title = cfg$labels[2])
    par(mfrow = c(1, 1))
    dev.off()
    cat("    13: signaling role scatter OK (side-by-side)\n")
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    # Fallback: individual files
    tryCatch({
      for (i in 1:2) {
        cc_i <- if (i == 1) cc1 else cc2
        png(file.path(out_dir, paste0("13_signaling_role_scatter_", cfg$conditions[i], ".png")),
            width = 8, height = 8, units = "in", res = 300)
        netAnalysis_signalingRole_scatter(cc_i, title = cfg$labels[i])
        dev.off()
      }
      cat("    13: signaling role scatter OK (individual)\n")
    }, error = function(e2) {
      try(dev.off(), silent = TRUE)
      cat("    13 FAILED:", conditionMessage(e2), "\n")
    })
  })

  # ---- 14-15: River/alluvial plots ----
  for (i in 1:2) {
    cc_i <- if (i == 1) cc1 else cc2
    cond <- cfg$conditions[i]

    for (pattern in c("outgoing", "incoming")) {
      idx <- ifelse(pattern == "outgoing", "14", "15")
      fname <- file.path(out_dir, paste0(idx, "_river_", pattern, "_", cond, ".png"))
      tryCatch({
        # identifyGlobalCommunicationPatterns first
        nPat <- min(3, length(cc_i@netP$pathways))
        if (nPat >= 2) {
          cc_i <- identifyCommunicationPatterns(cc_i, pattern = pattern, k = nPat)
          png(fname, width = 10, height = 6, units = "in", res = 200)
          netAnalysis_river(cc_i, pattern = pattern)
          dev.off()
        }
      }, error = function(e) {
        try(dev.off(), silent = TRUE)
      })
    }
  }
  n_river <- length(list.files(out_dir, pattern = "river"))
  cat("    14-15: river plots:", n_river, "generated\n")

  rm(cc1, cc2, cc_merged); gc(verbose = FALSE)
}

# Final summary
cat("\n", strrep("=", 60), "\n")
cat("Final plot counts:\n")
for (d in c("cellchat_L1_cohort_v2", "cellchat_L2_g3_v2")) {
  full_path <- file.path(BASE, "cci/plots", d)
  if (dir.exists(full_path)) {
    n_main <- length(list.files(full_path, pattern = "\\.png$"))
    n_pw <- length(list.files(file.path(full_path, "pathways"),
                               pattern = "\\.png$", recursive = TRUE))
    n_scatter <- length(list.files(file.path(full_path, "scatter"),
                                    pattern = "\\.png$", recursive = TRUE))
    cat("  ", d, ":", n_main, "main +", n_pw, "pathway +", n_scatter, "scatter PNG\n")
  }
}
cat("\nDone:", format(Sys.time()), "\n")

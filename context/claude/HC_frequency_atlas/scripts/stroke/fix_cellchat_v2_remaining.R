#!/usr/bin/env Rscript
# Fix remaining CellChat v2 plots: signaling patterns + scatter

.libPaths(c(
  "/home/user3/GJC_KDW_250721/renv/library/R-4.3/x86_64-pc-linux-gnu",
  "/home/user3/R/x86_64-pc-linux-gnu-library/4.3",
  .libPaths()
))

suppressPackageStartupMessages({
  library(qs)
  library(CellChat)
  library(ggplot2)
})

BASE <- "/data/user3/sobj/stroke_hc_v8_2"

for (layer in c("cellchat_L1_cohort_v2", "cellchat_L2_g3_v2")) {
  cat("\n=== Fixing:", layer, "===\n")
  out_dir <- file.path(BASE, "cci/plots", layer)
  cc_file <- file.path(out_dir, "cellchat_comparison.qs")
  if (!file.exists(cc_file)) {
    cat("  No comparison object, skipping\n")
    next
  }

  cellchat <- qread(cc_file)

  # Compute centrality for each dataset in the list
  cat("  Computing centrality...\n")
  for (i in 1:length(cellchat@netP)) {
    tryCatch({
      # Need to work on individual objects
      # netP is a list of lists, each with $pathways, $prob, $pval, etc.
    }, error = function(e) {
      cat("    Centrality calc error:", conditionMessage(e), "\n")
    })
  }

  # Load the individual condition objects and compute centrality
  conditions <- names(cellchat@netP)
  if (is.null(conditions)) conditions <- paste0("Dataset", 1:length(cellchat@netP))
  cat("  Conditions:", paste(conditions, collapse = ", "), "\n")

  cc_objects <- list()
  for (cond in conditions) {
    cc_file_cond <- file.path(out_dir, paste0("cellchat_", cond, ".qs"))
    if (file.exists(cc_file_cond)) {
      cc <- qread(cc_file_cond)
      cc <- netAnalysis_computeCentrality(cc, slot.name = "netP")
      cc_objects[[cond]] <- cc
      qsave(cc, cc_file_cond)
      cat("    ", cond, ": centrality computed\n")
    }
  }

  if (length(cc_objects) == 2) {
    # Re-merge with centrality
    all_ct <- sort(unique(unlist(lapply(cc_objects, function(x) levels(x@idents)))))
    cc_lifted <- lapply(cc_objects, function(x) liftCellChat(x, group.new = all_ct))
    cellchat2 <- mergeCellChat(cc_lifted, add.names = conditions)
    qsave(cellchat2, file.path(out_dir, "cellchat_comparison.qs"))

    # ---- 10-11: Signaling role heatmap ----
    for (pattern in c("outgoing", "incoming")) {
      idx <- ifelse(pattern == "outgoing", "10", "11")
      tryCatch({
        png(file.path(out_dir, paste0(idx, "_", pattern, "_pattern.png")),
            width = 10, height = 6, units = "in", res = 300)
        netAnalysis_signalingRole_heatmap(cellchat2, pattern = pattern,
                                          width = 10, height = 6)
        dev.off()
        cat("    ", idx, ": ", pattern, " pattern OK\n")
      }, error = function(e) cat("    ", idx, " FAILED:", conditionMessage(e), "\n"))
    }

    # ---- 12: Signaling changes scatter ----
    tryCatch({
      idents_use <- levels(cellchat2@idents[[1]])
      p_scatter <- netAnalysis_signalingChanges_scatter(cellchat2,
                                                        idents.use = idents_use)
      ggsave(file.path(out_dir, "12_signaling_scatter.png"), p_scatter,
             width = 8, height = 8, dpi = 300)
      cat("    12: scatter OK\n")
    }, error = function(e) cat("    12 scatter FAILED:", conditionMessage(e), "\n"))

    # ---- 13: River plot for top pathways ----
    tryCatch({
      # Find differentially expressed pathways
      pathways_1 <- cc_objects[[1]]@netP$pathways
      pathways_2 <- cc_objects[[2]]@netP$pathways
      common_pw <- intersect(pathways_1, pathways_2)

      if (length(common_pw) > 0) {
        for (cond_idx in 1:2) {
          cond <- conditions[cond_idx]
          tryCatch({
            png(file.path(out_dir, paste0("14_river_", cond, ".png")),
                width = 10, height = 6, units = "in", res = 200)
            netAnalysis_river(cc_objects[[cond_idx]], pattern = "outgoing")
            dev.off()
          }, error = function(e) {})

          tryCatch({
            png(file.path(out_dir, paste0("15_river_incoming_", cond, ".png")),
                width = 10, height = 6, units = "in", res = 200)
            netAnalysis_river(cc_objects[[cond_idx]], pattern = "incoming")
            dev.off()
          }, error = function(e) {})
        }
        cat("    14-15: river plots attempted\n")
      }
    }, error = function(e) cat("    14-15 river FAILED:", conditionMessage(e), "\n"))

    rm(cellchat2)
  }

  rm(cellchat, cc_objects); gc(verbose = FALSE)
}

# Count final results
cat("\n=== Final plot counts ===\n")
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

#!/usr/bin/env Rscript
# ==============================================================================
# CellChat Re-implementation: Proper Comparison Pipeline
# ==============================================================================
# Problem: Original merged objects have empty @data → mergeCellChat() fails
#          → Standard comparison functions (rankNet, diffInteraction, bubble) unavailable
#
# Solution: Load per-sample objects → liftCellChat() for cell type union →
#           merge per-condition → mergeCellChat for comparison
#
# Layers: L1 (HC vs Stroke), L2 (g3 Good vs Bad)
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
})

cat("=== CellChat Re-implementation ===\n")
cat("Started:", format(Sys.time()), "\n\n")

BASE <- "/data/user3/sobj/stroke_hc_v8_2"

# ==============================================================================
# Helper: Load and prepare per-sample CellChat objects
# ==============================================================================
load_sample_objects <- function(sample_dir) {
  dirs <- list.dirs(sample_dir, recursive = FALSE, full.names = TRUE)
  cat("  Loading", length(dirs), "sample objects...\n")

  cc_list <- list()
  failed <- c()

  for (d in dirs) {
    cc_file <- file.path(d, "cellchat.qs")
    if (!file.exists(cc_file)) next

    sample_name <- basename(d)
    tryCatch({
      cc <- qread(cc_file)
      # CellChat v2 compatibility
      cc <- updateCellChat(cc)
      cc_list[[sample_name]] <- cc
    }, error = function(e) {
      cat("    WARN: Failed to load", sample_name, ":", conditionMessage(e), "\n")
      failed <<- c(failed, sample_name)
    })
  }

  cat("  Loaded:", length(cc_list), "objects")
  if (length(failed) > 0) cat(" (", length(failed), "failed)")
  cat("\n")
  return(cc_list)
}

# ==============================================================================
# Helper: Get union of cell types and lift all objects
# ==============================================================================
lift_all <- function(cc_list) {
  all_ct <- sort(unique(unlist(lapply(cc_list, function(x) levels(x@idents)))))
  cat("  Cell type union:", paste(all_ct, collapse = ", "), "\n")

  lifted <- lapply(names(cc_list), function(nm) {
    tryCatch({
      liftCellChat(cc_list[[nm]], group.new = all_ct)
    }, error = function(e) {
      cat("    WARN: liftCellChat failed for", nm, ":", conditionMessage(e), "\n")
      return(NULL)
    })
  })
  names(lifted) <- names(cc_list)
  lifted <- lifted[!sapply(lifted, is.null)]
  cat("  Lifted:", length(lifted), "objects to", length(all_ct), "cell types\n")
  return(lifted)
}

# ==============================================================================
# Helper: Merge per-condition and create comparison object
# ==============================================================================
merge_and_compare <- function(cc_list, condition_map, layer_name, out_dir) {
  cat("\n", strrep("=", 60), "\n")
  cat("  Processing:", layer_name, "\n")
  cat(strrep("=", 60), "\n")

  # Split by condition
  conditions <- unique(condition_map$condition)
  cat("  Conditions:", paste(conditions, collapse = " vs "), "\n")

  cc_by_cond <- list()
  for (cond in conditions) {
    cond_samples <- condition_map$sample[condition_map$condition == cond]
    cond_cc <- cc_list[names(cc_list) %in% cond_samples]
    cat("  ", cond, ":", length(cond_cc), "samples\n")

    if (length(cond_cc) == 0) {
      cat("    ERROR: No samples for condition", cond, "\n")
      next
    }

    # Merge samples within condition
    tryCatch({
      if (length(cond_cc) == 1) {
        cc_by_cond[[cond]] <- cond_cc[[1]]
      } else {
        cc_by_cond[[cond]] <- mergeCellChat(
          cond_cc,
          add.names = names(cond_cc),
          cell.prefix = TRUE
        )
      }
      cat("    Merged OK:", ncol(cc_by_cond[[cond]]@data), "cells in @data\n")
    }, error = function(e) {
      cat("    ERROR merging", cond, ":", conditionMessage(e), "\n")
    })
  }

  if (length(cc_by_cond) < 2) {
    cat("  FATAL: Need 2 conditions for comparison, got", length(cc_by_cond), "\n")
    return(NULL)
  }

  # Create comparison object
  cat("\n  Creating comparison object...\n")
  cellchat <- tryCatch({
    mergeCellChat(cc_by_cond, add.names = conditions)
  }, error = function(e) {
    cat("    ERROR in final merge:", conditionMessage(e), "\n")
    return(NULL)
  })

  if (is.null(cellchat)) return(NULL)
  cat("  Comparison object created OK\n")

  # Save the properly merged object
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  qsave(cellchat, file.path(out_dir, "cellchat_comparison.qs"))
  cat("  Saved comparison object\n")

  # ==== Generate standard comparison plots ====
  cat("\n  Generating comparison plots...\n")

  # ---- 1. compareInteractions: count & weight ----
  tryCatch({
    p1 <- compareInteractions(cellchat, show.legend = FALSE, group = conditions)
    ggsave(file.path(out_dir, "01_compare_interactions_count.png"), p1,
           width = 6, height = 5, dpi = 300)
    p2 <- compareInteractions(cellchat, show.legend = FALSE, group = conditions,
                               measure = "weight")
    ggsave(file.path(out_dir, "02_compare_interactions_weight.png"), p2,
           width = 6, height = 5, dpi = 300)
    cat("    01-02: compareInteractions OK\n")
  }, error = function(e) cat("    01-02: compareInteractions FAILED:", conditionMessage(e), "\n"))

  # ---- 2. netVisual_diffInteraction: circle + heatmap ----
  tryCatch({
    png(file.path(out_dir, "03_diff_interaction_circle_count.png"),
        width = 8, height = 8, units = "in", res = 300)
    netVisual_diffInteraction(cellchat, weight.scale = TRUE)
    dev.off()

    png(file.path(out_dir, "04_diff_interaction_circle_weight.png"),
        width = 8, height = 8, units = "in", res = 300)
    netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "weight")
    dev.off()

    png(file.path(out_dir, "05_diff_interaction_heatmap_count.png"),
        width = 8, height = 8, units = "in", res = 300)
    netVisual_heatmap(cellchat)
    dev.off()

    png(file.path(out_dir, "06_diff_interaction_heatmap_weight.png"),
        width = 8, height = 8, units = "in", res = 300)
    netVisual_heatmap(cellchat, measure = "weight")
    dev.off()
    cat("    03-06: diffInteraction circle+heatmap OK\n")
  }, error = function(e) cat("    03-06: diffInteraction FAILED:", conditionMessage(e), "\n"))

  # ---- 3. rankNet: information flow ----
  tryCatch({
    p_rank <- rankNet(cellchat, mode = "comparison", stacked = TRUE,
                       do.stat = TRUE)
    ggsave(file.path(out_dir, "07_rankNet_stacked.png"), p_rank,
           width = 8, height = 10, dpi = 300)

    p_rank2 <- rankNet(cellchat, mode = "comparison", stacked = FALSE,
                        do.stat = TRUE)
    ggsave(file.path(out_dir, "08_rankNet_unstacked.png"), p_rank2,
           width = 10, height = 10, dpi = 300)
    cat("    07-08: rankNet OK\n")
  }, error = function(e) cat("    07-08: rankNet FAILED:", conditionMessage(e), "\n"))

  # ---- 4. Bubble plot ----
  tryCatch({
    p_bubble <- netVisual_bubble(cellchat, comparison = c(1, 2),
                                  remove.isolate = FALSE,
                                  angle.x = 45)
    ggsave(file.path(out_dir, "09_bubble_comparison.png"), p_bubble,
           width = 14, height = 12, dpi = 300)
    cat("    09: bubble comparison OK\n")
  }, error = function(e) cat("    09: bubble FAILED:", conditionMessage(e), "\n"))

  # ---- 5. Signaling role heatmap ----
  tryCatch({
    png(file.path(out_dir, "10_signaling_role_heatmap.png"),
        width = 12, height = 8, units = "in", res = 300)
    netAnalysis_signalingRole_heatmap(cellchat, pattern = "all", signaling = NULL,
                                      width = 12, height = 8)
    dev.off()
    cat("    10: signaling role heatmap OK\n")
  }, error = function(e) cat("    10: signaling role heatmap FAILED:", conditionMessage(e), "\n"))

  # ---- 6. Signaling changes scatter ----
  tryCatch({
    # Get common pathways
    pathways_1 <- cellchat@netP[[1]]$pathways
    pathways_2 <- cellchat@netP[[2]]$pathways
    common_pathways <- intersect(pathways_1, pathways_2)

    if (length(common_pathways) > 0) {
      # Overall signaling
      p_scatter <- netAnalysis_signalingChanges_scatter(cellchat,
                                                        idents.use = NULL)
      ggsave(file.path(out_dir, "11_signaling_scatter.png"), p_scatter,
             width = 8, height = 8, dpi = 300)
      cat("    11: signaling scatter OK\n")
    } else {
      cat("    11: No common pathways, scatter skipped\n")
    }
  }, error = function(e) cat("    11: signaling scatter FAILED:", conditionMessage(e), "\n"))

  # ---- 7. Per-cell-type outgoing/incoming ----
  tryCatch({
    ct_names <- levels(cellchat@idents[[1]])
    n_ct <- length(ct_names)

    if (n_ct > 0) {
      # Outgoing patterns
      png(file.path(out_dir, "12_outgoing_patterns.png"),
          width = 10, height = 6, units = "in", res = 300)
      netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",
                                        width = 10, height = 6)
      dev.off()

      # Incoming patterns
      png(file.path(out_dir, "13_incoming_patterns.png"),
          width = 10, height = 6, units = "in", res = 300)
      netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",
                                        width = 10, height = 6)
      dev.off()
      cat("    12-13: outgoing/incoming patterns OK\n")
    }
  }, error = function(e) cat("    12-13: patterns FAILED:", conditionMessage(e), "\n"))

  # ---- Summary ----
  n_plots <- length(list.files(out_dir, pattern = "\\.png$"))
  cat("\n  Total plots generated:", n_plots, "\n")
  cat("  Output directory:", out_dir, "\n")

  return(cellchat)
}


# ==============================================================================
# LAYER 1: HC vs Stroke
# ==============================================================================
cat(strrep("=", 60), "\n")
cat("LAYER 1: HC vs Stroke\n")
cat(strrep("=", 60), "\n\n")

# Load per-sample objects
sample_dir_l1 <- file.path(BASE, "cci/cellchat/L1_cohort_anno2/samples")
cc_list_l1 <- load_sample_objects(sample_dir_l1)

# Build condition map from directory names
# HC samples: h-prefix, Stroke: X-prefix (from make.names on numeric patient IDs)
l1_map <- data.frame(
  sample = names(cc_list_l1),
  condition = ifelse(grepl("^h", names(cc_list_l1)), "HC", "Stroke"),
  stringsAsFactors = FALSE
)
cat("  L1 mapping: HC =", sum(l1_map$condition == "HC"),
    ", Stroke =", sum(l1_map$condition == "Stroke"), "\n")

# Lift all objects to common cell type space
cc_lifted_l1 <- lift_all(cc_list_l1)

# Merge and compare
out_l1 <- file.path(BASE, "cci/plots/cellchat_L1_cohort_v2")
cellchat_l1 <- merge_and_compare(cc_lifted_l1, l1_map, "L1: HC vs Stroke", out_l1)

rm(cc_list_l1, cc_lifted_l1); gc(verbose = FALSE)

# ==============================================================================
# LAYER 2: g3 Good vs Bad
# ==============================================================================
cat("\n\n")
cat(strrep("=", 60), "\n")
cat("LAYER 2: g3 Good vs Bad\n")
cat(strrep("=", 60), "\n\n")

# Load per-sample objects
sample_dir_l2 <- file.path(BASE, "cci/cellchat/L2_g3_anno2/samples")
cc_list_l2 <- load_sample_objects(sample_dir_l2)

# Build condition map — need to read metadata
cat("  Loading metadata for g3 mapping...\n")
sobj_l2 <- qread(file.path(BASE, "5_2_is_g3.qs"))
g3_mapping <- unique(sobj_l2@meta.data[, c("patient_name", "g3")])
g3_mapping$sample <- make.names(g3_mapping$patient_name)
g3_mapping$condition <- ifelse(g3_mapping$g3 == 1, "g3_Good", "g3_Bad")
rm(sobj_l2); gc(verbose = FALSE)

l2_map <- data.frame(
  sample = names(cc_list_l2),
  stringsAsFactors = FALSE
)
l2_map <- merge(l2_map, g3_mapping[, c("sample", "condition")], by = "sample", all.x = TRUE)

# Remove any unmapped samples
unmapped <- sum(is.na(l2_map$condition))
if (unmapped > 0) {
  cat("  WARNING:", unmapped, "samples could not be mapped to g3\n")
  l2_map <- l2_map[!is.na(l2_map$condition), ]
}
cat("  L2 mapping: Good =", sum(l2_map$condition == "g3_Good"),
    ", Bad =", sum(l2_map$condition == "g3_Bad"), "\n")

# Lift all objects
cc_lifted_l2 <- lift_all(cc_list_l2)

# Merge and compare
out_l2 <- file.path(BASE, "cci/plots/cellchat_L2_g3_v2")
cellchat_l2 <- merge_and_compare(cc_lifted_l2, l2_map, "L2: g3 Good vs Bad", out_l2)

rm(cc_list_l2, cc_lifted_l2); gc(verbose = FALSE)

# ==============================================================================
# Summary
# ==============================================================================
cat("\n\n", strrep("=", 60), "\n")
cat("CellChat Re-implementation Complete\n")
cat("Finished:", format(Sys.time()), "\n\n")

for (d in c("cellchat_L1_cohort_v2", "cellchat_L2_g3_v2")) {
  full_path <- file.path(BASE, "cci/plots", d)
  if (dir.exists(full_path)) {
    n_png <- length(list.files(full_path, pattern = "\\.png$"))
    cat("  ", d, ":", n_png, "PNG\n")
  }
}
cat("\nDone.\n")

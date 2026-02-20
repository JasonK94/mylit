#!/usr/bin/env Rscript
# ==============================================================================
# Fix CellChat heatmap-based plots that fail with png() device
# Strategy: Use pdf() then convert, or use ragg::agg_png which handles grid better
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
  library(ComplexHeatmap)
})

BASE <- "/data/user3/sobj/stroke_hc_v8_2"

# Helper: save ComplexHeatmap-based plot to PNG via PDF intermediate
save_ht_png <- function(expr, fname, width = 10, height = 7, res = 300) {
  # CellChat heatmap functions use ComplexHeatmap which calls draw() internally
  # This can fail with png() device. Use pdf() as intermediate then convert.
  pdf_file <- sub("\\.png$", ".pdf", fname)
  pdf(pdf_file, width = width, height = height)
  result <- tryCatch(eval(expr), error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    NULL
  })
  dev.off()

  if (!is.null(result) || file.exists(pdf_file)) {
    # Convert PDF to PNG
    system2("convert", c("-density", as.character(res),
                          pdf_file, fname),
            stdout = FALSE, stderr = FALSE)
    if (file.exists(fname) && file.size(fname) > 0) {
      file.remove(pdf_file)
      return(TRUE)
    }
  }
  # Fallback: try direct png
  png(fname, width = width, height = height, units = "in", res = res)
  tryCatch({
    eval(expr)
    dev.off()
    return(TRUE)
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    return(FALSE)
  })
}

configs <- list(
  list(name = "cellchat_L1_cohort_v2",
       conditions = c("HC", "Stroke"),
       labels = c("HC", "Stroke")),
  list(name = "cellchat_L2_g3_v2",
       conditions = c("Good", "Bad"),
       labels = c("Good (g3=1)", "Bad (g3=2)"))
)

for (cfg in configs) {
  cat("\n=== Fixing heatmaps:", cfg$name, "===\n")
  out_dir <- file.path(BASE, "cci/plots", cfg$name)

  # Load individual objects
  cc1 <- qread(file.path(out_dir, paste0("cellchat_", cfg$conditions[1], ".qs")))
  cc2 <- qread(file.path(out_dir, paste0("cellchat_", cfg$conditions[2], ".qs")))
  cc_merged <- qread(file.path(out_dir, "cellchat_comparison.qs"))

  # ---- 05-06: Differential heatmaps (merged object) ----
  for (measure_type in c("count", "weight")) {
    idx <- ifelse(measure_type == "count", "05", "06")
    fname <- file.path(out_dir, paste0(idx, "_diff_heatmap_", measure_type, ".png"))
    if (file.exists(fname) && file.size(fname) > 1000) {
      cat("  ", idx, ": already exists (", file.size(fname), " bytes)\n")
      next
    }
    cat("  Generating", idx, "...\n")

    # Try direct png with explicit draw
    tryCatch({
      png(fname, width = 8, height = 7, units = "in", res = 300)
      if (measure_type == "count") {
        ht <- netVisual_heatmap(cc_merged)
      } else {
        ht <- netVisual_heatmap(cc_merged, measure = "weight")
      }
      # ComplexHeatmap returns the object; draw it explicitly
      if (inherits(ht, "Heatmap") || inherits(ht, "HeatmapList")) {
        ComplexHeatmap::draw(ht)
      }
      dev.off()
      cat("    ", idx, ": OK (", file.size(fname), " bytes)\n")
    }, error = function(e) {
      try(dev.off(), silent = TRUE)
      cat("    ", idx, ": FAILED:", conditionMessage(e), "\n")
    })
  }

  # ---- 10-11: Signaling role heatmaps (individual objects) ----
  for (pattern in c("outgoing", "incoming")) {
    idx <- ifelse(pattern == "outgoing", "10", "11")

    for (i in 1:2) {
      cc_i <- if (i == 1) cc1 else cc2
      fname <- file.path(out_dir, paste0(idx, "_", pattern, "_", cfg$conditions[i], ".png"))
      if (file.exists(fname) && file.size(fname) > 1000) {
        cat("  ", idx, " (", cfg$conditions[i], "): already exists\n")
        next
      }
      cat("  Generating", idx, "(", cfg$conditions[i], ")...\n")

      tryCatch({
        png(fname, width = 8, height = 6, units = "in", res = 300)
        ht <- netAnalysis_signalingRole_heatmap(
          cc_i, pattern = pattern,
          title = paste(cfg$labels[i], "-", pattern),
          width = 7, height = 5
        )
        if (inherits(ht, "Heatmap") || inherits(ht, "HeatmapList")) {
          ComplexHeatmap::draw(ht)
        }
        dev.off()
        if (file.exists(fname)) {
          cat("    OK (", file.size(fname), " bytes)\n")
        }
      }, error = function(e) {
        try(dev.off(), silent = TRUE)
        cat("    FAILED:", conditionMessage(e), "\n")
      })
    }
  }

  # ---- 13: Signaling role scatter (ggplot, per condition) ----
  for (i in 1:2) {
    cc_i <- if (i == 1) cc1 else cc2
    fname <- file.path(out_dir, paste0("13_signaling_role_scatter_", cfg$conditions[i], ".png"))
    if (file.exists(fname) && file.size(fname) > 1000) {
      cat("  13 (", cfg$conditions[i], "): already exists\n")
      next
    }
    tryCatch({
      p <- netAnalysis_signalingRole_scatter(cc_i, title = cfg$labels[i])
      if (is.ggplot(p)) {
        ggsave(fname, p, width = 8, height = 8, dpi = 300)
      } else {
        # Base graphics
        png(fname, width = 8, height = 8, units = "in", res = 300)
        netAnalysis_signalingRole_scatter(cc_i, title = cfg$labels[i])
        dev.off()
      }
      cat("  13 (", cfg$conditions[i], "): OK\n")
    }, error = function(e) {
      try(dev.off(), silent = TRUE)
      cat("  13 (", cfg$conditions[i], "): FAILED:", conditionMessage(e), "\n")
    })
  }

  rm(cc1, cc2, cc_merged); gc(verbose = FALSE)
}

# Final count
cat("\n=== Final counts ===\n")
for (d in c("cellchat_L1_cohort_v2", "cellchat_L2_g3_v2")) {
  full_path <- file.path(BASE, "cci/plots", d)
  all_png <- list.files(full_path, pattern = "\\.png$", recursive = TRUE)
  main_png <- list.files(full_path, pattern = "\\.png$", recursive = FALSE)
  pw_png <- list.files(file.path(full_path, "pathways"), pattern = "\\.png$")
  sc_png <- list.files(file.path(full_path, "scatter"), pattern = "\\.png$")
  cat("  ", d, ":", length(main_png), "main +", length(pw_png), "pathway +",
      length(sc_png), "scatter =", length(all_png), "total\n")
}
cat("\nDone:", format(Sys.time()), "\n")

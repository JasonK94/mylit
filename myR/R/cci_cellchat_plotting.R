#' Plot CellChat Comparison
#'
#' Generates comparison plots for merged CellChat objects.
#'
#' @param cc_merged A merged CellChat object (result of aggregate analysis).
#' @param output_dir Directory to save plots.
#' @param width Plot width.
#' @param height Plot height.
#' @export
plot_cellchat_comparison <- function(cc_merged, output_dir, width = 10, height = 10) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    # Ensure output_dir is absolute path for safety
    output_dir <- normalizePath(output_dir)

    message("Generating comparison plots in: ", output_dir)

    # Helper to save plot
    save_plot <- function(filename, plot_func, w = width, h = height) {
        file_path <- file.path(output_dir, filename)
        message("  Saving: ", filename)

        tryCatch(
            {
                pdf(file_path, width = w, height = h)
                plot_func()
                dev.off()

                if (file.size(file_path) < 4000) {
                    warning("  WARNING: ", filename, " seems empty (size < 4KB).")
                }
            },
            error = function(e) {
                warning("  ERROR generating ", filename, ": ", e$message)
                if (file.exists(file_path)) unlink(file_path) # Delete failed file
            }
        )
    }

    # 1. Differential Number & Strength (Bar plot)
    save_plot("diff_interaction_stats.pdf", function() {
        par(mfrow = c(1, 2))
        CellChat::netVisual_diffInteraction(cc_merged, weight.scale = T)
        CellChat::netVisual_diffInteraction(cc_merged, weight.scale = T, measure = "weight")
    }, w = 12, h = 6)

    # 2. Differential Interaction Network (Circle plot)
    # Red = Increased in second group, Blue = Decreased
    save_plot("diff_interaction_circle.pdf", function() {
        par(mfrow = c(1, 2))
        CellChat::netVisual_diffInteraction(cc_merged, weight.scale = T, measure = "count.merged", label.edge = T)
        CellChat::netVisual_diffInteraction(cc_merged, weight.scale = T, measure = "weight.merged", label.edge = T)
    }, w = 12, h = 12)

    # 3. Information Flow (Stacked Bar)
    # Significant pathways in each group
    save_plot("info_flow_pathways.pdf", function() {
        print(CellChat::rankNet(cc_merged, mode = "comparison", stacked = T, do.stat = TRUE))
    }, w = 10, h = 12)

    # 4. Bubble Plot Comparison
    # Compare ligand-receptor pairs
    save_plot("comparison_bubble.pdf", function() {
        print(CellChat::netVisual_bubble(cc_merged, remove.isolate = FALSE))
    }, w = 14, h = 10)

    message("Comparison plotting complete.")
}

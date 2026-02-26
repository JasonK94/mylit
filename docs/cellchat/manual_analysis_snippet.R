# RStudio Snippet for Manual CellChat Analysis
# Copy and paste this into an R script in RStudio.

# --- 0. Fix 'read failed' Error ---
# If you see "Error: read failed on .../CellChat.rdb":
# -> Please RESTART your R session (Session > Restart R) and try again.
#    This happens when the R session loses connection to the package files.

# 1. Setup Environment
renv_lib <- "/home/user3/GJC_KDW_250721/renv/library/R-4.3/x86_64-pc-linux-gnu"
if (dir.exists(renv_lib)) .libPaths(c(renv_lib, .libPaths()))

library(CellChat)
library(Seurat)
library(qs)
library(ggplot2)
library(NMF)
library(ggalluvial)

# 2. Load Data (Make sure to use the CORRECTED run)
# Note: permissive_run3 was likely empty due to a bug. Please use the new run.
message("Loading objects...")
# Update this path to your NEW run output (e.g., permissive_run4_fixed)
cc2 <- qs::qread("/data/user3/sobj/cci/cellchat/pooled_secreted_strict_c1/samples/X2/cellchat.qs")
cc1 <- qs::qread("/data/user3/sobj/cci/cellchat/pooled_secreted_strict_c1/samples/X1/cellchat.qs")
cc2 <- CellChat::updateCellChat(cc2)
cc1 <- CellChat::updateCellChat(cc1)
# 3. Check Data Quality (IMPORTANT)
# Q: Why is @data empty?
# A: mergeCellChat drops @data (expression) to save memory.
#    This is NORMAL. We analyze @net$prob (probabilities), not @data.

if (is.null(cc2@net$prob)) stop("Stroke object has no interactions!")
print(paste("Stroke Interactions:", length(cc2@net$prob[cc2@net$prob > 0])))
print(paste("Control Interactions:", length(cc1@net$prob[cc1@net$prob > 0])))

# Check for numeric cell identities
if (all(grepl("^[0-9]+$", levels(cc2@idents)))) {
    message("WARNING: Cell identities are numeric (e.g., '1', '2'). Plots will show '1->1'.")
    message("To fix, assign biological names: levels(cc2@idents) <- c('Neuron', 'Microglia', ...)")
}

# 4. Pattern Analysis (Identify Sender/Receiver functional groups)
# Visualize outgoing signaling patterns
# Note: selectK can take time.
# selectK(cc2, pattern = "outgoing")

# Run NMF with k=5 (Adjust based on selectK plot or biological intuition)
# If you get 'dim(X)' error here -> Your @net$prob is effectively empty (signal too weak upstream).
cc2 <- identifyCommunicationPatterns(cc2, pattern = "outgoing", k = 5)

# Plot
netAnalysis_river(cc2, pattern = "outgoing")
netAnalysis_dot(cc2, pattern = "outgoing")

# 5. Differential Analysis (Visualization)
# Merge for comparison
object.list <- list(Control = cc1, Stroke = cc2)
cc.merged <- mergeCellChat(object.list, add.names = names(object.list))

# Circle plot of differential interactions
# ★ Interpretation: Red = increased in Stroke (cc2), Blue = increased in Control (cc1)

# --- Option 1: Filter by specific cell types (RECOMMENDED) ---
# Select a subset of interesting cell types to reduce clutter
ids <- levels(cc.merged@idents$joint)
message("Available cell types (", length(ids), "): ", paste(head(ids, 10), collapse = ", "), "...")

# Select cell types of interest (adjust as needed)
cells.use <- ids[1:min(8, length(ids))] # First 8 cell types
# Or manually: cells.use <- c("Microglia", "Astrocyte", "Neuron", "Oligodendrocyte")

message("Plotting diff interaction for: ", paste(cells.use, collapse = ", "))

# --- Option 2: THRESHOLD-BASED FILTERING (to reduce clutter) ---
# Instead of netVisual_diffInteraction (which plots everything),
# we manually compute the difference matrix and filter by threshold.

# Compute difference matrices
count_ctrl <- cc.merged@net$Control$count
count_stroke <- cc.merged@net$Stroke$count
weight_ctrl <- cc.merged@net$Control$weight
weight_stroke <- cc.merged@net$Stroke$weight

# Ensure same dimensions (lift to common cell types)
common_cells <- intersect(rownames(count_ctrl), rownames(count_stroke))
count_ctrl <- count_ctrl[common_cells, common_cells]
count_stroke <- count_stroke[common_cells, common_cells]
weight_ctrl <- weight_ctrl[common_cells, common_cells]
weight_stroke <- weight_stroke[common_cells, common_cells]

# Calculate difference
count_diff <- count_stroke - count_ctrl
weight_diff <- weight_stroke - weight_ctrl

# ★ THRESHOLD FILTERING ★
# Set edges below threshold to 0 (they won't be plotted)
count_thresh <- 2 # Only show edges with |count difference| >= 2
weight_thresh <- 0.05 # Only show edges with |weight difference| >= 0.05

count_diff_filtered <- count_diff
count_diff_filtered[abs(count_diff) < count_thresh] <- 0

weight_diff_filtered <- weight_diff
weight_diff_filtered[abs(weight_diff) < weight_thresh] <- 0

# Subset to selected cell types
if (!is.null(cells.use)) {
    cells.use <- intersect(cells.use, common_cells)
    count_diff_filtered <- count_diff_filtered[cells.use, cells.use]
    weight_diff_filtered <- weight_diff_filtered[cells.use, cells.use]
}

message(
    "Edges remaining after count filter (thresh=", count_thresh, "): ",
    sum(count_diff_filtered != 0)
)
message(
    "Edges remaining after weight filter (thresh=", weight_thresh, "): ",
    sum(weight_diff_filtered != 0)
)

# Plot using netVisual_circle (base function)
par(mfrow = c(1, 2))
tryCatch(
    {
        # Positive values = Stroke up, Negative = Control up
        # We split into two colors by sign

        # Count Difference Plot
        netVisual_circle(count_diff_filtered,
            weight.scale = TRUE,
            edge.width.max = 10,
            title.name = paste0("Diff Count (thresh≥", count_thresh, ")\nRed=↑Stroke, Blue=↑Control")
        )

        # Weight Difference Plot
        netVisual_circle(weight_diff_filtered,
            weight.scale = TRUE,
            edge.width.max = 10,
            title.name = paste0("Diff Weight (thresh≥", weight_thresh, ")\nRed=↑Stroke, Blue=↑Control")
        )
    },
    error = function(e) message("Diff plot error: ", e$message)
)

# Information Flow (Pathway level changes)
# Shows which pathways are turned ON/OFF or increased/decreased
rankNet(cc.merged, mode = "comparison", stacked = T, do.stat = TRUE)

# Bubble Plot for specific pathways (Top significant ones)
if (length(cc2@netP$pathways) > 0) {
    stroke_pathways <- cc2@netP$pathways
    message("Attempting to find plottable pathways from top 5: ", paste(head(stroke_pathways, 5), collapse = ", "))

    # --- Filtering Options (User Config) ---
    # 1. Inspect P-values to understand why 'thresh' might not work
    #    (If all p-values are 0 or < 0.01, changing thresh won't filter much)
    message("P-value summary (Stroke):")
    print(summary(as.vector(cc.merged@net$Stroke$pval)))

    # 2. Filter by Cell Groups (Sender/Receiver) to reduce x-axis
    #    ★ KEY FIX: Select only a FEW cell types to avoid hairball plots ★
    ids <- levels(cc.merged@idents)
    message("Available cell types (", length(ids), "): ", paste(ids, collapse = ", "))

    # Select first 3 cell types (or manually specify)
    sources.use <- ids[1:min(3, length(ids))]
    targets.use <- ids[1:min(3, length(ids))]
    # Or: sources.use <- c("Microglia", "Astrocyte", "Neuron")

    message("Using sources: ", paste(sources.use, collapse = ", "))
    message("Using targets: ", paste(targets.use, collapse = ", "))

    # Test pathways one by one (simplified - no manual prob filtering)
    # NOTE: P-values are high (~1.0), so we set thresh=1.0 to show ALL interactions regardless of significance.
    #       This is for EXPLORATORY visualization. Interpret with caution.
    valid_plot_generated <- FALSE
    for (pw in head(stroke_pathways, 10)) {
        message("Testing pathway: ", pw)
        tryCatch(
            {
                p <- netVisual_bubble(cc.merged,
                    signaling = pw,
                    remove.isolate = TRUE,
                    comparison = c(1, 2),
                    angle.x = 45,
                    sources.use = sources.use,
                    targets.use = targets.use,
                    thresh = 0.001 # Set to 1.0 to bypass p-value filter (since all p~1)
                )
                print(p)
                valid_plot_generated <- TRUE
                message("  -> Success! Plotted ", pw)
            },
            error = function(e) {
                message("  -> Failed: ", e$message)
            }
        )

        if (valid_plot_generated) break
    }

    if (!valid_plot_generated) {
        message("Could not generate bubble plot for any of the top 10 pathways.")
        message("Tip: Your CellChat objects may have very weak signals. Check cc2@net$prob manually.")
    }
}

message("Analysis Ready.")

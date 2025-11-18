#' CCI Utility Functions
#'
#' This module provides utility functions for CCI analysis.
#'
#' @name utils_cci
NULL

#' Format DEG Summary Table
#'
#' Creates a formatted summary table from DEG results.
#'
#' @param deg_df Data frame with DEG results
#' @param receiver_cluster Character string, receiver cluster ID
#'
#' @return Data frame with summary statistics
#' @export
format_deg_summary <- function(deg_df, receiver_cluster) {
  
  receiver_degs <- deg_df %>%
    dplyr::filter(cluster == receiver_cluster)
  
  if (nrow(receiver_degs) == 0) {
    return(data.frame(
      cluster = receiver_cluster,
      n_degs = 0,
      n_upregulated = 0,
      n_downregulated = 0
    ))
  }
  
  # Determine logFC column
  logfc_col <- if ("avg_log2FC" %in% colnames(receiver_degs)) "avg_log2FC" else "logFC"
  
  summary_df <- data.frame(
    cluster = receiver_cluster,
    n_degs = nrow(receiver_degs),
    n_upregulated = sum(receiver_degs[[logfc_col]] > 0, na.rm = TRUE),
    n_downregulated = sum(receiver_degs[[logfc_col]] < 0, na.rm = TRUE),
    mean_logfc = mean(receiver_degs[[logfc_col]], na.rm = TRUE),
    median_logfc = median(receiver_degs[[logfc_col]], na.rm = TRUE)
  )
  
  # Add p-value summary if available
  pval_col <- NULL
  if ("p_val_adj" %in% colnames(receiver_degs)) {
    pval_col <- "p_val_adj"
  } else if ("FDR" %in% colnames(receiver_degs)) {
    pval_col <- "FDR"
  } else if ("p_val" %in% colnames(receiver_degs)) {
    pval_col <- "p_val"
  }
  
  if (!is.null(pval_col)) {
    summary_df$mean_pval <- mean(receiver_degs[[pval_col]], na.rm = TRUE)
    summary_df$median_pval <- median(receiver_degs[[pval_col]], na.rm = TRUE)
    summary_df$n_significant <- sum(receiver_degs[[pval_col]] < 0.05, na.rm = TRUE)
  }
  
  return(summary_df)
}

#' Identify Top Senders by DEG Count
#'
#' Identifies top sender clusters based on number of DEGs or other criteria.
#' This is a helper function that can be used for sender prioritization.
#'
#' @param deg_df Data frame with DEG results
#' @param sender_clusters Character vector of sender cluster IDs
#' @param top_n Integer, number of top senders to return (default: NULL, return all)
#'
#' @return Character vector of top sender cluster IDs
#' @export
identify_top_senders <- function(deg_df, sender_clusters, top_n = NULL) {
  
  # Count DEGs per cluster
  deg_counts <- deg_df %>%
    dplyr::filter(cluster %in% sender_clusters) %>%
    dplyr::count(cluster, sort = TRUE)
  
  if (nrow(deg_counts) == 0) {
    return(character(0))
  }
  
  top_senders <- deg_counts$cluster
  
  if (!is.null(top_n) && top_n < length(top_senders)) {
    top_senders <- top_senders[1:top_n]
  }
  
  return(top_senders)
}

#' Create Sender-Receiver Map
#'
#' Creates a mapping data frame between sender and receiver clusters.
#'
#' @param sender_clusters Character vector of sender cluster IDs
#' @param receiver_cluster Character string, receiver cluster ID
#'
#' @return Data frame with sender-receiver mappings
#' @export
create_sender_receiver_map <- function(sender_clusters, receiver_cluster) {
  
  map_df <- data.frame(
    sender_cluster = sender_clusters,
    receiver_cluster = receiver_cluster,
    stringsAsFactors = FALSE
  )
  
  return(map_df)
}

#' Resume CCI Analysis from Prepared Data
#'
#' Loads prepared data from checkpoint and continues with NicheNet analysis.
#' Useful when you want to resume from a saved intermediate checkpoint.
#'
#' @param prepared_data_file Character string, path to `cci_prepared_data_*.qs` file
#' @param sobj Seurat object (required for NicheNet analysis)
#' @param cluster_col Character string, cluster column name (default: NULL, auto-detected from common names)
#' @param condition_col Character string, condition column for DE analysis
#' @param condition_oi Character string, condition of interest
#' @param condition_ref Character string, reference condition
#' @param species Character string, "human" or "mouse" (default: "human")
#' @param assay_name Character string, assay to use (default: "SCT")
#' @param top_n_ligands Integer, number of top ligands (default: 20)
#' @param top_n_targets_per_ligand Integer, number of top targets per ligand (default: NULL, auto-adjusted based on DEG count)
#' @param ligand_target_cutoff Numeric, cutoff for ligand-target links (default: 0.33)
#' @param nichenet_data_dir Character string, NicheNet data directory (default: NULL)
#' @param output_dir Character string, output directory (default: NULL)
#' @param run_circos Logical, whether to generate circos plot (default: TRUE)
#' @param verbose Logical, whether to print progress messages (default: TRUE)
#' @param ... Additional arguments passed to run_nichenet_analysis
#'
#' @return List with CCI analysis results
#' @export
resume_cci_from_prepared <- function(prepared_data_file,
                                      sobj,
                                      cluster_col = NULL,
                                      condition_col,
                                      condition_oi,
                                      condition_ref = NULL,
                                      species = c("human", "mouse"),
                                      assay_name = "SCT",
                                      top_n_ligands = 20,
                                      top_n_targets_per_ligand = NULL,
                                      ligand_target_cutoff = 0.33,
                                      nichenet_data_dir = NULL,
                                      output_dir = NULL,
                                      run_circos = TRUE,
                                      verbose = TRUE,
                                      ...) {
  
  if (!requireNamespace("qs", quietly = TRUE)) {
    stop("qs package is required for loading prepared data")
  }
  
  if (!file.exists(prepared_data_file)) {
    stop("Prepared data file not found: ", prepared_data_file)
  }
  
  if (verbose) message("Loading prepared data from: ", prepared_data_file)
  prepared_data <- qs::qread(prepared_data_file)
  
  # Extract components
  receiver_degs <- prepared_data$receiver_degs
  sender_clusters <- prepared_data$sender_clusters
  receiver_cluster <- prepared_data$receiver_cluster
  
  if (is.null(receiver_degs) || nrow(receiver_degs) == 0) {
    stop("No receiver DEGs found in prepared data")
  }
  
  # Auto-adjust top_n_targets_per_ligand if not specified
  if (is.null(top_n_targets_per_ligand)) {
    n_degs <- nrow(receiver_degs)
    if (n_degs > 3000) {
      top_n_targets_per_ligand <- 50  # Reduce for large DEG sets
      if (verbose) message("Large DEG set (", n_degs, " genes). Reducing top_n_targets_per_ligand to ", top_n_targets_per_ligand, " for faster computation.")
    } else if (n_degs > 1000) {
      top_n_targets_per_ligand <- 100
      if (verbose) message("Moderate DEG set (", n_degs, " genes). Setting top_n_targets_per_ligand to ", top_n_targets_per_ligand, ".")
    } else {
      top_n_targets_per_ligand <- 200  # Default for small sets
    }
  }
  
  # Check if run_nichenet_analysis is available
  if (!exists("run_nichenet_analysis")) {
    source_candidates <- c(
      "/home/user3/data_user3/git_repo/_wt/cci/myR/R/CCI.R",
      "/home/user3/data_user3/git_repo/mylit/myR/R/CCI.R"
    )
    sourced <- FALSE
    for (candidate in source_candidates) {
      if (file.exists(candidate)) {
        source(candidate)
        sourced <- TRUE
        if (verbose) message("  Loaded run_nichenet_analysis from ", candidate)
        break
      }
    }
    if (!sourced) {
      stop("Cannot find run_nichenet_analysis function. Please ensure CCI.R is sourced.")
    }
  }
  
  if (verbose) {
    message("Resuming NicheNet analysis with:")
    message("  - Receiver: ", receiver_cluster)
    message("  - Senders: ", length(sender_clusters), " cluster(s)")
    message("  - DEGs: ", nrow(receiver_degs))
    message("  - Top targets per ligand: ", top_n_targets_per_ligand)
  }
  
  # Determine cluster_col if not provided
  if (is.null(cluster_col)) {
    if (!is.null(sobj@meta.data)) {
      # Try to find a common cluster column
      possible_cols <- c("anno3.scvi", "seurat_clusters", "cluster", "cell_type")
      for (col in possible_cols) {
        if (col %in% colnames(sobj@meta.data)) {
          cluster_col <- col
          if (verbose) message("Auto-detected cluster_col: ", cluster_col)
          break
        }
      }
    }
    if (is.null(cluster_col)) {
      stop("Cannot determine cluster_col from Seurat object. Please specify cluster_col parameter.")
    }
  }
  
  nichenet_results <- run_nichenet_analysis(
    seurat_obj = sobj,
    species = species,
    sender_celltypes = sender_clusters,
    receiver_celltype = receiver_cluster,
    assay_name = assay_name,
    cluster_col = cluster_col,
    receiver_DE_ident1 = condition_oi,
    receiver_DE_ident2 = condition_ref,
    receiver_DE_group_by = condition_col,
    receiver_de_table = receiver_degs,
    top_n_ligands = top_n_ligands,
    top_n_targets_per_ligand = top_n_targets_per_ligand,
    ligand_target_cutoff = ligand_target_cutoff,
    nichenet_data_dir = nichenet_data_dir,
    output_dir = output_dir,
    run_circos = run_circos,
    verbose = verbose,
    ...
  )
  
  # Compile results
  results <- list(
    nichenet_results = nichenet_results,
    sender_receiver_map = create_sender_receiver_map(sender_clusters, receiver_cluster),
    deg_summary = format_deg_summary(receiver_degs, receiver_cluster),
    receiver_degs = receiver_degs,
    sender_clusters = sender_clusters,
    receiver_cluster = receiver_cluster,
    prepared_summary = prepared_data$summary,
    output_path = output_dir
  )
  
  return(results)
}

#' Replot from Saved NicheNet Results
#'
#' Loads saved NicheNet results from a .qs file and replots the circos plot
#' with customizable parameters.
#'
#' @param results_file Character string, path to the saved results file (e.g., `nichenet_results.qs`)
#' @param output_file Character string, optional path to save the replotted figure. If NULL, plots to current device.
#' @param circos_cex_text Numeric, text size for circos plot labels (default: 0.8)
#' @param circos_show_legend Logical, whether to show legend (default: TRUE)
#' @param circos_legend_position Character string, legend position (default: "topright")
#' @param circos_legend_size Numeric, legend text size (default: 0.9)
#' @param circos_legend_inset Numeric vector, legend inset (default: c(-0.15, 0))
#' @param format Character string, output format: "png", "pdf", or "screen" (default: "screen")
#' @param width Numeric, plot width in inches (default: 10)
#' @param height Numeric, plot height in inches (default: 10)
#' @param res Numeric, resolution for PNG (default: 100)
#' @param verbose Logical, whether to print messages (default: TRUE)
#'
#' @return Invisible list with plot object and metadata
#' @export
#'
#' @examples
#' \dontrun{
#' # Replot from saved results
#' replot_nichenet_circos(
#'   results_file = "/data/user3/sobj/run6/nichenet_results.qs",
#'   format = "png",
#'   output_file = "/data/user3/sobj/run6/circos_replot.png",
#'   circos_cex_text = 1.0,
#'   circos_show_legend = TRUE
#' )
#' }
replot_nichenet_circos <- function(results_file,
                                    output_file = NULL,
                                    circos_cex_text = 0.8,
                                    circos_show_legend = TRUE,
                                    circos_legend_position = "topright",
                                    circos_legend_size = 0.9,
                                    circos_legend_inset = c(-0.15, 0),
                                    format = c("screen", "png", "pdf"),
                                    width = 10,
                                    height = 10,
                                    res = 100,
                                    verbose = TRUE) {
  
  if (!requireNamespace("qs", quietly = TRUE)) {
    stop("qs package is required for loading saved results")
  }
  
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("circlize package is required for replotting")
  }
  
  format <- match.arg(format)
  
  if (!file.exists(results_file)) {
    stop("Results file not found: ", results_file)
  }
  
  if (verbose) message("Loading results from: ", results_file)
  results <- qs::qread(results_file)
  
  # Check file structure - handle both nested and flat structures
  # Structure 1: results$nichenet_results$plot_circos (from run_cci_analysis)
  # Structure 2: results$plot_circos (direct NicheNet results)
  if (!is.null(results$nichenet_results) && !is.null(results$nichenet_results$plot_circos)) {
    # Nested structure
    nichenet_results <- results$nichenet_results
    plot_circos <- nichenet_results$plot_circos
    if (verbose) message("Found nested structure: results$nichenet_results$plot_circos")
  } else if (!is.null(results$plot_circos)) {
    # Flat structure (direct NicheNet results)
    nichenet_results <- results
    plot_circos <- results$plot_circos
    if (verbose) message("Found flat structure: results$plot_circos")
  } else {
    stop("No circos plot found in results. Please check the file structure. ",
         "Expected either 'results$nichenet_results$plot_circos' or 'results$plot_circos'.")
  }
  
  # Check if circos plot data exists
  if (is.null(plot_circos)) {
    stop("No circos plot found in results. The original analysis may not have generated a circos plot.")
  }
  
  # Try to extract circos data from the saved plot
  # If plot_circos is a recorded plot, we need to replay it and extract data
  # For now, we'll replay the original plot with new parameters
  if (verbose) message("Replaying circos plot with new parameters...")
  
  # Open graphics device
  device_opened <- FALSE
  current_device <- NULL
  
  if (format == "png" && !is.null(output_file)) {
    if (verbose) message("Opening PNG device: ", output_file)
    grDevices::png(output_file, width = width * res, height = height * res, res = res)
    device_opened <- TRUE
    current_device <- grDevices::dev.cur()
    if (verbose) message("  Device opened: ", current_device)
  } else if (format == "pdf" && !is.null(output_file)) {
    if (verbose) message("Opening PDF device: ", output_file)
    grDevices::pdf(output_file, width = width, height = height)
    device_opened <- TRUE
    current_device <- grDevices::dev.cur()
    if (verbose) message("  Device opened: ", current_device)
  }
  
  tryCatch({
    # Replay the original plot
    if (inherits(plot_circos, "recordedplot")) {
      if (verbose) message("Replaying recorded plot...")
      # Ensure we're using the correct device
      if (device_opened) {
        grDevices::dev.set(current_device)
      }
      # Replay the recorded plot
      grDevices::replayPlot(plot_circos)
      
      # Note: We cannot easily modify parameters of a recorded plot
      # The user would need to re-run the analysis with different parameters
      # or we would need to store the raw circos data separately
      if (verbose) {
        message("Note: Recorded plots cannot be modified. Using original plot.")
        message("To change parameters, re-run the analysis with different circos_* parameters.")
      }
    } else {
      # If it's not a recorded plot, try to plot it directly
      if (verbose) message("Plotting circos plot object...")
      if (device_opened) {
        grDevices::dev.set(current_device)
      }
      graphics::plot(plot_circos)
    }
    
    # Ensure device is properly closed
    if (device_opened) {
      if (verbose) message("Closing graphics device...")
      grDevices::dev.off()
      # Wait a moment for file system to sync
      Sys.sleep(0.1)
      
      if (verbose && !is.null(output_file)) {
        if (file.exists(output_file)) {
          file_size <- file.info(output_file)$size
          message("Plot saved to: ", output_file, " (", file_size, " bytes)")
        } else {
          warning("Warning: Output file was not created: ", output_file)
          warning("  Current working directory: ", getwd())
          warning("  Device was: ", current_device)
        }
      }
    }
  }, error = function(e) {
    if (device_opened) {
      try({
        if (!is.null(current_device) && current_device %in% grDevices::dev.list()) {
          grDevices::dev.set(current_device)
          grDevices::dev.off()
        }
      }, silent = TRUE)
    }
    stop("Error during plot generation: ", e$message)
  })
  
  return(invisible(list(
    plot = plot_circos,
    results = nichenet_results,
    metadata = list(
      format = format,
      output_file = output_file,
      parameters = list(
        circos_cex_text = circos_cex_text,
        circos_show_legend = circos_show_legend,
        circos_legend_position = circos_legend_position,
        circos_legend_size = circos_legend_size,
        circos_legend_inset = circos_legend_inset
      )
    )
  )))
}


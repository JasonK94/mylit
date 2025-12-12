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
      top_n_targets_per_ligand <- 50 # Reduce for large DEG sets
      if (verbose) message("Large DEG set (", n_degs, " genes). Reducing top_n_targets_per_ligand to ", top_n_targets_per_ligand, " for faster computation.")
    } else if (n_degs > 1000) {
      top_n_targets_per_ligand <- 100
      if (verbose) message("Moderate DEG set (", n_degs, " genes). Setting top_n_targets_per_ligand to ", top_n_targets_per_ligand, ".")
    } else {
      top_n_targets_per_ligand <- 200 # Default for small sets
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
#' Replots the circos plot from NicheNet results, either from a saved .qs file
#' or from a results object already in memory.
#'
#' @param results_file Character string, path to the saved results file (e.g., `nichenet_results.qs`).
#'                    If NULL, `results` parameter must be provided.
#' @param results List, NicheNet results object already in memory. If provided, `results_file` is ignored.
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
#' # Replot from saved file
#' replot_nichenet_circos(
#'   results_file = "/data/user3/sobj/run6/nichenet_results.qs",
#'   format = "png",
#'   output_file = "/data/user3/sobj/run6/circos_replot.png"
#' )
#'
#' # Replot from results object in memory (faster for testing)
#' results <- run_cci_analysis(...)
#' replot_nichenet_circos(
#'   results = results$nichenet_results,
#'   format = "png",
#'   output_file = "circos_replot.png"
#' )
#' }
replot_nichenet_circos <- function(results_file = NULL,
                                   results = NULL,
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
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("circlize package is required for replotting")
  }

  format <- match.arg(format)

  # Load results from file or use provided results object
  if (!is.null(results)) {
    # Use provided results object (faster, no file I/O)
    if (verbose) message("Using results object from memory")
    results_obj <- results
  } else if (!is.null(results_file)) {
    # Load from file
    if (!requireNamespace("qs", quietly = TRUE)) {
      stop("qs package is required for loading saved results")
    }
    if (!file.exists(results_file)) {
      stop("Results file not found: ", results_file)
    }
    if (verbose) message("Loading results from: ", results_file)
    results_obj <- qs::qread(results_file)
  } else {
    stop("Either 'results_file' or 'results' must be provided")
  }

  # Check file structure - handle both nested and flat structures
  # Structure 1: results$nichenet_results$plot_circos (from run_cci_analysis)
  # Structure 2: results$plot_circos (direct NicheNet results)
  if (!is.null(results_obj$nichenet_results) && !is.null(results_obj$nichenet_results$plot_circos)) {
    # Nested structure
    nichenet_results <- results_obj$nichenet_results
    plot_circos <- nichenet_results$plot_circos
    if (verbose) message("Found nested structure: results$nichenet_results$plot_circos")
  } else if (!is.null(results_obj$plot_circos)) {
    # Flat structure (direct NicheNet results)
    nichenet_results <- results_obj
    plot_circos <- results_obj$plot_circos
    if (verbose) message("Found flat structure: results$plot_circos")
  } else {
    stop(
      "No circos plot found in results. Please check the structure. ",
      "Expected either 'results$nichenet_results$plot_circos' or 'results$plot_circos'."
    )
  }

  # Check if circos plot data exists
  if (is.null(plot_circos)) {
    stop("No circos plot found in results. The original analysis may not have generated a circos plot.")
  }

  # Try to extract circos data from the saved plot
  # If plot_circos is a recorded plot, we need to replay it and extract data
  # For now, we'll replay the original plot with new parameters
  if (verbose) message("Replaying circos plot with new parameters...")

  # Check for existing devices and close non-null devices if needed
  existing_devices <- grDevices::dev.list()
  if (verbose && length(existing_devices) > 0) {
    message("  Existing devices before opening: ", paste(existing_devices, collapse = ", "))
  }

  # Open graphics device
  device_opened <- FALSE
  current_device <- NULL

  if (format == "png" && !is.null(output_file)) {
    if (verbose) message("Opening PNG device: ", output_file)
    # Ensure absolute path
    abs_output_file <- normalizePath(output_file, mustWork = FALSE)
    if (!file.exists(dirname(abs_output_file))) {
      dir.create(dirname(abs_output_file), recursive = TRUE, showWarnings = FALSE)
    }
    grDevices::png(abs_output_file, width = width * res, height = height * res, res = res)
    device_opened <- TRUE
    current_device <- grDevices::dev.cur()
    if (verbose) message("  Device opened: ", current_device)
  } else if (format == "pdf" && !is.null(output_file)) {
    if (verbose) message("Opening PDF device: ", output_file)
    # Ensure absolute path
    abs_output_file <- normalizePath(output_file, mustWork = FALSE)
    if (!file.exists(dirname(abs_output_file))) {
      dir.create(dirname(abs_output_file), recursive = TRUE, showWarnings = FALSE)
    }
    grDevices::pdf(abs_output_file, width = width, height = height)
    device_opened <- TRUE
    current_device <- grDevices::dev.cur()
    if (verbose) message("  Device opened: ", current_device)
  }

  tryCatch(
    {
      # Replay the original plot
      if (inherits(plot_circos, "recordedplot")) {
        if (verbose) message("Replaying recorded plot...")
        # Ensure we're using the correct device
        if (device_opened) {
          active_dev <- grDevices::dev.set(current_device)
          if (verbose) message("  Active device set to: ", active_dev)
          if (active_dev != current_device) {
            warning("Warning: Device mismatch. Expected ", current_device, " but got ", active_dev)
          }
        }

        # Check if plot is valid before replaying
        if (is.null(plot_circos)) {
          stop("Plot object is NULL. Cannot replay plot.")
        }

        if (inherits(plot_circos, "recordedplot")) {
          if (length(plot_circos) == 0) {
            stop(
              "Plot object is empty (length 0). Cannot replay plot. ",
              "This may indicate the plot was not properly recorded during analysis."
            )
          }
          # Check if plot has data: either [[1]] (plot commands) or [[2]] (raw graphics data) should exist
          # Note: [[1]] can be empty (length 0) but [[2]] may contain raw graphics data
          has_plot_data <- FALSE
          if (length(plot_circos) >= 1 && !is.null(plot_circos[[1]]) && length(plot_circos[[1]]) > 0) {
            has_plot_data <- TRUE
          } else if (length(plot_circos) >= 2 && !is.null(plot_circos[[2]]) && length(plot_circos[[2]]) > 0) {
            # [[2]] contains raw graphics data (can be raw bytes)
            has_plot_data <- TRUE
          }
          if (!has_plot_data) {
            stop(
              "Plot object contains no plot data. Cannot replay plot. ",
              "This may indicate the plot was not properly recorded during analysis. ",
              "Please re-run the analysis to regenerate the plot."
            )
          }
        }

        # Replay the recorded plot
        # Use tryCatch to catch any errors during replay
        replay_success <- tryCatch(
          {
            grDevices::replayPlot(plot_circos)
            TRUE
          },
          error = function(e) {
            if (verbose) warning("Error during replayPlot: ", e$message)
            FALSE
          }
        )

        if (!replay_success) {
          stop("Failed to replay plot. The recorded plot may be empty or corrupted.")
        }

        # Force flush the device to ensure plot is written
        if (device_opened) {
          grDevices::dev.flush()
          # Verify device is still active
          if (grDevices::dev.cur() != current_device) {
            grDevices::dev.set(current_device)
            if (verbose) message("  Device reset to: ", current_device)
          }
          if (verbose) message("  Device flushed and verified")
        }

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
          active_dev <- grDevices::dev.set(current_device)
          if (verbose) message("  Active device set to: ", active_dev)
        }
        graphics::plot(plot_circos)

        # Force flush the device
        if (device_opened) {
          grDevices::dev.flush()
          if (verbose) message("  Device flushed")
        }
      }

      # Ensure device is properly closed
      if (device_opened) {
        if (verbose) {
          message("Closing graphics device...")
          message("  Current device before close: ", grDevices::dev.cur())
        }

        # Verify we're closing the correct device
        if (grDevices::dev.cur() != current_device) {
          if (verbose) warning("Warning: Current device (", grDevices::dev.cur(), ") differs from expected (", current_device, ")")
          grDevices::dev.set(current_device)
        }

        # Force flush before closing to ensure all data is written
        grDevices::dev.flush()

        # Close the device
        grDevices::dev.off()

        # Wait a moment for file system to sync (increased from 0.3 to 0.5)
        Sys.sleep(0.5)

        if (verbose && !is.null(output_file)) {
          # Check if file exists (use absolute path)
          abs_output_file <- normalizePath(output_file, mustWork = FALSE)
          if (file.exists(abs_output_file)) {
            file_size <- file.info(abs_output_file)$size
            message("Plot saved to: ", abs_output_file, " (", file_size, " bytes)")
          } else if (file.exists(output_file)) {
            file_size <- file.info(output_file)$size
            message("Plot saved to: ", output_file, " (", file_size, " bytes)")
          } else {
            warning("Warning: Output file was not created: ", output_file)
            warning("  Absolute path: ", abs_output_file)
            warning("  Current working directory: ", getwd())
            warning("  Device was: ", current_device)
            warning("  Device list after close: ", paste(grDevices::dev.list(), collapse = ", "))
            # Try to check if file exists in parent directory
            parent_dir <- dirname(abs_output_file)
            if (dir.exists(parent_dir)) {
              files_in_dir <- list.files(parent_dir, pattern = basename(output_file))
              if (length(files_in_dir) > 0) {
                warning("  Found similar files in directory: ", paste(files_in_dir, collapse = ", "))
              }
            }
          }
        }
      }
    },
    error = function(e) {
      if (device_opened) {
        try(
          {
            # Try to close the device we opened
            if (!is.null(current_device)) {
              if (current_device %in% grDevices::dev.list()) {
                grDevices::dev.set(current_device)
                grDevices::dev.off()
              }
            }
            # Also try to close any remaining non-null devices
            remaining_devices <- grDevices::dev.list()
            if (length(remaining_devices) > 0 && !is.null(remaining_devices)) {
              for (dev in remaining_devices) {
                if (dev != grDevices::dev.cur()) {
                  try(grDevices::dev.off(dev), silent = TRUE)
                }
              }
            }
          },
          silent = TRUE
        )
      }
      stop("Error during plot generation: ", e$message)
    }
  )

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

#' Load NicheNet Reference Data
#'
#' Loads LR network and Ligand-Target matrix from shared directory or downloads them.
#' Handles .rds/.qs conversion.
#'
#' @param species Character string, "human" or "mouse"
#' @param data_dir Character string, custom data directory (optional)
#' @param verbose Logical, whether to print messages (default: TRUE)
#'
#' @return List containing 'lr_network' and 'ligand_target_matrix'
#' @export
load_nichenet_reference <- function(species = "human", data_dir = NULL, verbose = TRUE) {
  if (!requireNamespace("qs", quietly = TRUE)) {
    stop("qs package is required.")
  }

  # Default shared directory
  shared_base <- "/data/user3/git_repo"
  shared_dir <- file.path(shared_base, species)
  user_data_path <- "/home/user3/GJC_KDW_250721/nichenet_data_human"

  if (is.null(data_dir)) {
    if (species == "human" && dir.exists(user_data_path)) {
      data_dir <- user_data_path
      if (verbose) message("Using user-specific NicheNet data directory: ", data_dir)
    } else if (dir.exists(shared_dir)) {
      data_dir <- shared_dir
      if (verbose) message("Using shared NicheNet data directory: ", data_dir)
    } else {
      data_dir <- file.path(getwd(), paste0("nichenet_data_", species))
      if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)
      if (verbose) message("Using local NicheNet data directory: ", data_dir)
    }
  }

  base_url <- "https://zenodo.org/records/7074291/files/"

  files_spec <- list(
    human = list(
      lr_network = "lr_network_human_21122021",
      ligand_target_matrix = "ligand_target_matrix_nsga2r_final"
    ),
    mouse = list(
      lr_network = "lr_network_mouse_21122021",
      ligand_target_matrix = "ligand_target_matrix_nsga2r_final_mouse"
    )
  )

  if (!species %in% names(files_spec)) stop("Unsupported species: ", species)

  req_files <- files_spec[[species]]

  results <- list()

  for (key in names(req_files)) {
    f_base <- req_files[[key]]
    qs_path <- file.path(data_dir, paste0(f_base, ".qs"))
    rds_path <- file.path(data_dir, paste0(f_base, ".rds"))

    if (file.exists(qs_path)) {
      if (verbose) message("Loading ", basename(qs_path), "...")
      results[[key]] <- qs::qread(qs_path)
    } else if (file.exists(rds_path)) {
      if (verbose) message("Loading ", basename(rds_path), " (converting to .qs)...")
      obj <- readRDS(rds_path)
      qs::qsave(obj, qs_path) # Convert for future use
      if (verbose) message("  Saved .qs version for faster loading.")
      results[[key]] <- obj
    } else {
      # Download RDS
      if (verbose) message("Downloading ", f_base, ".rds ...")
      download.file(paste0(base_url, f_base, ".rds"), rds_path, mode = "wb", quiet = !verbose)
      obj <- readRDS(rds_path)
      qs::qsave(obj, qs_path)
      results[[key]] <- obj
    }
  }

  if ("lr_network" %in% names(results)) {
    results$lr_network <- results$lr_network %>% dplyr::distinct(from, to)
  }

  return(results)
}

#' Save Parameter Log
#'
#' Saves a list of parameters to a log file.
#'
#' @param params List of parameters
#' @param output_dir Directory to save the log file
#' @param prefix Prefix for the log file name (default: "analysis_params")
#' @export
save_parameters_log <- function(params, output_dir, prefix = "analysis_params") {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  log_file <- file.path(output_dir, paste0(prefix, ".txt"))

  sink(log_file)
  cat("Analysis Parameters\n")
  cat("===================\n")
  cat("Date: ", as.character(Sys.time()), "\n\n")

  for (name in names(params)) {
    val <- params[[name]]
    if (length(val) > 1) {
      val_str <- paste(val, collapse = ", ")
    } else if (is.null(val)) {
      val_str <- "NULL"
    } else {
      val_str <- as.character(val)
    }
    cat(sprintf("%-20s: %s\n", name, val_str))
  }
  sink()
  message("Parameters saved to: ", log_file)
}

#' Get Safe Filename
#'
#' Generates a unique filename by appending an incrementing suffix if the file already exists.
#'
#' @param filepath Full path to the desired file
#' @return A unique filepath
#' @export
get_safe_filename <- function(filepath) {
  if (!file.exists(filepath)) {
    return(filepath)
  }

  dir <- dirname(filepath)
  fname <- basename(filepath)
  ext <- tools::file_ext(fname)
  name <- tools::file_path_sans_ext(fname)

  i <- 1
  new_filepath <- file.path(dir, paste0(name, "_", i, ".", ext))
  while (file.exists(new_filepath)) {
    i <- i + 1
    new_filepath <- file.path(dir, paste0(name, "_", i, ".", ext))
  }
  return(new_filepath)
}

#' Save Plot in Multiple Formats
#'
#' Saves the current plot or a plot object to PDF and PNG.
#'
#' @param plot_obj Detailed plot object (for ggplot usually) or NULL (for base plot recorded with recordPlot) or "base" (to execute code)
#' @param filename Base filename without extension
#' @param outdir Output directory
#' @param width Width in inches
#' @param height Height in inches
#' @param plot_code Expression or function to draw the base plot (optional, for base graphics)
#' @export
save_plot_formats <- function(filename, outdir, width = 10, height = 8, plot_obj = NULL, plot_code = NULL) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  # Safe filenames
  pdf_path <- get_safe_filename(file.path(outdir, paste0(filename, ".pdf")))
  png_path <- get_safe_filename(file.path(outdir, paste0(filename, ".png")))

  # PDF
  pdf(pdf_path, width = width, height = height)
  if (!is.null(plot_obj) && (inherits(plot_obj, "ggplot") || inherits(plot_obj, "grob"))) {
    print(plot_obj)
  } else if (!is.null(plot_code)) {
    if (is.function(plot_code)) plot_code() else eval(plot_code)
  } else if (!is.null(plot_obj) && inherits(plot_obj, "recordedplot")) {
    replayPlot(plot_obj)
  }
  dev.off()

  # PNG
  png(png_path, width = width * 100, height = height * 100, res = 100)
  if (!is.null(plot_obj) && (inherits(plot_obj, "ggplot") || inherits(plot_obj, "grob"))) {
    print(plot_obj)
  } else if (!is.null(plot_code)) {
    if (is.function(plot_code)) plot_code() else eval(plot_code)
  } else if (!is.null(plot_obj) && inherits(plot_obj, "recordedplot")) {
    replayPlot(plot_obj)
  }
  dev.off()

  message("Saved plots: ", pdf_path, " / ", png_path)
  return(c(pdf = pdf_path, png = png_path))
}

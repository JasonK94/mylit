#' Run CellChat Analysis with Proper Aggregation
#'
#' This function implements the correct CellChat workflow:
#' 1. Split by individual samples (e.g., patients)
#' 2. Run CellChat on each sample independently
#' 3. Merge by condition (e.g., disease status) using mergeCellChat()
#' 4. Compare conditions using rankNet/compareInteractions
#'
#' This is analogous to pseudobulk DEG analysis - more rigorous than pooling all cells.
#'
#' @param input_data A Seurat object or path to .qs/.rds file.
#' @param group.by Character string. Cell type grouping column.
#' @param split.by Character string. Sample-level split (e.g., hos_no for patients). Each will get independent CellChat analysis.
#' @param aggregate.by Character string. Condition for aggregation/comparison (e.g., g3). CellChat objects will be merged by this. Default: NULL (no aggregation).
#' @param subset.aggregate Character vector. Specific conditions to analyze. Default: NULL (all conditions). (e.g. c("2", "1"))
#' @param species Character string. "human" or "mouse". Default: "human".
#' @param db.use Character vector. DB types to analyze separately. If multiple, runs analysis for each type.
#'        Options: "Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact". Default: NULL (uses all in one analysis).
#' @param assay_name Character string. Default: "RNA".
#' @param min.cells.group Integer. Minimum cells per cell type. Default: 10.
#' @param prob.threshold Numeric. Probability threshold for computeCommunProb. Default: 0.05.
#' @param population.size Logical. Scale by population size. Default: FALSE.
#' @param output_dir Character string. Default: NULL.
#' @param n_cores Integer. Default: 16.
#' @param do_parallel Logical. Default: TRUE.
#' @param verbose Logical. Default: TRUE.
#' @param run_name Character string. Default: NULL (auto-generated).
#'
#' @return If aggregate.by is NULL: list of CellChat objects per sample.
#'         If aggregate.by is provided: list of merged CellChat objects per condition.
#' @export
run_cellchat_analysis <- function(input_data,
                                  group.by,
                                  split.by = NULL,
                                  aggregate.by = NULL,
                                  subset.aggregate = NULL,
                                  species = c("human", "mouse"),
                                  db.use = NULL,
                                  assay_name = "RNA",
                                  min.cells.group = 10,
                                  prob.threshold = 0.05,
                                  population.size = FALSE,
                                  output_dir = NULL,
                                  n_cores = 16,
                                  do_parallel = TRUE,
                                  verbose = TRUE,
                                  run_name = NULL) {
  species <- match.arg(species)

  # --- Setup Logging ---
  log_root <- file.path(getwd(), "logs", "cc")
  if (!dir.exists(log_root)) dir.create(log_root, recursive = TRUE)

  if (is.null(run_name)) {
    existing_runs <- list.dirs(log_root, full.names = FALSE, recursive = FALSE)
    run_nums <- as.numeric(sub("run", "", grep("^run[0-9]+$", existing_runs, value = TRUE)))
    next_num <- if (length(run_nums) == 0) 1 else max(run_nums) + 1
    run_name <- paste0("run", next_num)
  }

  run_log_dir <- file.path(log_root, run_name)
  if (!dir.exists(run_log_dir)) dir.create(run_log_dir, recursive = TRUE)

  log_file <- file.path(run_log_dir, "cellchat_analysis.log")

  log_msg <- function(msg) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    formatted_msg <- paste0("[", timestamp, "] ", msg)
    if (verbose) message(formatted_msg)
    cat(formatted_msg, "\n", file = log_file, append = TRUE)
  }

  log_msg(paste0("=== CellChat Analysis (Proper Aggregation Method) ==="))
  log_msg(paste0("Run: ", run_name))

  # --- Load Data ---
  sobj <- NULL
  if (is.character(input_data)) {
    log_msg(paste0("Loading Seurat object from: ", input_data))
    if (!file.exists(input_data)) stop("Input file not found: ", input_data)

    if (grepl("\\.qs$", input_data)) {
      if (!requireNamespace("qs", quietly = TRUE)) stop("Package 'qs' needed")
      sobj <- qs::qread(input_data)
    } else if (grepl("\\.rds$", input_data)) {
      sobj <- readRDS(input_data)
    } else {
      stop("Unsupported format")
    }
  } else if (inherits(input_data, "Seurat")) {
    log_msg("Using provided Seurat object")
    sobj <- input_data
  } else {
    stop("input_data must be Seurat object or path")
  }

  # --- Setup Output Directory ---
  if (is.null(output_dir)) {
    output_dir <- file.path(getwd(), "docs", "cellchat", run_name)
  }
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  log_msg(paste0("Output directory: ", output_dir))

  # --- Handle Multiple DB Types ---
  if (!is.null(db.use) && length(db.use) > 1) {
    log_msg(paste0("Multiple DB types specified: ", paste(db.use, collapse = ", ")))
    log_msg("Running separate analysis for each DB type...")

    results_by_db <- list()
    for (db_type in db.use) {
      log_msg(paste0("\n==== Processing DB type: ", db_type, " ===="))

      db_run_name <- paste0(run_name, "_", gsub("[^A-Za-z0-9]", "", db_type))
      db_output_dir <- file.path(output_dir, make.names(db_type))

      result <- run_cellchat_analysis(
        input_data = sobj,
        group.by = group.by,
        split.by = split.by,
        aggregate.by = aggregate.by,
        subset.aggregate = subset.aggregate,
        species = species,
        db.use = db_type, # Single DB type
        assay_name = assay_name,
        min.cells.group = min.cells.group,
        prob.threshold = prob.threshold,
        population.size = population.size,
        output_dir = db_output_dir,
        n_cores = n_cores,
        do_parallel = do_parallel,
        verbose = verbose,
        run_name = db_run_name
      )

      results_by_db[[db_type]] <- result
    }

    log_msg(paste0("\n=== Completed analysis for ", length(db.use), " DB types ==="))
    return(results_by_db)
  }

  # --- Main Analysis Logic ---

  # If no split.by, do single analysis
  if (is.null(split.by)) {
    log_msg("No split.by specified - running single CellChat analysis")
    log_msg("WARNING: This pools all samples together. Consider using split.by for proper sample-level analysis.")

    result <- run_single_cellchat(
      sobj = sobj,
      group.by = group.by,
      species = species,
      db.use = db.use,
      assay_name = assay_name,
      min.cells.group = min.cells.group,
      prob.threshold = prob.threshold,
      output_dir = output_dir,
      n_cores = n_cores,
      do_parallel = do_parallel,
      log_msg = log_msg,
      run_name = run_name
    )

    return(result)
  }

  # --- Split-based Analysis (PROPER METHOD) ---
  if (!split.by %in% colnames(sobj@meta.data)) {
    stop("split.by column '", split.by, "' not found")
  }

  split_values <- unique(sobj@meta.data[[split.by]])
  split_values <- split_values[!is.na(split_values)]

  log_msg(paste0("split.by: ", split.by, " (", length(split_values), " samples)"))

  # Get aggregate mapping if provided
  aggregate_mapping <- NULL
  if (!is.null(aggregate.by)) {
    if (!aggregate.by %in% colnames(sobj@meta.data)) {
      stop("aggregate.by column '", aggregate.by, "' not found")
    }

    # Create mapping: split_value -> aggregate_value
    aggregate_mapping <- sobj@meta.data %>%
      dplyr::select(all_of(c(split.by, aggregate.by))) %>%
      dplyr::distinct() %>%
      tibble::as_tibble()

    colnames(aggregate_mapping) <- c("split_val", "aggregate_val")

    # Filter by subset.aggregate if provided
    if (!is.null(subset.aggregate)) {
      aggregate_mapping <- aggregate_mapping %>%
        dplyr::filter(aggregate_val %in% subset.aggregate)

      if (nrow(aggregate_mapping) == 0) {
        stop("No samples found for subset.aggregate values")
      }
    }

    split_values <- aggregate_mapping$split_val

    log_msg(paste0("aggregate.by: ", aggregate.by))
    log_msg(paste0("Conditions: ", paste(unique(aggregate_mapping$aggregate_val), collapse = ", ")))
    log_msg(paste0("This is the PROPER method (like pseudobulk DEG analysis)"))
  }

  log_msg(paste0("Analyzing ", length(split_values), " samples independently..."))

  # Run CellChat for each sample
  cellchat_list <- list()
  for (i in seq_along(split_values)) {
    split_val <- split_values[i]
    log_msg(paste0("\n[", i, "/", length(split_values), "] Processing: ", split_val))

    # Subset data
    sobj_subset <- sobj[, sobj@meta.data[[split.by]] == split_val]
    n_cells <- ncol(sobj_subset)

    log_msg(paste0("  Cells: ", n_cells))

    if (n_cells < 100) {
      log_msg(paste0("  SKIP: Too few cells (\u003c100)"))
      next
    }

    # Run analysis
    sample_output_dir <- file.path(output_dir, "samples", make.names(split_val))

    cc <- tryCatch(
      {
        run_single_cellchat(
          sobj = sobj_subset,
          group.by = group.by,
          species = species,
          db.use = db.use,
          assay_name = assay_name,
          min.cells.group = min.cells.group,
          prob.threshold = prob.threshold,
          output_dir = sample_output_dir,
          n_cores = n_cores,
          do_parallel = do_parallel,
          log_msg = log_msg,
          run_name = paste0(run_name, "_", make.names(split_val)),
          skip_plots = TRUE # Skip individual plots for now
        )
      },
      error = function(e) {
        log_msg(paste0("  ERROR: ", e$message))
        return(NULL)
      }
    )

    if (!is.null(cc)) {
      cellchat_list[[as.character(split_val)]] <- cc
    }
  }

  log_msg(paste0("\nCompleted ", length(cellchat_list), "/", length(split_values), " samples"))

  # --- Aggregate by Condition ---
  if (!is.null(aggregate.by) && length(cellchat_list) > 0) {
    log_msg("\n=== Aggregating by condition (mergeCellChat) ===")

    # Group by aggregate_val
    cellchat_by_condition <- list()
    for (cond in unique(aggregate_mapping$aggregate_val)) {
      samples_in_cond <- aggregate_mapping %>%
        dplyr::filter(aggregate_val == cond) %>%
        dplyr::pull(split_val) %>%
        as.character()

      # Get CellChat objects for this condition
      cc_objects <- cellchat_list[samples_in_cond]
      cc_objects <- cc_objects[!sapply(cc_objects, is.null)]

      if (length(cc_objects) == 0) {
        log_msg(paste0("  ", cond, ": No valid samples, skipping"))
        next
      }

      log_msg(paste0("  ", cond, ": Merging ", length(cc_objects), " samples"))

      # Ensure objects are updated for v2 before merging
      cc_objects <- lapply(cc_objects, CellChat::updateCellChat)

      # Merge using CellChat::mergeCellChat
      cc_merged <- CellChat::mergeCellChat(
        cc_objects,
        add.names = names(cc_objects),
        cell.prefix = TRUE
      )

      cellchat_by_condition[[cond]] <- cc_merged

      # Save merged object
      cond_output_dir <- file.path(output_dir, "merged", make.names(cond))
      if (!dir.exists(cond_output_dir)) dir.create(cond_output_dir, recursive = TRUE)

      qs_file <- file.path(cond_output_dir, "cellchat_merged.qs")
      qs::qsave(cc_merged, qs_file)
      log_msg(paste0("    Saved: ", qs_file))
    }

    log_msg("\n=== Aggregation Complete ===")
    return(list(
      by_sample = cellchat_list,
      by_condition = cellchat_by_condition,
      mapping = aggregate_mapping
    ))
  }

  # Return sample-level results if no aggregation
  return(cellchat_list)
}


#' Run Single CellChat Analysis
#'
#' Internal function for running CellChat on a single dataset
#'
#' @keywords internal
run_single_cellchat <- function(sobj, group.by, species, db.use, assay_name,
                                min.cells.group, prob.threshold, output_dir,
                                n_cores, do_parallel, log_msg, run_name,
                                skip_plots = FALSE) {
  # Check packages
  if (!requireNamespace("CellChat", quietly = TRUE)) {
    stop("CellChat package required")
  }

  # Prepare
  if (!group.by %in% colnames(sobj@meta.data)) {
    stop("group.by '", group.by, "' not found")
  }

  # Checkpoint directory
  checkpoint_dir <- file.path(output_dir, "checkpoints")
  if (!dir.exists(checkpoint_dir)) dir.create(checkpoint_dir, recursive = TRUE)

  # Helper functions
  save_ckpt <- function(obj, name) {
    file <- file.path(checkpoint_dir, paste0("cc_", name, ".qs"))
    qs::qsave(obj, file)
  }

  load_ckpt <- function(name) {
    file <- file.path(checkpoint_dir, paste0("cc_", name, ".qs"))
    if (file.exists(file) && file.size(file) > 0) {
      log_msg(paste0("    Found checkpoint: ", name))
      return(qs::qread(file))
    }
    return(NULL)
  }

  # Step 1: Create CellChat object
  cc <- load_ckpt("step1")
  if (is.null(cc)) {
    log_msg("  Creating CellChat object...")
    cc <- CellChat::createCellChat(sobj, group.by = group.by, assay = assay_name)

    # Set DB
    if (species == "human") {
      CellChatDB <- CellChat::CellChatDB.human
    } else {
      CellChatDB <- CellChat::CellChatDB.mouse
    }

    if (!is.null(db.use)) {
      if (length(db.use) == 1) {
        CellChatDB <- CellChat::subsetDB(CellChatDB, search = db.use)
        log_msg(paste0("    DB: ", db.use))
      }
    }
    cc@DB <- CellChatDB

    # Preprocess
    cc <- CellChat::subsetData(cc)

    # Dynamic Parallel Configuration
    n_cells <- ncol(cc@data.signaling)
    optimal_workers <- 1

    if (do_parallel && n_cells >= 1000) {
      if (n_cells < 5000) {
        optimal_workers <- min(4, n_cores)
      } else if (n_cells < 10000) {
        optimal_workers <- min(8, n_cores)
      } else {
        optimal_workers <- n_cores # Use user limit for large data
      }
    }

    if (optimal_workers > 1) {
      # Choose backend: multicore (Linux/Mac CLI) vs multisession (RStudio/Windows)
      is_rstudio <- Sys.getenv("RSTUDIO") == "1"
      is_windows <- .Platform$OS.type == "windows"

      backend <- "multisession"
      if (!is_rstudio && !is_windows) {
        backend <- "multicore"
      }

      log_msg(paste0("    Parallel: ", backend, " with ", optimal_workers, " workers (Cells: ", n_cells, ")"))
      future::plan(backend, workers = optimal_workers)
    } else {
      log_msg(paste0("    Parallel: OFF (Too few cells: ", n_cells, ")"))
      future::plan("sequential")
    }

    # Corrected variable usage: prob.threshold
    cc <- CellChat::identifyOverExpressedGenes(cc, thresh.p = prob.threshold)

    cc <- CellChat::identifyOverExpressedInteractions(cc)

    save_ckpt(cc, "step1")
  } else {
    log_msg("  Loaded checkpoint: step1")
    # Restore parallel plan for subsequent steps if needed
    if (do_parallel) future::plan("multisession", workers = n_cores)
    # Ensure v2 compatibility for loaded objects
    cc <- CellChat::updateCellChat(cc)
  }

  # Step 2: Communication probability
  cc_step2 <- load_ckpt("step2")
  if (is.null(cc_step2)) {
    log_msg(paste0("  Computing probabilities..."))
    cc <- CellChat::computeCommunProb(cc, type = "triMean")
    cc <- CellChat::filterCommunication(cc, min.cells = min.cells.group)
    save_ckpt(cc, "step2")
  } else {
    cc <- cc_step2
    log_msg("  Loaded checkpoint: step2")
    cc <- CellChat::updateCellChat(cc)
  }

  # Step 3: Pathways
  cc_step3 <- load_ckpt("step3")
  if (is.null(cc_step3)) {
    log_msg("  Computing pathways...")
    cc <- CellChat::computeCommunProbPathway(cc)

    # Check if aggregateNet is still valid in v2 or needs args
    cc <- CellChat::aggregateNet(cc)
    save_ckpt(cc, "step3")
  } else {
    cc <- cc_step3
    log_msg("  Loaded checkpoint: step3")
    cc <- CellChat::updateCellChat(cc)
  }

  # Save final
  if (!skip_plots) {
    log_msg("  Generating plots...")

    # Update for plotting if needed
    cc <- CellChat::updateCellChat(cc)

    # Circle plot
    pdf(file.path(output_dir, "net_circle.pdf"), width = 10, height = 10)
    groupSize <- as.numeric(table(cc@idents))
    CellChat::netVisual_circle(cc@net$count,
      vertex.weight = groupSize,
      weight.scale = TRUE, label.edge = FALSE,
      title.name = "Number of interactions"
    )
    dev.off()

    # Bubble plot
    if (length(cc@netP$pathways) > 0) {
      pdf(file.path(output_dir, "net_bubble.pdf"), width = 12, height = 8)
      tryCatch(
        {
          print(CellChat::netVisual_bubble(cc, remove.isolate = FALSE))
        },
        error = function(e) log_msg(paste0("    Bubble plot error: ", e$message))
      )
      dev.off()
    }
  }

  # Save object
  qs::qsave(cc, file.path(output_dir, "cellchat.qs"))

  return(cc)
}

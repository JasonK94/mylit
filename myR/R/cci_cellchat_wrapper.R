#' Run CellChat Analysis
#'
#' This function wraps the CellChat analysis pipeline.
#' It takes a Seurat object (or path) and runs the standard CellChat workflow.
#'
#' @param input_data A Seurat object or a character string path to a .qs/.rds file containing a Seurat object.
#' @param group.by Character string. The metadata column to use for cell labels (e.g., cell types).
#' @param species Character string. "human" or "mouse". Default: "human".
#' @param db.use Character vector. Subset of CellChatDB to use.
#'        Options: "Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact".
#'        If NULL, uses the entire database. Default: NULL.
#' @param assay_name Character string. Assay to use. Default: "RNA".
#' @param split.by Character string. Metadata column to split analysis by (e.g., patient ID, condition). If provided, runs separate CellChat for each group. Default: NULL.
#' @param subset.split Character vector. Specific groups from split.by to analyze. If NULL, analyzes all groups. Default: NULL.
#' @param min.cells.group Integer. Minimum number of cells per group.by category within each split group. Default: 10.
#' @param population.size Logical. Whether to scale communication probability by cell population sizes. Default: FALSE.
#' @param output_dir Character string. Path to save results. Default: NULL.
#' @param n_cores Integer. Number of cores for parallel processing. Default: 16.
#' @param do_parallel Logical. Whether to use parallel processing. Default: TRUE.
#' @param verbose Logical. Default: TRUE.
#' @param run_name Character string. Name of the run for logging purposes. If NULL, auto-generated.
#'
#' @return A CellChat object with results.
#' @export
#' @importFrom CellChat createCellChat setIdent subsetData identifyOverExpressedGenes identifyOverExpressedInteractions computeCommunProb filterCommunication computeCommunProbPathway aggregateNet netVisual_circle
#' @importFrom Seurat DefaultAssay
#' @importFrom future plan multiprocess
#' @importFrom qs qread
run_cellchat_analysis <- function(input_data,
                                  group.by,
                                  species = c("human", "mouse"),
                                  db.use = NULL,
                                  assay_name = "RNA",
                                  split.by = NULL,
                                  subset.split = NULL,
                                  min.cells.group = 10,
                                  population.size = FALSE,
                                  output_dir = NULL,
                                  n_cores = 16,
                                  do_parallel = TRUE,
                                  verbose = TRUE,
                                  run_name = NULL) {
  species <- match.arg(species)

  # --- Setup Logging and Directories ---
  # Determine run name and log directory
  log_root <- file.path(getwd(), "logs", "cc")
  if (!dir.exists(log_root)) dir.create(log_root, recursive = TRUE)

  if (is.null(run_name)) {
    # Auto-generate run name: run1, run2, ...
    existing_runs <- list.dirs(log_root, full.names = FALSE, recursive = FALSE)
    run_nums <- as.numeric(sub("run", "", grep("^run[0-9]+$", existing_runs, value = TRUE)))
    next_num <- if (length(run_nums) == 0) 1 else max(run_nums) + 1
    run_name <- paste0("run", next_num)
  }

  run_log_dir <- file.path(log_root, run_name)
  if (!dir.exists(run_log_dir)) dir.create(run_log_dir, recursive = TRUE)

  # Setup log file
  log_file <- file.path(run_log_dir, "cellchat_analysis.log")

  # Helper for logging
  log_msg <- function(msg) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    formatted_msg <- paste0("[", timestamp, "] ", msg)
    if (verbose) message(formatted_msg)
    cat(formatted_msg, "\n", file = log_file, append = TRUE)
  }

  log_msg(paste0("=== Starting CellChat Analysis (Run: ", run_name, ") ==="))
  log_msg(paste0("Log file: ", log_file))

  # --- Load Data ---
  sobj <- NULL
  if (is.character(input_data)) {
    log_msg(paste0("Loading Seurat object from: ", input_data))
    if (!file.exists(input_data)) stop("Input file not found: ", input_data)

    if (grepl("\\.qs$", input_data)) {
      if (!requireNamespace("qs", quietly = TRUE)) stop("Package 'qs' needed for .qs files.")
      sobj <- qs::qread(input_data)
    } else if (grepl("\\.rds$", input_data)) {
      sobj <- readRDS(input_data)
    } else {
      stop("Unsupported file format. Use .qs or .rds")
    }
  } else if (inherits(input_data, "Seurat")) {
    log_msg("Using provided Seurat object.")
    sobj <- input_data
  } else {
    stop("input_data must be a Seurat object or a path to one.")
  }

  # --- Handle Sample Splitting ---
  if (!is.null(split.by)) {
    if (!split.by %in% colnames(sobj@meta.data)) {
      stop("split.by column '", split.by, "' not found in metadata.")
    }

    split_groups <- unique(sobj@meta.data[[split.by]])
    if (!is.null(subset.split)) {
      split_groups <- intersect(split_groups, subset.split)
      if (length(split_groups) == 0) {
        stop("No valid split groups found. Check subset.split values.")
      }
    }

    log_msg(paste0("Splitting analysis by '", split.by, "' into ", length(split_groups), " groups: ", paste(split_groups, collapse = ", ")))

    # Run analysis for each split group
    results_list <- list()
    for (split_val in split_groups) {
      log_msg(paste0("\n=== Processing split group: ", split_val, " ==="))

      # Subset data
      sobj_subset <- sobj[, sobj@meta.data[[split.by]] == split_val]
      log_msg(paste0("Subset contains ", ncol(sobj_subset), " cells"))

      # Run analysis with modified run_name and output_dir
      split_run_name <- paste0(run_name, "_", make.names(split_val))
      split_output_dir <- if (!is.null(output_dir)) file.path(output_dir, make.names(split_val)) else NULL

      # Recursively call with split.by = NULL to prevent infinite recursion
      result <- run_cellchat_analysis(
        input_data = sobj_subset,
        group.by = group.by,
        species = species,
        db.use = db.use,
        assay_name = assay_name,
        split.by = NULL, # Important: prevent recursion
        min.cells.group = min.cells.group,
        population.size = population.size,
        output_dir = split_output_dir,
        n_cores = n_cores,
        do_parallel = do_parallel,
        verbose = verbose,
        run_name = split_run_name
      )

      results_list[[split_val]] <- result
    }

    log_msg(paste0("\n=== Split analysis complete. Processed ", length(results_list), " groups ==="))
    return(results_list)
  }

  # --- Setup Output Directory ---
  if (is.null(output_dir)) {
    output_dir <- file.path(getwd(), "docs", "cellchat", run_name)
  }
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  log_msg(paste0("Output directory: ", output_dir))

  # --- Checkpoint Helper ---
  # Saves intermediate result to output_dir/checkpoints/
  checkpoint_dir <- file.path(output_dir, "checkpoints")
  if (!dir.exists(checkpoint_dir)) dir.create(checkpoint_dir, recursive = TRUE)

  save_checkpoint <- function(obj, step_name) {
    file_path_qs <- file.path(checkpoint_dir, paste0("cellchat_", step_name, ".qs"))
    log_msg(paste0("Saving checkpoint: ", step_name, " to ", file_path_qs))
    if (requireNamespace("qs", quietly = TRUE)) {
      qs::qsave(obj, file_path_qs)
    } else {
      file_path_rds <- file.path(checkpoint_dir, paste0("cellchat_", step_name, ".rds"))
      log_msg(paste0("qs package not found. Saving checkpoint to ", file_path_rds))
      saveRDS(obj, file_path_rds)
    }
  }

  load_checkpoint <- function(step_name) {
    file_path_qs <- file.path(checkpoint_dir, paste0("cellchat_", step_name, ".qs"))
    file_path_rds <- file.path(checkpoint_dir, paste0("cellchat_", step_name, ".rds"))

    if (file.exists(file_path_qs)) {
      log_msg(paste0("Loading checkpoint: ", step_name, " from ", file_path_qs))
      if (requireNamespace("qs", quietly = TRUE)) {
        return(qs::qread(file_path_qs))
      } else {
        log_msg("qs package not found but checkpoint is .qs. Cannot load.")
        return(NULL)
      }
    } else if (file.exists(file_path_rds)) {
      log_msg(paste0("Loading checkpoint: ", step_name, " from ", file_path_rds))
      return(readRDS(file_path_rds))
    }
    return(NULL)
  }

  # --- Analysis Steps ---

  # Check packages
  if (!requireNamespace("CellChat", quietly = TRUE)) {
    stop("Package 'CellChat' is required but not installed.")
  }

  # Prepare Seurat object
  log_msg("Preparing Seurat object...")
  if (!group.by %in% colnames(sobj@meta.data)) {
    stop("Column '", group.by, "' not found in metadata.")
  }

  # Step 1: Create CellChat Object
  # Check if we can skip this
  cellchat <- load_checkpoint("step1_created")

  if (is.null(cellchat)) {
    log_msg("Creating CellChat object...")
    cellchat <- CellChat::createCellChat(object = sobj, group.by = group.by, assay = assay_name)

    # Set database
    log_msg(paste0("Setting CellChatDB for ", species, "..."))
    if (species == "human") {
      CellChatDB <- CellChat::CellChatDB.human
    } else {
      CellChatDB <- CellChat::CellChatDB.mouse
    }

    if (!is.null(db.use)) {
      CellChatDB.use <- CellChat::subsetDB(CellChatDB, search = db.use)
      log_msg(paste0("Using subset of DB: ", paste(db.use, collapse = ", ")))
    } else {
      CellChatDB.use <- CellChatDB
      log_msg("Using full DB")
    }
    cellchat@DB <- CellChatDB.use

    # Pre-processing
    log_msg("Preprocessing expression data (subsetData)...")
    cellchat <- CellChat::subsetData(cellchat) # Subset to genes in DB

    if (do_parallel) {
      log_msg(paste0("Setting up parallel processing with ", n_cores, " cores..."))
      future::plan("multisession", workers = n_cores)
    }

    log_msg("Identifying over-expressed genes...")
    cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
    log_msg("Identifying over-expressed interactions...")
    cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)

    save_checkpoint(cellchat, "step1_created")
  } else {
    log_msg("Skipped Step 1 (Loaded from checkpoint).")
    # Re-init parallel if needed (checkpoints don't save future plan)
    if (do_parallel) {
      log_msg(paste0("Setting up parallel processing with ", n_cores, " cores..."))
      future::plan("multisession", workers = n_cores)
    }
  }

  # Step 2: Compute Communication Probability
  cellchat_step2 <- load_checkpoint("step2_prob")
  if (is.null(cellchat_step2)) {
    log_msg("Computing communication probabilities (computeCommunProb)... This may take a while.")
    cellchat <- CellChat::computeCommunProb(cellchat, type = "triMean")

    log_msg(paste0("Filtering communication (filterCommunication) with min.cells = ", min.cells.group, "..."))
    cellchat <- CellChat::filterCommunication(cellchat, min.cells = min.cells.group)

    save_checkpoint(cellchat, "step2_prob")
  } else {
    cellchat <- cellchat_step2
    log_msg("Skipped Step 2 (Loaded from checkpoint).")
  }

  # Step 3: Pathway Probabilities & Aggregation
  cellchat_step3 <- load_checkpoint("step3_pathway")
  if (is.null(cellchat_step3)) {
    log_msg("Computing pathway probabilities (computeCommunProbPathway)...")
    cellchat <- CellChat::computeCommunProbPathway(cellchat)

    log_msg("Aggregating network (aggregateNet)...")
    cellchat <- CellChat::aggregateNet(cellchat)

    save_checkpoint(cellchat, "step3_pathway")
  } else {
    cellchat <- cellchat_step3
    log_msg("Skipped Step 3 (Loaded from checkpoint).")
  }

  # --- Visualization and Saving ---
  log_msg("Generating and saving plots...")

  # 1. Circle plot
  pdf_file <- file.path(output_dir, "net_visual_circle.pdf")
  log_msg(paste0("Saving circle plot to ", pdf_file))
  pdf(pdf_file, width = 10, height = 10)
  groupSize <- as.numeric(table(cellchat@idents))
  CellChat::netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = "Number of interactions")
  CellChat::netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = "Interaction weights/strength")
  dev.off()

  # 2. Bubble plot
  if (length(cellchat@netP$pathways) > 0) {
    pdf_file_bubble <- file.path(output_dir, "net_visual_bubble.pdf")
    log_msg(paste0("Saving bubble plot to ", pdf_file_bubble))
    pdf(pdf_file_bubble, width = 12, height = 8)
    tryCatch(
      {
        print(CellChat::netVisual_bubble(cellchat, sources.use = NULL, targets.use = NULL, remove.isolate = FALSE))
      },
      error = function(e) log_msg(paste0("Error in bubble plot: ", e$message))
    )
    dev.off()
  } else {
    log_msg("No pathways found, skipping bubble plot.")
  }

  # Save Final Object
  qs_file <- file.path(output_dir, "cellchat_object.qs")
  log_msg(paste0("Saving final CellChat object to ", qs_file))
  if (requireNamespace("qs", quietly = TRUE)) {
    qs::qsave(cellchat, qs_file)
  } else {
    rds_file <- file.path(output_dir, "cellchat_object.rds")
    log_msg(paste0("qs package not found. Saving to ", rds_file))
    saveRDS(cellchat, rds_file)
  }

  log_msg(paste0("=== CellChat Analysis Complete (Run: ", run_name, ") ==="))
  return(cellchat)
}

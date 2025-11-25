
#' @export
ligand_to_target=function(ligand,target,NN_data=NicheNetData,ligand_tf_matrix=NULL,weighted_networks=NULL){
  if(is.null(ligand_tf_matrix)){
    ligand_tf_matrix=NicheNetData$ligand_tf_matrix
  }
  if(is.null(weighted_networks)){
    weighted_networks=NicheNetData$weighted_networks
  }
  
  
  active_signaling_network = get_ligand_signaling_path (
    ligand_tf_matrix = ligand_tf_matrix, ligands_all = ligand, 
    targets_all = target, weighted_networks = weighted_networks)
  
  # For better visualization of edge weigths: normalize edge weights to 
  # make them comparable between signaling and gene regulatory interactions
  active_signaling_network_min_max = active_signaling_network
  
  active_signaling_network_min_max$sig = active_signaling_network_min_max$sig %>% 
    mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
  
  active_signaling_network_min_max$gr = active_signaling_network_min_max$gr %>% 
    mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
  
  graph_min_max = diagrammer_format_signaling_graph (
    signaling_graph_list = active_signaling_network_min_max, 
    ligands_all = ligand, targets_all = target, 
    sig_color = "indianred", gr_color = "steelblue")
  
  # To render the graph: uncomment following line of code
  plot=DiagrammeR::render_graph(graph_min_max, layout = "tree")
  
  return(plot)
}


#' Perform NicheNet Cell-Cell Communication Analysis
#'
#' @description
#' This function wraps the NicheNet analysis pipeline to identify key ligands,
#' their target genes, and receptors involved in cell-cell communication between
#' defined sender and receiver cell populations. It handles loading of NicheNet's
#' prior interaction data, with options for local caching, global environment storage,
#' and downloading from Zenodo if necessary.
#'
#' @param seurat_obj A Seurat object containing the single-cell expression data.
#' @param species Character string, either "human" or "mouse", to load the appropriate NicheNet databases.
#' @param sender_celltypes Character vector of cell types acting as senders.
#' @param receiver_celltype Character string of the cell type acting as the receiver.
#' @param assay_name Character string, the name of the assay in the Seurat object to use (e.g., "SCT", "RNA"). Default: "SCT".
#' @param cluster_col Character string, the name of the metadata column in `seurat_obj` that contains cell type or cluster annotations.
#'        This will be used to set `Idents(seurat_obj)`.
#' @param receiver_DE_ident1 Character string, the identity for `ident.1` in `FindMarkers` for DE analysis within the `receiver_celltype`.
#' @param receiver_DE_ident2 Character string, the identity for `ident.2` in `FindMarkers` (e.g., control group). Can be NULL if comparing against all other cells in `receiver_DE_group_by`.
#' @param receiver_DE_group_by Character string, the metadata column used for `group.by` in `FindMarkers` for DE analysis within the `receiver_celltype`.
#' @param min_pct_expressed Numeric, minimum fraction of cells in which a gene must be expressed in sender/receiver groups for initial filtering. Default: 0.10.
#' @param p_val_adj_cutoff Numeric, adjusted p-value cutoff for DEGs. Default: 0.05.
#' @param logfc_cutoff Numeric, log fold-change cutoff for DEGs (positive values for upregulation in ident1). Default: 0.25.
#' @param top_n_ligands Integer, number of top ligands to consider for detailed analysis and visualization. Default: 20.
#' @param top_n_targets_per_ligand Integer, number of top target genes to show per ligand in visualizations. Default: 200.
#' @param ligand_target_cutoff Numeric, cutoff for visualizing ligand-target links (e.g., based on regulatory potential). Default: 0.33.
#' @param nichenet_data_dir Character string, path to the directory for local caching of NicheNet's RDS data files.
#'        Used if data is not found in the global environment object specified by `nichenet_data_name`.
#'        If NULL, a species-specific default path within the working directory will be used (e.g., "./nichenet_data_human/").
#' @param nichenet_data_name Character string or NULL. The name of the R object in the Global Environment
#'        where NicheNet datasets are stored or will be stored. If NULL, data is not stored/retrieved globally by this function.
#'        Default: "NicheNetData".
#' @param output_dir Character string, path to a main directory where run-specific subfolders for outputs will be created. If NULL, plots are returned but not saved to files. Default: NULL.
#' @param run_circos Logical, whether to attempt to generate and return a circos plot of ligand-receptor interactions. Requires `circlize`, `RColorBrewer`, and `evaluate` packages. Default: TRUE.
#' @param run_signaling_path_inference Logical, whether to infer and visualize ligand-to-target signaling paths for the top 1 ligand. Requires `DiagrammeR`, `DiagrammeRsvg`, `rsvg`. Default: FALSE.
#' @param verbose Logical, whether to print progress messages and summary of warnings. Default: TRUE.
#'
#' @return A list containing:
#'  \itemize{
#'    \item `ligand_activities`: DataFrame of ligand activities.
#'    \item `best_upstream_ligands`: Character vector of top N ligands.
#'    \item `active_ligand_target_links_df`: DataFrame of active ligand-target links.
#'    \item `plot_ligand_target_network`: ggplot object of the ligand-target heatmap.
#'    \item `ligand_receptor_network_df`: DataFrame of prioritized ligand-receptor interactions.
#'    \item `plot_ligand_receptor_network`: ggplot object of the ligand-receptor heatmap.
#'    \item `plot_ligand_activity_hist`: ggplot object of ligand activity histogram.
#'    \item `plot_ligand_aupr_heatmap`: ggplot object of top ligand AUPR heatmap.
#'    \item `DE_table_receiver`: DataFrame of differentially expressed genes in the receiver.
#'    \item `plot_circos`: A recorded circos plot object (replay with `plot()`) or a status message.
#'    \item `circos_notes`: Character vector of "Note:" messages captured during circos plot generation.
#'    \item `signaling_path_graph`: A DiagrammeR graph object or NULL.
#'    \item `nichenet_data_loaded_info`: List describing how NicheNet data was loaded.
#'    \item `output_path_run`: Path to the run-specific output directory if `output_dir` was provided.
#'    \item `warnings_occurred`: A character vector of warning messages collected during execution, if any.
#'  }
#' @details
#' To visualize ligand expression (e.g., as a DotPlot or Heatmap), you can use the returned
#' `best_upstream_ligands` and your original `seurat_obj`. For example:
#' `Seurat::DotPlot(seurat_obj, features = results$best_upstream_ligands, group.by = "your_celltype_column")`
#' For a more integrated graphical legend with the Circos plot, consider using the `ComplexHeatmap` package's `Legend` function externally.
#'
#' @importFrom Seurat DefaultAssay Idents FindMarkers GetAssayData Cells
#' @importFrom dplyr %>% filter select distinct arrange desc top_n pull inner_join mutate group_by summarise across ends_with bind_rows rename
#' @importFrom tidyr spread drop_na
#' @importFrom tibble column_to_rownames rownames_to_column as_tibble tibble
#' @importFrom magrittr set_rownames
#' @importFrom nichenetr get_expressed_genes predict_ligand_activities get_weighted_ligand_target_links prepare_ligand_target_visualization make_heatmap_ggplot get_ligand_signaling_path diagrammer_format_signaling_graph
#' @importFrom ggplot2 ggplot aes geom_histogram geom_vline labs theme_classic theme element_text scale_fill_gradient2 ggsave
#' @importFrom utils download.file capture.output
#' @importFrom stats dist hclust
#' @importFrom grDevices recordPlot dev.off pdf png
#' @importFrom DiagrammeRsvg export_svg
#' @importFrom rsvg rsvg_png
#' @importFrom RColorBrewer brewer.pal
#' @importFrom circlize circos.clear chordDiagram circos.track circos.text uy mm_h CELL_META
#' @importFrom graphics legend
#' @importFrom evaluate evaluate
#' @export
run_nichenet_analysis <- function(seurat_obj,
                                  species = c("human", "mouse"),
                                  sender_celltypes,
                                  receiver_celltype,
                                  assay_name = "SCT",
                                  cluster_col,
                                  receiver_DE_ident1,
                                  receiver_DE_ident2 = NULL,
                                  receiver_DE_group_by,
                                  receiver_de_table = NULL,
                                  receiver_gene_col = "gene",
                                  receiver_logfc_col = NULL,
                                  receiver_pval_col = NULL,
                                  min_pct_expressed = 0.10,
                                  p_val_adj_cutoff = 0.05,
                                  logfc_cutoff = 0.25,
                                  top_n_ligands = 20,
                                  top_n_targets_per_ligand = 200,
                                  ligand_target_cutoff = 0.33,
                                  nichenet_data_dir = NULL,
                                  nichenet_data_name = "NicheNetData",
                                  output_dir = NULL,
                                  run_circos = TRUE,
                                  ...,
                                  verbose = TRUE) {
  args_option=list(...)
  # --- 0. Initialize Warning Collection & Argument Setup ---
  warnings_collected <- character(0)
  circos_notes_collected <- character(0) # For circos specific notes
  
  # Capture function arguments for saving (before withCallingHandlers)
  function_args_list <- as.list(match.call(expand.dots = FALSE))[-1]  # Remove function name
  # Get actual argument values from current environment
  captured_args <- list()
  for (arg_name in names(function_args_list)) {
    if (arg_name == "...") {
      captured_args[["..."]] <- args_option
    } else {
      # Try to evaluate the argument
      tryCatch({
        captured_args[[arg_name]] <- eval(function_args_list[[arg_name]], envir = parent.frame())
      }, error = function(e) {
        # If evaluation fails, store the expression as string
        captured_args[[arg_name]] <- deparse(function_args_list[[arg_name]])
      })
    }
  }
  
  analysis_results <- withCallingHandlers({
    # Capture call early for argument saving
    call_obj <- match.call(expand.dots = TRUE)
    species <- match.arg(species)
    
    required_pkgs <- c("nichenetr", "Seurat", "dplyr", "tidyr", "tibble", "ggplot2", "rlang")
    missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
      stop("Required package(s) not found: ", paste(missing_pkgs, collapse=", "), ". Please install them.", call. = FALSE)
    }
    
    current_run_output_dir <- NULL
    if (!is.null(output_dir)) {
      if (!dir.exists(output_dir)) {
        if (verbose) message("Main output directory '", output_dir, "' does not exist. Creating it.")
        tryCatch(dir.create(output_dir, recursive = TRUE, showWarnings = FALSE),
                 error = function(e) stop("Failed to create main output directory: ", output_dir, ". Error: ", e$message))
      }
      existing_runs <- list.dirs(output_dir, full.names = FALSE, recursive = FALSE)
      run_numbers <- grep("^run[0-9]+$", existing_runs, value = TRUE) %>%
        sub("run", "", .) %>%
        as.numeric()
      next_run_number <- if (length(run_numbers) == 0 || all(is.na(run_numbers))) 1 else max(run_numbers, na.rm = TRUE) + 1
      run_folder_name <- paste0("run", next_run_number)
      current_run_output_dir <- file.path(output_dir, run_folder_name)
      if (!dir.exists(current_run_output_dir)) {
        if (verbose) message("Creating run-specific output directory: ", current_run_output_dir)
        tryCatch(dir.create(current_run_output_dir, recursive = TRUE, showWarnings = FALSE),
                 error = function(e) stop("Failed to create run-specific output directory: ", current_run_output_dir, ". Error: ", e$message))
      }
      args_file_path <- file.path(current_run_output_dir, paste0(run_folder_name, "_args.txt"))
      tryCatch({
        # Use captured arguments from before withCallingHandlers
        args_summary <- paste0(
          "# NicheNet Analysis Arguments\n",
          "# Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
          "# Function: run_nichenet_analysis\n\n",
          "# Arguments:\n"
        )
        # Format each captured argument
        for (arg_name in names(captured_args)) {
          if (arg_name == "...") {
            if (length(captured_args[["..."]]) > 0) {
              args_summary <- paste0(args_summary, "\n# Additional arguments (...):\n")
              for (opt_name in names(captured_args[["..."]])) {
                arg_val <- captured_args[["..."]][[opt_name]]
                if (is.character(arg_val) && length(arg_val) == 1) {
                  args_summary <- paste0(args_summary, opt_name, " = \"", arg_val, "\"\n")
                } else if (is.numeric(arg_val) && length(arg_val) == 1) {
                  args_summary <- paste0(args_summary, opt_name, " = ", arg_val, "\n")
                } else {
                  args_summary <- paste0(args_summary, opt_name, " = ", paste(utils::capture.output(print(arg_val)), collapse = " "), "\n")
                }
              }
            }
            next
          }
          arg_val <- captured_args[[arg_name]]
          
          # Format based on type
          if (is.character(arg_val) && length(arg_val) == 1) {
            args_summary <- paste0(args_summary, arg_name, " = \"", arg_val, "\"\n")
          } else if (is.character(arg_val) && length(arg_val) > 1) {
            args_summary <- paste0(args_summary, arg_name, " = c(\"", paste(arg_val, collapse = "\", \""), "\")\n")
          } else if (is.numeric(arg_val) && length(arg_val) == 1) {
            args_summary <- paste0(args_summary, arg_name, " = ", arg_val, "\n")
          } else if (is.numeric(arg_val) && length(arg_val) > 1) {
            args_summary <- paste0(args_summary, arg_name, " = c(", paste(arg_val, collapse = ", "), ")\n")
          } else if (is.logical(arg_val) && length(arg_val) == 1) {
            args_summary <- paste0(args_summary, arg_name, " = ", arg_val, "\n")
          } else if (is.null(arg_val)) {
            args_summary <- paste0(args_summary, arg_name, " = NULL\n")
          } else if (inherits(arg_val, "Seurat")) {
            args_summary <- paste0(args_summary, arg_name, " = <Seurat object: ", ncol(arg_val), " cells, ", nrow(arg_val), " features>\n")
          } else if (is.data.frame(arg_val)) {
            args_summary <- paste0(args_summary, arg_name, " = <data.frame: ", nrow(arg_val), " rows, ", ncol(arg_val), " cols>\n")
          } else if (is.character(arg_val) && length(arg_val) == 1 && grepl("^<", arg_val)) {
            # Already a string representation
            args_summary <- paste0(args_summary, arg_name, " = ", arg_val, "\n")
          } else {
            # Fallback: try to capture output
            arg_str <- tryCatch({
              paste(utils::capture.output(print(arg_val)), collapse = " ")
            }, error = function(e) {
              paste0("<", class(arg_val)[1], ">")
            })
            if (nchar(arg_str) > 200) {
              arg_str <- paste0(substr(arg_str, 1, 197), "...")
            }
            args_summary <- paste0(args_summary, arg_name, " = ", arg_str, "\n")
          }
        }
        writeLines(args_summary, args_file_path)
      }, error = function(e) {
        warning("NN_Warning: Failed to save arguments to file: ", args_file_path, ". Error: ", e$message)
        # Fallback: try simple deparse of call
        tryCatch({
          call_str <- deparse(call_obj, width.cutoff = 500)
          writeLines(c("# NicheNet Analysis Arguments (fallback)", 
                       paste("# Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
                       "",
                       paste(call_str, collapse = "\n")), 
                     args_file_path)
        }, error = function(e2) {
          warning("NN_Warning: Fallback save also failed: ", e2$message)
        })
      })
      if (verbose && file.exists(args_file_path)) {
        file_size <- file.info(args_file_path)$size
        if (file_size > 0) {
          message("Arguments saved to: ", args_file_path, " (", file_size, " bytes)")
        } else {
          warning("NN_Warning: Arguments file created but is empty: ", args_file_path)
        }
      }
    } else {
      if (verbose) message("`output_dir` is NULL. Plots will be returned as objects but not saved to files.")
    }
    
    # --- I. NicheNet Data Loading ---
    if (verbose) message("Attempting to load NicheNet prior model data for ", species, "...")
    base_url <- "https://zenodo.org/records/7074291/files/"
    files_to_load_spec <- list(
      human = list(lr_network = "lr_network_human_21122021.rds", sig_network = "signaling_network_human_21122021.rds", gr_network = "gr_network_human_21122021.rds", weighted_networks = "weighted_networks_nsga2r_final.rds", ligand_target_matrix = "ligand_target_matrix_nsga2r_final.rds", ligand_tf_matrix = "ligand_tf_matrix_nsga2r_final.rds"),
      mouse = list(lr_network = "lr_network_mouse_21122021.rds", sig_network = "signaling_network_mouse_21122021.rds", gr_network = "gr_network_mouse_21122021.rds", weighted_networks = "weighted_networks_nsga2r_final_mouse.rds", ligand_target_matrix = "ligand_target_matrix_nsga2r_final_mouse.rds", ligand_tf_matrix = "ligand_tf_matrix_nsga2r_final_mouse.rds")
    )
    required_files_map <- files_to_load_spec[[species]]
    nichenet_data_loaded_from_global <- FALSE
    temp_data_env <- new.env()
    
    if (!is.null(nichenet_data_name) && nzchar(nichenet_data_name)) {
      if (exists(nichenet_data_name, envir = .GlobalEnv)) {
        if (verbose) message("Found '", nichenet_data_name, "' in global environment. Validating...")
        global_data_env <- get(nichenet_data_name, envir = .GlobalEnv)
        if (is.list(global_data_env) && all(names(required_files_map) %in% names(global_data_env))) {
          are_all_components_valid <- all(sapply(names(required_files_map), function(comp_name) !is.null(global_data_env[[comp_name]])))
          has_wn_lr <- ("weighted_networks_lr" %in% names(global_data_env) && !is.null(global_data_env[["weighted_networks_lr"]])) ||
            (!is.null(global_data_env[["weighted_networks"]]) && !is.null(global_data_env[["lr_network"]]))
          if (are_all_components_valid && has_wn_lr) {
            if (verbose) message("Global data '", nichenet_data_name, "' is valid and complete. Using it.")
            for (comp_name in names(required_files_map)) {
              assign(comp_name, global_data_env[[comp_name]], envir = temp_data_env)
            }
            if ("weighted_networks_lr" %in% names(global_data_env) && !is.null(global_data_env[["weighted_networks_lr"]])) {
              assign("weighted_networks_lr", global_data_env[["weighted_networks_lr"]], envir = temp_data_env)
            } else {
              temp_lr_network_global <- global_data_env[["lr_network"]] %>% dplyr::distinct(from, to)
              assign("weighted_networks_lr", global_data_env[["weighted_networks"]]$lr_sig %>%
                       dplyr::inner_join(temp_lr_network_global, by = c("from", "to")), envir = temp_data_env)
            }
            nichenet_data_loaded_from_global <- TRUE
          } else {
            if (verbose) message("Global data '", nichenet_data_name, "' is incomplete. Will load from files/download.")
          }
        } else {
          if (verbose) message("Global data '", nichenet_data_name, "' is not a valid list or missing key components. Will load from files/download.")
        }
      } else {
        if (verbose) message("'", nichenet_data_name, "' not found in global environment. Will load from files/download.")
      }
    } else {
      if (verbose) message("`nichenet_data_name` not provided or empty. Data will not be managed via global environment for loading.")
    }
    
    if (!nichenet_data_loaded_from_global) {
      data_dir_to_use <- nichenet_data_dir
      if (is.null(data_dir_to_use)) {
        # Try default paths: first check /data/user3/git_repo/human or /data/user3/git_repo/mouse
        default_paths <- c(
          file.path("/data/user3/git_repo", species),
          file.path(getwd(), paste0("nichenet_data_", species))
        )
        data_dir_to_use <- NULL
        for (path in default_paths) {
          if (dir.exists(path)) {
            # Check if at least one required file exists
            test_file <- file.path(path, required_files_map[[1]])
            if (file.exists(test_file)) {
              data_dir_to_use <- path
              if (verbose) message("Found NicheNet data directory: ", data_dir_to_use)
              break
            }
          }
        }
        if (is.null(data_dir_to_use)) {
          data_dir_to_use <- file.path(getwd(), paste0("nichenet_data_", species))
          if (verbose) message("`nichenet_data_dir` is NULL, using default local cache: ", data_dir_to_use)
        }
      }
      if (!dir.exists(data_dir_to_use)) {
        if (verbose) message("Creating NicheNet data directory: ", data_dir_to_use)
        tryCatch(dir.create(data_dir_to_use, recursive = TRUE, showWarnings = FALSE),
                 error = function(e) stop("Failed to create NicheNet data directory: ", data_dir_to_use, ". Error: ", e$message))
      }
      loaded_data_for_global_assignment <- list()
      n_files <- length(required_files_map)
      if (verbose) message("Loading ", n_files, " NicheNet data files from: ", data_dir_to_use)
      for (i in seq_along(required_files_map)) {
        simple_name <- names(required_files_map)[i]
        rds_filename <- required_files_map[[simple_name]]
        file_path <- file.path(data_dir_to_use, rds_filename)
        qs_file_path <- file.path(data_dir_to_use, sub("\\.rds$", ".qs", rds_filename))
        
        # Progress indicator
        if (verbose) {
          pct <- round(i / n_files * 100, 1)
          bar_length <- 30
          filled <- round(bar_length * i / n_files)
          bar <- paste0(rep("=", filled), collapse = "")
          spaces <- paste0(rep(" ", bar_length - filled), collapse = "")
          message("  [", bar, spaces, "] ", pct, "% (", i, "/", n_files, ") Loading ", rds_filename)
        }
        
        # Try .qs file first (faster), then .rds
        data_loaded <- FALSE
        if (file.exists(qs_file_path) && requireNamespace("qs", quietly = TRUE)) {
          tryCatch({
            assign(simple_name, qs::qread(qs_file_path), envir = temp_data_env)
            loaded_data_for_global_assignment[[simple_name]] <- get(simple_name, envir = temp_data_env)
            data_loaded <- TRUE
            if (verbose) message("    Loaded from .qs file (faster)")
          }, error = function(e) {
            if (verbose) message("    Failed to load .qs, trying .rds: ", e$message)
          })
        }
        
        if (!data_loaded) {
          if (!file.exists(file_path)) {
            if (verbose) message("    File ", rds_filename, " not found. Downloading from Zenodo...")
            download_url <- paste0(base_url, rds_filename, "?download=1")
            tryCatch({
              utils::download.file(url = download_url, destfile = file_path, mode = "wb", quiet = !verbose)
              if (verbose) message("    Downloaded ", rds_filename, " successfully")
            }, error = function(e) {
              stop("Failed to download ", rds_filename, ". Error: ", e$message, call. = FALSE)
            })
          }
          assign(simple_name, readRDS(file_path), envir = temp_data_env)
          loaded_data_for_global_assignment[[simple_name]] <- get(simple_name, envir = temp_data_env)
          
          # Optionally save as .qs for faster future loading
          if (requireNamespace("qs", quietly = TRUE) && !file.exists(qs_file_path)) {
            tryCatch({
              qs::qsave(get(simple_name, envir = temp_data_env), qs_file_path)
              if (verbose) message("    Saved as .qs for faster future loading")
            }, error = function(e) {
              # Silently fail - .qs conversion is optional
            })
          }
        }
      }
      temp_data_env$lr_network <- temp_data_env$lr_network %>% dplyr::distinct(from, to)
      loaded_data_for_global_assignment[["lr_network"]] <- temp_data_env$lr_network
      temp_data_env$weighted_networks_lr <- temp_data_env$weighted_networks$lr_sig %>%
        dplyr::inner_join(temp_data_env$lr_network, by = c("from", "to"))
      loaded_data_for_global_assignment[["weighted_networks_lr"]] <- temp_data_env$weighted_networks_lr
      if (!is.null(nichenet_data_name) && nzchar(nichenet_data_name)) {
        assign(nichenet_data_name, loaded_data_for_global_assignment, envir = .GlobalEnv)
        if (verbose) message("NicheNet data loaded from files/download and assigned to '", nichenet_data_name, "' in global environment.")
      } else {
        if (verbose) message("NicheNet data loaded from files/download. Not assigned to global environment.")
      }
    }
    
    lr_network             <- get("lr_network", envir = temp_data_env)
    sig_network            <- get("sig_network", envir = temp_data_env)
    gr_network             <- get("gr_network", envir = temp_data_env)
    weighted_networks      <- get("weighted_networks", envir = temp_data_env)
    ligand_target_matrix   <- get("ligand_target_matrix", envir = temp_data_env)
    ligand_tf_matrix       <- get("ligand_tf_matrix", envir = temp_data_env)
    weighted_networks_lr   <- get("weighted_networks_lr", envir = temp_data_env)
    rm(temp_data_env)
    
    nichenet_data_loaded_info <- list(
      source = if(nichenet_data_loaded_from_global) paste("Global Env:", nichenet_data_name) else paste("Files/Download from:", data_dir_to_use %||% file.path(getwd(), paste0("nichenet_data_", species))), # ensure data_dir_to_use is defined
      species = species, files = names(required_files_map)
    )
    if (verbose) message("NicheNet data successfully prepared for analysis.")
    
    # --- II. Seurat Object Setup & Expressed Genes ---
    if (verbose) message("Setting up Seurat object and identifying expressed genes...")
    start_time_genes <- Sys.time()
    Seurat::DefaultAssay(seurat_obj) <- assay_name
    Seurat::Idents(seurat_obj) <- cluster_col
    if(!all(sender_celltypes %in% levels(Seurat::Idents(seurat_obj)))) stop("One or more sender_celltypes not found in Seurat Idents (column '", cluster_col, "').", call. = FALSE)
    if(!receiver_celltype %in% levels(Seurat::Idents(seurat_obj))) stop("Receiver_celltype '", receiver_celltype, "' not found in Seurat Idents (column '", cluster_col, "').", call. = FALSE)
    
    if (verbose) {
      message("  Checking expressed genes in ", length(sender_celltypes), " sender cell type(s)...")
      message("    [This may take a while for large datasets]")
    }
    # Process with progress logging for large sender lists
    n_senders <- length(sender_celltypes)
    list_expressed_genes_sender <- list()
    for (i in seq_along(sender_celltypes)) {
      s_type <- sender_celltypes[i]
      if (verbose) {
        pct <- round(i / n_senders * 100, 1)
        bar_length <- 30
        filled <- round(bar_length * i / n_senders)
        bar <- paste0(rep("=", filled), collapse = "")
        spaces <- paste0(rep(" ", bar_length - filled), collapse = "")
        message("    [", bar, spaces, "] ", pct, "% (", i, "/", n_senders, ") Processing sender: ", s_type)
      }
      list_expressed_genes_sender[[i]] <- nichenetr::get_expressed_genes(s_type, seurat_obj, min_pct_expressed)
    }
    expressed_genes_sender <- unique(unlist(list_expressed_genes_sender))
    if (verbose) message("  Checking expressed genes in receiver cell type...")
    expressed_genes_receiver <- nichenetr::get_expressed_genes(receiver_celltype, seurat_obj, min_pct_expressed)
    elapsed_genes <- difftime(Sys.time(), start_time_genes, units = "secs")
    if (verbose) message("  Expressed genes identified (", round(elapsed_genes, 1), " seconds)")
    background_expressed_genes <- expressed_genes_receiver[expressed_genes_receiver %in% rownames(ligand_target_matrix)]
    if(length(background_expressed_genes) < 20) warning("NN_Warning: Very few (<20) background expressed genes found in receiver (", length(background_expressed_genes), "). NicheNet might not perform well.")
    
    # --- III. Define Gene Set of Interest (DEGs in receiver) ---
    if (!is.null(receiver_de_table)) {
      if (verbose) message("Using precomputed receiver DE table (", nrow(receiver_de_table), " genes)")
      DE_table_receiver <- receiver_de_table
    } else {
      if (verbose) message("Performing DE analysis for receiver: ", receiver_celltype, " (", receiver_DE_ident1, " vs ", if(!is.null(receiver_DE_ident2)) receiver_DE_ident2 else "all others", " in group ", receiver_DE_group_by,")")
      if (!receiver_DE_group_by %in% colnames(seurat_obj@meta.data)) stop("`receiver_DE_group_by` column '", receiver_DE_group_by, "' not found in Seurat object metadata.", call. = FALSE)
      start_time_de <- Sys.time()
      DE_table_receiver <- Seurat::FindMarkers(
        seurat_obj,
        ident.1 = receiver_DE_ident1,
        ident.2 = receiver_DE_ident2,
        group.by = receiver_DE_group_by,
        subset.ident = receiver_celltype,
        assay = assay_name,
        min.pct = min_pct_expressed,
        logfc.threshold = 0
      ) %>%
        tibble::rownames_to_column("gene")
      elapsed_de <- difftime(Sys.time(), start_time_de, units = "secs")
      if (verbose) message("  DE analysis completed (", round(elapsed_de, 1), " seconds)")
    }

    DE_table_receiver <- tibble::as_tibble(DE_table_receiver)
    if (receiver_gene_col %in% colnames(DE_table_receiver) && receiver_gene_col != "gene") {
      DE_table_receiver <- DE_table_receiver %>% dplyr::rename(gene = !!rlang::sym(receiver_gene_col))
    } else if (!"gene" %in% colnames(DE_table_receiver)) {
      if (!is.null(rownames(DE_table_receiver))) {
        DE_table_receiver <- DE_table_receiver %>% tibble::rownames_to_column("gene")
      } else {
        stop("Receiver DE table must contain a 'gene' column or gene rownames.", call. = FALSE)
      }
    }

    detect_logfc_col <- function(df) {
      if ("avg_log2FC" %in% colnames(df)) return("avg_log2FC")
      if ("avg_logFC" %in% colnames(df)) return("avg_logFC")
      if ("logFC" %in% colnames(df)) return("logFC")
      NULL
    }
    detect_pval_col <- function(df) {
      if ("p_val_adj" %in% colnames(df)) return("p_val_adj")
      if ("FDR" %in% colnames(df)) return("FDR")
      if ("p_val" %in% colnames(df)) return("p_val")
      NULL
    }

    logfc_col <- if (!is.null(receiver_logfc_col)) receiver_logfc_col else detect_logfc_col(DE_table_receiver)
    pval_col <- if (!is.null(receiver_pval_col)) receiver_pval_col else detect_pval_col(DE_table_receiver)

    # If top_n_targets_per_ligand is specified, skip cutoff filtering
    # (DEGs were already filtered by top_n in extract_receiver_degs)
    # Otherwise, apply cutoff filters
    if (is.null(top_n_targets_per_ligand) || top_n_targets_per_ligand == 0) {
      if (!is.null(logfc_col) && logfc_col %in% colnames(DE_table_receiver)) {
        DE_table_receiver <- DE_table_receiver %>%
          dplyr::filter(!!rlang::sym(logfc_col) > logfc_cutoff)
      } else if (is.null(receiver_de_table)) {
        warning("NN_Warning: Receiver DE table lacks a logFC column. Skipping logFC filtering.")
      }

      if (!is.null(pval_col) && pval_col %in% colnames(DE_table_receiver)) {
        # Only apply p-value filter if cutoff is meaningful (< 1.0)
        if (p_val_adj_cutoff < 1.0) {
          DE_table_receiver <- DE_table_receiver %>%
            dplyr::filter(!!rlang::sym(pval_col) < p_val_adj_cutoff)
        }
      } else if (is.null(receiver_de_table)) {
        warning("NN_Warning: Receiver DE table lacks a p-value column. Skipping p-value filtering.")
      }
    } else {
      if (verbose) message("  top_n_targets_per_ligand is specified (", top_n_targets_per_ligand, "). Skipping cutoff filtering (DEGs already filtered by top_n).")
    }

    if (verbose) message("  Receiver DEGs after column filtering: ", nrow(DE_table_receiver))

    geneset_oi <- DE_table_receiver$gene[DE_table_receiver$gene %in% rownames(ligand_target_matrix)]
    if(length(geneset_oi) == 0) {
      if (verbose) {
        message("  DEGs in DE_table_receiver: ", nrow(DE_table_receiver))
        message("  Genes in ligand_target_matrix: ", nrow(ligand_target_matrix))
        overlap_before <- base::intersect(DE_table_receiver$gene, rownames(ligand_target_matrix))
        message("  Overlap: ", length(overlap_before), " genes")
      }
      stop("No DEGs found in the receiver cell type that are also in the NicheNet ligand_target_matrix. Number of DEGs before NicheNet filter: ", nrow(DE_table_receiver), call. = FALSE)
    }
    if (verbose) {
      message(length(geneset_oi), " DEGs (upregulated in ident.1) from receiver '", receiver_celltype, "' identified for NicheNet analysis.")
      message("  Sample genes in geneset_oi: ", paste(head(geneset_oi, 5), collapse = ", "), if(length(geneset_oi) > 5) "..." else "")
    }
    
    # --- IV. Define Potential Ligands ---
    if (verbose) message("Defining potential ligands...")
    ligands_db <- lr_network %>% dplyr::pull(from) %>% unique()
    receptors_db <- lr_network %>% dplyr::pull(to) %>% unique()
    expressed_ligands <- base::intersect(ligands_db, expressed_genes_sender)
    expressed_receptors <- base::intersect(receptors_db, expressed_genes_receiver)
    potential_ligands <- lr_network %>% dplyr::filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% dplyr::pull(from) %>% unique()
    if(length(potential_ligands) == 0) stop("No potential ligands found based on expression in sender/receiver and L-R database.", call. = FALSE)
    if (verbose) message(length(potential_ligands), " potential ligands identified.")
    
    # --- V. Perform NicheNet Ligand Activity Analysis ---
    if (verbose) {
      message("Predicting ligand activities for ", length(potential_ligands), " potential ligands...")
      message("  [Processing ", length(geneset_oi), " DEGs - this may take several minutes for large DEG sets]")
      message("  [", paste0(rep("=", 30), collapse = ""), "] Starting ligand activity prediction...")
    }
    start_time_ligand <- Sys.time()
    ligand_activities <- nichenetr::predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands) %>% dplyr::arrange(dplyr::desc(aupr_corrected))
    elapsed_ligand <- difftime(Sys.time(), start_time_ligand, units = "secs")
    if (verbose) {
      message("  [", paste0(rep("=", 30), collapse = ""), "] 100% Ligand activity prediction completed (", round(elapsed_ligand, 1), " seconds)")
      if (elapsed_ligand > 60) {
        message("  [Note: Large DEG sets take longer. Consider filtering DEGs if runtime is a concern.]")
      }
    }
    best_upstream_ligands <- ligand_activities %>% dplyr::top_n(min(top_n_ligands, nrow(ligand_activities)), aupr_corrected) %>% dplyr::pull(test_ligand) %>% unique()
    if (verbose) {
      message("Top ", length(best_upstream_ligands), " ligands selected based on AUPR.")
      target_est <- max(1, round(length(best_upstream_ligands) * top_n_targets_per_ligand / 800, 1))
      message("  (Ligand-target inference ~", target_est, " min heuristic)")
    }
    
    # --- VI. Infer Active Target Genes of Top Ligands ---
    if (verbose) {
      message("Inferring active target genes for top ", length(best_upstream_ligands), " ligands...")
      message("  [Each ligand will infer up to ", top_n_targets_per_ligand, " target genes]")
    }
    start_time_targets <- Sys.time()
    # Process ligands with progress logging
    active_ligand_target_links_list <- list()
    n_ligands <- length(best_upstream_ligands)
    for (i in seq_along(best_upstream_ligands)) {
      lig <- best_upstream_ligands[i]
      if (verbose) {
        pct <- round(i / n_ligands * 100, 1)
        bar_length <- 30
        filled <- round(bar_length * i / n_ligands)
        bar <- paste0(rep("=", filled), collapse = "")
        spaces <- paste0(rep(" ", bar_length - filled), collapse = "")
        message("  [", bar, spaces, "] ", pct, "% (", i, "/", n_ligands, ") Processing ligand: ", lig)
      }
      active_ligand_target_links_list[[i]] <- tryCatch({
        links_result <- nichenetr::get_weighted_ligand_target_links(
          ligand = lig,  # Note: parameter is 'ligand' (singular), not 'ligands'
          geneset = geneset_oi,
          ligand_target_matrix = ligand_target_matrix,
          n = top_n_targets_per_ligand
        )
        # Debug: check if links were found
        if (verbose && nrow(links_result) == 0 && i == 1) {
          message("    Note: No links found for first ligand '", lig, "'. Checking geneset overlap...")
          if (length(geneset_oi) > 0) {
            message("      geneset_oi size: ", length(geneset_oi))
            message("      Ligand in matrix: ", lig %in% colnames(ligand_target_matrix))
            if (lig %in% colnames(ligand_target_matrix)) {
              ligand_targets <- rownames(ligand_target_matrix)[ligand_target_matrix[, lig] > 0]
              overlap <- base::intersect(geneset_oi, ligand_targets)
              message("      Overlap with geneset_oi: ", length(overlap), " genes")
            }
          }
        }
        links_result
      }, error = function(e) {
        warning("NN_Warning: Failed to get links for ligand ", lig, ": ", e$message)
        return(tibble::tibble())
      })
    }
    active_ligand_target_links_df <- dplyr::bind_rows(active_ligand_target_links_list) %>% tidyr::drop_na()
    elapsed_targets <- difftime(Sys.time(), start_time_targets, units = "secs")
    if (verbose) {
      message("  Ligand-target inference completed (", round(elapsed_targets, 1), " seconds, ", nrow(active_ligand_target_links_df), " links found)")
      if (nrow(active_ligand_target_links_df) == 0) {
        message("  WARNING: No ligand-target links found. This may indicate:")
        message("    - DEGs are not in NicheNet ligand_target_matrix")
        message("    - Ligands have no predicted targets")
        message("    - Consider checking DEG relevance or lowering ligand_target_cutoff")
      }
    }
    p_ligand_target_network <- NULL
    if(nrow(active_ligand_target_links_df) > 0){
      if (verbose) {
        message("  Preparing ligand-target visualization matrix...")
        message("    Input links: ", nrow(active_ligand_target_links_df), " links")
        message("    Using cutoff: ", ligand_target_cutoff)
      }
      start_time_viz <- Sys.time()
      active_ligand_target_links_matrix <- nichenetr::prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = ligand_target_cutoff)
      elapsed_viz <- difftime(Sys.time(), start_time_viz, units = "secs")
      if (verbose) {
        if (!is.null(active_ligand_target_links_matrix)) {
          message("  Visualization matrix prepared (", round(elapsed_viz, 1), " seconds, ", nrow(active_ligand_target_links_matrix), " targets x ", ncol(active_ligand_target_links_matrix), " ligands)")
        } else {
          message("  Visualization matrix is NULL after preparation")
        }
      }
      # Auto-adjust cutoff if matrix is empty
      current_cutoff <- ligand_target_cutoff
      max_attempts <- 5
      attempt <- 1
      
      while ((is.null(active_ligand_target_links_matrix) || nrow(active_ligand_target_links_matrix) == 0 || ncol(active_ligand_target_links_matrix) == 0) && attempt <= max_attempts) {
        if (attempt > 1) {
          # Lower cutoff progressively
          current_cutoff <- max(0.1, current_cutoff - 0.05)
          if (verbose) {
            message("  Attempt ", attempt, ": Retrying with lower cutoff (", current_cutoff, ")...")
          }
          active_ligand_target_links_matrix <- nichenetr::prepare_ligand_target_visualization(
            ligand_target_df = active_ligand_target_links_df, 
            ligand_target_matrix = ligand_target_matrix, 
            cutoff = current_cutoff
          )
        }
        
        if (is.null(active_ligand_target_links_matrix) || nrow(active_ligand_target_links_matrix) == 0 || ncol(active_ligand_target_links_matrix) == 0) {
          attempt <- attempt + 1
        } else {
          if (verbose && attempt > 1) {
            message("  Successfully created matrix with cutoff ", current_cutoff, " (", nrow(active_ligand_target_links_matrix), " targets x ", ncol(active_ligand_target_links_matrix), " ligands)")
          }
          break
        }
      }
      
      if (is.null(active_ligand_target_links_matrix) || nrow(active_ligand_target_links_matrix) == 0 || ncol(active_ligand_target_links_matrix) == 0) {
        warning("NN_Warning: Ligand-target matrix is empty/small after ", max_attempts, " attempts with cutoffs down to ", current_cutoff, ". ",
                "Input had ", nrow(active_ligand_target_links_df), " links. ",
                "Consider checking if DEGs are relevant for NicheNet analysis.")
        p_ligand_target_network <- ggplot2::ggplot() + ggplot2::labs(title=paste0("No L-T links for heatmap after preparation (cutoff=", current_cutoff, ")"))
      } else {
        ligands_for_ordering <- base::intersect(best_upstream_ligands, colnames(active_ligand_target_links_matrix)) %>% rev()
        targets_for_ordering <- active_ligand_target_links_df$target %>% unique() %>% base::intersect(rownames(active_ligand_target_links_matrix))
        rownames(active_ligand_target_links_matrix) <- make.names(rownames(active_ligand_target_links_matrix))
        colnames(active_ligand_target_links_matrix) <- make.names(colnames(active_ligand_target_links_matrix))
        order_ligands_final <- make.names(ligands_for_ordering)
        order_targets_final <- make.names(targets_for_ordering)
        order_ligands_final_intersect <- base::intersect(order_ligands_final, colnames(active_ligand_target_links_matrix))
        order_targets_final_intersect <- base::intersect(order_targets_final, rownames(active_ligand_target_links_matrix))
        if (length(order_targets_final_intersect) == 0 || length(order_ligands_final_intersect) == 0) {
          warning("NN_Warning: No common targets or ligands for L-T heatmap after make.names. Skipping.")
          p_ligand_target_network <- ggplot2::ggplot() + ggplot2::labs(title="No common targets/ligands for L-T heatmap.")
        } else {
          active_ligand_target_links_vis <- active_ligand_target_links_matrix[order_targets_final_intersect, order_ligands_final_intersect, drop = FALSE]
          p_ligand_target_network <- nichenetr::make_heatmap_ggplot(active_ligand_target_links_vis %>% t(), "Prioritized ligands", "Predicted target genes", color = "purple", legend_position = "top", x_axis_position = "top", legend_title = "Regulatory potential") + ggplot2::theme(axis.text.x = ggplot2::element_text(face = "italic"))
        }
      }
    } else {
      warning("NN_Warning: active_ligand_target_links_df is empty. Skipping L-T heatmap.")
      p_ligand_target_network <- ggplot2::ggplot() + ggplot2::labs(title="No active L-T links found.")
    }
    
    # --- VII. Infer Receptors of Top-Ranked Ligands ---
    if (verbose) message("Inferring receptors for top ligands...")
    p_ligand_receptor_network <- NULL
    lr_network_top_df_large <- dplyr::tibble()
    lr_network_top <- lr_network %>% dplyr::filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% dplyr::distinct(from, to)
    best_upstream_receptors <- lr_network_top %>% dplyr::pull(to) %>% unique()
    if(length(best_upstream_ligands) > 0 && length(best_upstream_receptors) > 0){
      lr_network_top_df_large <- weighted_networks_lr %>% dplyr::filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
      valid_ligands_in_df <- base::intersect(best_upstream_ligands, unique(lr_network_top_df_large$from))
      valid_receptors_in_df <- base::intersect(best_upstream_receptors, unique(lr_network_top_df_large$to))
      if(length(valid_ligands_in_df) > 0 && length(valid_receptors_in_df) > 0) {
        lr_network_top_df_filtered <- lr_network_top_df_large %>% dplyr::filter(from %in% valid_ligands_in_df, to %in% valid_receptors_in_df)
        if (nrow(lr_network_top_df_filtered) > 0) {
          lr_network_top_df <- lr_network_top_df_filtered %>% tidyr::spread(key = "from", value = "weight", fill = 0)
          lr_network_top_matrix <- lr_network_top_df %>% dplyr::select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)
          if(nrow(lr_network_top_matrix) > 1 && ncol(lr_network_top_matrix) > 1) {
            dist_receptors <- stats::dist(lr_network_top_matrix, method = "binary"); hclust_receptors <- stats::hclust(dist_receptors, method = "ward.D2"); order_receptors <- hclust_receptors$labels[hclust_receptors$order]
            dist_ligands <- stats::dist(t(lr_network_top_matrix), method = "binary"); hclust_ligands <- stats::hclust(dist_ligands, method = "ward.D2"); order_ligands_receptor <- hclust_ligands$labels[hclust_ligands$order]
          } else { order_receptors <- rownames(lr_network_top_matrix); order_ligands_receptor <- colnames(lr_network_top_matrix) }
          vis_ligand_receptor_network <- lr_network_top_matrix[order_receptors, order_ligands_receptor, drop=FALSE]
          rownames(vis_ligand_receptor_network) <- make.names(rownames(vis_ligand_receptor_network)); colnames(vis_ligand_receptor_network) <- make.names(colnames(vis_ligand_receptor_network))
          p_ligand_receptor_network <- nichenetr::make_heatmap_ggplot(vis_ligand_receptor_network %>% t(), "Prioritized Ligands", "Receptors", color = "mediumvioletred", x_axis_position = "top", legend_title = "Prior interaction potential")
        } else { warning("NN_Warning: Not enough valid L-R pairs for heatmap after filtering."); p_ligand_receptor_network <- ggplot2::ggplot() + ggplot2::labs(title="Not enough valid L-R pairs for heatmap.") }
      } else { warning("NN_Warning: Not enough valid L-R pairs for heatmap after initial filtering."); p_ligand_receptor_network <- ggplot2::ggplot() + ggplot2::labs(title="Not enough valid L-R pairs for heatmap.") }
    } else { warning("NN_Warning: No L-R interactions to display (no best ligands/receptors)."); p_ligand_receptor_network <- ggplot2::ggplot() + ggplot2::labs(title="No L-R interactions to display.") }
    
    # --- VIII. Additional Plots ---
    p_hist_lig_activity <- ggplot2::ggplot(ligand_activities, ggplot2::aes(x = aupr_corrected)) + ggplot2::geom_histogram(color = "black", fill = "darkorange", bins=30) + ggplot2::geom_vline(ggplot2::aes(xintercept = min(ligand_activities %>% dplyr::top_n(min(top_n_ligands, nrow(ligand_activities)), aupr_corrected) %>% dplyr::pull(aupr_corrected))), color = "red", linetype = "dashed", linewidth = 1) + ggplot2::labs(x = "Ligand Activity (AUPR corrected)", y = "# Ligands") + ggplot2::theme_classic()
    vis_ligand_aupr <- ligand_activities %>% dplyr::filter(test_ligand %in% best_upstream_ligands) %>% tibble::column_to_rownames("test_ligand") %>% dplyr::select(aupr_corrected) %>% dplyr::arrange(aupr_corrected) %>% as.matrix()
    p_ligand_aupr_heatmap <- nichenetr::make_heatmap_ggplot(vis_ligand_aupr, "Prioritized ligands", "Ligand activity", legend_title = "AUPR", color = "darkorange") + ggplot2::theme(axis.text.x.top = ggplot2::element_blank())
    
    # --- IX. (Optional) Signaling Path Inference ---
    # THIS ENTIRE SECTION IS REMOVED
    signaling_results <- NULL
    
    # --- X. (Optional) Circos Plot ---
    recorded_circos_plot <- NULL
    # temp_circos_notes <- character(0) # This will be handled by circos_notes_collected
    circos_status_message <- "Circos plot not generated: run_circos was FALSE or required packages missing."
    
    if(run_circos){
      # Ensure Seurat is available for GetAssayData, Idents, Cells
      if(requireNamespace("circlize", quietly = TRUE) &&
         requireNamespace("RColorBrewer", quietly = TRUE) &&
         requireNamespace("Seurat", quietly = TRUE)){ # Seurat needed for Cells, Idents, GetAssayData
        
        if (verbose) message("Generating Circos plot for ligand-receptor interactions with sender information...")
        
        if(exists("lr_network_top_df_large") && is.data.frame(lr_network_top_df_large) && nrow(lr_network_top_df_large) > 0 && all(c("from", "to", "weight") %in% colnames(lr_network_top_df_large)) ){
          circos_lr_data_raw <- lr_network_top_df_large %>% dplyr::filter(from %in% best_upstream_ligands)
          
          if(nrow(circos_lr_data_raw) > 0) {
            # Data Preparation for Circos (ensure all pkg calls are prefixed if not in @importFrom)
            unique_ligands_in_circos <- unique(circos_lr_data_raw$from)
            n_ligands_circos <- length(unique_ligands_in_circos)
            n_senders_circos <- length(sender_celltypes)
            if (verbose) {
              message("  Mapping ", n_ligands_circos, " ligands to ", n_senders_circos, " sender cell types...")
              message("  [", paste0(rep("=", 30), collapse = ""), "] Starting ligand-sender mapping...")
            }
            start_time_circos_map <- Sys.time()
            # Pre-compute cells per sender type to avoid repeated calls
            if (verbose) message("    Step 1/3: Pre-computing cells per sender type...")
            cells_per_sender <- lapply(sender_celltypes, function(s_type) {
              Seurat::Cells(seurat_obj)[Seurat::Idents(seurat_obj) == s_type]
            })
            names(cells_per_sender) <- sender_celltypes
            
            # Optimize: Only load expression data for ligands that are actually in the assay
            # Check which ligands are available in the assay first
            if (verbose) message("    Step 2/3: Loading expression data and filtering available ligands...")
            assay_data_available <- Seurat::GetAssayData(seurat_obj, assay = assay_name, slot = "data")
            available_ligands <- base::intersect(unique_ligands_in_circos, rownames(assay_data_available))
            if (verbose && length(available_ligands) < length(unique_ligands_in_circos)) {
              message("    Note: ", length(unique_ligands_in_circos) - length(available_ligands), " ligands not found in assay, will skip mapping for those.")
            }
            
            # Pre-compute expression percentages per sender for available ligands only
            # This is much faster than checking each ligand individually
            if (verbose) message("    Step 3/3: Computing expression percentages per sender...")
            if (length(available_ligands) > 0 && length(sender_celltypes) > 0) {
              # Use sparse matrix operations for efficiency
              ligand_sender_map <- character(length(unique_ligands_in_circos))
              names(ligand_sender_map) <- unique_ligands_in_circos
              
              # For each sender, compute expression percentage for all ligands at once
              for (i in seq_along(sender_celltypes)) {
                s_type <- sender_celltypes[i]
                if (verbose && n_senders_circos > 5) {
                  pct <- round(i / n_senders_circos * 100, 1)
                  bar_length <- 20
                  filled <- round(bar_length * i / n_senders_circos)
                  bar <- paste0(rep("=", filled), collapse = "")
                  spaces <- paste0(rep(" ", bar_length - filled), collapse = "")
                  message("      [", bar, spaces, "] ", pct, "% (", i, "/", n_senders_circos, ") Processing sender: ", s_type)
                }
                cells_in_sender <- cells_per_sender[[s_type]]
                if (length(cells_in_sender) > 0) {
                  # Get expression for all available ligands in this sender
                  expr_subset <- assay_data_available[available_ligands, cells_in_sender, drop = FALSE]
                  # Compute percentage of cells expressing each ligand
                  pct_expressing <- Matrix::rowSums(expr_subset > 0) / length(cells_in_sender)
                  # Find ligands that meet threshold and haven't been assigned yet
                  ligands_meeting_threshold <- names(pct_expressing)[pct_expressing >= min_pct_expressed]
                  for (lig in ligands_meeting_threshold) {
                    if (lig %in% names(ligand_sender_map) && ligand_sender_map[lig] == "") {
                      ligand_sender_map[lig] <- s_type
                    }
                  }
                }
              }
              # Assign "UnknownSender" to ligands that weren't assigned
              ligand_sender_map[ligand_sender_map == ""] <- "UnknownSender"
            } else {
              ligand_sender_map <- setNames(rep("UnknownSender", length(unique_ligands_in_circos)), unique_ligands_in_circos)
            }
            
            elapsed_circos_map <- difftime(Sys.time(), start_time_circos_map, units = "secs")
            if (verbose) {
              message("  [", paste0(rep("=", 30), collapse = ""), "] 100% Ligand-sender mapping completed (", round(elapsed_circos_map, 1), " seconds)")
              message("  Generating Circos plot visualization...")
            }
            
            ligand_df_circos <- tibble::tibble(
              item = unique_ligands_in_circos,
              color_group = ligand_sender_map[item]
            )
            receptor_df_circos <- tibble::tibble(
              item = unique(circos_lr_data_raw$to),
              color_group = "Receptor"
            )
            links_df_circos <- circos_lr_data_raw %>%
              dplyr::select(from, to, weight) %>% # Use 'weight' for link value if desired by chordDiagram aesthetic
              dplyr::rename(origin = from, target = to)
            
            # Color palette generation
            unique_sender_color_groups <- unique(ligand_df_circos$color_group)
            n_actual_senders <- length(base::setdiff(unique_sender_color_groups, "UnknownSender"))
            sender_palette_values <- character(0)
            if (n_actual_senders > 0) {
              pal_name <- "Set1"
              max_cols_for_palette <- RColorBrewer::brewer.pal.info[pal_name, "maxcolors"]
              num_cols_to_request <- min(n_actual_senders, max_cols_for_palette)
              sender_palette_values <- RColorBrewer::brewer.pal(n = num_cols_to_request, name = pal_name)
              if (n_actual_senders > length(sender_palette_values)) {
                sender_palette_values <- rep(sender_palette_values, length.out = n_actual_senders)
              }
            }
            sender_palette <- stats::setNames(sender_palette_values, base::setdiff(unique_sender_color_groups, "UnknownSender"))
            if ("UnknownSender" %in% unique_sender_color_groups) {
              sender_palette["UnknownSender"] <- "grey50"
            }
            receptor_color <- RColorBrewer::brewer.pal(n = 3, name = "Set2")[1]
            
            # Sector order and color vector
            ordered_ligands <- ligand_df_circos %>% dplyr::arrange(color_group, item) %>% dplyr::pull(item)
            ordered_receptors <- receptor_df_circos %>% dplyr::arrange(item) %>% dplyr::pull(item)
            sector_order <- c(ordered_ligands, ordered_receptors)
            
            grid_cols_ordered <- stats::setNames(character(length(sector_order)), sector_order)
            for(lig in ordered_ligands) {
              sender_grp <- ligand_sender_map[lig] # Get the specific sender group for this ligand
              grid_cols_ordered[lig] <- sender_palette[sender_grp] %||% "grey70" # Fallback
            }
            for(rec in ordered_receptors) {
              grid_cols_ordered[rec] <- receptor_color
            }
            
            
            # Define defaults if not provided by the user
            circos_cex_text <- args_option$circos_cex_text %||% 0.8  # Increased from 0.4 for better visibility
            circos_show_legend <- args_option$circos_show_legend %||% TRUE
            circos_legend_position <- args_option$circos_legend_position %||% "topright"
            circos_legend_size <- args_option$circos_legend_size %||% 0.9
            circos_legend_inset <- args_option$circos_legend_inset %||% c(-0.15, 0)
            circos_png_width <- args_option$circos_png_width %||% 1000
            circos_png_height <- args_option$circos_png_height %||% 1000
            circos_png_res <- args_option$circos_png_res %||% 100
            circos_pdf_width <- args_option$circos_pdf_width %||% 10
            circos_pdf_height <- args_option$circos_pdf_height %||% 10
            # Add any other parameters you want to control for circlize::chordDiagram, etc.
            # e.g., link_transparency <- args_option$link_transparency %||% 0.25
            
            # Prepare legend data
            legend_labels <- character(0)
            legend_colors <- character(0)
            if (circos_show_legend) {
              # Add sender groups to legend
              for (s_group in names(sender_palette)) {
                legend_labels <- c(legend_labels, paste0(s_group, " (Ligands)"))
                legend_colors <- c(legend_colors, sender_palette[s_group])
              }
              # Add receptor to legend
              legend_labels <- c(legend_labels, "Receptor")
              legend_colors <- c(legend_colors, receptor_color)
            }
            
            # Define the plotting function (without legend)
            do_circos_plotting <- function(show_legend = FALSE) {
              circlize::circos.clear()
              circlize::chordDiagram(links_df_circos,
                                     order = sector_order,
                                     grid.col = grid_cols_ordered,
                                     transparency = 0.25,
                                     annotationTrack = "grid",
                                     preAllocateTracks = list(track.height = circlize::mm_h(4)),
                                     directional = 1, direction.type = "arrows",
                                     link.arr.type = "big.arrow", link.arr.length = 0.05, link.arr.width = 0.05)
              circlize::circos.track(track.index = 1, panel.fun = function(x, y) {
                circlize::circos.text(circlize::CELL_META$xcenter, circlize::CELL_META$ylim[1] + circlize::uy(5, "mm"),
                                      circlize::CELL_META$sector.index,
                                      facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),
                                      cex=circos_cex_text)
              }, bg.border = NA)
              
              # Add legend if requested - MUST be after chordDiagram
              if (show_legend && length(legend_labels) > 0) {
                # Set xpd=TRUE to allow plotting outside plot region
                # Use inset = c(0, 0) to place legend at the edge
                graphics::par(xpd = TRUE)
                graphics::legend(
                  x = circos_legend_position,
                  legend = legend_labels,
                  fill = legend_colors,
                  cex = circos_legend_size,
                  inset = c(0, 0),  # No inset - place at edge
                  bg = "white",
                  box.lty = 1,
                  box.lwd = 1,
                  title = "Cell Types"
                )
              }
            }
            
            # Function to plot legend only
            do_legend_only <- function() {
              graphics::par(xpd = TRUE, mar = c(0, 0, 0, 0))
              graphics::plot.new()
              graphics::legend(
                x = "center",
                legend = legend_labels,
                fill = legend_colors,
                cex = circos_legend_size,
                bg = "white",
                box.lty = 1,
                box.lwd = 1,
                title = "Cell Types"
              )
            }
            
            # Capture messages (notes) from circlize
            # Plot to a dummy PDF to ensure full execution for message capture
            temp_pdf_for_capture_notes <- tempfile(fileext = ".pdf")
            grDevices::pdf(temp_pdf_for_capture_notes)
            # Note: The actual capture happens on the expression inside capture.output
            captured_circlize_messages <- utils::capture.output({
              do_circos_plotting() # This call's messages will be captured
            }, type = "message")
            grDevices::dev.off() # Close the dummy PDF device
            if(file.exists(temp_pdf_for_capture_notes)) unlink(temp_pdf_for_capture_notes) # Clean up
            
            # Filter for "Note:" lines and limit display
            temp_notes <- grep("^Note:", captured_circlize_messages, value = TRUE)
            if (length(temp_notes) > 0) {
              # Store all notes but only show first 5 in verbose output
              circos_notes_collected <<- c(circos_notes_collected, temp_notes)
              if (verbose && length(temp_notes) > 5) {
                message("  Note: ", length(temp_notes), " circlize notes captured (showing first 5):")
                for (i in 1:min(5, length(temp_notes))) {
                  message("    ", temp_notes[i])
                }
                message("    ... and ", length(temp_notes) - 5, " more (suppressed)")
              } else if (verbose && length(temp_notes) <= 5) {
                for (note in temp_notes) {
                  message("  ", note)
                }
              }
            }
            
            
            tryCatch({
              if(!is.null(current_run_output_dir)){
                # Save plot without legend
                grDevices::pdf(file.path(current_run_output_dir, "NicheNet_Circos_LR_no_legend.pdf"), width=circos_pdf_width, height=circos_pdf_height); do_circos_plotting(show_legend = FALSE); grDevices::dev.off()
                grDevices::png(file.path(current_run_output_dir, "NicheNet_Circos_LR_no_legend.png"), width=circos_png_width, height=circos_png_height, res=circos_png_res); do_circos_plotting(show_legend = FALSE); grDevices::dev.off()
                
                # Save plot with legend
                if (circos_show_legend && length(legend_labels) > 0) {
                  grDevices::pdf(file.path(current_run_output_dir, "NicheNet_Circos_LR_with_legend.pdf"), width=circos_pdf_width, height=circos_pdf_height); do_circos_plotting(show_legend = TRUE); grDevices::dev.off()
                  grDevices::png(file.path(current_run_output_dir, "NicheNet_Circos_LR_with_legend.png"), width=circos_png_width, height=circos_png_height, res=circos_png_res); do_circos_plotting(show_legend = TRUE); grDevices::dev.off()
                  
                  # Save legend only
                  grDevices::pdf(file.path(current_run_output_dir, "NicheNet_Circos_LR_legend_only.pdf"), width=4, height=6); do_legend_only(); grDevices::dev.off()
                  grDevices::png(file.path(current_run_output_dir, "NicheNet_Circos_LR_legend_only.png"), width=400, height=600, res=circos_png_res); do_legend_only(); grDevices::dev.off()
                }
                
                # Save legend text file
                legend_text <- "Circos Legend (Sender Cell Type for Ligands):\n"
                for(s_group in names(sender_palette)){ # Iterate over actual palette used
                  legend_text <- paste0(legend_text, s_group, " (Ligands): Color ", sender_palette[s_group], "\n")
                }
                legend_text <- paste0(legend_text, "Receptor: Color ", receptor_color, "\n")
                cat(legend_text, file=file.path(current_run_output_dir, "NicheNet_Circos_LR_legend.txt"))
                if(verbose) message("Circos plots & legend saved to: ", current_run_output_dir)
              }
              
              # Record plot with legend if requested
              # Open a temporary device to record the plot
              temp_pdf_for_record <- tempfile(fileext = ".pdf")
              grDevices::pdf(temp_pdf_for_record, width = circos_pdf_width, height = circos_pdf_height)
              circlize::circos.clear() # Clear before plotting
              do_circos_plotting(show_legend = circos_show_legend && length(legend_labels) > 0)
              # Force flush to ensure plot is fully rendered before recording
              grDevices::dev.flush()
              # Record the plot while device is still open
              recorded_circos_plot <- grDevices::recordPlot()
              # Verify plot was recorded
              if (is.null(recorded_circos_plot) || (inherits(recorded_circos_plot, "recordedplot") && length(recorded_circos_plot) > 0 && length(recorded_circos_plot[[1]]) == 0)) {
                if (verbose) warning("NN_Warning: recordPlot() produced empty plot. Attempting to re-record...")
                # Try recording again
                grDevices::dev.flush()
                recorded_circos_plot <- grDevices::recordPlot()
              }
              grDevices::dev.off() # Close the temporary device
              if(file.exists(temp_pdf_for_record)) unlink(temp_pdf_for_record) # Clean up
              circlize::circos.clear() # Clear after recordPlot
              
              circos_status_message <- "Circos plot generated and recorded (sender-grouped)."
            }, error = function(e) {
              if (verbose) message("Could not generate or save Sender-grouped Circos plot: ", e$message)
              circos_status_message <- paste("Sender-grouped Circos plot failed:", e$message)
              try(circlize::circos.clear(), silent = TRUE) # Attempt to clean up circlize state
            })
          } else {
            if (verbose) message("No data for Circos plot after filtering (circos_lr_data_raw empty).");
            circos_status_message <- "No data for Circos plot after filtering."
          }
        } else {
          if (verbose) message("Skipping Circos plot: lr_network_top_df_large data unavailable/malformed.");
          circos_status_message <- "Circos plot skipped: lr_network_top_df_large data."
        }
      } else {
        circos_status_message <- "Circos plot skipped: 'circlize', 'RColorBrewer', or 'Seurat' package not fully available."
        if (verbose) message(circos_status_message)
      }
    }
    
    # --- XI. Collate and Return Results ---
    if (verbose) message("NicheNet analysis complete. Collating results...")
    results_list <- list(
      ligand_activities = ligand_activities, best_upstream_ligands = best_upstream_ligands,
      active_ligand_target_links_df = if(exists("active_ligand_target_links_df")) active_ligand_target_links_df else dplyr::tibble(),
      plot_ligand_target_network = p_ligand_target_network,
      ligand_receptor_network_df = lr_network_top_df_large,
      plot_ligand_receptor_network = p_ligand_receptor_network,
      plot_ligand_activity_hist = p_hist_lig_activity, plot_ligand_aupr_heatmap = p_ligand_aupr_heatmap,
      DE_table_receiver = DE_table_receiver,
      plot_circos = if(!is.null(recorded_circos_plot)) recorded_circos_plot else circos_status_message,
      circos_notes = if(length(circos_notes_collected) > 0) circos_notes_collected else "No notes from circos plot generation.", # ADDED
      # signaling_path_graph = signaling_results, # REMOVED
      nichenet_data_loaded_info = nichenet_data_loaded_info,
      output_path_run = if(!is.null(current_run_output_dir)) current_run_output_dir else "Output not saved to files."
    )

    if (!is.null(current_run_output_dir)) {
      checkpoint_path <- file.path(current_run_output_dir, "nichenet_results.qs")
      if (requireNamespace("qs", quietly = TRUE)) {
        tryCatch({
          qs::qsave(results_list, checkpoint_path)
          if (verbose) message("NicheNet checkpoint saved to: ", checkpoint_path)
        }, error = function(e) warning("NN_Warning: Failed to save NicheNet checkpoint: ", e$message))
      }
    }

    
    if (!is.null(current_run_output_dir)) {
      plot_list_gg <- list(NicheNet_Ligand_Target_Heatmap = results_list$plot_ligand_target_network, NicheNet_Ligand_Receptor_Heatmap = results_list$plot_ligand_receptor_network, NicheNet_Ligand_Activity_Histogram = results_list$plot_ligand_activity_hist, NicheNet_Ligand_AUPR_Heatmap = results_list$plot_ligand_aupr_heatmap)
      plot_dims <- list(NicheNet_Ligand_Target_Heatmap = list(w=10,h=8), NicheNet_Ligand_Receptor_Heatmap = list(w=8,h=6), NicheNet_Ligand_Activity_Histogram = list(w=6,h=5), NicheNet_Ligand_AUPR_Heatmap = list(w=6,h=8))
      for(plot_name in names(plot_list_gg)){ 
        p_obj <- plot_list_gg[[plot_name]]
        if(!is.null(p_obj) && inherits(p_obj, "ggplot")){ 
          # Check if plot has actual data (not just empty plot with title)
          plot_file <- file.path(current_run_output_dir, paste0(plot_name, ".png"))
          tryCatch({
            # Use print() to ensure plot is rendered before saving
            ggplot2::ggsave(filename = plot_file, plot = p_obj, width = plot_dims[[plot_name]]$w, height = plot_dims[[plot_name]]$h)
            # Verify file was created and has reasonable size (>1KB)
            if (file.exists(plot_file)) {
              file_size <- file.info(plot_file)$size
              if (file_size < 1024) {
                warning("NN_Warning: Plot ", plot_name, " file is very small (", file_size, " bytes). May be empty.")
              } else if (verbose) {
                message("  Saved ", plot_name, " (", round(file_size/1024, 1), " KB)")
              }
            } else {
              warning("NN_Warning: Plot file was not created: ", plot_file)
            }
          }, error = function(e) {
            warning("NN_Warning: Failed to save plot ", plot_name, ": ", e$message)
          })
        } else if (verbose && !is.null(p_obj)) {
          message("  Skipping ", plot_name, " (not a ggplot object, class: ", class(p_obj)[1], ")")
        } else if (verbose && is.null(p_obj)) {
          message("  Skipping ", plot_name, " (NULL)")
        }
      }
    }

    objects_to_clear <- c("ligand_activities", "active_ligand_target_links_df", "active_ligand_target_links_matrix",
                          "lr_network_top_df_large", "list_expressed_genes_sender", "expressed_genes_sender",
                          "expressed_genes_receiver", "ligand_df_circos", "receptor_df_circos",
                          "links_df_circos", "captured_circlize_messages")
    rm(list = base::intersect(objects_to_clear, ls()), envir = environment())
    gc(verbose = FALSE)

    return(results_list)
  },
  warning = function(w) {
    if (startsWith(conditionMessage(w), "NN_Warning:")) {
      msg <- sub("^NN_Warning: ", "", conditionMessage(w))
      warnings_collected <<- c(warnings_collected, paste("NicheNet Function:", msg))
    } else {
      warnings_collected <<- c(warnings_collected, paste("R Warning:", conditionMessage(w)))
    }
    invokeRestart("muffleWarning")
  })
  
  # Add collected warnings to the final results
  analysis_results$warnings_occurred <- if(length(warnings_collected) > 0) warnings_collected else "No warnings reported."
  
  # Ensure circos_notes is part of the final returned list, even if an error occurred before it was assigned in the main block
  if(!"circos_notes" %in% names(analysis_results) || is.null(analysis_results$circos_notes)){
    analysis_results$circos_notes <- if(length(circos_notes_collected) > 0) circos_notes_collected else "No notes from circos plot generation (or circos plot did not run)."
  }
  
  
  if (verbose && length(warnings_collected) > 0) {
    message("--------------------------------------------------------------------")
    message("NicheNet Analysis completed with ", length(warnings_collected), " warning(s). See 'warnings_occurred' in output list for details.")
    message("Summary of warnings:")
    for(i in seq_along(warnings_collected)){
      message(paste0(i, ". ", warnings_collected[i]))
      if(i >= 5 && length(warnings_collected) > 5) {
        message("... (and ", length(warnings_collected) - 5, " more warnings)")
        break
      }
    }
    message("--------------------------------------------------------------------")
  } else if (verbose) {
    message("NicheNet Analysis completed successfully with no warnings reported.")
  }
  # ADDED: Summary of circos_notes if verbose
  if (verbose && length(analysis_results$circos_notes) > 0 && !all(analysis_results$circos_notes %in% c("No notes from circos plot generation.", "No notes from circos plot generation (or circos plot did not run).")) ) {
    message("--------------------------------------------------------------------")
    message("Notes from Circos plot generation (see 'circos_notes' in output list):")
    for(i in seq_along(analysis_results$circos_notes)){
      note_to_print <- analysis_results$circos_notes[i]
      # Ensure it's a single string for message()
      if(is.list(note_to_print)) note_to_print <- paste(unlist(note_to_print), collapse=" ")
      message(paste0("- ", note_to_print)) # Make sure note_to_print is a single string
      if(i >= 3 && length(analysis_results$circos_notes) > 3) { # Print max 3 notes in summary
        message("... (and ", length(analysis_results$circos_notes) - 3, " more notes)")
        break
      }
    }
    message("--------------------------------------------------------------------")
  }
  
  
  return(analysis_results)
}

# Helper for path construction (from nichenetr utils, or base R `||` if appropriate)
"%||%" <- function(a, b) {
  if (!is.null(a) && !is.logical(a) && (is.character(a) || is.numeric(a)) && length(a) == 1 && !is.na(a) && nzchar(as.character(a)) ) a else b
}



#' Prepare Data for NicheNet Circos Plot
#'
#' @description
#' This function prepares the necessary data structures for generating a NicheNet
#' ligand-receptor circos plot. It identifies expressing sender cell types for
#' ligands and prepares link dataframes and color/order information.
#'
#' @param lr_network_top_df_large A data frame of ligand-receptor interactions,
#'        typically from NicheNet analysis (e.g., `results$ligand_receptor_network_df`).
#'        Must contain 'from' (ligand), 'to' (receptor), and 'weight' columns.
#' @param best_upstream_ligands A character vector of top N ligands to include.
#' @param seurat_obj A Seurat object (used if `ligand_sender_source` is "seurat_obj").
#' @param sender_celltypes Character vector of sender cell types (used if `ligand_sender_source` is "seurat_obj").
#' @param assay_name Name of the assay in Seurat object (used if `ligand_sender_source` is "seurat_obj").
#' @param min_pct_expressed Minimum percentage of cells in a sender type expressing a ligand
#'        for it to be attributed to that sender (used if `ligand_sender_source` is "seurat_obj").
#' @param ligand_sender_map Optional. A named character vector where names are ligands
#'        and values are their primary sender cell types. If provided, `seurat_obj`,
#'        `sender_celltypes`, etc., are ignored for sender mapping.
#' @param circos_args A list of additional arguments for customizing colors or palettes, e.g.,
#'        `list(sender_palette_name = "Set1", receptor_color_name = "Set2", receptor_color_index = 1)`.
#'
#' @return A list containing:
#'   \itemize{
#'     \item `links_df`: Data frame for `chordDiagram` (origin, target, weight).
#'     \item `sector_order`: Character vector defining the order of items in the circos plot.
#'     \item `grid_cols_ordered`: Named character vector of colors for each sector.
#'     \item `sender_palette_legend`: Named character vector for legend (sender type -> color).
#'     \item `receptor_color_legend`: Named character vector for legend (receptor group -> color).
#'   }
#' Returns NULL if essential input is missing or no valid links can be formed.
#'
#' @importFrom dplyr %>% filter select rename arrange pull distinct
#' @importFrom tibble tibble
#' @importFrom Seurat Cells Idents GetAssayData
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
#' @importFrom stats setNames
#' @export
prepare_nichenet_circos_data <- function(lr_network_top_df_large,
                                         best_upstream_ligands,
                                         seurat_obj = NULL,
                                         sender_celltypes = NULL,
                                         assay_name = "RNA", # Default, ensure it matches your usage
                                         min_pct_expressed = 0.10,
                                         ligand_sender_map = NULL,
                                         circos_args = list()) {
  
  if (is.null(lr_network_top_df_large) || nrow(lr_network_top_df_large) == 0 ||
      !all(c("from", "to", "weight") %in% colnames(lr_network_top_df_large))) {
    warning("`lr_network_top_df_large` is missing, empty, or lacks required columns (from, to, weight).")
    return(NULL)
  }
  if (is.null(best_upstream_ligands) || length(best_upstream_ligands) == 0) {
    warning("`best_upstream_ligands` is empty.")
    return(NULL)
  }
  
  circos_lr_data_raw <- lr_network_top_df_large %>%
    dplyr::filter(from %in% best_upstream_ligands)
  
  if (nrow(circos_lr_data_raw) == 0) {
    warning("No ligand-receptor interactions found for the `best_upstream_ligands`.")
    return(NULL)
  }
  
  unique_ligands_in_circos <- unique(circos_lr_data_raw$from)
  
  if (is.null(ligand_sender_map)) {
    if (is.null(seurat_obj) || is.null(sender_celltypes)) {
      warning("If `ligand_sender_map` is not provided, `seurat_obj` and `sender_celltypes` must be supplied.")
      return(NULL)
    }
    # Determine ligand_sender_map from Seurat object
    ligand_sender_map <- sapply(unique_ligands_in_circos, function(lig) {
      expressing_senders <- c()
      for (s_type in sender_celltypes) {
        # Ensure Idents are set on the Seurat object if not already
        # current_idents <- Seurat::Idents(seurat_obj) # store current idents
        # Seurat::Idents(seurat_obj) <- seurat_obj@meta.data[[cluster_col_name]] # If cluster_col_name needed
        
        cells_in_sender <- Seurat::Cells(seurat_obj)[Seurat::Idents(seurat_obj) == s_type]
        if (length(cells_in_sender) > 0) {
          if (lig %in% rownames(Seurat::GetAssayData(seurat_obj, assay = assay_name, slot = "data"))) {
            gene_expr_in_sender <- Seurat::GetAssayData(seurat_obj, assay = assay_name, slot = "data")[lig, cells_in_sender, drop=FALSE]
            if (sum(gene_expr_in_sender > 0) / length(cells_in_sender) >= min_pct_expressed) {
              expressing_senders <- c(expressing_senders, s_type)
            }
          }
        }
        # Seurat::Idents(seurat_obj) <- current_idents # restore idents
      }
      if (length(expressing_senders) > 0) return(expressing_senders[1])
      return("UnknownSender")
    }, USE.NAMES = TRUE)
  }
  
  
  ligand_df_circos <- tibble::tibble(
    item = unique_ligands_in_circos,
    color_group = ligand_sender_map[item] # item is the ligand name
  )
  receptor_df_circos <- tibble::tibble(
    item = unique(circos_lr_data_raw$to),
    color_group = "Receptor"
  )
  links_df <- circos_lr_data_raw %>%
    dplyr::select(from, to, weight) %>%
    dplyr::rename(origin = from, target = to)
  
  # Color palette generation
  sender_palette_name <- circos_args$sender_palette_name %||% "Set1"
  receptor_color_name <- circos_args$receptor_color_name %||% "Set2"
  receptor_color_idx <- circos_args$receptor_color_idx %||% 1
  
  
  unique_sender_color_groups <- unique(ligand_df_circos$color_group)
  n_actual_senders <- length(base::setdiff(unique_sender_color_groups, "UnknownSender"))
  sender_palette_values <- character(0)
  
  if (n_actual_senders > 0) {
    max_cols_for_palette <- RColorBrewer::brewer.pal.info[sender_palette_name, "maxcolors"]
    num_cols_to_request <- min(n_actual_senders, max_cols_for_palette)
    sender_palette_values <- RColorBrewer::brewer.pal(n = num_cols_to_request, name = sender_palette_name)
    if (n_actual_senders > length(sender_palette_values)) {
      sender_palette_values <- rep(sender_palette_values, length.out = n_actual_senders)
    }
  }
  sender_palette_legend <- stats::setNames(sender_palette_values, base::setdiff(unique_sender_color_groups, "UnknownSender"))
  if ("UnknownSender" %in% unique_sender_color_groups) {
    sender_palette_legend["UnknownSender"] <- "grey50"
  }
  receptor_color_legend <- stats::setNames(
    RColorBrewer::brewer.pal(n = max(3, receptor_color_idx), name = receptor_color_name)[receptor_color_idx],
    "Receptor"
  )
  
  
  # Sector order and color vector
  ordered_ligands <- ligand_df_circos %>% dplyr::arrange(color_group, item) %>% dplyr::pull(item)
  ordered_receptors <- receptor_df_circos %>% dplyr::arrange(item) %>% dplyr::pull(item)
  sector_order <- c(ordered_ligands, ordered_receptors)
  
  grid_cols_ordered <- stats::setNames(character(length(sector_order)), sector_order)
  for(lig in ordered_ligands) {
    sender_grp_for_lig <- ligand_sender_map[lig]
    grid_cols_ordered[lig] <- sender_palette_legend[sender_grp_for_lig] %||% "grey70" # Fallback
  }
  for(rec in ordered_receptors) {
    grid_cols_ordered[rec] <- receptor_color_legend["Receptor"]
  }
  
  return(list(
    links_df = links_df,
    sector_order = sector_order,
    grid_cols_ordered = grid_cols_ordered,
    sender_palette_legend = sender_palette_legend,
    receptor_color_legend = receptor_color_legend
  ))
}


#' Draw NicheNet Ligand-Receptor Circos Plot
#'
#' @description
#' This function takes pre-formatted data and generates a circos plot
#' visualizing ligand-receptor interactions.
#'
#' @param circos_data A list containing `links_df`, `sector_order`, and
#'        `grid_cols_ordered`, typically the output of `prepare_nichenet_circos_data()`.
#' @param cex_text Numeric, expansion factor for text labels. Default: 0.4.
#' @param link_transparency Numeric, transparency for links (0-1). Default: 0.25.
#' @param plot_title Optional character string for the plot title.
#' @param show_legend Logical, whether to attempt to draw a legend on the plot (basic). Default: FALSE.
#' @param sender_palette_legend Named character vector (sender type -> color) for legend.
#' @param receptor_color_legend Named character vector (receptor group -> color) for legend.
#'
#' @return A recorded plot object if successful, otherwise NULL.
#'         The plot is also drawn to the current graphics device.
#'
#' @importFrom circlize circos.clear chordDiagram circos.track circos.text uy mm_h CELL_META
#' @importFrom graphics par legend title
#' @importFrom grDevices recordPlot
#' @export
draw_nichenet_circos_plot <- function(circos_data,
                                      cex_text = 0.4,
                                      link_transparency = 0.25,
                                      plot_title = NULL,
                                      show_legend = FALSE,
                                      sender_palette_legend = NULL, # From circos_data
                                      receptor_color_legend = NULL, # From circos_data
                                      legend_size=0.9,
                                      legend_inset=c(-0.15,0)
) {
  if (is.null(circos_data) ||
      !all(c("links_df", "sector_order", "grid_cols_ordered") %in% names(circos_data))) {
    warning("Invalid `circos_data` input.")
    return(NULL)
  }
  
  links_df <- circos_data$links_df
  sector_order <- circos_data$sector_order
  grid_cols_ordered <- circos_data$grid_cols_ordered
  
  # Use legend info from circos_data if not passed directly
  if(is.null(sender_palette_legend) && !is.null(circos_data$sender_palette_legend)){
    sender_palette_legend <- circos_data$sender_palette_legend
  }
  if(is.null(receptor_color_legend) && !is.null(circos_data$receptor_color_legend)){
    receptor_color_legend <- circos_data$receptor_color_legend
  }
  
  
  recorded_plot <- NULL
  captured_notes_list <- list() # To store notes if any
  
  tryCatch({
    # Capture notes
    temp_pdf_for_capture_notes <- tempfile(fileext = ".pdf")
    grDevices::pdf(temp_pdf_for_capture_notes)
    
    current_notes <- utils::capture.output({
      circlize::circos.clear()
      circlize::chordDiagram(links_df,
                             order = sector_order,
                             grid.col = grid_cols_ordered,
                             transparency = link_transparency,
                             annotationTrack = "grid",
                             preAllocateTracks = list(track.height = circlize::mm_h(4)),
                             directional = 1, direction.type = "arrows",
                             link.arr.type = "big.arrow", link.arr.length = 0.05, link.arr.width = 0.05)
      circlize::circos.track(track.index = 1, panel.fun = function(x, y) {
        circlize::circos.text(circlize::CELL_META$xcenter,
                              circlize::CELL_META$ylim[1] + circlize::uy(5, "mm"),
                              circlize::CELL_META$sector.index,
                              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = cex_text)
      }, bg.border = NA)
    }, type = "message")
    grDevices::dev.off()
    if(file.exists(temp_pdf_for_capture_notes)) unlink(temp_pdf_for_capture_notes)
    
    captured_notes_list$notes <- grep("^Note:", current_notes, value = TRUE)
    
    
    # Actual plotting for recording
    circlize::circos.clear()
    circlize::chordDiagram(links_df,
                           order = sector_order,
                           grid.col = grid_cols_ordered,
                           transparency = link_transparency,
                           annotationTrack = "grid",
                           preAllocateTracks = list(track.height = circlize::mm_h(4)),
                           directional = 1, direction.type = "arrows",
                           link.arr.type = "big.arrow", link.arr.length = 0.05, link.arr.width = 0.05)
    circlize::circos.track(track.index = 1, panel.fun = function(x, y) {
      circlize::circos.text(circlize::CELL_META$xcenter,
                            circlize::CELL_META$ylim[1] + circlize::uy(5, "mm"),
                            circlize::CELL_META$sector.index,
                            facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = cex_text)
    }, bg.border = NA)
    
    if (!is.null(plot_title)) {
      graphics::title(main = plot_title)
    }
    
    if (show_legend && !is.null(sender_palette_legend) && !is.null(receptor_color_legend)) {
      # Combine legends
      all_legend_items <- c(sender_palette_legend, receptor_color_legend)
      if (length(all_legend_items) > 0) {
        # Need to plot legend outside the circos plot area, adjust plot margins
        # This is a basic legend, ComplexHeatmap::Legend is better for complex legends
        current_par <- graphics::par(no.readonly = TRUE)
        on.exit(graphics::par(current_par)) # Reset par on exit
        graphics::par(mar = c(1, 1, 1, 7) + 0.1, xpd=TRUE) # Adjust right margin for legend
        
        legend_names <- names(all_legend_items)
        legend_colors <- all_legend_items
        
        # Filter out "UnknownSender" if its color is grey50 to avoid redundant grey in legend
        # Or, make sure it's named distinctively if "grey50" is a general fallback.
        # For now, include all named entries.
        
        graphics::legend("topright",
                         inset = legend_inset, # Adjust inset to place it outside plot
                         legend = legend_names,
                         fill = legend_colors,
                         cex = legend_size, # Adjust legend text size
                         title = "Color Legend",
                         bty = "n") # No box
      }
    }
    
    recorded_plot <- grDevices::recordPlot()
    circlize::circos.clear()
    
  }, error = function(e) {
    warning("Error during circos plot generation: ", e$message)
    try(circlize::circos.clear(), silent = TRUE)
    return(NULL)
  })
  
  # Add notes to the recorded_plot object as an attribute if you want to carry them along
  if (!is.null(recorded_plot) && length(captured_notes_list$notes) > 0) {
    attr(recorded_plot, "circos_notes") <- captured_notes_list$notes
  }
  
  return(recorded_plot)
}
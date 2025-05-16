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
#'        Default: "NicheNet_Global_Dataset".
#' @param output_dir Character string, path to a main directory where run-specific subfolders for outputs will be created. If NULL, plots are returned but not saved to files. Default: NULL.
#' @param run_circos Logical, whether to attempt to generate and return a circos plot of ligand-receptor interactions. Requires `circlize` and `RColorBrewer` packages. Default: TRUE.
#' @param run_signaling_path_inference Logical, whether to infer and visualize ligand-to-target signaling paths for the top 1 ligand. Requires `DiagrammeR`, `DiagrammeRsvg`, `rsvg`. Default: FALSE.
#' @param verbose Logical, whether to print progress messages and summary of warnings. Default: TRUE.
#'
#' @return A list containing:
#'   \itemize{
#'     \item `ligand_activities`: DataFrame of ligand activities.
#'     \item `best_upstream_ligands`: Character vector of top N ligands.
#'     \item `active_ligand_target_links_df`: DataFrame of active ligand-target links.
#'     \item `plot_ligand_target_network`: ggplot object of the ligand-target heatmap.
#'     \item `ligand_receptor_network_df`: DataFrame of prioritized ligand-receptor interactions.
#'     \item `plot_ligand_receptor_network`: ggplot object of the ligand-receptor heatmap.
#'     \item `plot_ligand_activity_hist`: ggplot object of ligand activity histogram.
#'     \item `plot_ligand_aupr_heatmap`: ggplot object of top ligand AUPR heatmap.
#'     \item `DE_table_receiver`: DataFrame of differentially expressed genes in the receiver.
#'     \item `plot_circos`: A recorded circos plot object (replay with `plot()`) or a status message.
#'     \item `signaling_path_graph`: A DiagrammeR graph object or NULL.
#'     \item `nichenet_data_loaded_info`: List describing how NicheNet data was loaded.
#'     \item `output_path_run`: Path to the run-specific output directory if `output_dir` was provided.
#'     \item `warnings_occurred`: A character vector of warning messages collected during execution, if any.
#'   }
#' @details
#' To visualize ligand expression (e.g., as a DotPlot or Heatmap), you can use the returned
#' `best_upstream_ligands` and your original `seurat_obj`. For example:
#' `Seurat::DotPlot(seurat_obj, features = results$best_upstream_ligands, group.by = "your_celltype_column")`
#'
#' @importFrom Seurat DefaultAssay Idents FindMarkers GetAssayData Cells
#' @importFrom dplyr %>% filter select distinct arrange desc top_n pull inner_join mutate group_by summarise across ends_with
#' @importFrom tidyr spread drop_na
#' @importFrom tibble column_to_rownames rownames_to_column as_tibble
#' @importFrom magrittr set_rownames
#' @importFrom nichenetr get_expressed_genes predict_ligand_activities get_weighted_ligand_target_links prepare_ligand_target_visualization make_heatmap_ggplot get_ligand_signaling_path diagrammer_format_signaling_graph
#' @importFrom ggplot2 ggplot aes geom_histogram geom_vline labs theme_classic theme element_text scale_fill_gradient2 RotatedAxis ggsave
#' @importFrom utils download.file capture.output
#' @importFrom stats dist hclust
#' @importFrom grDevices recordPlot dev.off pdf png
#' @importFrom DiagrammeRsvg export_svg
#' @importFrom rsvg rsvg_png
#' @export
run_nichenet_analysis_legacy <- function(seurat_obj,
                                  species = c("human", "mouse"),
                                  sender_celltypes,
                                  receiver_celltype,
                                  assay_name = "SCT",
                                  cluster_col,
                                  receiver_DE_ident1,
                                  receiver_DE_ident2 = NULL,
                                  receiver_DE_group_by,
                                  min_pct_expressed = 0.10,
                                  p_val_adj_cutoff = 0.05,
                                  logfc_cutoff = 0.25,
                                  top_n_ligands = 20,
                                  top_n_targets_per_ligand = 200,
                                  ligand_target_cutoff = 0.33,
                                  nichenet_data_dir = NULL,
                                  nichenet_data_name = "NicheNet_Global_Dataset",
                                  output_dir = NULL,
                                  run_circos = TRUE,
                                  run_signaling_path_inference = FALSE,
                                  verbose = TRUE) {
  
  # --- 0. Initialize Warning Collection & Argument Setup ---
  warnings_collected <- character(0)
  
  # Capture warnings using withCallingHandlers
  # The main analysis block will be wrapped in this
  # Note: `<<-` is used to modify `warnings_collected` in the parent environment of the handler
  
  analysis_results <- withCallingHandlers({
    
    call_obj <- match.call()
    species <- match.arg(species)
    
    # Package checks
    required_pkgs <- c("nichenetr", "Seurat", "dplyr", "tidyr", "tibble", "ggplot2")
    missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
      stop("Required package(s) not found: ", paste(missing_pkgs, collapse=", "), ". Please install them.", call. = FALSE)
    }
    
    # --- Output Directory Management ---
    current_run_output_dir <- NULL
    if (!is.null(output_dir)) {
      # (Logic from previous response for creating run-specific output dir and saving args)
      # This ensures current_run_output_dir is set if output_dir is provided.
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
        con_args <- file(args_file_path, "w")
        cat(paste(utils::capture.output(print(call_obj)), collapse = "\n"), file = con_args)
        close(con_args)
      }, error = function(e) warning("NN_Warning: Failed to save arguments to file: ", args_file_path, ". Error: ", e$message))
      if (verbose && file.exists(args_file_path)) message("Arguments saved to: ", args_file_path)
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
    # Create a temporary environment for loading data to avoid polluting function's env directly before checks
    temp_data_env <- new.env()
    
    if (!is.null(nichenet_data_name) && nzchar(nichenet_data_name)) {
      if (exists(nichenet_data_name, envir = .GlobalEnv)) {
        if (verbose) message("Found '", nichenet_data_name, "' in global environment. Validating...")
        global_data_env <- get(nichenet_data_name, envir = .GlobalEnv)
        if (is.list(global_data_env) && all(names(required_files_map) %in% names(global_data_env))) {
          are_all_components_valid <- all(sapply(names(required_files_map), function(comp_name) !is.null(global_data_env[[comp_name]])))
          # Check for weighted_networks_lr presence or reconstructability
          has_wn_lr <- ("weighted_networks_lr" %in% names(global_data_env) && !is.null(global_data_env[["weighted_networks_lr"]])) ||
            (!is.null(global_data_env[["weighted_networks"]]) && !is.null(global_data_env[["lr_network"]]))
          
          if (are_all_components_valid && has_wn_lr) {
            if (verbose) message("Global data '", nichenet_data_name, "' is valid and complete. Using it.")
            for (comp_name in names(required_files_map)) {
              assign(comp_name, global_data_env[[comp_name]], envir = temp_data_env)
            }
            if ("weighted_networks_lr" %in% names(global_data_env) && !is.null(global_data_env[["weighted_networks_lr"]])) {
              assign("weighted_networks_lr", global_data_env[["weighted_networks_lr"]], envir = temp_data_env)
            } else { # Reconstruct if not directly present but base components are
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
        data_dir_to_use <- file.path(getwd(), paste0("nichenet_data_", species))
        if (verbose) message("`nichenet_data_dir` is NULL, using default local cache: ", data_dir_to_use)
      }
      if (!dir.exists(data_dir_to_use)) {
        if (verbose) message("Creating NicheNet data directory: ", data_dir_to_use)
        tryCatch(dir.create(data_dir_to_use, recursive = TRUE, showWarnings = FALSE),
                 error = function(e) stop("Failed to create NicheNet data directory: ", data_dir_to_use, ". Error: ", e$message))
      }
      
      loaded_data_for_global_assignment <- list()
      for (simple_name in names(required_files_map)) {
        rds_filename <- required_files_map[[simple_name]]
        file_path <- file.path(data_dir_to_use, rds_filename)
        if (!file.exists(file_path)) {
          if (verbose) message("File ", rds_filename, " not found in '", data_dir_to_use, "'. Downloading from Zenodo...")
          download_url <- paste0(base_url, rds_filename, "?download=1")
          tryCatch({
            utils::download.file(url = download_url, destfile = file_path, mode = "wb", quiet = !verbose)
            if (verbose) message("Downloaded ", rds_filename, " successfully to ", data_dir_to_use)
          }, error = function(e) {
            stop("Failed to download ", rds_filename, ". Error: ", e$message, call. = FALSE)
          })
        }
        assign(simple_name, readRDS(file_path), envir = temp_data_env)
        loaded_data_for_global_assignment[[simple_name]] <- get(simple_name, envir = temp_data_env)
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
    
    # Make data from temp_data_env available in function's execution environment
    # These are the key NicheNet data objects needed by the rest of the function.
    lr_network             <- get("lr_network", envir = temp_data_env)
    sig_network            <- get("sig_network", envir = temp_data_env)
    gr_network             <- get("gr_network", envir = temp_data_env)
    weighted_networks      <- get("weighted_networks", envir = temp_data_env)
    ligand_target_matrix   <- get("ligand_target_matrix", envir = temp_data_env)
    ligand_tf_matrix       <- get("ligand_tf_matrix", envir = temp_data_env)
    weighted_networks_lr   <- get("weighted_networks_lr", envir = temp_data_env)
    rm(temp_data_env) # Clean up temp env
    
    nichenet_data_loaded_info <- list(
      source = if(nichenet_data_loaded_from_global) paste("Global Env:", nichenet_data_name) else paste("Files/Download from:", nichenet_data_dir %||% "default path"),
      species = species, files = names(required_files_map)
    )
    if (verbose) message("NicheNet data successfully prepared for analysis.")
    
    # --- II. Seurat Object Setup & Expressed Genes ---
    if (verbose) message("Setting up Seurat object and identifying expressed genes...")
    Seurat::DefaultAssay(seurat_obj) <- assay_name
    Seurat::Idents(seurat_obj) <- cluster_col
    if(!all(sender_celltypes %in% levels(Seurat::Idents(seurat_obj)))) stop("One or more sender_celltypes not found in Seurat Idents (column '", cluster_col, "').", call. = FALSE)
    if(!receiver_celltype %in% levels(Seurat::Idents(seurat_obj))) stop("Receiver_celltype '", receiver_celltype, "' not found in Seurat Idents (column '", cluster_col, "').", call. = FALSE)
    
    list_expressed_genes_sender <- lapply(sender_celltypes, nichenetr::get_expressed_genes, seurat_obj, min_pct_expressed)
    expressed_genes_sender <- unique(unlist(list_expressed_genes_sender))
    expressed_genes_receiver <- nichenetr::get_expressed_genes(receiver_celltype, seurat_obj, min_pct_expressed)
    background_expressed_genes <- expressed_genes_receiver[expressed_genes_receiver %in% rownames(ligand_target_matrix)]
    if(length(background_expressed_genes) < 20) warning("NN_Warning: Very few (<20) background expressed genes found in receiver (", length(background_expressed_genes), "). NicheNet might not perform well.")
    
    # --- III. Define Gene Set of Interest (DEGs in receiver) ---
    if (verbose) message("Performing DE analysis for receiver: ", receiver_celltype, " (", receiver_DE_ident1, " vs ", if(!is.null(receiver_DE_ident2)) receiver_DE_ident2 else "all others", " in group ", receiver_DE_group_by,")")
    if (!receiver_DE_group_by %in% colnames(seurat_obj@meta.data)) stop("`receiver_DE_group_by` column '", receiver_DE_group_by, "' not found in Seurat object metadata.", call. = FALSE)
    DE_table_receiver <- Seurat::FindMarkers(seurat_obj, ident.1 = receiver_DE_ident1, ident.2 = receiver_DE_ident2, group.by = receiver_DE_group_by, subset.ident = receiver_celltype, assay = assay_name, min.pct = 0.10, logfc.threshold = 0 ) %>% # Use 0.10 for min.pct in FindMarkers too for consistency, or make it another param
      tibble::rownames_to_column("gene") %>%
      dplyr::filter(p_val_adj < p_val_adj_cutoff, avg_log2FC > logfc_cutoff)
    geneset_oi <- DE_table_receiver$gene[DE_table_receiver$gene %in% rownames(ligand_target_matrix)]
    if(length(geneset_oi) == 0) stop("No DEGs found in the receiver cell type that are also in the NicheNet ligand_target_matrix. Number of DEGs before NicheNet filter: ", nrow(DE_table_receiver), call. = FALSE)
    if (verbose) message(length(geneset_oi), " DEGs (upregulated in ident.1) from receiver '", receiver_celltype, "' identified for NicheNet analysis.")
    
    # --- IV. Define Potential Ligands ---
    if (verbose) message("Defining potential ligands...")
    ligands_db <- lr_network %>% dplyr::pull(from) %>% unique()
    receptors_db <- lr_network %>% dplyr::pull(to) %>% unique()
    expressed_ligands <- intersect(ligands_db, expressed_genes_sender)
    expressed_receptors <- intersect(receptors_db, expressed_genes_receiver)
    potential_ligands <- lr_network %>% dplyr::filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% dplyr::pull(from) %>% unique()
    if(length(potential_ligands) == 0) stop("No potential ligands found based on expression in sender/receiver and L-R database.", call. = FALSE)
    if (verbose) message(length(potential_ligands), " potential ligands identified.")
    
    # --- V. Perform NicheNet Ligand Activity Analysis ---
    if (verbose) message("Predicting ligand activities...")
    ligand_activities <- nichenetr::predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands) %>% dplyr::arrange(dplyr::desc(aupr_corrected))
    best_upstream_ligands <- ligand_activities %>% dplyr::top_n(min(top_n_ligands, nrow(ligand_activities)), aupr_corrected) %>% dplyr::pull(test_ligand) %>% unique()
    if (verbose) message("Top ", length(best_upstream_ligands), " ligands selected based on AUPR.")
    
    # --- VI. Infer Active Target Genes of Top Ligands ---
    # (No changes needed to the corrected block from previous interactions for this part's core logic)
    if (verbose) message("Inferring active target genes for top ligands...")
    active_ligand_target_links_df <- best_upstream_ligands %>%
      lapply(nichenetr::get_weighted_ligand_target_links, geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = top_n_targets_per_ligand) %>%
      dplyr::bind_rows() %>% tidyr::drop_na()
    p_ligand_target_network <- NULL
    if(nrow(active_ligand_target_links_df) > 0){
      active_ligand_target_links_matrix <- nichenetr::prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = ligand_target_cutoff)
      if (is.null(active_ligand_target_links_matrix) || nrow(active_ligand_target_links_matrix) == 0 || ncol(active_ligand_target_links_matrix) == 0) {
        warning("NN_Warning: Ligand-target matrix is empty/small after preparation. Skipping L-T heatmap.")
        p_ligand_target_network <- ggplot2::ggplot() + ggplot2::labs(title="No L-T links for heatmap after preparation.")
      } else { # (rest of the make.names and plotting logic for L-T heatmap - unchanged)
        ligands_for_ordering <- intersect(best_upstream_ligands, colnames(active_ligand_target_links_matrix)) %>% rev()
        targets_for_ordering <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links_matrix))
        rownames(active_ligand_target_links_matrix) <- make.names(rownames(active_ligand_target_links_matrix))
        colnames(active_ligand_target_links_matrix) <- make.names(colnames(active_ligand_target_links_matrix))
        order_ligands_final <- make.names(ligands_for_ordering)
        order_targets_final <- make.names(targets_for_ordering)
        order_ligands_final_intersect <- intersect(order_ligands_final, colnames(active_ligand_target_links_matrix))
        order_targets_final_intersect <- intersect(order_targets_final, rownames(active_ligand_target_links_matrix))
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
    # (No major changes to core logic)
    if (verbose) message("Inferring receptors for top ligands...")
    p_ligand_receptor_network <- NULL 
    lr_network_top_df_large <- dplyr::tibble() 
    lr_network_top <- lr_network %>% dplyr::filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% dplyr::distinct(from, to)
    best_upstream_receptors <- lr_network_top %>% dplyr::pull(to) %>% unique()
    if(length(best_upstream_ligands) > 0 && length(best_upstream_receptors) > 0){
      lr_network_top_df_large <- weighted_networks_lr %>% dplyr::filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
      valid_ligands_in_df <- intersect(best_upstream_ligands, unique(lr_network_top_df_large$from))
      valid_receptors_in_df <- intersect(best_upstream_receptors, unique(lr_network_top_df_large$to))
      if(length(valid_ligands_in_df) > 0 && length(valid_receptors_in_df) > 0) {
        lr_network_top_df_filtered <- lr_network_top_df_large %>% dplyr::filter(from %in% valid_ligands_in_df, to %in% valid_receptors_in_df)
        if (nrow(lr_network_top_df_filtered) > 0) { 
          lr_network_top_df <- lr_network_top_df_filtered %>% tidyr::spread(key = "from", value = "weight", fill = 0) 
          lr_network_top_matrix <- lr_network_top_df %>% dplyr::select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)
          if(nrow(lr_network_top_matrix) > 1 && ncol(lr_network_top_matrix) > 1) { # (Clustering logic unchanged)
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
    signaling_results <- NULL
    if(run_signaling_path_inference){
      if(requireNamespace("DiagrammeR", quietly = TRUE) && requireNamespace("DiagrammeRsvg", quietly = TRUE) && requireNamespace("rsvg", quietly = TRUE) ){
        if(length(best_upstream_ligands) > 0){ # (Logic for selecting targets_for_path and path inference unchanged but uses verbose for messages now)
          if (verbose) message("Inferring signaling paths for top ligand(s)...")
          ligands_for_path <- best_upstream_ligands[1]; targets_for_path <- character(0)
          if(exists("active_ligand_target_links_df") && is.data.frame(active_ligand_target_links_df) && nrow(active_ligand_target_links_df) > 0 && all(c("test_ligand", "target", "weight") %in% colnames(active_ligand_target_links_df))) {
            if (ligands_for_path %in% active_ligand_target_links_df$test_ligand) { targets_for_path <- active_ligand_target_links_df %>% dplyr::filter(test_ligand == ligands_for_path) %>% dplyr::top_n(5, weight) %>% dplyr::pull(target) %>% unique()
            } else { if (verbose) message("Top ligand '", ligands_for_path, "' not found in L-T links. Cannot select targets for signaling path.") }
          } else { if (verbose) message("L-T links df is empty/malformed. Cannot select targets for signaling path.") }
          if(length(targets_for_path) > 0 && !is.null(ligand_tf_matrix) && !is.null(ligand_target_matrix) && ligands_for_path %in% rownames(ligand_tf_matrix) && all(targets_for_path %in% colnames(ligand_target_matrix))) {
            active_signaling_network <- nichenetr::get_ligand_signaling_path(ligand_tf_matrix = ligand_tf_matrix, ligands_all = ligands_for_path, targets_all = targets_for_path, weighted_networks = weighted_networks )
            if(nrow(active_signaling_network$sig) > 0 || nrow(active_signaling_network$gr) > 0) { # (Normalization logic unchanged)
              if(nrow(active_signaling_network$sig) > 0 && sum(!is.na(active_signaling_network$sig$weight)) > 0) { min_w_sig <- min(active_signaling_network$sig$weight, na.rm=TRUE); max_w_sig <- max(active_signaling_network$sig$weight, na.rm=TRUE); if (min_w_sig < max_w_sig) { active_signaling_network$sig <- active_signaling_network$sig %>% dplyr::mutate(weight = ((weight - min_w_sig) / (max_w_sig - min_w_sig)) + 0.75) } else { active_signaling_network$sig <- active_signaling_network$sig %>% dplyr::mutate(weight = 1.0) }}
              if(nrow(active_signaling_network$gr) > 0 && sum(!is.na(active_signaling_network$gr$weight)) > 0) { min_w_gr <- min(active_signaling_network$gr$weight, na.rm=TRUE); max_w_gr <- max(active_signaling_network$gr$weight, na.rm=TRUE); if (min_w_gr < max_w_gr) { active_signaling_network$gr <- active_signaling_network$gr %>% dplyr::mutate(weight = ((weight - min_w_gr) / (max_w_gr - min_w_gr)) + 0.75) } else { active_signaling_network$gr <- active_signaling_network$gr %>% dplyr::mutate(weight = 1.0) }}
              signaling_results <- nichenetr::diagrammer_format_signaling_graph(signaling_graph_list = active_signaling_network, ligands_all = ligands_for_path, targets_all = targets_for_path, sig_color = "indianred", gr_color = "steelblue")
              if (verbose) message("Signaling path inference complete for ", ligands_for_path)
            } else { if(verbose) message("No signaling/GR interactions found for selected L-T pairs in signaling path inference."); signaling_results <- NULL }
          } else { if (verbose) message("Could not run signaling path: Conditions not met."); signaling_results <- NULL }
        } else { if (verbose) message("Skipping signaling path inference: No best_upstream_ligands."); signaling_results <- NULL }
      } else { if (verbose) message("Skipping signaling path inference: 'DiagrammeR', 'DiagrammeRsvg', or 'rsvg' package not found.") }
    }
    
    # --- X. (Optional) Circos Plot ---
    recorded_circos_plot <- NULL 
    circos_status_message <- "Circos plot not generated: run_circos was FALSE or required packages missing."
    if(run_circos){
      if(requireNamespace("circlize", quietly = TRUE) && requireNamespace("RColorBrewer", quietly = TRUE) && requireNamespace("Seurat", quietly=TRUE)){ # Added Seurat for GetAssayData
        library(circlize) 
        library(RColorBrewer)
        if (verbose) message("Generating Circos plot for ligand-receptor interactions with sender information...")
        
        if(exists("lr_network_top_df_large") && is.data.frame(lr_network_top_df_large) && nrow(lr_network_top_df_large) > 0 && all(c("from", "to", "weight") %in% colnames(lr_network_top_df_large)) ){
          circos_lr_data_raw <- lr_network_top_df_large %>% dplyr::filter(from %in% best_upstream_ligands) 
          
          if(nrow(circos_lr_data_raw) > 0) {
            # Determine expressing sender for each ligand in circos_lr_data_raw$from
            unique_ligands_in_circos <- unique(circos_lr_data_raw$from)
            ligand_sender_map <- sapply(unique_ligands_in_circos, function(lig) {
              expressing_senders <- c()
              for (s_type in sender_celltypes) {
                # Check expression in the specific sender cell type
                # Subsetting Seurat object for specific cell type
                # This requires DefaultAssay and Idents to be set correctly on seurat_obj
                cells_in_sender <- Seurat::Cells(seurat_obj)[Seurat::Idents(seurat_obj) == s_type]
                if (length(cells_in_sender) > 0) {
                  # Check if gene is in the assay for this object
                  if (lig %in% rownames(Seurat::GetAssayData(seurat_obj, assay = assay_name, slot = "data"))) {
                    gene_expr_in_sender <- Seurat::GetAssayData(seurat_obj, assay = assay_name, slot = "data")[lig, cells_in_sender, drop=FALSE]
                    if (sum(gene_expr_in_sender > 0) / length(cells_in_sender) >= min_pct_expressed) {
                      expressing_senders <- c(expressing_senders, s_type)
                    }
                  } # else: ligand not in assay data, treat as not expressed by this sender
                }
              }
              if (length(expressing_senders) > 0) return(expressing_senders[1]) # Take first for simplicity
              return("UnknownSender") 
            }, USE.NAMES = TRUE)
            
            # Create a data frame for circos: gene, group, display_name
            ligand_df_circos <- tibble::tibble(
              item = unique_ligands_in_circos,
              group = paste0(ligand_sender_map[item], "_L"), # e.g., CD4_T_Cell_L
              color_group = ligand_sender_map[item] # For coloring by sender
            )
            receptor_df_circos <- tibble::tibble(
              item = unique(circos_lr_data_raw$to),
              group = "Receptor_R", # Single group for all receptors
              color_group = "Receptor"
            )
            all_items_df_circos <- dplyr::bind_rows(ligand_df_circos, receptor_df_circos) %>% dplyr::distinct(item, .keep_all = TRUE)
            
            # Prepare links (from ligand to receptor)
            links_df_circos <- circos_lr_data_raw %>% 
              dplyr::select(from, to, weight) %>%
              dplyr::rename(origin = from, target = to) # For clarity
            
            # Define colors for sender groups and receptor group
            unique_sender_groups_for_color <- unique(ligand_df_circos$color_group)
            sender_palette <- RColorBrewer::brewer.pal(n = min(length(unique_sender_groups_for_color), 9), name = "Set1")
            if(length(unique_sender_groups_for_color) > length(sender_palette) && length(sender_palette) > 0) sender_palette <- rep(sender_palette, length.out = length(unique_sender_groups_for_color))
            names(sender_palette) <- unique_sender_groups_for_color
            
            receptor_color <- RColorBrewer::brewer.pal(n = 3, name = "Set2")[1] # Pick one color for receptors
            
            # Full color map for all items based on their 'color_group'
            item_colors <- sapply(all_items_df_circos$item, function(it) {
              cg <- all_items_df_circos$color_group[all_items_df_circos$item == it]
              if (cg == "Receptor") return(receptor_color)
              return(sender_palette[cg] %||% "grey50") # Fallback color
            }, USE.NAMES = TRUE)
            
            plot_circos_func <- function() {
              circos.clear()
              # Gap between sender groups and receptor group
              # Define order of groups: Senders first, then Receptor
              group_order <- c(paste0(unique(ligand_sender_map[unique_ligands_in_circos]), "_L"), "Receptor_R")
              
              # Create factors for proper ordering of sectors if needed
              # all_items_df_circos$group_factor <- factor(all_items_df_circos$group, levels = group_order)
              # For chordDiagram, the 'group' argument can be a named character vector for items
              # Or, if items_df_circos provides the grouping directly by having a 'group' column with items.
              
              # Prepare data for chordDiagram in the required format: a data frame of links,
              # and optionally a data frame or named vector for sector order/grouping.
              # The items themselves (ligands, receptors) will form the sectors.
              # We use 'group' argument of chordDiagram to color sectors by a higher-level group (Sender/Receptor).
              # Let's try to provide item_colors directly to grid.col
              
              # Create the order for items (sectors)
              # Order ligands by their sender group, then alphabetically
              # Order receptors alphabetically
              ordered_ligands <- ligand_df_circos %>% dplyr::arrange(group, item) %>% dplyr::pull(item)
              ordered_receptors <- receptor_df_circos %>% dplyr::arrange(item) %>% dplyr::pull(item)
              sector_order <- c(ordered_ligands, ordered_receptors)
              
              # Make sure item_colors are in this sector_order
              final_item_colors_ordered <- item_colors[sector_order]
              
              chordDiagram(links_df_circos,
                           order = sector_order, # Order of sectors
                           grid.col = final_item_colors_ordered, # Colors for each sector (item)
                           transparency = 0.25,
                           annotationTrack = "grid",
                           preAllocateTracks = list(track.height = circlize::mm_h(4)),
                           directional = 1, 
                           direction.type = "arrows",
                           link.arr.type = "big.arrow",
                           link.arr.length = 0.05, link.arr.width = 0.05 
              )
              circos.track(track.index = 1, panel.fun = function(x, y) {
                circos.text(CELL_META$xcenter, CELL_META$ylim[1] + circlize::uy(5, "mm"), 
                            CELL_META$sector.index, # This is the item name
                            facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex=0.4) # Reduced cex further
              }, bg.border = NA) 
              
              # Add legend manually (ComplexHeatmap::Legend is good for this but adds dependency)
              # Simple text legend for now if verbose and outputting to file
              if(verbose && !is.null(current_run_output_dir)) {
                legend_text <- "Circos Legend:\n"
                for(s_group in names(sender_palette)){
                  legend_text <- paste0(legend_text, s_group, " (Ligands): ", sender_palette[s_group], "\n")
                }
                legend_text <- paste0(legend_text, "Receptor: ", receptor_color, "\n")
                cat(legend_text, file=file.path(current_run_output_dir, "NicheNet_Circos_LR_legend.txt"))
              }
            }
            
            tryCatch({ # (Saving and recording logic unchanged from previous response)
              if(!is.null(current_run_output_dir)){
                pdf(file.path(current_run_output_dir, "NicheNet_Circos_LR.pdf"), width=10, height=10); plot_circos_func(); dev.off()
                png(file.path(current_run_output_dir, "NicheNet_Circos_LR.png"), width=1000, height=1000, res=100); plot_circos_func(); dev.off()
                if(verbose) message("Circos plot & legend info saved.")
              }
              plot_circos_func(); recorded_circos_plot <- recordPlot(); circos.clear()
              circos_status_message <- "Circos plot generated and recorded (sender-grouped)."
            }, error = function(e) {
              if (verbose) message("Could not generate or save Sender-grouped Circos plot: ", e$message)
              circos_status_message <- paste("Sender-grouped Circos plot failed:", e$message)
              if(grDevices::dev.cur() != 1) grDevices::dev.off()
            })
          } else { if (verbose) message("No data for Circos plot after filtering."); circos_status_message <- "No data for Circos plot after filtering." }
        } else { if (verbose) message("Skipping Circos plot: lr_network_top_df_large data unavailable/malformed."); circos_status_message <- "Circos plot skipped: lr_network_top_df_large data." }
      } else { circos_status_message <- "Circos plot skipped: 'circlize', 'RColorBrewer', or 'Seurat' package not fully available."; if (verbose) message(circos_status_message) }
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
      signaling_path_graph = signaling_results,
      nichenet_data_loaded_info = nichenet_data_loaded_info,
      output_path_run = if(!is.null(current_run_output_dir)) current_run_output_dir else "Output not saved to files."
      # Warnings will be added outside withCallingHandlers
    )
    
    # Save ggplot objects
    if (!is.null(current_run_output_dir)) {
      # (Plot saving logic unchanged from previous response)
      plot_list_gg <- list(NicheNet_Ligand_Target_Heatmap = results_list$plot_ligand_target_network, NicheNet_Ligand_Receptor_Heatmap = results_list$plot_ligand_receptor_network, NicheNet_Ligand_Activity_Histogram = results_list$plot_ligand_activity_hist, NicheNet_Ligand_AUPR_Heatmap = results_list$plot_ligand_aupr_heatmap)
      plot_dims <- list(NicheNet_Ligand_Target_Heatmap = list(w=10,h=8), NicheNet_Ligand_Receptor_Heatmap = list(w=8,h=6), NicheNet_Ligand_Activity_Histogram = list(w=6,h=5), NicheNet_Ligand_AUPR_Heatmap = list(w=6,h=8))
      for(plot_name in names(plot_list_gg)){ p_obj <- plot_list_gg[[plot_name]]; if(!is.null(p_obj) && inherits(p_obj, "ggplot")){ tryCatch(ggplot2::ggsave(filename = file.path(current_run_output_dir, paste0(plot_name, ".png")), plot = p_obj, width = plot_dims[[plot_name]]$w, height = plot_dims[[plot_name]]$h), error = function(e) warning("NN_Warning: Failed to save plot ", plot_name, ": ", e$message))}}
      if(run_signaling_path_inference && !is.null(results_list$signaling_path_graph) && requireNamespace("DiagrammeRsvg", quietly = TRUE) && requireNamespace("rsvg", quietly = TRUE)){ tryCatch({ DiagrammeRsvg::export_svg(results_list$signaling_path_graph) %>% charToRaw() %>% rsvg::rsvg_png(file.path(current_run_output_dir, "NicheNet_Signaling_Path.png")); if (verbose) message("Signaling path graph saved.") }, error = function(e) warning("NN_Warning: Failed to save signaling path graph: ", e$message)) }
    }
    return(results_list) # Return results from within withCallingHandlers
  }, 
  warning = function(w) {
    # Check if it's one of my custom warnings or a general R warning
    if (startsWith(conditionMessage(w), "NN_Warning:")) {
      msg <- sub("^NN_Warning: ", "", conditionMessage(w))
      warnings_collected <<- c(warnings_collected, paste("NicheNet Function:", msg))
    } else {
      # For general R warnings, capture the message as is
      warnings_collected <<- c(warnings_collected, paste("R Warning:", conditionMessage(w)))
    }
    invokeRestart("muffleWarning") # Suppress original warning from printing
  }) # End of withCallingHandlers
  
  # Add collected warnings to the final results
  analysis_results$warnings_occurred <- if(length(warnings_collected) > 0) warnings_collected else "No warnings reported."
  
  if (verbose && length(warnings_collected) > 0) {
    message("--------------------------------------------------------------------")
    message("NicheNet Analysis completed with ", length(warnings_collected), " warning(s). See 'warnings_occurred' in output list for details.")
    message("Summary of warnings:")
    for(i in seq_along(warnings_collected)){
      message(paste0(i, ". ", warnings_collected[i]))
      if(i >= 5 && length(warnings_collected) > 5) { # Print max 5 warnings in summary
        message("... (and ", length(warnings_collected) - 5, " more warnings)")
        break
      }
    }
    message("--------------------------------------------------------------------")
  } else if (verbose) {
    message("NicheNet Analysis completed successfully with no warnings reported.")
  }
  
  return(analysis_results)
}

# Helper for path construction (from nichenetr utils, or base R `||` if appropriate)
"%||%" <- function(a, b) {
  if (!is.null(a) && !is.logical(a) && (is.character(a) || is.numeric(a)) && length(a) == 1 && !is.na(a) && nzchar(as.character(a)) ) a else b
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
#'        Default: "NicheNet_Global_Dataset".
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
#' @importFrom ggplot2 ggplot aes geom_histogram geom_vline labs theme_classic theme element_text scale_fill_gradient2 RotatedAxis ggsave
#' @importFrom utils download.file capture.output
#' @importFrom stats dist hclust
#' @importFrom grDevices recordPlot dev.off pdf png
#' @importFrom DiagrammeRsvg export_svg
#' @importFrom rsvg rsvg_png
#' @importFrom RColorBrewer brewer.pal
#' @importFrom circlize circos.clear chordDiagram circos.track circos.text uy mm_h CELL_META
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
                                  min_pct_expressed = 0.10,
                                  p_val_adj_cutoff = 0.05,
                                  logfc_cutoff = 0.25,
                                  top_n_ligands = 20,
                                  top_n_targets_per_ligand = 200,
                                  ligand_target_cutoff = 0.33,
                                  nichenet_data_dir = NULL,
                                  nichenet_data_name = "NicheNet_Global_Dataset",
                                  output_dir = NULL,
                                  run_circos = TRUE,
                                  run_signaling_path_inference = FALSE,
                                  verbose = TRUE) {
  
  # --- 0. Initialize Warning Collection & Argument Setup ---
  warnings_collected <- character(0)
  circos_notes_collected <- character(0) # For circos specific notes
  
  analysis_results <- withCallingHandlers({
    call_obj <- match.call()
    species <- match.arg(species)
    
    required_pkgs <- c("nichenetr", "Seurat", "dplyr", "tidyr", "tibble", "ggplot2")
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
        con_args <- file(args_file_path, "w")
        cat(paste(utils::capture.output(print(call_obj)), collapse = "\n"), file = con_args)
        close(con_args)
      }, error = function(e) warning("NN_Warning: Failed to save arguments to file: ", args_file_path, ". Error: ", e$message))
      if (verbose && file.exists(args_file_path)) message("Arguments saved to: ", args_file_path)
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
        data_dir_to_use <- file.path(getwd(), paste0("nichenet_data_", species))
        if (verbose) message("`nichenet_data_dir` is NULL, using default local cache: ", data_dir_to_use)
      }
      if (!dir.exists(data_dir_to_use)) {
        if (verbose) message("Creating NicheNet data directory: ", data_dir_to_use)
        tryCatch(dir.create(data_dir_to_use, recursive = TRUE, showWarnings = FALSE),
                 error = function(e) stop("Failed to create NicheNet data directory: ", data_dir_to_use, ". Error: ", e$message))
      }
      loaded_data_for_global_assignment <- list()
      for (simple_name in names(required_files_map)) {
        rds_filename <- required_files_map[[simple_name]]
        file_path <- file.path(data_dir_to_use, rds_filename)
        if (!file.exists(file_path)) {
          if (verbose) message("File ", rds_filename, " not found in '", data_dir_to_use, "'. Downloading from Zenodo...")
          download_url <- paste0(base_url, rds_filename, "?download=1")
          tryCatch({
            utils::download.file(url = download_url, destfile = file_path, mode = "wb", quiet = !verbose)
            if (verbose) message("Downloaded ", rds_filename, " successfully to ", data_dir_to_use)
          }, error = function(e) {
            stop("Failed to download ", rds_filename, ". Error: ", e$message, call. = FALSE)
          })
        }
        assign(simple_name, readRDS(file_path), envir = temp_data_env)
        loaded_data_for_global_assignment[[simple_name]] <- get(simple_name, envir = temp_data_env)
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
    Seurat::DefaultAssay(seurat_obj) <- assay_name
    Seurat::Idents(seurat_obj) <- cluster_col
    if(!all(sender_celltypes %in% levels(Seurat::Idents(seurat_obj)))) stop("One or more sender_celltypes not found in Seurat Idents (column '", cluster_col, "').", call. = FALSE)
    if(!receiver_celltype %in% levels(Seurat::Idents(seurat_obj))) stop("Receiver_celltype '", receiver_celltype, "' not found in Seurat Idents (column '", cluster_col, "').", call. = FALSE)
    
    list_expressed_genes_sender <- lapply(sender_celltypes, nichenetr::get_expressed_genes, seurat_obj, min_pct_expressed)
    expressed_genes_sender <- unique(unlist(list_expressed_genes_sender))
    expressed_genes_receiver <- nichenetr::get_expressed_genes(receiver_celltype, seurat_obj, min_pct_expressed)
    background_expressed_genes <- expressed_genes_receiver[expressed_genes_receiver %in% rownames(ligand_target_matrix)]
    if(length(background_expressed_genes) < 20) warning("NN_Warning: Very few (<20) background expressed genes found in receiver (", length(background_expressed_genes), "). NicheNet might not perform well.")
    
    # --- III. Define Gene Set of Interest (DEGs in receiver) ---
    if (verbose) message("Performing DE analysis for receiver: ", receiver_celltype, " (", receiver_DE_ident1, " vs ", if(!is.null(receiver_DE_ident2)) receiver_DE_ident2 else "all others", " in group ", receiver_DE_group_by,")")
    if (!receiver_DE_group_by %in% colnames(seurat_obj@meta.data)) stop("`receiver_DE_group_by` column '", receiver_DE_group_by, "' not found in Seurat object metadata.", call. = FALSE)
    DE_table_receiver <- Seurat::FindMarkers(seurat_obj, ident.1 = receiver_DE_ident1, ident.2 = receiver_DE_ident2, group.by = receiver_DE_group_by, subset.ident = receiver_celltype, assay = assay_name, min.pct = min_pct_expressed, logfc.threshold = 0 ) %>%
      tibble::rownames_to_column("gene") %>%
      dplyr::filter(p_val_adj < p_val_adj_cutoff, avg_log2FC > logfc_cutoff)
    geneset_oi <- DE_table_receiver$gene[DE_table_receiver$gene %in% rownames(ligand_target_matrix)]
    if(length(geneset_oi) == 0) stop("No DEGs found in the receiver cell type that are also in the NicheNet ligand_target_matrix. Number of DEGs before NicheNet filter: ", nrow(DE_table_receiver), call. = FALSE)
    if (verbose) message(length(geneset_oi), " DEGs (upregulated in ident.1) from receiver '", receiver_celltype, "' identified for NicheNet analysis.")
    
    # --- IV. Define Potential Ligands ---
    if (verbose) message("Defining potential ligands...")
    ligands_db <- lr_network %>% dplyr::pull(from) %>% unique()
    receptors_db <- lr_network %>% dplyr::pull(to) %>% unique()
    expressed_ligands <- intersect(ligands_db, expressed_genes_sender)
    expressed_receptors <- intersect(receptors_db, expressed_genes_receiver)
    potential_ligands <- lr_network %>% dplyr::filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% dplyr::pull(from) %>% unique()
    if(length(potential_ligands) == 0) stop("No potential ligands found based on expression in sender/receiver and L-R database.", call. = FALSE)
    if (verbose) message(length(potential_ligands), " potential ligands identified.")
    
    # --- V. Perform NicheNet Ligand Activity Analysis ---
    if (verbose) message("Predicting ligand activities...")
    ligand_activities <- nichenetr::predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands) %>% dplyr::arrange(dplyr::desc(aupr_corrected))
    best_upstream_ligands <- ligand_activities %>% dplyr::top_n(min(top_n_ligands, nrow(ligand_activities)), aupr_corrected) %>% dplyr::pull(test_ligand) %>% unique()
    if (verbose) message("Top ", length(best_upstream_ligands), " ligands selected based on AUPR.")
    
    # --- VI. Infer Active Target Genes of Top Ligands ---
    if (verbose) message("Inferring active target genes for top ligands...")
    active_ligand_target_links_df <- best_upstream_ligands %>%
      lapply(nichenetr::get_weighted_ligand_target_links, geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = top_n_targets_per_ligand) %>%
      dplyr::bind_rows() %>% tidyr::drop_na()
    p_ligand_target_network <- NULL
    if(nrow(active_ligand_target_links_df) > 0){
      active_ligand_target_links_matrix <- nichenetr::prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = ligand_target_cutoff)
      if (is.null(active_ligand_target_links_matrix) || nrow(active_ligand_target_links_matrix) == 0 || ncol(active_ligand_target_links_matrix) == 0) {
        warning("NN_Warning: Ligand-target matrix is empty/small after preparation. Skipping L-T heatmap.")
        p_ligand_target_network <- ggplot2::ggplot() + ggplot2::labs(title="No L-T links for heatmap after preparation.")
      } else {
        ligands_for_ordering <- intersect(best_upstream_ligands, colnames(active_ligand_target_links_matrix)) %>% rev()
        targets_for_ordering <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links_matrix))
        rownames(active_ligand_target_links_matrix) <- make.names(rownames(active_ligand_target_links_matrix))
        colnames(active_ligand_target_links_matrix) <- make.names(colnames(active_ligand_target_links_matrix))
        order_ligands_final <- make.names(ligands_for_ordering)
        order_targets_final <- make.names(targets_for_ordering)
        order_ligands_final_intersect <- intersect(order_ligands_final, colnames(active_ligand_target_links_matrix))
        order_targets_final_intersect <- intersect(order_targets_final, rownames(active_ligand_target_links_matrix))
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
      valid_ligands_in_df <- intersect(best_upstream_ligands, unique(lr_network_top_df_large$from))
      valid_receptors_in_df <- intersect(best_upstream_receptors, unique(lr_network_top_df_large$to))
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
    signaling_results <- NULL
    if(run_signaling_path_inference){
      if(requireNamespace("DiagrammeR", quietly = TRUE) && requireNamespace("DiagrammeRsvg", quietly = TRUE) && requireNamespace("rsvg", quietly = TRUE) ){
        if(length(best_upstream_ligands) > 0){
          if (verbose) message("Inferring signaling paths for top ligand(s)...")
          ligands_for_path <- best_upstream_ligands[1]; targets_for_path <- character(0)
          if(exists("active_ligand_target_links_df") && is.data.frame(active_ligand_target_links_df) && nrow(active_ligand_target_links_df) > 0 && all(c("test_ligand", "target", "weight") %in% colnames(active_ligand_target_links_df))) {
            if (ligands_for_path %in% active_ligand_target_links_df$test_ligand) { targets_for_path <- active_ligand_target_links_df %>% dplyr::filter(test_ligand == ligands_for_path) %>% dplyr::top_n(5, weight) %>% dplyr::pull(target) %>% unique()
            } else { if (verbose) message("Top ligand '", ligands_for_path, "' not found in L-T links. Cannot select targets for signaling path.") }
          } else { if (verbose) message("L-T links df is empty/malformed. Cannot select targets for signaling path.") }
          if(length(targets_for_path) > 0 && !is.null(ligand_tf_matrix) && !is.null(ligand_target_matrix) && ligands_for_path %in% rownames(ligand_tf_matrix) && all(targets_for_path %in% colnames(ligand_target_matrix))) {
            active_signaling_network <- nichenetr::get_ligand_signaling_path(ligand_tf_matrix = ligand_tf_matrix, ligands_all = ligands_for_path, targets_all = targets_for_path, weighted_networks = weighted_networks )
            if(nrow(active_signaling_network$sig) > 0 || nrow(active_signaling_network$gr) > 0) {
              if(nrow(active_signaling_network$sig) > 0 && sum(!is.na(active_signaling_network$sig$weight)) > 0) { min_w_sig <- min(active_signaling_network$sig$weight, na.rm=TRUE); max_w_sig <- max(active_signaling_network$sig$weight, na.rm=TRUE); if (min_w_sig < max_w_sig) { active_signaling_network$sig <- active_signaling_network$sig %>% dplyr::mutate(weight = ((weight - min_w_sig) / (max_w_sig - min_w_sig)) + 0.75) } else { active_signaling_network$sig <- active_signaling_network$sig %>% dplyr::mutate(weight = 1.0) }}
              if(nrow(active_signaling_network$gr) > 0 && sum(!is.na(active_signaling_network$gr$weight)) > 0) { min_w_gr <- min(active_signaling_network$gr$weight, na.rm=TRUE); max_w_gr <- max(active_signaling_network$gr$weight, na.rm=TRUE); if (min_w_gr < max_w_gr) { active_signaling_network$gr <- active_signaling_network$gr %>% dplyr::mutate(weight = ((weight - min_w_gr) / (max_w_gr - min_w_gr)) + 0.75) } else { active_signaling_network$gr <- active_signaling_network$gr %>% dplyr::mutate(weight = 1.0) }}
              signaling_results <- nichenetr::diagrammer_format_signaling_graph(signaling_graph_list = active_signaling_network, ligands_all = ligands_for_path, targets_all = targets_for_path, sig_color = "indianred", gr_color = "steelblue")
              if (verbose) message("Signaling path inference complete for ", ligands_for_path)
            } else { if(verbose) message("No signaling/GR interactions found for selected L-T pairs in signaling path inference."); signaling_results <- NULL }
          } else { if (verbose) message("Could not run signaling path: Conditions not met."); signaling_results <- NULL }
        } else { if (verbose) message("Skipping signaling path inference: No best_upstream_ligands."); signaling_results <- NULL }
      } else { if (verbose) message("Skipping signaling path inference: 'DiagrammeR', 'DiagrammeRsvg', or 'rsvg' package not found.") }
    }
    
    # --- X. (Optional) Circos Plot ---
    recorded_circos_plot <- NULL
    temp_circos_notes <- character(0) # Local temporary storage for notes
    circos_status_message <- "Circos plot not generated: run_circos was FALSE or required packages missing."
    
    if(run_circos){
      if(requireNamespace("circlize", quietly = TRUE) &&
         requireNamespace("RColorBrewer", quietly = TRUE) &&
         requireNamespace("Seurat", quietly=TRUE) && # Already imported but good to check
         requireNamespace("evaluate", quietly=TRUE) ){ # For capturing notes
        
        if (verbose) message("Generating Circos plot for ligand-receptor interactions with sender information...")
        
        if(exists("lr_network_top_df_large") && is.data.frame(lr_network_top_df_large) && nrow(lr_network_top_df_large) > 0 && all(c("from", "to", "weight") %in% colnames(lr_network_top_df_large)) ){
          circos_lr_data_raw <- lr_network_top_df_large %>% dplyr::filter(from %in% best_upstream_ligands)
          
          if(nrow(circos_lr_data_raw) > 0) {
            unique_ligands_in_circos <- unique(circos_lr_data_raw$from)
            ligand_sender_map <- sapply(unique_ligands_in_circos, function(lig) {
              expressing_senders <- c()
              for (s_type in sender_celltypes) {
                cells_in_sender <- Seurat::Cells(seurat_obj)[Seurat::Idents(seurat_obj) == s_type]
                if (length(cells_in_sender) > 0) {
                  if (lig %in% rownames(Seurat::GetAssayData(seurat_obj, assay = assay_name, slot = "data"))) {
                    gene_expr_in_sender <- Seurat::GetAssayData(seurat_obj, assay = assay_name, slot = "data")[lig, cells_in_sender, drop=FALSE]
                    if (sum(gene_expr_in_sender > 0) / length(cells_in_sender) >= min_pct_expressed) {
                      expressing_senders <- c(expressing_senders, s_type)
                    }
                  }
                }
              }
              if (length(expressing_senders) > 0) return(expressing_senders[1])
              return("UnknownSender")
            }, USE.NAMES = TRUE)
            
            ligand_df_circos <- tibble::tibble(
              item = unique_ligands_in_circos,
              group = paste0(ligand_sender_map[item], "_L"),
              color_group = ligand_sender_map[item]
            )
            receptor_df_circos <- tibble::tibble(
              item = unique(circos_lr_data_raw$to),
              group = "Receptor_R",
              color_group = "Receptor"
            )
            all_items_df_circos <- dplyr::bind_rows(ligand_df_circos, receptor_df_circos) %>% dplyr::distinct(item, .keep_all = TRUE)
            links_df_circos <- circos_lr_data_raw %>%
              dplyr::select(from, to, weight) %>%
              dplyr::rename(origin = from, target = to)
            
            unique_sender_groups_for_color <- unique(ligand_df_circos$color_group)
            # Ensure "UnknownSender" gets a color if present
            if ("UnknownSender" %in% unique_sender_groups_for_color && length(unique_sender_groups_for_color) == 1) {
              sender_palette_values <- "grey50" # Special case for only UnknownSender
            } else {
              n_senders_for_palette <- length(unique_sender_groups_for_color[unique_sender_groups_for_color != "UnknownSender"])
              if (n_senders_for_palette > 0) {
                base_palette <- RColorBrewer::brewer.pal(n = min(n_senders_for_palette, 9), name = "Set1")
                if(n_senders_for_palette > length(base_palette)) base_palette <- rep(base_palette, length.out = n_senders_for_palette)
                sender_palette_values <- base_palette
              } else {
                sender_palette_values <- character(0) # No non-unknown senders
              }
            }
            
            sender_palette <- stats::setNames(sender_palette_values, unique_sender_groups_for_color[unique_sender_groups_for_color != "UnknownSender"])
            if ("UnknownSender" %in% unique_sender_groups_for_color) {
              sender_palette["UnknownSender"] <- "grey50"
            }
            
            
            receptor_color <- RColorBrewer::brewer.pal(n = 3, name = "Set2")[1]
            item_colors <- sapply(all_items_df_circos$item, function(it) {
              cg <- all_items_df_circos$color_group[all_items_df_circos$item == it]
              if (cg == "Receptor") return(receptor_color)
              return(sender_palette[cg] %||% "grey70") # Fallback for any unmapped sender
            }, USE.NAMES = TRUE)
            
            plot_circos_func_env <- new.env() # Create a new environment for the function
            # Assign necessary variables to this environment
            plot_circos_func_env$links_df_circos <- links_df_circos
            plot_circos_func_env$ligand_df_circos <- ligand_df_circos
            plot_circos_func_env$receptor_df_circos <- receptor_df_circos
            plot_circos_func_env$item_colors <- item_colors
            plot_circos_func_env$verbose <- verbose # Pass verbose
            plot_circos_func_env$current_run_output_dir <- current_run_output_dir # Pass output dir
            plot_circos_func_env$sender_palette <- sender_palette
            plot_circos_func_env$receptor_color <- receptor_color
            
            
            plot_circos_func <- function() {
              # This function will now be executed within plot_circos_func_env by evaluate
              circlize::circos.clear()
              ordered_ligands <- get("ligand_df_circos", envir=plot_circos_func_env) %>% dplyr::arrange(group, item) %>% dplyr::pull(item)
              ordered_receptors <- get("receptor_df_circos", envir=plot_circos_func_env) %>% dplyr::arrange(item) %>% dplyr::pull(item)
              sector_order <- c(ordered_ligands, ordered_receptors)
              final_item_colors_ordered <- get("item_colors", envir=plot_circos_func_env)[sector_order]
              
              circlize::chordDiagram(get("links_df_circos", envir=plot_circos_func_env),
                                     order = sector_order,
                                     grid.col = final_item_colors_ordered,
                                     transparency = 0.25,
                                     annotationTrack = "grid",
                                     preAllocateTracks = list(track.height = circlize::mm_h(4)),
                                     directional = 1,
                                     direction.type = "arrows",
                                     link.arr.type = "big.arrow",
                                     link.arr.length = 0.05, link.arr.width = 0.05
              )
              circlize::circos.track(track.index = 1, panel.fun = function(x, y) {
                circlize::circos.text(circlize::CELL_META$xcenter, circlize::CELL_META$ylim[1] + circlize::uy(5, "mm"),
                                      circlize::CELL_META$sector.index,
                                      facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex=0.4)
              }, bg.border = NA)
              
              # Legend generation (text file)
              temp_current_run_output_dir <- get("current_run_output_dir", envir=plot_circos_func_env)
              if(get("verbose", envir=plot_circos_func_env) && !is.null(temp_current_run_output_dir)) {
                legend_text <- "Circos Legend (Sender Cell Type for Ligands):\n"
                temp_sender_palette <- get("sender_palette", envir=plot_circos_func_env)
                for(s_group in names(temp_sender_palette)){ # Iterate over named palette
                  if (s_group != "UnknownSender") { # Exclude UnknownSender from primary list if desired, or handle separately
                    legend_text <- paste0(legend_text, s_group, " (Ligands): Color ", temp_sender_palette[s_group], "\n")
                  }
                }
                if ("UnknownSender" %in% names(temp_sender_palette)) {
                  legend_text <- paste0(legend_text, "UnknownSender (Ligands): Color ", temp_sender_palette["UnknownSender"], "\n")
                }
                legend_text <- paste0(legend_text, "Receptor: Color ", get("receptor_color", envir=plot_circos_func_env), "\n")
                tryCatch(
                  cat(legend_text, file=file.path(temp_current_run_output_dir, "NicheNet_Circos_LR_legend.txt")),
                  error = function(e_cat) warning("NN_Warning: Failed to write circos legend file: ", e_cat$message)
                )
                
              }
            }
            # End of plot_circos_func definition
            
            tryCatch({
              # Capture all output from plot_circos_func using evaluate
              # The plot itself will render to the active graphics device if any
              # Suppress plot rendering during evaluate by setting a temporary null device if needed,
              # then replay for recordPlot. For now, direct call for simplicity if interactive.
              # If running non-interactively, pdf/png will capture it.
              
              # Use a temporary file for capturing output of evaluate
              tmp_out_file <- tempfile()
              # Execute and capture output
              captured_eval <- evaluate::evaluate(
                quote(plot_circos_func()), # plot_circos_func will be found in plot_circos_func_env
                envir = plot_circos_func_env,
                output_handler = evaluate::new_output_handler(
                  value = function(x, visible=TRUE) x, # Keep the return value of the function (if any)
                  text = function(x) cat(x, file = tmp_out_file, append = TRUE, sep = "\n"),
                  message = function(m) cat(paste("Message:", conditionMessage(m)), file = tmp_out_file, append = TRUE, sep = "\n"),
                  warning = function(w) cat(paste("Warning:", conditionMessage(w)), file = tmp_out_file, append = TRUE, sep = "\n"),
                  error = function(e) cat(paste("Error:", conditionMessage(e)), file = tmp_out_file, append = TRUE, sep = "\n")
                )
              )
              
              all_captured_output <- if(file.exists(tmp_out_file)) readLines(tmp_out_file, warn = FALSE) else character(0)
              if(file.exists(tmp_out_file)) unlink(tmp_out_file)
              
              # Filter for "Note:" messages specifically from circlize
              temp_circos_notes <- grep("^Note:", all_captured_output, value = TRUE)
              
              
              if(!is.null(current_run_output_dir)){
                grDevices::pdf(file.path(current_run_output_dir, "NicheNet_Circos_LR.pdf"), width=10, height=10); plot_circos_func_env$plot_circos_func(); grDevices::dev.off()
                grDevices::png(file.path(current_run_output_dir, "NicheNet_Circos_LR.png"), width=1000, height=1000, res=100); plot_circos_func_env$plot_circos_func(); grDevices::dev.off()
                if(verbose) message("Circos plot & legend info saved.")
              }
              
              # Plot to active device to record it
              plot_circos_func_env$plot_circos_func()
              recorded_circos_plot <- grDevices::recordPlot()
              circlize::circos.clear() # Crucial after recordPlot
              
              circos_status_message <- "Circos plot generated and recorded (sender-grouped)."
            }, error = function(e) {
              if (verbose) message("Could not generate or save Sender-grouped Circos plot: ", e$message)
              circos_status_message <- paste("Sender-grouped Circos plot failed:", e$message)
              if(grDevices::dev.cur() != 1) grDevices::dev.off() # Close device if plotting failed mid-way
              # Also capture error message into notes if it's a circlize specific one.
              # temp_circos_notes <<- c(temp_circos_notes, paste("Error during circos generation:", e$message)) # Use <<- if outside handler
            })
            # Assign collected notes to the parent scope variable
            circos_notes_collected <<- c(circos_notes_collected, temp_circos_notes)
            
          } else { if (verbose) message("No data for Circos plot after filtering."); circos_status_message <- "No data for Circos plot after filtering." }
        } else { if (verbose) message("Skipping Circos plot: lr_network_top_df_large data unavailable/malformed."); circos_status_message <- "Circos plot skipped: lr_network_top_df_large data." }
      } else {
        circos_status_message <- "Circos plot skipped: 'circlize', 'RColorBrewer', 'Seurat', or 'evaluate' package not fully available."
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
      circos_notes = if(length(circos_notes_collected) > 0) circos_notes_collected else "No notes from circos plot generation.", # Add collected notes
      signaling_path_graph = signaling_results,
      nichenet_data_loaded_info = nichenet_data_loaded_info,
      output_path_run = if(!is.null(current_run_output_dir)) current_run_output_dir else "Output not saved to files."
    )
    
    if (!is.null(current_run_output_dir)) {
      plot_list_gg <- list(NicheNet_Ligand_Target_Heatmap = results_list$plot_ligand_target_network, NicheNet_Ligand_Receptor_Heatmap = results_list$plot_ligand_receptor_network, NicheNet_Ligand_Activity_Histogram = results_list$plot_ligand_activity_hist, NicheNet_Ligand_AUPR_Heatmap = results_list$plot_ligand_aupr_heatmap)
      plot_dims <- list(NicheNet_Ligand_Target_Heatmap = list(w=10,h=8), NicheNet_Ligand_Receptor_Heatmap = list(w=8,h=6), NicheNet_Ligand_Activity_Histogram = list(w=6,h=5), NicheNet_Ligand_AUPR_Heatmap = list(w=6,h=8))
      for(plot_name in names(plot_list_gg)){ p_obj <- plot_list_gg[[plot_name]]; if(!is.null(p_obj) && inherits(p_obj, "ggplot")){ tryCatch(ggplot2::ggsave(filename = file.path(current_run_output_dir, paste0(plot_name, ".png")), plot = p_obj, width = plot_dims[[plot_name]]$w, height = plot_dims[[plot_name]]$h), error = function(e) warning("NN_Warning: Failed to save plot ", plot_name, ": ", e$message))}}
      if(run_signaling_path_inference && !is.null(results_list$signaling_path_graph) && requireNamespace("DiagrammeRsvg", quietly = TRUE) && requireNamespace("rsvg", quietly = TRUE)){ tryCatch({ DiagrammeRsvg::export_svg(results_list$signaling_path_graph) %>% charToRaw() %>% rsvg::rsvg_png(file.path(current_run_output_dir, "NicheNet_Signaling_Path.png")); if (verbose) message("Signaling path graph saved.") }, error = function(e) warning("NN_Warning: Failed to save signaling path graph: ", e$message)) }
    }
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
  
  analysis_results$warnings_occurred <- if(length(warnings_collected) > 0) warnings_collected else "No warnings reported."
  # Ensure circos_notes is part of the final returned list, even if an error occurred before it was assigned
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
  if (verbose && length(analysis_results$circos_notes) > 0 && !all(analysis_results$circos_notes %in% c("No notes from circos plot generation.", "No notes from circos plot generation (or circos plot did not run).")) ) {
    message("Notes from Circos plot generation (see 'circos_notes' in output list):")
    for(i in seq_along(analysis_results$circos_notes)){
      message(paste0("- ", analysis_results$circos_notes[i]))
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
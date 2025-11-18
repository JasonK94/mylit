#' Run CCI Analysis
#'
#' Main function for Cell-to-Cell Interaction analysis using NicheNet.
#' Takes a Seurat object, cluster information, and DEG list as inputs.
#'
#' @param sobj Seurat object
#' @param cluster_col Character string, name of metadata column with cluster annotations
#' @param deg_df Data frame with DEG results. Must contain columns: gene, cluster, and either logFC/avg_log2FC, p_val_adj
#' @param receiver_cluster Character string, receiver cluster ID
#' @param sender_clusters Character vector, sender cluster IDs (optional, NULL for auto-identification)
#' @param condition_col Character string, condition column for DE analysis (optional, if DEG already computed)
#' @param condition_oi Character string, condition of interest (optional)
#' @param condition_ref Character string, reference condition (optional)
#' @param species Character string, "human" or "mouse" (default: "human")
#' @param assay_name Character string, assay to use (default: "SCT")
#' @param min_pct_expressed Numeric, minimum fraction of cells expressing gene (default: 0.10)
#' @param p_val_adj_cutoff Numeric, adjusted p-value cutoff (default: 0.05)
#' @param logfc_cutoff Numeric, log fold-change cutoff (default: 0.25)
#' @param top_n_ligands Integer, number of top ligands (default: 20)
#' @param top_n_targets_per_ligand Integer, number of top targets per ligand (default: NULL, auto-adjusted: 50 for >3000 DEGs, 100 for >1000 DEGs, 200 otherwise)
#' @param ligand_target_cutoff Numeric, cutoff for ligand-target links (default: 0.33)
#' @param nichenet_data_dir Character string, NicheNet data directory (default: NULL)
#' @param nichenet_data_name Character string, NicheNet data name in global env (default: "NicheNetData")
#' @param output_dir Character string, output directory (default: NULL)
#' @param run_circos Logical, whether to generate circos plot (default: TRUE)
#' @param circos_title Character string, title for circos plot (default: NULL, auto-generated)
#' @param circos_show_legend Logical, whether to show legend in circos plot (default: TRUE)
#' @param circos_legend_position Character string, legend position (default: "topright")
#' @param circos_legend_size Numeric, legend text size (default: 0.9)
#' @param circos_legend_inset Numeric vector, legend inset (default: c(-0.15, 0))
#' @param auto_save Logical, whether to automatically save final results (default: TRUE)
#' @param save_prepared_data Logical, whether to save intermediate prepared data (default: FALSE)
#' @param verbose Logical, whether to print progress messages (default: TRUE)
#' @param ... Additional arguments passed to run_nichenet_analysis (circos params excluded)
#'
#' @return List with CCI analysis results
#' @export
run_cci_analysis <- function(sobj,
                              cluster_col,
                              deg_df,
                              receiver_cluster,
                              sender_clusters = NULL,
                              condition_col = NULL,
                              condition_oi = NULL,
                              condition_ref = NULL,
                              species = c("human", "mouse"),
                              assay_name = "SCT",
                              min_pct_expressed = 0.10,
                              p_val_adj_cutoff = 0.05,
                              logfc_cutoff = 0.25,
                              top_n_ligands = 20,
                              top_n_targets_per_ligand = NULL,  # Auto-adjusted based on DEG count if NULL
                              ligand_target_cutoff = 0.33,
                              nichenet_data_dir = NULL,
                              nichenet_data_name = "NicheNetData",
                              output_dir = NULL,
                              run_circos = TRUE,
                              circos_title = NULL,
                              circos_show_legend = TRUE,
                              circos_legend_position = "topright",
                              circos_legend_size = 0.9,
                              circos_legend_inset = c(-0.15, 0),
                              auto_save = TRUE,
                              save_prepared_data = FALSE,
                              verbose = TRUE,
                              ...) {
  
  species <- match.arg(species)
  
  if (verbose) {
    message("=== Starting CCI Analysis ===")
    message("Total steps: 7 (Validation → DEG extraction → Sender ID → Expressed genes → Summary → NicheNet → Compile)")
  }
  
  # Step 1: Validate inputs
  if (verbose) message("Step 1/7 (14%): Validating inputs...")
  validation <- validate_cci_inputs(sobj, cluster_col, deg_df, receiver_cluster, sender_clusters)
  
  # Step 2: Extract receiver DEGs
  if (verbose) message("Step 2/7 (29%): Extracting receiver DEGs...")
  receiver_degs <- extract_receiver_degs(
    deg_df, 
    receiver_cluster,
    p_val_adj_cutoff = p_val_adj_cutoff,
    logfc_cutoff = logfc_cutoff,
    only_upregulated = TRUE
  )
  
  if (nrow(receiver_degs) == 0) {
    stop("No DEGs found for receiver cluster '", receiver_cluster, "' after filtering")
  }
  
  if (verbose) {
    message("  Found ", nrow(receiver_degs), " DEGs for receiver cluster after filtering")
  }
  
  # Auto-adjust top_n_targets_per_ligand if not explicitly set and DEG count is large
  dots_list <- list(...)
  if (is.null(top_n_targets_per_ligand) && !"top_n_targets_per_ligand" %in% names(dots_list)) {
    n_degs <- nrow(receiver_degs)
    if (n_degs > 3000) {
      top_n_targets_per_ligand <- 50
      if (verbose) message("  Large DEG set detected (", n_degs, " genes). Auto-adjusting top_n_targets_per_ligand to ", top_n_targets_per_ligand, " for faster computation.")
    } else if (n_degs > 1000) {
      top_n_targets_per_ligand <- 100
      if (verbose) message("  Moderate DEG set (", n_degs, " genes). Auto-adjusting top_n_targets_per_ligand to ", top_n_targets_per_ligand, ".")
    } else {
      top_n_targets_per_ligand <- 200  # Keep default for small sets
    }
  }
  
  # Step 3: Identify sender clusters
  if (verbose) message("Step 3/7 (43%): Identifying sender clusters...")
  sender_clusters_final <- identify_sender_clusters(
    sobj, 
    cluster_col, 
    receiver_cluster, 
    sender_clusters
  )
  
  if (verbose) message("  Using ", length(sender_clusters_final), " sender cluster(s)")
  
  # Step 4: Get expressed genes
  if (verbose) message("Step 4/7 (57%): Identifying expressed genes...")
  
  # Set Idents for nichenetr::get_expressed_genes
  Seurat::DefaultAssay(sobj) <- assay_name
  Seurat::Idents(sobj) <- sobj@meta.data[[cluster_col]]
  
  expressed_genes_sender <- filter_expressed_genes(
    sobj,
    sender_clusters_final,
    min_pct_expressed = min_pct_expressed,
    assay_name = assay_name,
    cluster_col = cluster_col
  )
  
  expressed_genes_receiver <- filter_expressed_genes(
    sobj,
    receiver_cluster,
    min_pct_expressed = min_pct_expressed,
    assay_name = assay_name,
    cluster_col = cluster_col
  )
  
  if (verbose) message("  Sender: ", length(expressed_genes_sender), " expressed genes")
  if (verbose) message("  Receiver: ", length(expressed_genes_receiver), " expressed genes")
  
  # Step 5: Prepare data summary
  if (verbose) message("Step 5/7 (71%): Preparing data summary...")
  prepared_summary <- prepare_cci_summary(
    receiver_degs,
    sender_clusters_final,
    receiver_cluster,
    expressed_genes_sender,
    expressed_genes_receiver
  )
  
  # Save intermediate results if save_prepared_data is TRUE
  if (save_prepared_data) {
    prepared_data <- list(
      receiver_degs = receiver_degs,
      sender_clusters = sender_clusters_final,
      receiver_cluster = receiver_cluster,
      expressed_genes_sender = expressed_genes_sender,
      expressed_genes_receiver = expressed_genes_receiver,
      summary = prepared_summary
    )
    # Use output_dir if provided, otherwise use default cci directory
    # save_cci_intermediate will automatically detect run folders and save to run*/prepared_data/
    prep_output_dir <- if (!is.null(output_dir)) output_dir else "/data/user3/sobj/cci"
    save_cci_intermediate(prepared_data, output_dir = prep_output_dir)
  }
  
  # Step 6: Run NicheNet analysis
  # We need to use the existing run_nichenet_analysis function
  # But it requires DE analysis. We'll create a workaround:
  # Create a temporary condition column that will allow us to use FindMarkers
  # with our pre-computed DEGs, or we'll modify the approach
  
  if (verbose) {
    message("Step 6/7 (86%): Running NicheNet analysis...")
    message("  [", paste0(rep("=", 30), collapse = ""), "] Starting NicheNet pipeline...")
    message("  This step includes:")
    message("    - NicheNet data loading")
    message("    - Expressed genes identification")
    message("    - Ligand activity prediction")
    message("    - Ligand-target inference")
    message("    - Visualization generation")
    if (run_circos) {
      message("    - Circos plot generation")
    }
  }
  
  # The challenge: run_nichenet_analysis does DE analysis internally
  # We have pre-computed DEGs, so we need to work around this
  
  # Option 1: Use run_nichenet_analysis with a dummy condition
  # Option 2: Call NicheNet functions directly (more complex but cleaner)
  
  # For now, we'll use Option 1 with a workaround:
  # Create a temporary condition column that will make FindMarkers return our genes
  
  # Actually, a better approach: we'll modify the Seurat object temporarily
  # to have a condition that matches our DEGs, then use run_nichenet_analysis
  
  # But this is complex. Let's use a simpler approach:
  # We'll call run_nichenet_analysis but we need to ensure it uses our DEGs
  
  # For now, let's use the existing function but we'll need to handle the DE analysis part
  # We'll create a modified version that accepts pre-computed DEGs
  
  # Since modifying run_nichenet_analysis is complex, we'll use it as-is
  # but we'll need to provide condition_col, condition_oi, condition_ref
  # and ensure FindMarkers returns genes that match our DEG list
  
  # Actually, the simplest approach: if condition_col is provided, use run_nichenet_analysis
  # Otherwise, we'll need to use a different strategy
  
  # Check if run_nichenet_analysis function is available
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
        if (isTRUE(verbose)) message("  Loaded run_nichenet_analysis from ", candidate)
        break
      }
    }
    if (!sourced) {
      stop(
        "Cannot find run_nichenet_analysis function. Please ensure CCI.R is sourced. ",
        "Tried: ", paste(source_candidates, collapse = ", ")
      )
    }
  }
  
  if (!is.null(condition_col) && !is.null(condition_oi)) {
    # Use run_nichenet_analysis with conditions
    # Note: run_nichenet_analysis will perform its own DE analysis using FindMarkers
    # The deg_df provided is used for validation/preparation, but the actual DE analysis
    # is done by run_nichenet_analysis internally
    
    if (verbose) {
      message("  Preparing NicheNet run with ", length(sender_clusters_final), " sender cluster(s) and ", nrow(receiver_degs), " receiver DEGs.")
      message("  Precomputed receiver DE tables will be reused when available to avoid redundant FindMarkers calls.")
      # More realistic estimate: DEG count significantly affects ligand activity prediction time
      # Base time for small DEG sets, then scale
      base_time_minutes <- 1
      if (nrow(receiver_degs) > 3000) {
        est_minutes <- max(10, round(nrow(receiver_degs) / 200, 1))  # Conservative for large sets
        message("  Estimated runtime: ~", est_minutes, " min (conservative estimate for large DEG set).")
        message("  [Note: Ligand activity prediction scales with DEG count - may take 10-20+ minutes]")
      } else if (nrow(receiver_degs) > 1000) {
        est_minutes <- max(3, round(nrow(receiver_degs) / 300, 1))
        message("  Estimated runtime: ~", est_minutes, " min (heuristic).")
      } else {
        est_minutes <- base_time_minutes
        message("  Estimated runtime: ~", est_minutes, " min (heuristic).")
      }
    }
    
    # Extract circos-related parameters from function arguments and ...
    # Priority: function arguments > ... arguments
    dots_list <- list(...)
    circos_params <- list(
      circos_title = if(!is.null(circos_title)) circos_title else dots_list$circos_title,
      circos_show_legend = if(!is.null(circos_show_legend)) circos_show_legend else dots_list$circos_show_legend,
      circos_legend_position = if(!is.null(circos_legend_position)) circos_legend_position else dots_list$circos_legend_position,
      circos_legend_size = if(!is.null(circos_legend_size)) circos_legend_size else dots_list$circos_legend_size,
      circos_legend_inset = if(!is.null(circos_legend_inset)) circos_legend_inset else dots_list$circos_legend_inset
    )
    
    # Remove circos params from dots_list to avoid duplication
    dots_list_clean <- dots_list
    for (param in c("circos_title", "circos_show_legend", "circos_legend_position", 
                    "circos_legend_size", "circos_legend_inset")) {
      dots_list_clean[[param]] <- NULL
    }
    
    # Generate circos title if not provided
    if (run_circos && is.null(circos_params$circos_title)) {
      # Create informative title with condition and cell types
      condition_str <- paste0(condition_col, ": ", condition_oi)
      if (!is.null(condition_ref)) {
        condition_str <- paste0(condition_str, " vs ", condition_ref)
      } else {
        condition_str <- paste0(condition_str, " vs others")
      }
      sender_str <- if(length(sender_clusters_final) <= 3) {
        paste(sender_clusters_final, collapse = ", ")
      } else {
        paste0(paste(head(sender_clusters_final, 2), collapse = ", "), " (+", length(sender_clusters_final) - 2, " more)")
      }
      circos_title_auto <- paste0("Ligand-Receptor Interactions\n", 
                                  "Receiver: ", receiver_cluster, " | ",
                                  "Senders: ", sender_str, "\n",
                                  "Condition: ", condition_str)
    } else {
      circos_title_auto <- circos_params$circos_title
    }
    
    # Call run_nichenet_analysis
    # Note: We pass more lenient parameters to FindMarkers to ensure DEGs are found
    # The actual filtering will be done by p_val_adj_cutoff and logfc_cutoff
    nichenet_results <- do.call(run_nichenet_analysis, c(
      list(
        seurat_obj = sobj,
        species = species,
        sender_celltypes = sender_clusters_final,
        receiver_celltype = receiver_cluster,
        assay_name = assay_name,
        cluster_col = cluster_col,
        receiver_DE_ident1 = condition_oi,
        receiver_DE_ident2 = condition_ref,
        receiver_DE_group_by = condition_col,
        receiver_de_table = receiver_degs,
        min_pct_expressed = max(0.05, min_pct_expressed),  # More lenient min_pct
        p_val_adj_cutoff = p_val_adj_cutoff,
        logfc_cutoff = max(0.1, logfc_cutoff),  # More lenient logfc threshold for FindMarkers
        top_n_ligands = top_n_ligands,
        top_n_targets_per_ligand = top_n_targets_per_ligand,
        ligand_target_cutoff = ligand_target_cutoff,
        nichenet_data_dir = nichenet_data_dir,
        nichenet_data_name = nichenet_data_name,
        output_dir = output_dir,
        run_circos = run_circos,
        verbose = verbose
      ),
      # Add circos parameters (with defaults)
      list(
        circos_title = circos_title_auto,
        circos_show_legend = if(is.null(circos_params$circos_show_legend)) TRUE else circos_params$circos_show_legend,
        circos_legend_position = if(is.null(circos_params$circos_legend_position)) "topright" else circos_params$circos_legend_position,
        circos_legend_size = circos_params$circos_legend_size,
        circos_legend_inset = circos_params$circos_legend_inset
      ),
      # Add remaining ... parameters (excluding circos params)
      dots_list_clean
    ))
  } else {
    # No condition provided - we need to use pre-computed DEGs
    # This requires a different approach since run_nichenet_analysis does DE analysis
    
    stop("When using pre-computed DEG list, you must provide condition_col, condition_oi (and optionally condition_ref) to use run_nichenet_analysis. ",
         "Alternatively, the DEG list should match what FindMarkers would return for the specified conditions.")
    
    # TODO: Implement direct NicheNet analysis with pre-computed DEGs
    # This would require calling NicheNet functions directly instead of using run_nichenet_analysis
  }
  
  # Step 7: Create sender-receiver map
  sender_receiver_map <- create_sender_receiver_map(sender_clusters_final, receiver_cluster)
  
  # Step 8: Format DEG summary
  deg_summary <- format_deg_summary(deg_df, receiver_cluster)
  
  # Step 9: Compile results
  if (verbose) message("Step 7/7 (100%): Compiling results...")
  results <- list(
    nichenet_results = nichenet_results,
    sender_receiver_map = sender_receiver_map,
    deg_summary = deg_summary,
    receiver_degs = receiver_degs,
    sender_clusters = sender_clusters_final,
    receiver_cluster = receiver_cluster,
    prepared_summary = prepared_summary,
    output_path = if (!is.null(output_dir)) output_dir else NULL
  )
  
  # Step 10: Save final results if auto_save is TRUE
  # Note: run_nichenet_analysis already saves results to output_dir/run*/nichenet_results.qs
  # This saves the full CCI results (including wrapper info) to the same location
  if (auto_save) {
    saved_path <- save_cci_final(results, output_dir = output_dir)
    if (!is.null(saved_path)) {
      results$saved_path <- saved_path
      if (verbose) message("Final CCI results saved to: ", saved_path)
    }
  }
  
  if (verbose) message("=== CCI Analysis Complete ===")
  
  return(results)
}


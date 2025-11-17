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
#' @param top_n_targets_per_ligand Integer, number of top targets per ligand (default: 200)
#' @param ligand_target_cutoff Numeric, cutoff for ligand-target links (default: 0.33)
#' @param nichenet_data_dir Character string, NicheNet data directory (default: NULL)
#' @param nichenet_data_name Character string, NicheNet data name in global env (default: "NicheNetData")
#' @param output_dir Character string, output directory (default: NULL)
#' @param run_circos Logical, whether to generate circos plot (default: TRUE)
#' @param auto_save Logical, whether to automatically save results (default: TRUE)
#' @param verbose Logical, whether to print progress messages (default: TRUE)
#' @param ... Additional arguments passed to run_nichenet_analysis
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
                              top_n_targets_per_ligand = 200,
                              ligand_target_cutoff = 0.33,
                              nichenet_data_dir = NULL,
                              nichenet_data_name = "NicheNetData",
                              output_dir = NULL,
                              run_circos = TRUE,
                              auto_save = TRUE,
                              verbose = TRUE,
                              ...) {
  
  species <- match.arg(species)
  
  if (verbose) message("=== Starting CCI Analysis ===")
  
  # Step 1: Validate inputs
  if (verbose) message("Step 1: Validating inputs...")
  validation <- validate_cci_inputs(sobj, cluster_col, deg_df, receiver_cluster, sender_clusters)
  
  # Step 2: Extract receiver DEGs
  if (verbose) message("Step 2: Extracting receiver DEGs...")
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
  
  if (verbose) message("  Found ", nrow(receiver_degs), " DEGs for receiver cluster")
  
  # Step 3: Identify sender clusters
  if (verbose) message("Step 3: Identifying sender clusters...")
  sender_clusters_final <- identify_sender_clusters(
    sobj, 
    cluster_col, 
    receiver_cluster, 
    sender_clusters
  )
  
  if (verbose) message("  Using ", length(sender_clusters_final), " sender cluster(s)")
  
  # Step 4: Get expressed genes
  if (verbose) message("Step 4: Identifying expressed genes...")
  
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
  if (verbose) message("Step 5: Preparing data summary...")
  prepared_summary <- prepare_cci_summary(
    receiver_degs,
    sender_clusters_final,
    receiver_cluster,
    expressed_genes_sender,
    expressed_genes_receiver
  )
  
  # Save intermediate results if auto_save is TRUE
  if (auto_save) {
    prepared_data <- list(
      receiver_degs = receiver_degs,
      sender_clusters = sender_clusters_final,
      receiver_cluster = receiver_cluster,
      expressed_genes_sender = expressed_genes_sender,
      expressed_genes_receiver = expressed_genes_receiver,
      summary = prepared_summary
    )
    save_cci_intermediate(prepared_data)
  }
  
  # Step 6: Run NicheNet analysis
  # We need to use the existing run_nichenet_analysis function
  # But it requires DE analysis. We'll create a workaround:
  # Create a temporary condition column that will allow us to use FindMarkers
  # with our pre-computed DEGs, or we'll modify the approach
  
  if (verbose) message("Step 6: Running NicheNet analysis...")
  
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
    # Try to source it
    source_file <- "/home/user3/data_user3/git_repo/mylit/myR/R/CCI.R"
    if (file.exists(source_file)) {
      source(source_file)
    } else {
      # Try alternative path
      source_file2 <- "/home/user3/data_user3/git_repo/_wt/cci/myR/R/CCI.R"
      if (file.exists(source_file2)) {
        source(source_file2)
      } else {
        stop("Cannot find run_nichenet_analysis function. Please ensure CCI.R is sourced. ",
             "Expected locations: ", source_file, " or ", source_file2)
      }
    }
  }
  
  if (!is.null(condition_col) && !is.null(condition_oi)) {
    # Use run_nichenet_analysis with conditions
    # Note: run_nichenet_analysis will perform its own DE analysis using FindMarkers
    # The deg_df provided is used for validation/preparation, but the actual DE analysis
    # is done by run_nichenet_analysis internally
    
    if (verbose) {
      message("  Note: run_nichenet_analysis will perform DE analysis using FindMarkers.")
      message("  The provided deg_df is used for validation, but internal DE analysis will be used.")
    }
    
    # Generate circos title if not provided
    if (run_circos && is.null(list(...)$circos_title)) {
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
      circos_title_auto <- list(...)$circos_title
    }
    
    # Call run_nichenet_analysis
    # Note: We pass more lenient parameters to FindMarkers to ensure DEGs are found
    # The actual filtering will be done by p_val_adj_cutoff and logfc_cutoff
    nichenet_results <- run_nichenet_analysis(
      seurat_obj = sobj,
      species = species,
      sender_celltypes = sender_clusters_final,
      receiver_celltype = receiver_cluster,
      assay_name = assay_name,
      cluster_col = cluster_col,
      receiver_DE_ident1 = condition_oi,
      receiver_DE_ident2 = condition_ref,
      receiver_DE_group_by = condition_col,
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
      circos_title = circos_title_auto,  # Pass auto-generated or user-provided title
      circos_show_legend = if(is.null(list(...)$circos_show_legend)) TRUE else list(...)$circos_show_legend,
      circos_legend_position = if(is.null(list(...)$circos_legend_position)) "topright" else list(...)$circos_legend_position,
      verbose = verbose,
      ...
    )
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
  if (verbose) message("Step 7: Compiling results...")
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
  
  # Step 10: Save results if auto_save is TRUE
  if (auto_save) {
    saved_path <- save_cci_final(results)
    results$saved_path <- saved_path
  }
  
  if (verbose) message("=== CCI Analysis Complete ===")
  
  return(results)
}


#' NicheNet Cell-Cell Communication Analysis
#'
#' This module provides functions for performing NicheNet-based cell-cell
#' communication analysis, including ligand activity prediction, receptor
#' inference, and visualization.
#'
#' @name nichenet_analysis
NULL

#' Get Ligand-Target Links from NicheNet
#'
#' Retrieves ligand-target gene regulatory potential from NicheNet databases.
#'
#' @param ligand Character vector of ligand genes
#' @param target Character vector of target genes
#' @param lr_network Ligand-receptor network (from NicheNet)
#' @param sig_network Signaling network (from NicheNet)
#' @param gr_network Gene regulatory network (from NicheNet)
#' @param top_n_targets Number of top targets to return per ligand (default: 250)
#'
#' @return Data frame with ligand-target regulatory potential scores
#'
#' @export
ligand_to_target <- function(ligand,
                             target,
                             lr_network,
                             sig_network,
                             gr_network,
                             top_n_targets = 250) {
  
  # Get weighted networks
  weighted_networks <- nichenetr::construct_weighted_networks(
    lr_network = lr_network,
    sig_network = sig_network,
    gr_network = gr_network,
    source_weights_df = NULL
  )
  
  # Get ligand-target matrix
  ligand_target_matrix <- nichenetr::construct_ligand_target_matrix(
    weighted_networks = weighted_networks,
    ligands = ligand,
    targets = target,
    ltf_cutoff = 0.99,
    algorithm = "PPR",
    damping_factor = 0.5,
    secondary_targets = FALSE
  )
  
  # Convert to data frame
  result <- as.data.frame(as.table(ligand_target_matrix))
  colnames(result) <- c("ligand", "target", "regulatory_potential")
  
  # Filter and sort
  result <- result[result$regulatory_potential > 0, ]
  result <- result[order(result$ligand, -result$regulatory_potential), ]
  
  # Keep top N targets per ligand
  result <- do.call(rbind, lapply(split(result, result$ligand), function(x) {
    head(x, top_n_targets)
  }))
  
  return(result)
}

#' Run Complete NicheNet Analysis
#'
#' Performs a comprehensive NicheNet cell-cell communication analysis workflow,
#' including ligand activity prediction, receptor inference, and visualization.
#'
#' @param seurat_obj Seurat object
#' @param sender_cells Vector of sender cell identities or logical vector
#' @param receiver_cells Vector of receiver cell identities or logical vector
#' @param condition_oi Condition of interest for receiver DE analysis
#' @param condition_ref Reference condition for receiver DE analysis
#' @param ident_col Identity column in metadata (default: "seurat_clusters")
#' @param condition_col Condition column in metadata (required for DE)
#' @param lr_network Ligand-receptor network (from NicheNet)
#' @param sig_network Signaling network (from NicheNet)
#' @param gr_network Gene regulatory network (from NicheNet)
#' @param ligand_target_matrix Pre-computed ligand-target matrix (optional)
#' @param expressed_pct Expression percentage threshold (default: 0.10)
#' @param top_n_ligands Number of top ligands to analyze (default: 20)
#' @param top_n_targets Number of top targets per ligand (default: 200)
#' @param plot_circos Generate Circos plot (default: FALSE)
#' @param output_dir Directory for saving plots (default: NULL)
#'
#' @return List containing:
#'   \item{ligand_activities}{Predicted ligand activities}
#'   \item{best_ligands}{Top prioritized ligands}
#'   \item{ligand_target_links}{Ligand-target regulatory links}
#'   \item{ligand_receptor_links}{Ligand-receptor interactions}
#'   \item{plots}{Generated plots}
#'   \item{warnings}{Any warnings encountered}
#'
#' @examples
#' \dontrun{
#' # Load NicheNet networks
#' lr_network <- readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
#' sig_network <- readRDS(url("https://zenodo.org/record/3260758/files/signaling_network.rds"))
#' gr_network <- readRDS(url("https://zenodo.org/record/3260758/files/gr_network.rds"))
#' 
#' # Run analysis
#' results <- run_nichenet_analysis(
#'   seurat_obj,
#'   sender_cells = c("0", "1", "2"),
#'   receiver_cells = "3",
#'   condition_oi = "Stimulated",
#'   condition_ref = "Control",
#'   condition_col = "condition",
#'   lr_network = lr_network,
#'   sig_network = sig_network,
#'   gr_network = gr_network
#' )
#' }
#'
#' @export
run_nichenet_analysis <- function(seurat_obj,
                                  sender_cells,
                                  receiver_cells,
                                  condition_oi = NULL,
                                  condition_ref = NULL,
                                  ident_col = "seurat_clusters",
                                  condition_col = NULL,
                                  lr_network,
                                  sig_network,
                                  gr_network,
                                  ligand_target_matrix = NULL,
                                  expressed_pct = 0.10,
                                  top_n_ligands = 20,
                                  top_n_targets = 200,
                                  plot_circos = FALSE,
                                  output_dir = NULL) {
  
  warning_messages <- character(0)
  
  # Validate inputs
  if (!ident_col %in% colnames(seurat_obj@meta.data)) {
    stop("ident_col '", ident_col, "' not found in metadata")
  }
  
  if (!is.null(condition_col) && !condition_col %in% colnames(seurat_obj@meta.data)) {
    stop("condition_col '", condition_col, "' not found in metadata")
  }
  
  # Define expressed genes
  message("1. Identifying expressed genes...")
  expressed_genes <- rownames(seurat_obj)[
    Matrix::rowMeans(Seurat::GetAssayData(seurat_obj, slot = "data") > 0) > expressed_pct
  ]
  
  message("   ", length(expressed_genes), " genes expressed in > ", 
          expressed_pct * 100, "% of cells")
  
  # Get expressed genes in sender and receiver cells
  sender_logical <- if (is.logical(sender_cells)) {
    sender_cells
  } else {
    seurat_obj@meta.data[[ident_col]] %in% sender_cells
  }
  
  receiver_logical <- if (is.logical(receiver_cells)) {
    receiver_cells
  } else {
    seurat_obj@meta.data[[ident_col]] %in% receiver_cells
  }
  
  expressed_genes_sender <- unique(rownames(seurat_obj)[
    Matrix::rowMeans(Seurat::GetAssayData(seurat_obj[, sender_logical], slot = "data") > 0) > expressed_pct
  ])
  
  expressed_genes_receiver <- unique(rownames(seurat_obj)[
    Matrix::rowMeans(Seurat::GetAssayData(seurat_obj[, receiver_logical], slot = "data") > 0) > expressed_pct
  ])
  
  # Get DE genes in receiver cells (if conditions specified)
  if (!is.null(condition_oi) && !is.null(condition_ref) && !is.null(condition_col)) {
    message("2. Finding DE genes in receiver cells...")
    
    Seurat::Idents(seurat_obj) <- seurat_obj@meta.data[[condition_col]]
    
    receiver_de <- tryCatch({
      Seurat::FindMarkers(
        seurat_obj[, receiver_logical],
        ident.1 = condition_oi,
        ident.2 = condition_ref,
        min.pct = 0.10,
        logfc.threshold = 0.25
      )
    }, error = function(e) {
      warning_msg <- paste("Error in FindMarkers:", e$message)
      warning_messages <<- c(warning_messages, warning_msg)
      return(NULL)
    })
    
    if (is.null(receiver_de) || nrow(receiver_de) == 0) {
      stop("No DE genes found in receiver cells")
    }
    
    genes_oi <- rownames(receiver_de)[receiver_de$p_val_adj < 0.05 & abs(receiver_de$avg_log2FC) > 0.25]
    message("   ", length(genes_oi), " DE genes (padj < 0.05, |logFC| > 0.25)")
    
  } else {
    # Use all expressed genes in receiver
    message("2. Using all expressed genes as targets (no DE analysis)...")
    genes_oi <- expressed_genes_receiver
  }
  
  # Define background genes
  background_genes <- expressed_genes
  
  # Get potential ligands
  message("3. Identifying potential ligands...")
  ligands <- lr_network %>% dplyr::pull(from) %>% unique()
  expressed_ligands <- intersect(ligands, expressed_genes_sender)
  
  receptors <- lr_network %>% dplyr::pull(to) %>% unique()
  expressed_receptors <- intersect(receptors, expressed_genes_receiver)
  
  potential_ligands <- lr_network %>%
    dplyr::filter(from %in% expressed_ligands & to %in% expressed_receptors) %>%
    dplyr::pull(from) %>%
    unique()
  
  message("   ", length(potential_ligands), " potential ligands")
  
  if (length(potential_ligands) == 0) {
    stop("No potential ligands found")
  }
  
  # Construct or use ligand-target matrix
  if (is.null(ligand_target_matrix)) {
    message("4. Constructing ligand-target matrix...")
    weighted_networks <- nichenetr::construct_weighted_networks(
      lr_network = lr_network,
      sig_network = sig_network,
      gr_network = gr_network
    )
    
    ligand_target_matrix <- nichenetr::construct_ligand_target_matrix(
      weighted_networks = weighted_networks,
      ligands = potential_ligands,
      targets = background_genes
    )
  } else {
    message("4. Using provided ligand-target matrix...")
  }
  
  # Predict ligand activities
  message("5. Predicting ligand activities...")
  ligand_activities <- nichenetr::predict_ligand_activities(
    geneset = genes_oi,
    background_expressed_genes = background_genes,
    ligand_target_matrix = ligand_target_matrix,
    potential_ligands = potential_ligands
  )
  
  ligand_activities <- ligand_activities %>%
    dplyr::arrange(dplyr::desc(pearson))
  
  best_ligands <- ligand_activities %>%
    dplyr::top_n(top_n_ligands, pearson) %>%
    dplyr::pull(test_ligand)
  
  message("   Top ligand: ", best_ligands[1], 
          " (Pearson = ", round(ligand_activities$pearson[1], 3), ")")
  
  # Get ligand-target links
  message("6. Extracting ligand-target links...")
  active_ligand_target_links <- ligand_target_matrix[best_ligands, genes_oi] %>%
    as.matrix() %>%
    reshape2::melt() %>%
    dplyr::rename(ligand = Var1, target = Var2, regulatory_potential = value) %>%
    dplyr::filter(regulatory_potential > 0) %>%
    dplyr::arrange(dplyr::desc(regulatory_potential))
  
  # Limit targets per ligand
  active_ligand_target_links <- active_ligand_target_links %>%
    dplyr::group_by(ligand) %>%
    dplyr::top_n(top_n_targets, regulatory_potential) %>%
    dplyr::ungroup()
  
  # Get ligand-receptor links
  message("7. Identifying ligand-receptor pairs...")
  ligand_receptor_links <- lr_network %>%
    dplyr::filter(from %in% best_ligands & to %in% expressed_receptors) %>%
    dplyr::rename(ligand = from, receptor = to)
  
  # Create visualizations
  message("8. Creating visualizations...")
  plots <- list()
  
  # Ligand activity heatmap
  tryCatch({
    ligand_pearson_matrix <- ligand_activities %>%
      dplyr::top_n(top_n_ligands, pearson) %>%
      dplyr::select(test_ligand, pearson) %>%
      tibble::column_to_rownames("test_ligand") %>%
      as.matrix()
    
    plots$ligand_activity <- pheatmap::pheatmap(
      ligand_pearson_matrix,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      main = "Ligand Activity",
      color = viridisLite::viridis(100),
      border_color = NA
    )
  }, error = function(e) {
    warning_msg <- paste("Failed to create ligand activity heatmap:", e$message)
    warning_messages <<- c(warning_messages, warning_msg)
  })
  
  # Circos plot (optional)
  if (plot_circos && requireNamespace("circlize", quietly = TRUE)) {
    message("9. Creating Circos plot...")
    tryCatch({
      circos_data <- prepare_nichenet_circos_data(
        ligand_receptor_links,
        active_ligand_target_links,
        ligand_activities,
        top_n_ligands = min(10, top_n_ligands),
        top_n_targets = min(20, top_n_targets)
      )
      
      if (!is.null(output_dir)) {
        circos_file <- file.path(output_dir, "nichenet_circos.pdf")
        pdf(circos_file, width = 10, height = 10)
      }
      
      draw_nichenet_circos_plot(circos_data)
      
      if (!is.null(output_dir)) {
        dev.off()
        message("   Circos plot saved: ", circos_file)
      }
      
      plots$circos <- "Circos plot generated"
      
    }, error = function(e) {
      warning_msg <- paste("Failed to create Circos plot:", e$message)
      warning_messages <<- c(warning_messages, warning_msg)
    })
  }
  
  message("\n=== NicheNet Analysis Complete ===")
  
  # Return results
  return(list(
    ligand_activities = ligand_activities,
    best_ligands = best_ligands,
    ligand_target_links = active_ligand_target_links,
    ligand_receptor_links = ligand_receptor_links,
    plots = plots,
    warnings = warning_messages
  ))
}

#' Prepare Data for NicheNet Circos Plot
#'
#' Prepares ligand-receptor and ligand-target data for Circos visualization.
#'
#' @param ligand_receptor_links Data frame with ligand-receptor pairs
#' @param ligand_target_links Data frame with ligand-target regulatory links
#' @param ligand_activities Data frame with ligand activity scores
#' @param top_n_ligands Number of top ligands to include (default: 10)
#' @param top_n_targets Number of top targets per ligand (default: 20)
#'
#' @return List with formatted data for Circos plot
#'
#' @export
prepare_nichenet_circos_data <- function(ligand_receptor_links,
                                         ligand_target_links,
                                         ligand_activities,
                                         top_n_ligands = 10,
                                         top_n_targets = 20) {
  
  # Select top ligands
  top_ligands <- ligand_activities %>%
    dplyr::top_n(top_n_ligands, pearson) %>%
    dplyr::pull(test_ligand)
  
  # Filter links
  lr_links_subset <- ligand_receptor_links %>%
    dplyr::filter(ligand %in% top_ligands)
  
  lt_links_subset <- ligand_target_links %>%
    dplyr::filter(ligand %in% top_ligands) %>%
    dplyr::group_by(ligand) %>%
    dplyr::top_n(top_n_targets, regulatory_potential) %>%
    dplyr::ungroup()
  
  # Get unique nodes
  receptors <- unique(lr_links_subset$receptor)
  targets <- unique(lt_links_subset$target)
  
  return(list(
    ligands = top_ligands,
    receptors = receptors,
    targets = targets,
    lr_links = lr_links_subset,
    lt_links = lt_links_subset,
    ligand_activities = ligand_activities %>% dplyr::filter(test_ligand %in% top_ligands)
  ))
}

#' Draw NicheNet Circos Plot
#'
#' Creates a Circos plot showing ligand-receptor and ligand-target interactions.
#'
#' @param circos_data List from prepare_nichenet_circos_data
#' @param ligand_color Color for ligands (default: "#E41A1C")
#' @param receptor_color Color for receptors (default: "#377EB8")
#' @param target_color Color for targets (default: "#4DAF4A")
#'
#' @return NULL (plot is drawn to current device)
#'
#' @export
draw_nichenet_circos_plot <- function(circos_data,
                                      ligand_color = "#E41A1C",
                                      receptor_color = "#377EB8",
                                      target_color = "#4DAF4A") {
  
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("circlize package is required for Circos plots")
  }
  
  # Combine all nodes
  all_nodes <- c(circos_data$ligands, circos_data$receptors, circos_data$targets)
  
  # Assign colors
  node_colors <- c(
    rep(ligand_color, length(circos_data$ligands)),
    rep(receptor_color, length(circos_data$receptors)),
    rep(target_color, length(circos_data$targets))
  )
  names(node_colors) <- all_nodes
  
  # Combine links
  lr_links_plot <- circos_data$lr_links %>%
    dplyr::select(ligand, receptor) %>%
    dplyr::rename(from = ligand, to = receptor)
  
  lt_links_plot <- circos_data$lt_links %>%
    dplyr::select(ligand, target) %>%
    dplyr::rename(from = ligand, to = target)
  
  all_links <- rbind(lr_links_plot, lt_links_plot)
  
  # Initialize Circos
  circlize::circos.clear()
  circlize::circos.par(gap.after = c(
    rep(1, length(circos_data$ligands) - 1), 10,
    rep(1, length(circos_data$receptors) - 1), 10,
    rep(1, length(circos_data$targets) - 1), 10
  ))
  
  # Create plot
  circlize::chordDiagram(
    all_links,
    grid.col = node_colors,
    transparency = 0.5,
    directional = 1,
    direction.type = "arrows",
    link.arr.type = "big.arrow"
  )
  
  # Add title
  title("NicheNet: Ligand-Receptor-Target Network")
  
  circlize::circos.clear()
}


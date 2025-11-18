#' Prepare Data for CCI Analysis
#'
#' This module provides functions for preparing and validating data for
#' Cell-to-Cell Interaction (CCI) analysis.
#'
#' @name prepare_cci_data
NULL

#' Validate CCI Inputs
#'
#' Validates inputs for CCI analysis including Seurat object, cluster column,
#' DEG dataframe, and receiver cluster.
#'
#' @param sobj Seurat object
#' @param cluster_col Character string, name of the metadata column containing cluster annotations
#' @param deg_df Data frame with DEG results. Must contain columns: gene, cluster, and either logFC/avg_log2FC, p_val_adj
#' @param receiver_cluster Character string, the receiver cluster ID
#' @param sender_clusters Character vector, sender cluster IDs (optional, can be NULL)
#'
#' @return List with validation results and processed data
#' @export
validate_cci_inputs <- function(sobj, cluster_col, deg_df, receiver_cluster, sender_clusters = NULL) {
  
  # Validate Seurat object
  if (!inherits(sobj, "Seurat")) {
    stop("`sobj` must be a Seurat object")
  }
  
  # Validate cluster column
  if (!cluster_col %in% colnames(sobj@meta.data)) {
    stop("`cluster_col` '", cluster_col, "' not found in Seurat object metadata")
  }
  
  # Validate DEG dataframe
  if (!is.data.frame(deg_df) || nrow(deg_df) == 0) {
    stop("`deg_df` must be a non-empty data frame")
  }
  
  # Check required columns in DEG dataframe
  required_cols <- c("gene", "cluster")
  missing_cols <- base::setdiff(required_cols, colnames(deg_df))
  if (length(missing_cols) > 0) {
    stop("`deg_df` missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Check for logFC column (either logFC or avg_log2FC)
  logfc_col <- NULL
  if ("logFC" %in% colnames(deg_df)) {
    logfc_col <- "logFC"
  } else if ("avg_log2FC" %in% colnames(deg_df)) {
    logfc_col <- "avg_log2FC"
  } else {
    stop("`deg_df` must contain either 'logFC' or 'avg_log2FC' column")
  }
  
  # Check for p-value column
  pval_col <- NULL
  if ("p_val_adj" %in% colnames(deg_df)) {
    pval_col <- "p_val_adj"
  } else if ("FDR" %in% colnames(deg_df)) {
    pval_col <- "FDR"
  } else if ("p_val" %in% colnames(deg_df)) {
    pval_col <- "p_val"
  } else {
    warning("`deg_df` does not contain p-value column. Will proceed without filtering.")
  }
  
  # Validate receiver cluster exists in DEG dataframe
  if (!receiver_cluster %in% deg_df$cluster) {
    stop("`receiver_cluster` '", receiver_cluster, "' not found in `deg_df$cluster`")
  }
  
  # Validate receiver cluster exists in Seurat object
  if (!receiver_cluster %in% sobj@meta.data[[cluster_col]]) {
    stop("`receiver_cluster` '", receiver_cluster, "' not found in Seurat object metadata column '", cluster_col, "'")
  }
  
  # Validate sender clusters if provided
  if (!is.null(sender_clusters)) {
    if (!all(sender_clusters %in% sobj@meta.data[[cluster_col]])) {
      missing_senders <- base::setdiff(sender_clusters, sobj@meta.data[[cluster_col]])
      stop("Some sender clusters not found in Seurat object: ", paste(missing_senders, collapse = ", "))
    }
  }
  
  return(list(
    valid = TRUE,
    logfc_col = logfc_col,
    pval_col = pval_col
  ))
}

#' Extract Receiver DEGs
#'
#' Extracts DEGs for the receiver cluster from a DEG dataframe.
#'
#' @param deg_df Data frame with DEG results
#' @param receiver_cluster Character string, the receiver cluster ID
#' @param p_val_adj_cutoff Numeric, adjusted p-value cutoff (default: 0.05)
#' @param logfc_cutoff Numeric, log fold-change cutoff (default: 0.25)
#' @param only_upregulated Logical, whether to only include upregulated genes (default: TRUE)
#'
#' @return Data frame with filtered DEGs for the receiver cluster
#' @export
extract_receiver_degs <- function(deg_df, receiver_cluster, 
                                   p_val_adj_cutoff = 0.05, 
                                   logfc_cutoff = 0.25,
                                   only_upregulated = TRUE) {
  
  # Filter for receiver cluster
  receiver_degs <- deg_df %>%
    dplyr::filter(cluster == receiver_cluster)
  
  if (nrow(receiver_degs) == 0) {
    stop("No DEGs found for receiver cluster '", receiver_cluster, "'")
  }
  
  # Determine logFC column
  logfc_col <- if ("avg_log2FC" %in% colnames(receiver_degs)) "avg_log2FC" else "logFC"
  
  # Filter by logFC
  if (only_upregulated) {
    receiver_degs <- receiver_degs %>%
      dplyr::filter(!!dplyr::sym(logfc_col) > logfc_cutoff)
  } else {
    receiver_degs <- receiver_degs %>%
      dplyr::filter(abs(!!dplyr::sym(logfc_col)) > logfc_cutoff)
  }
  
  # Filter by p-value if available
  pval_col <- NULL
  if ("p_val_adj" %in% colnames(receiver_degs)) {
    pval_col <- "p_val_adj"
  } else if ("FDR" %in% colnames(receiver_degs)) {
    pval_col <- "FDR"
  } else if ("p_val" %in% colnames(receiver_degs)) {
    pval_col <- "p_val"
  }
  
  if (!is.null(pval_col)) {
    receiver_degs <- receiver_degs %>%
      dplyr::filter(!!dplyr::sym(pval_col) < p_val_adj_cutoff)
  }
  
  if (nrow(receiver_degs) == 0) {
    warning("No DEGs passed filtering criteria for receiver cluster '", receiver_cluster, "'")
  }
  
  return(receiver_degs)
}

#' Identify Sender Clusters
#'
#' Identifies potential sender clusters based on expressed ligands.
#' If sender_clusters is provided, validates and returns them.
#' Otherwise, identifies all clusters except the receiver as potential senders.
#'
#' @param sobj Seurat object
#' @param cluster_col Character string, name of the metadata column containing cluster annotations
#' @param receiver_cluster Character string, the receiver cluster ID
#' @param sender_clusters Character vector, sender cluster IDs (optional, can be NULL for auto-identification)
#'
#' @return Character vector of sender cluster IDs
#' @export
identify_sender_clusters <- function(sobj, cluster_col, receiver_cluster, sender_clusters = NULL) {
  
  # Get all clusters
  all_clusters <- unique(sobj@meta.data[[cluster_col]])
  all_clusters <- all_clusters[!is.na(all_clusters)]
  
  # Remove receiver from potential senders
  potential_senders <- base::setdiff(all_clusters, receiver_cluster)
  
  if (is.null(sender_clusters)) {
    # Auto-identify: use all clusters except receiver
    message("Auto-identifying sender clusters: using all clusters except receiver (", 
            length(potential_senders), " clusters)")
    return(potential_senders)
  } else {
    # Validate provided sender clusters
    invalid_senders <- base::setdiff(sender_clusters, potential_senders)
    if (length(invalid_senders) > 0) {
      warning("Some sender clusters are invalid or same as receiver: ", 
              paste(invalid_senders, collapse = ", "))
      sender_clusters <- base::setdiff(sender_clusters, invalid_senders)
    }
    
    if (length(sender_clusters) == 0) {
      stop("No valid sender clusters provided")
    }
    
    message("Using ", length(sender_clusters), " specified sender cluster(s)")
    return(sender_clusters)
  }
}

#' Filter Expressed Genes
#'
#' Filters genes based on expression percentage in specified cell types.
#'
#' @param sobj Seurat object
#' @param cell_types Character vector, cell type IDs to check expression
#' @param min_pct_expressed Numeric, minimum percentage of cells expressing the gene (default: 0.10)
#' @param assay_name Character string, assay to use (default: "SCT")
#' @param cluster_col Character string, name of metadata column with cluster annotations (default: NULL, uses Idents)
#'
#' @return Character vector of expressed gene names
#' @export
filter_expressed_genes <- function(sobj, cell_types, min_pct_expressed = 0.10, assay_name = "SCT", cluster_col = NULL) {
  
  if (!requireNamespace("nichenetr", quietly = TRUE)) {
    stop("nichenetr package is required for filtering expressed genes")
  }
  
  Seurat::DefaultAssay(sobj) <- assay_name
  
  # Set Idents if cluster_col is provided
  if (!is.null(cluster_col)) {
    Seurat::Idents(sobj) <- sobj@meta.data[[cluster_col]]
  }
  
  # Get expressed genes for each cell type with progress logging
  n_celltypes <- length(cell_types)
  expressed_genes_list <- list()
  for (i in seq_along(cell_types)) {
    ct <- cell_types[i]
    if (n_celltypes > 5 && i %% max(1, floor(n_celltypes / 10)) == 0) {
      message("    Progress: ", round(i / n_celltypes * 100, 0), "% (", i, "/", n_celltypes, " cell types)")
    }
    tryCatch({
      # nichenetr::get_expressed_genes expects the cell type to be in Idents
      # It uses Seurat::GetAssayData and checks expression percentage
      expressed_genes_list[[i]] <- nichenetr::get_expressed_genes(ct, sobj, min_pct_expressed)
    }, error = function(e) {
      warning("Error getting expressed genes for ", ct, ": ", e$message)
      return(character(0))
    })
  }
  
  # Combine and get unique genes
  expressed_genes <- unique(unlist(expressed_genes_list))
  
  return(expressed_genes)
}

#' Prepare CCI Data Summary
#'
#' Creates a summary of prepared data for CCI analysis.
#'
#' @param receiver_degs Data frame with receiver DEGs
#' @param sender_clusters Character vector of sender cluster IDs
#' @param receiver_cluster Character string, receiver cluster ID
#' @param expressed_genes_sender Character vector of expressed genes in senders
#' @param expressed_genes_receiver Character vector of expressed genes in receiver
#'
#' @return List with summary information
#' @export
prepare_cci_summary <- function(receiver_degs, sender_clusters, receiver_cluster,
                                 expressed_genes_sender, expressed_genes_receiver) {
  
  summary_list <- list(
    receiver_cluster = receiver_cluster,
    sender_clusters = sender_clusters,
    n_sender_clusters = length(sender_clusters),
    n_receiver_degs = nrow(receiver_degs),
    n_expressed_genes_sender = length(expressed_genes_sender),
    n_expressed_genes_receiver = length(expressed_genes_receiver),
    receiver_deg_genes = unique(receiver_degs$gene)
  )
  
  return(summary_list)
}


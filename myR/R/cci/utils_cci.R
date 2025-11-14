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


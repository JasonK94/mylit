# ============================================================================
# Aggregate Cluster-level DEG Consensus Results
# ============================================================================
# 여러 클러스터의 DEG consensus 결과를 gene level로 통합
# ============================================================================

#' Aggregate Cluster-level DEG Consensus Results
#'
#' @description
#' Aggregates DEG consensus results across multiple clusters to generate
#' gene-level meta-statistics.
#'
#' @param deg_list Named list of consensus_scores data frames from each cluster.
#'   Each element should be a data frame with columns: gene, mean_beta, meta_p,
#'   meta_p_adj, avg_log2FC, p_val, p_val_adj, cluster (optional).
#'   List names should be cluster names.
#' @param meta_p_method Method for meta-analysis of p-values across clusters:
#'   "stouffer" (default), "fisher", "inverse_variance"
#' @param mean_beta_method Method for combining mean_beta across clusters:
#'   "simple_mean" (default), "weighted_mean", "inverse_variance_weighted"
#' @param se_matrix Optional standard error matrix (genes × clusters) for
#'   inverse_variance methods
#' @param cluster_weights Optional weights for each cluster (named vector or list)
#'
#' @details
#' This function combines DEG results across clusters by:
#' 1. Identifying all unique genes across clusters
#' 2. For each gene, computing:
#'    - n_clusters: number of clusters where gene is present
#'    - mean_mean_beta: mean of mean_beta across clusters
#'    - sd_mean_beta: standard deviation of mean_beta
#'    - meta_meta_p: meta-analysis of meta_p across clusters
#'    - sd_meta_p: standard deviation of meta_p
#'    - concordance: proportion of clusters where mean_beta has same sign as mean_mean_beta
#'    - cluster_presence: binary indicators for each cluster (0/1 or TRUE/FALSE)
#' 3. Using specified meta-analysis methods for p-values and effect sizes
#'
#' @return Data frame with columns:
#'   - gene: gene identifier
#'   - n_clusters: number of clusters where gene is present
#'   - mean_mean_beta: mean of mean_beta across clusters
#'   - sd_mean_beta: standard deviation of mean_beta
#'   - meta_meta_p: meta-analysis p-value across clusters
#'   - meta_meta_p_adj: BH-adjusted meta_meta_p
#'   - sd_meta_p: standard deviation of meta_p
#'   - concordance: proportion of clusters with consistent sign
#'   - cluster_*: binary indicators (0/1) for each cluster presence
#'
#' @export
aggregate_cluster_deg_consensus <- function(
  deg_list,
  meta_p_method = c("stouffer", "fisher", "inverse_variance"),
  mean_beta_method = c("simple_mean", "weighted_mean", "inverse_variance_weighted"),
  se_matrix = NULL,
  cluster_weights = NULL
) {
  
  if (length(deg_list) == 0) {
    stop("deg_list가 비어있습니다.")
  }
  
  meta_p_method <- match.arg(meta_p_method)
  mean_beta_method <- match.arg(mean_beta_method)
  
  # Extract cluster names
  cluster_names <- names(deg_list)
  if (is.null(cluster_names)) {
    cluster_names <- paste0("cluster_", seq_along(deg_list))
    names(deg_list) <- cluster_names
  }
  
  # Standardize column names and ensure required columns exist
  deg_list_std <- lapply(seq_along(deg_list), function(i) {
    df <- deg_list[[i]]
    cluster_name <- cluster_names[i]
    
    # Ensure cluster column exists
    if (!"cluster" %in% colnames(df)) {
      df$cluster <- cluster_name
    }
    
    # Standardize column names
    if ("p_val" %in% colnames(df) && !"meta_p" %in% colnames(df)) {
      df$meta_p <- df$p_val
    }
    if ("p_val_adj" %in% colnames(df) && !"meta_p_adj" %in% colnames(df)) {
      df$meta_p_adj <- df$p_val_adj
    }
    if ("avg_log2FC" %in% colnames(df) && !"mean_beta" %in% colnames(df)) {
      df$mean_beta <- df$avg_log2FC
    }
    
    # Ensure required columns
    if (!"gene" %in% colnames(df)) {
      stop(sprintf("Cluster '%s'의 데이터에 'gene' 컬럼이 없습니다.", cluster_name))
    }
    if (!"mean_beta" %in% colnames(df)) {
      warning(sprintf("Cluster '%s'의 데이터에 'mean_beta' 컬럼이 없습니다. NA로 처리합니다.", cluster_name))
      df$mean_beta <- NA_real_
    }
    if (!"meta_p" %in% colnames(df)) {
      warning(sprintf("Cluster '%s'의 데이터에 'meta_p' 컬럼이 없습니다. NA로 처리합니다.", cluster_name))
      df$meta_p <- NA_real_
    }
    
    return(df)
  })
  names(deg_list_std) <- cluster_names
  
  # Collect all unique genes
  all_genes <- unique(unlist(lapply(deg_list_std, function(x) x$gene)))
  all_genes <- sort(all_genes[!is.na(all_genes)])
  
  if (length(all_genes) == 0) {
    stop("유전자가 발견되지 않았습니다.")
  }
  
  # Initialize result data frame
  result <- data.frame(
    gene = all_genes,
    n_clusters = 0L,
    mean_mean_beta = NA_real_,
    sd_mean_beta = NA_real_,
    meta_meta_p = NA_real_,
    meta_meta_p_adj = NA_real_,
    sd_meta_p = NA_real_,
    concordance = NA_real_,
    stringsAsFactors = FALSE
  )
  
  # Add cluster presence columns
  for (clust in cluster_names) {
    result[[paste0("cluster_", make.names(clust))]] <- 0L
  }
  
  # Cluster weights
  if (is.null(cluster_weights)) {
    cluster_weights <- rep(1, length(cluster_names))
    names(cluster_weights) <- cluster_names
  } else {
    if (is.null(names(cluster_weights))) {
      if (length(cluster_weights) == length(cluster_names)) {
        names(cluster_weights) <- cluster_names
      } else {
        warning("cluster_weights 길이가 cluster 수와 일치하지 않습니다. 동일 가중치를 사용합니다.")
        cluster_weights <- rep(1, length(cluster_names))
        names(cluster_weights) <- cluster_names
      }
    }
    # Fill missing weights with 1
    missing_weights <- setdiff(cluster_names, names(cluster_weights))
    if (length(missing_weights) > 0) {
      cluster_weights[missing_weights] <- 1
    }
  }
  
  # Process each gene
  for (i in seq_along(all_genes)) {
    gene <- all_genes[i]
    
    # Collect values across clusters
    beta_values <- numeric()
    p_values <- numeric()
    cluster_present <- character()
    beta_signs <- integer()
    
    for (clust in cluster_names) {
      df <- deg_list_std[[clust]]
      gene_idx <- which(df$gene == gene)
      
      if (length(gene_idx) > 0) {
        gene_idx <- gene_idx[1]  # Take first if duplicate
        clust_col <- paste0("cluster_", make.names(clust))
        result[[clust_col]][i] <- 1L
        
        beta_val <- df$mean_beta[gene_idx]
        p_val <- df$meta_p[gene_idx]
        
        if (!is.na(beta_val)) {
          beta_values <- c(beta_values, beta_val)
          beta_signs <- c(beta_signs, sign(beta_val))
          cluster_present <- c(cluster_present, clust)
        }
        
        if (!is.na(p_val) && p_val > 0 && p_val < 1) {
          p_values <- c(p_values, p_val)
        }
      }
    }
    
    # n_clusters
    result$n_clusters[i] <- length(beta_values)
    
    if (result$n_clusters[i] == 0) next
    
    # mean_mean_beta and sd_mean_beta (based on mean_beta_method)
    if (length(beta_values) > 0) {
      # Calculate mean_mean_beta based on method
      if (mean_beta_method == "simple_mean") {
        result$mean_mean_beta[i] <- mean(beta_values)
      } else if (mean_beta_method == "weighted_mean") {
        w_valid <- cluster_weights[cluster_present]
        w_valid[is.na(w_valid)] <- 1
        result$mean_mean_beta[i] <- sum(beta_values * w_valid) / sum(w_valid)
      } else if (mean_beta_method == "inverse_variance_weighted") {
        if (is.null(se_matrix)) {
          # Fallback to simple mean
          warning(sprintf("se_matrix가 없어 gene '%s'에 대해 simple_mean을 사용합니다.", gene))
          result$mean_mean_beta[i] <- mean(beta_values)
        } else {
          # Get SE values for this gene across clusters
          se_valid <- numeric(length(cluster_present))
          for (j in seq_along(cluster_present)) {
            clust_j <- cluster_present[j]
            if (gene %in% rownames(se_matrix) && clust_j %in% colnames(se_matrix)) {
              se_valid[j] <- se_matrix[gene, clust_j]
            } else {
              se_valid[j] <- NA_real_
            }
          }
          se_valid[is.na(se_valid) | se_valid <= 0] <- Inf
          weights_iv <- 1 / (se_valid^2)
          weights_iv[is.infinite(weights_iv)] <- 0
          if (sum(weights_iv) == 0) {
            # Fallback to simple mean
            result$mean_mean_beta[i] <- mean(beta_values)
          } else {
            result$mean_mean_beta[i] <- sum(beta_values * weights_iv) / sum(weights_iv)
          }
        }
      }
      
      # sd_mean_beta (standard deviation across clusters)
      if (length(beta_values) > 1) {
        result$sd_mean_beta[i] <- stats::sd(beta_values)
      } else {
        result$sd_mean_beta[i] <- NA_real_
      }
      
      # Concordance: proportion with same sign as mean_mean_beta
      mean_sign <- sign(result$mean_mean_beta[i])
      if (mean_sign == 0) {
        # If mean is exactly 0, count all non-zero values as concordant
        concordant <- sum(beta_signs != 0)
        result$concordance[i] <- if (length(beta_signs) > 0) concordant / length(beta_signs) else NA_real_
      } else {
        concordant <- sum(beta_signs == mean_sign | beta_signs == 0)
        result$concordance[i] <- concordant / length(beta_signs)
      }
    }
    
    # sd_meta_p
    if (length(p_values) > 1) {
      result$sd_meta_p[i] <- stats::sd(p_values)
    } else {
      result$sd_meta_p[i] <- NA_real_
    }
    
    # meta_meta_p: meta-analysis of p-values across clusters
    if (length(p_values) > 0) {
      # Collect p-values with cluster information for weighting
      p_valid <- p_values
      beta_for_p <- beta_values[match(cluster_present, cluster_names)]
      valid_for_p <- !is.na(p_valid) & p_valid > 0 & p_valid < 1 & !is.na(beta_for_p)
      
      if (sum(valid_for_p) > 0) {
        p_valid <- p_valid[valid_for_p]
        beta_valid <- beta_for_p[valid_for_p]
        clust_valid <- cluster_present[valid_for_p]
        
        # Clip p-values
        p_valid <- pmin(pmax(p_valid, 1e-300), 1 - 1e-16)
        
        if (meta_p_method == "stouffer") {
          # Stouffer's Z-score
          z_scores <- sign(beta_valid) * stats::qnorm(1 - p_valid / 2)
          w_valid <- cluster_weights[clust_valid]
          w_valid[is.na(w_valid)] <- 1
          z_comb <- sum(w_valid * z_scores) / sqrt(sum(w_valid^2))
          result$meta_meta_p[i] <- 2 * (1 - stats::pnorm(abs(z_comb)))
          
        } else if (meta_p_method == "fisher") {
          # Fisher's combined p-value
          chi2 <- -2 * sum(log(p_valid))
          df_fisher <- 2 * length(p_valid)
          p_combined <- stats::pchisq(chi2, df_fisher, lower.tail = FALSE)
          beta_sign_consistent <- all(sign(beta_valid) >= 0) || all(sign(beta_valid) <= 0)
          if (beta_sign_consistent) {
            result$meta_meta_p[i] <- p_combined / 2
          } else {
            result$meta_meta_p[i] <- p_combined
          }
          
        } else if (meta_p_method == "inverse_variance") {
          # Inverse variance weighting (requires se)
          if (is.null(se_matrix)) {
            warning("se_matrix가 없어 inverse_variance 방법을 사용할 수 없습니다. Stouffer로 대체합니다.")
            z_scores <- sign(beta_valid) * stats::qnorm(1 - p_valid / 2)
            z_comb <- mean(z_scores)
            result$meta_meta_p[i] <- 2 * (1 - stats::pnorm(abs(z_comb)))
          } else {
            # Use se_matrix for weighting
            # se_matrix should be genes × clusters matrix or data.frame
            se_valid <- numeric(length(clust_valid))
            for (j in seq_along(clust_valid)) {
              clust_j <- clust_valid[j]
              if (gene %in% rownames(se_matrix) && clust_j %in% colnames(se_matrix)) {
                se_valid[j] <- se_matrix[gene, clust_j]
              } else {
                se_valid[j] <- NA_real_
              }
            }
            se_valid[is.na(se_valid) | se_valid <= 0] <- Inf
            
            weights_iv <- 1 / (se_valid^2)
            weights_iv[is.infinite(weights_iv)] <- 0
            if (sum(weights_iv) == 0) {
              # Fallback
              z_scores <- sign(beta_valid) * stats::qnorm(1 - p_valid / 2)
              z_comb <- mean(z_scores)
              result$meta_meta_p[i] <- 2 * (1 - stats::pnorm(abs(z_comb)))
            } else {
              beta_weighted <- sum(beta_valid * weights_iv) / sum(weights_iv)
              se_combined <- sqrt(1 / sum(weights_iv))
              z_comb <- beta_weighted / se_combined
              result$meta_meta_p[i] <- 2 * (1 - stats::pnorm(abs(z_comb)))
            }
          }
        }
      }
    }
  }
  
  # BH adjustment for meta_meta_p
  valid_meta <- !is.na(result$meta_meta_p)
  if (any(valid_meta)) {
    result$meta_meta_p_adj[valid_meta] <- stats::p.adjust(result$meta_meta_p[valid_meta], method = "BH")
  }
  
  return(result)
}


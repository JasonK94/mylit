# ============================================================================
# Consensus Analysis (Phase 4-5)
# ============================================================================
# 방법론 수준 클러스터링 및 Consensus DEG signature 생성
# ============================================================================

#' Compute Agreement Scores
#'
#' @description
#' Computes agreement scores for each gene across methods.
#'
#' @param significance_matrix Significance indicator matrix (genes × methods)
#'
#' @return Named vector of agreement scores (one per gene)
#'
#' @export
compute_agreement_scores <- function(significance_matrix) {
  if (is.null(significance_matrix) || nrow(significance_matrix) == 0) {
    stop("significance_matrix가 비어있습니다.")
  }
  
  # 각 유전자에 대한 방법론 간 일치도 계산
  agreement_scores <- rowMeans(significance_matrix, na.rm = TRUE)
  names(agreement_scores) <- rownames(significance_matrix)
  
  return(agreement_scores)
}

#' Perform PCA on DEG Methods
#'
#' @description
#' Performs PCA on the method space using beta, logp, and significance matrices.
#'
#' @param deg_matrices List of DEG matrices from build_deg_matrices()
#' @param n_components Number of principal components (default: min(10, n_methods-1))
#'
#' @return List with PCA results and coordinates
#'
#' @export
perform_deg_pca <- function(
  deg_matrices,
  n_components = NULL
) {
  
  if (!requireNamespace("stats", quietly = TRUE)) {
    stop("stats 패키지가 필요합니다.")
  }
  
  beta_matrix <- deg_matrices$beta
  logp_matrix <- deg_matrices$logp
  significance_matrix <- deg_matrices$significance
  
  # 결합된 특성 행렬 생성 (각 방법론에 대해)
  # 방법론별로 평균 beta, 평균 logp, 평균 significance를 사용
  method_features <- cbind(
    colMeans(abs(beta_matrix), na.rm = TRUE),  # 절댓값 사용
    colMeans(logp_matrix, na.rm = TRUE),
    colMeans(significance_matrix, na.rm = TRUE),
    colSums(significance_matrix, na.rm = TRUE)  # 유의한 유전자 수
  )
  colnames(method_features) <- c("mean_abs_beta", "mean_logp", "mean_significance", "n_significant")
  rownames(method_features) <- colnames(beta_matrix)
  
  # PCA 수행
  if (is.null(n_components)) {
    n_components <- min(10, ncol(method_features) - 1, nrow(method_features) - 1)
  }
  
  pca_result <- stats::prcomp(method_features, scale. = TRUE, center = TRUE)
  
  return(list(
    pca = pca_result,
    coordinates = pca_result$x[, 1:min(n_components, ncol(pca_result$x)), drop = FALSE],
    variance_explained = summary(pca_result)$importance[2, ],
    method_features = method_features
  ))
}

#' Cluster DEG Methods
#'
#' @description
#' Performs clustering on DEG methods using hierarchical clustering.
#'
#' @param deg_matrices List of DEG matrices from build_deg_matrices()
#' @param method Clustering method: "hierarchical", "kmeans" (default: "hierarchical")
#' @param k Number of clusters for kmeans (default: NULL, auto-detect)
#' @param distance Distance metric for hierarchical clustering (default: "euclidean")
#'
#' @return List with clustering results
#'
#' @export
cluster_deg_methods <- function(
  deg_matrices,
  method = c("hierarchical", "kmeans"),
  k = NULL,
  distance = "euclidean"
) {
  
  method <- match.arg(method)
  
  beta_matrix <- deg_matrices$beta
  logp_matrix <- deg_matrices$logp
  significance_matrix <- deg_matrices$significance
  
  # 방법론별 특성 추출
  method_features <- cbind(
    colMeans(beta_matrix, na.rm = TRUE),
    colMeans(logp_matrix, na.rm = TRUE),
    colMeans(significance_matrix, na.rm = TRUE),
    colSums(significance_matrix, na.rm = TRUE)  # 유의한 유전자 수
  )
  colnames(method_features) <- c("mean_beta", "mean_logp", "mean_significance", "n_significant")
  rownames(method_features) <- colnames(beta_matrix)
  
  # 정규화
  method_features_scaled <- scale(method_features)
  
  if (method == "hierarchical") {
    # 거리 행렬 계산
    dist_matrix <- dist(t(method_features_scaled), method = distance)
    
    # 계층적 클러스터링
    hclust_result <- stats::hclust(dist_matrix, method = "ward.D2")
    
    # 클러스터 할당 (k가 지정되지 않으면 자동 선택)
    if (is.null(k)) {
      # 간단한 휴리스틱: 방법론 수의 절반
      k <- max(2, floor(nrow(method_features) / 2))
    }
    
    clusters <- stats::cutree(hclust_result, k = k)
    
    return(list(
      hclust = hclust_result,
      clusters = clusters,
      method_features = method_features,
      k = k
    ))
    
  } else if (method == "kmeans") {
    if (is.null(k)) {
      k <- max(2, floor(nrow(method_features) / 2))
    }
    
    kmeans_result <- stats::kmeans(method_features_scaled, centers = k, nstart = 10)
    
    return(list(
      kmeans = kmeans_result,
      clusters = kmeans_result$cluster,
      method_features = method_features,
      k = k
    ))
  }
}

#' Compute Consensus Scores
#'
#' @description
#' Computes consensus scores for each gene based on agreement and clustering.
#'
#' @param deg_matrices List of DEG matrices from build_deg_matrices()
#' @param agreement_scores Agreement scores from compute_agreement_scores()
#' @param clustering_results Clustering results from cluster_deg_methods()
#' @param weights Optional weights for each method (default: NULL, equal weights)
#' 
#' @details
#' In addition to simple averages, this function performs a meta-analysis of
#' raw p-values across methods using Stouffer's Z-score method (with sign taken
#' from logFC). It returns per-gene meta-analysis statistics:
#'   - \code{meta_z}: combined Z-score
#'   - \code{meta_p}: two-sided meta-analysis p-value
#'   - \code{meta_p_adj}: BH-adjusted meta-analysis p-value (across genes)
#' The final \code{consensus_score} combines effect size, agreement, and
#' meta-analysis evidence as:
#'   \deqn{consensus\_score = (agreement × |weighted\_beta|) × -log10(meta\_p)}
#'
#' @return Data frame with consensus scores
#'
#' @export
compute_consensus_scores <- function(
  deg_matrices,
  agreement_scores,
  clustering_results = NULL,
  weights = NULL
) {
  
  beta_matrix <- deg_matrices$beta
  pvalue_matrix <- deg_matrices$pvalue
  significance_matrix <- deg_matrices$significance
  methods <- colnames(beta_matrix)
  
  # 가중치 설정
  if (is.null(weights)) {
    weights <- rep(1, length(methods))
    names(weights) <- methods
  }
  
  # 각 유전자에 대한 consensus score 계산
  genes <- rownames(beta_matrix)
  
  consensus_scores <- data.frame(
    gene = genes,
    mean_beta = rowMeans(beta_matrix, na.rm = TRUE),
    mean_logp = rowMeans(deg_matrices$logp, na.rm = TRUE),
    agreement = agreement_scores[genes],
    n_significant = rowSums(significance_matrix, na.rm = TRUE),
    n_methods = length(methods),
    stringsAsFactors = FALSE
  )
  
  # 가중 평균 beta
  consensus_scores$weighted_beta <- apply(beta_matrix, 1, function(x) {
    valid <- !is.na(x)
    if (sum(valid) == 0) return(NA_real_)
    method_names <- colnames(beta_matrix)[valid]
    x_valid <- x[valid]
    w_valid <- weights[method_names]
    # weights가 없는 경우 1로 처리
    w_valid[is.na(w_valid)] <- 1
    sum(x_valid * w_valid) / sum(w_valid)
  })
  
  # --- Meta-analysis of p-values across methods (Stouffer's Z) ---
  meta_z <- rep(NA_real_, length(genes))
  meta_p <- rep(NA_real_, length(genes))
  
  for (i in seq_along(genes)) {
    p_row <- pvalue_matrix[i, ]
    beta_row <- beta_matrix[i, ]
    valid <- !is.na(p_row) & !is.na(beta_row) & p_row > 0 & p_row < 1
    if (!any(valid)) next
    
    p_valid <- p_row[valid]
    beta_valid <- beta_row[valid]
    method_names <- methods[valid]
    
    # clip p-values for numerical stability
    p_valid <- pmin(pmax(p_valid, 1e-300), 1 - 1e-16)
    # two-sided Z-score with sign from beta
    z_scores <- sign(beta_valid) * stats::qnorm(1 - p_valid / 2)
    
    w_valid <- weights[method_names]
    w_valid[is.na(w_valid)] <- 1
    
    z_comb <- sum(w_valid * z_scores) / sqrt(sum(w_valid^2))
    meta_z[i] <- z_comb
    meta_p[i] <- 2 * (1 - stats::pnorm(abs(z_comb)))
  }
  
  consensus_scores$meta_z <- meta_z
  consensus_scores$meta_p <- meta_p
  
  # BH-adjusted meta p-values (across genes)
  meta_p_adj <- rep(NA_real_, length(meta_p))
  valid_meta <- !is.na(meta_p)
  if (any(valid_meta)) {
    meta_p_adj[valid_meta] <- stats::p.adjust(meta_p[valid_meta], method = "BH")
  }
  consensus_scores$meta_p_adj <- meta_p_adj
  
  # Signal/evidence components and final consensus score
  eps <- 1e-300
  consensus_scores$score_signal <- consensus_scores$agreement * abs(consensus_scores$weighted_beta)
  consensus_scores$score_evidence <- -log10(pmax(consensus_scores$meta_p, eps))
  consensus_scores$consensus_score <- consensus_scores$score_signal * consensus_scores$score_evidence
  
  # 클러스터 정보 추가 (있는 경우)
  if (!is.null(clustering_results) && "clusters" %in% names(clustering_results)) {
    # 클러스터별 가중치 계산 가능
    # 일단 기본 구현만
  }
  
  return(consensus_scores)
}

#' Generate Consensus DEG List
#'
#' @description
#' Generates final consensus DEG list based on consensus scores.
#'
#' @param consensus_scores Consensus scores from compute_consensus_scores()
#' @param fdr_threshold FDR threshold (default: 0.05)
#' @param agreement_threshold Minimum agreement score (default: 0.5)
#' @param min_methods Minimum number of methods showing significance (default: NULL)
#'
#' @return Data frame with consensus DEG list
#'
#' @export
generate_consensus_deg_list <- function(
  consensus_scores,
  fdr_threshold = 0.05,
  agreement_threshold = 0.5,
  min_methods = NULL
) {
  
  # 필터링 조건
  keep <- rep(TRUE, nrow(consensus_scores))
  
  # Agreement threshold
  if (!is.null(agreement_threshold)) {
    keep <- keep & (!is.na(consensus_scores$agreement)) & 
                   (consensus_scores$agreement >= agreement_threshold)
  }
  
  # Minimum methods
  if (!is.null(min_methods)) {
    keep <- keep & (!is.na(consensus_scores$n_significant)) & 
                   (consensus_scores$n_significant >= min_methods)
  }
  
  # Meta-analysis FDR threshold (if meta_p_adj available)
  if (!is.null(fdr_threshold) && "meta_p_adj" %in% colnames(consensus_scores)) {
    keep <- keep & (!is.na(consensus_scores$meta_p_adj)) &
                   (consensus_scores$meta_p_adj <= fdr_threshold)
  }
  
  consensus_deg <- consensus_scores[keep, ]
  
  if (nrow(consensus_deg) == 0) {
    warning("필터링 후 유전자가 없습니다. threshold를 낮춰보세요.")
    return(consensus_deg)
  }
  
  # 정렬 (consensus_score 기준, NA는 마지막으로)
  consensus_deg <- consensus_deg[order(-consensus_deg$consensus_score, na.last = TRUE), ]
  
  return(consensus_deg)
}


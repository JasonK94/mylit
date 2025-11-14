# ============================================================================
# DEG Results Standardization (Phase 3)
# ============================================================================
# 여러 방법론의 결과를 표준화된 형식으로 통합
# ============================================================================

#' Standardize DEG Results from Different Methods
#'
#' @description
#' Converts DEG results from different methods into a standardized format
#' with consistent column names and structure.
#'
#' @param deg_result DEG result from a method (data frame or list)
#' @param method_name Name of the method (e.g., "muscat-edgeR", "limma-voom")
#' @param source_function Name of the source function (optional, for debugging)
#'
#' @return Standardized data frame with columns:
#'   - gene: Gene identifier
#'   - logFC: Log fold change (effect size)
#'   - pvalue: P-value
#'   - pvalue_adj: Adjusted p-value (FDR)
#'   - statistic: Test statistic
#'   - se: Standard error (if available)
#'   - method: Method name
#'   - cluster_id: Cluster ID (if available)
#'
#' @export
standardize_deg_results <- function(
  deg_result,
  method_name,
  source_function = NULL
) {
  
  if (is.null(deg_result)) {
    stop("deg_result가 NULL입니다.")
  }
  
  # 결과 형식 확인 및 변환
  if (is.data.frame(deg_result)) {
    df <- deg_result
  } else if (is.list(deg_result)) {
    # 리스트인 경우
    if ("data.frame" %in% class(deg_result)) {
      # 이미 데이터프레임인 리스트 (예: tibble)
      df <- as.data.frame(deg_result)
    } else if (length(deg_result) > 0 && is.data.frame(deg_result[[1]])) {
      # 리스트의 첫 번째 요소가 데이터프레임
      df <- deg_result[[1]]
      warning(sprintf("리스트의 첫 번째 요소를 사용합니다 (method: %s)", method_name))
    } else {
      stop(sprintf("지원하지 않는 결과 형식입니다 (method: %s, class: %s)", 
                   method_name, paste(class(deg_result), collapse=", ")))
    }
  } else {
    stop(sprintf("지원하지 않는 결과 형식입니다 (method: %s, class: %s)", 
                 method_name, paste(class(deg_result), collapse=", ")))
  }
  
  # 결과가 비어있는 경우
  if (nrow(df) == 0) {
    warning(sprintf("방법론 '%s'의 결과가 비어있습니다.", method_name))
    return(data.frame(
      gene = character(0),
      logFC = numeric(0),
      pvalue = numeric(0),
      pvalue_adj = numeric(0),
      statistic = numeric(0),
      method = character(0),
      stringsAsFactors = FALSE
    ))
  }
  
  # 표준화된 결과 데이터프레임 생성
  result_std <- data.frame(
    gene = character(nrow(df)),
    logFC = numeric(nrow(df)),
    pvalue = numeric(nrow(df)),
    pvalue_adj = numeric(nrow(df)),
    statistic = numeric(nrow(df)),
    method = rep(method_name, nrow(df)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # se와 cluster_id는 선택적이므로 나중에 추가
  
  # gene 컬럼 찾기
  gene_cols <- c("gene", "Gene", "gene_id", "GeneID", "primerid", "gene_name", "symbol")
  gene_col <- NULL
  for (col in gene_cols) {
    if (col %in% colnames(df)) {
      gene_col <- col
      break
    }
  }
  if (is.null(gene_col) && "gene" %in% colnames(df)) {
    gene_col <- "gene"
  } else if (is.null(gene_col)) {
    # rownames 사용
    if (!is.null(rownames(df)) && length(rownames(df)) > 0) {
      result_std$gene <- rownames(df)
    } else {
      stop(sprintf("유전자 ID를 찾을 수 없습니다 (method: %s)", method_name))
    }
  } else {
    result_std$gene <- as.character(df[[gene_col]])
  }
  
  # logFC 컬럼 찾기
  logfc_cols <- c("logFC", "log2FoldChange", "log_fold_change", "avg_log2FC", "beta", "coef", "coefficient")
  logfc_col <- NULL
  for (col in logfc_cols) {
    if (col %in% colnames(df)) {
      logfc_col <- col
      break
    }
  }
  if (!is.null(logfc_col)) {
    result_std$logFC <- as.numeric(df[[logfc_col]])
  } else {
    result_std$logFC <- NA_real_
  }
  
  # pvalue 컬럼 찾기
  pval_cols <- c("pvalue", "p_value", "PValue", "P.Value", "pval", "p", "p_val")
  pval_col <- NULL
  for (col in pval_cols) {
    if (col %in% colnames(df)) {
      pval_col <- col
      break
    }
  }
  if (!is.null(pval_col)) {
    result_std$pvalue <- as.numeric(df[[pval_col]])
  } else {
    result_std$pvalue <- NA_real_
  }
  
  # pvalue_adj 컬럼 찾기 (muscat의 p_adj.loc, p_adj.glb 포함)
  # 우선순위: local FDR (클러스터별) > global FDR > 기타
  # Note: muscat의 경우 local FDR이 더 많은 유의한 유전자를 제공
  padj_cols <- c("pvalue_adj", "p_val_adj", "adj_p_val", "FDR", "padj", "adj.P.Val", "fdr", "qvalue", "q",
                 "p_adj.loc", "p_adj_loc",   # local FDR (클러스터별, 우선)
                 "p_adj.glb", "p_adj_glb")   # global FDR (차선)
  padj_col <- NULL
  for (col in padj_cols) {
    if (col %in% colnames(df)) {
      padj_col <- col
      break
    }
  }
  if (!is.null(padj_col)) {
    result_std$pvalue_adj <- as.numeric(df[[padj_col]])
  } else {
    # pvalue로부터 계산 시도
    if (!is.null(pval_col) && !all(is.na(result_std$pvalue))) {
      result_std$pvalue_adj <- stats::p.adjust(result_std$pvalue, method = "BH")
    } else {
      result_std$pvalue_adj <- NA_real_
    }
  }
  
  # statistic 컬럼 찾기
  stat_cols <- c("statistic", "stat", "t", "F", "LR", "z", "Wald", "LRT")
  stat_col <- NULL
  for (col in stat_cols) {
    if (col %in% colnames(df)) {
      stat_col <- col
      break
    }
  }
  if (!is.null(stat_col)) {
    result_std$statistic <- as.numeric(df[[stat_col]])
  } else {
    result_std$statistic <- NA_real_
  }
  
  # se 컬럼 찾기 (선택적) - 나중에 추가
  se_cols <- c("se", "SE", "std_err", "standard_error", "lfcSE")
  se_col <- NULL
  for (col in se_cols) {
    if (col %in% colnames(df)) {
      se_col <- col
      break
    }
  }
  if (!is.null(se_col)) {
    result_std$se <- as.numeric(df[[se_col]])
  }
  
  # cluster_id 컬럼 찾기 (선택적)
  cluster_cols <- c("cluster_id", "cluster", "cluster_label", "cell_type")
  cluster_col <- NULL
  for (col in cluster_cols) {
    if (col %in% colnames(df)) {
      cluster_col <- col
      break
    }
  }
  if (!is.null(cluster_col)) {
    result_std$cluster_id <- as.character(df[[cluster_col]])
  } else {
    # cluster_id가 없으면 NA로 채움
    result_std$cluster_id <- NA_character_
  }
  
  # NA 제거 (gene은 필수)
  result_std <- result_std[!is.na(result_std$gene) & result_std$gene != "", ]
  
  # 컬럼 순서 정리
  base_cols <- c("gene", "logFC", "pvalue", "pvalue_adj", "statistic", "method")
  optional_cols <- c("se", "cluster_id")
  final_cols <- c(base_cols, optional_cols[optional_cols %in% colnames(result_std)])
  result_std <- result_std[, final_cols, drop = FALSE]
  
  return(result_std)
}

#' Build DEG Matrices from Standardized Results
#'
#' @description
#' Constructs gene × method matrices from standardized DEG results.
#'
#' @param standardized_results_list Named list of standardized DEG results
#' @param genes Optional vector of gene IDs. If NULL, uses intersection of all methods.
#' @param fdr_threshold FDR threshold for significance (default: 0.05)
#'
#' @return List with matrices:
#'   - beta: Effect size matrix (genes × methods)
#'   - pvalue: P-value matrix
#'   - logp: -log10(pvalue) matrix
#'   - significance: Significance indicator matrix (FDR < threshold)
#'   - statistic: Test statistic matrix
#'   - se: Standard error matrix (if available)
#'
#' @export
build_deg_matrices <- function(
  standardized_results_list,
  genes = NULL,
  fdr_threshold = 0.05
) {
  
  if (length(standardized_results_list) == 0) {
    stop("standardized_results_list가 비어있습니다.")
  }
  
  # 모든 방법론에서 유전자 수집
  all_genes <- unique(unlist(lapply(standardized_results_list, function(x) x$gene)))
  
  # genes가 지정되지 않으면 교집합 사용
  if (is.null(genes)) {
    genes <- Reduce(intersect, lapply(standardized_results_list, function(x) x$gene))
    if (length(genes) == 0) {
      warning("모든 방법론에서 공통 유전자가 없습니다. 합집합을 사용합니다.")
      genes <- all_genes
    }
  }
  
  genes <- sort(unique(genes))
  methods <- names(standardized_results_list)
  
  # 행렬 초기화
  beta_matrix <- matrix(NA_real_, nrow = length(genes), ncol = length(methods),
                       dimnames = list(genes, methods))
  pvalue_matrix <- matrix(NA_real_, nrow = length(genes), ncol = length(methods),
                          dimnames = list(genes, methods))
  logp_matrix <- matrix(NA_real_, nrow = length(genes), ncol = length(methods),
                       dimnames = list(genes, methods))
  significance_matrix <- matrix(FALSE, nrow = length(genes), ncol = length(methods),
                               dimnames = list(genes, methods))
  statistic_matrix <- matrix(NA_real_, nrow = length(genes), ncol = length(methods),
                            dimnames = list(genes, methods))
  
  # se 행렬 (선택적)
  has_se <- any(vapply(standardized_results_list, function(x) "se" %in% colnames(x), logical(1)))
  if (has_se) {
    se_matrix <- matrix(NA_real_, nrow = length(genes), ncol = length(methods),
                       dimnames = list(genes, methods))
  } else {
    se_matrix <- NULL
  }
  
  # 각 방법론별로 행렬 채우기
  for (method in methods) {
    result <- standardized_results_list[[method]]
    
    # 클러스터별 결과가 있는 경우 (같은 유전자가 여러 클러스터에 나타날 수 있음)
    # 각 유전자에 대해 가장 유의한 클러스터의 결과를 선택
    if ("cluster_id" %in% colnames(result) && any(duplicated(result$gene))) {
      # 유전자별로 가장 유의한 결과 선택 (최소 pvalue_adj 또는 최대 |logFC|)
      result_agg <- result[order(result$gene, 
                                  if ("pvalue_adj" %in% colnames(result)) {
                                    ifelse(is.na(result$pvalue_adj), Inf, result$pvalue_adj)
                                  } else {
                                    ifelse(is.na(result$pvalue), Inf, result$pvalue)
                                  }), ]
      result_agg <- result_agg[!duplicated(result_agg$gene), ]
      result <- result_agg
    }
    
    # gene을 인덱스로 사용
    match_idx <- match(result$gene, genes)
    valid_idx <- !is.na(match_idx)
    
    if (sum(valid_idx) == 0) {
      warning(sprintf("방법론 '%s'에서 유전자 매칭 실패", method))
      next
    }
    
    idx <- match_idx[valid_idx]
    
    # beta (logFC)
    if ("logFC" %in% colnames(result)) {
      beta_matrix[idx, method] <- result$logFC[valid_idx]
    }
    
    # pvalue
    if ("pvalue" %in% colnames(result)) {
      pvalue_matrix[idx, method] <- result$pvalue[valid_idx]
      logp_matrix[idx, method] <- -log10(pmax(result$pvalue[valid_idx], 1e-300))
    }
    
    # significance (FDR < threshold)
    if ("pvalue_adj" %in% colnames(result)) {
      significance_matrix[idx, method] <- !is.na(result$pvalue_adj[valid_idx]) & 
                                          result$pvalue_adj[valid_idx] < fdr_threshold
    }
    
    # statistic
    if ("statistic" %in% colnames(result)) {
      statistic_matrix[idx, method] <- result$statistic[valid_idx]
    }
    
    # se
    if (has_se && "se" %in% colnames(result)) {
      se_matrix[idx, method] <- result$se[valid_idx]
    }
  }
  
  result_list <- list(
    beta = beta_matrix,
    pvalue = pvalue_matrix,
    logp = logp_matrix,
    significance = significance_matrix,
    statistic = statistic_matrix,
    genes = genes,
    methods = methods
  )
  
  if (!is.null(se_matrix)) {
    result_list$se <- se_matrix
  }
  
  return(result_list)
}


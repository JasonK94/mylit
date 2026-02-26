#' Plot Cross-Cluster GSEA/GO Heatmap
#'
#' Creates a heatmap visualization comparing pathway enrichment results across
#' multiple clusters. Pathways are filtered by significance and frequency across
#' clusters, then displayed as a heatmap with hierarchical clustering.
#'
#' @param merged_df Data frame containing merged pathway analysis results from
#'   multiple clusters. Must contain columns specified by \code{x_var}, \code{y_var},
#'   \code{fill_var}, and \code{p_var}.
#' @param top_n_per_cluster Numeric. Number of top pathways to select per cluster
#'   (default: 10). Pathways are ranked by absolute value of \code{fill_var}.
#' @param n_clusters Numeric vector or NULL. Filter pathways that appear in
#'   exactly this many clusters (e.g., \code{c(2,3)} for pathways in 2 or 3 clusters).
#'   If NULL, no filtering is applied (default: NULL).
#' @param x_var Character string. Column name for x-axis grouping (default: "cluster").
#' @param y_var Character string. Column name for y-axis pathway labels (default: "description").
#' @param fill_var Character string. Column name for heatmap fill values (default: "effect_size").
#'   For GSEA results, use "ES" or "NES"; for GO/KEGG, use "effect_size".
#' @param color_label Character string or NULL. Label for the color scale legend.
#'   If NULL, uses default label (default: NULL).
#' @param p_var Character string. Column name for p-value/adjusted p-value
#'   (default: "adj_p_value"). Used for significance filtering and star annotations.
#' @param wrap_width Numeric or NULL. Width for text wrapping of pathway names.
#'   If NULL, no wrapping is applied (default: 50).
#' @return A ggplot2 object showing the cross-cluster heatmap, or NULL if no
#'   pathways match the filtering criteria.
#' @examples
#' \dontrun{
#' # Basic usage with GSEA results
#' heatmap_plot <- plot_compare_gsea_heatmap(
#'   all_gsea_data,
#'   top_n_per_cluster = 50,
#'   fill_var = "ES",
#'   wrap_width = 40
#' )
#' 
#' # Filter pathways appearing in 2 clusters
#' heatmap_plot <- plot_compare_gsea_heatmap(
#'   all_go_data,
#'   top_n_per_cluster = 50,
#'   n_clusters = 2,
#'   color_label = "Effect Size",
#'   wrap_width = 40
#' )
#' }
#' @importFrom dplyr group_by arrange slice_head ungroup filter pull mutate select
#' @importFrom dplyr n_distinct summarise case_when
#' @importFrom tidyr pivot_wider
#' @importFrom ggplot2 ggplot aes geom_tile geom_text scale_fill_gradient2 labs
#' @importFrom ggplot2 theme_bw theme element_text coord_fixed
#' @importFrom rlang .data sym
#' @importFrom stringr str_wrap
#' @importFrom stats dist hclust
#' @export
plot_compare_gsea_heatmap <- function(merged_df, 
                                      top_n_per_cluster = 10, 
                                      n_clusters = NULL,
                                      x_var = "cluster", 
                                      y_var = "description",
                                      fill_var = "effect_size",
                                      color_label = NULL,
                                      p_var = "adj_p_value",
                                      wrap_width = 50) {
  
  # [Step 0] 중복 제거 (GO 분석 대응)
  deduplicated_df <- merged_df %>%
    dplyr::group_by(.data[[x_var]], .data[[y_var]]) %>%
    dplyr::arrange(.data[[p_var]]) %>% 
    dplyr::slice_head(n = 1) %>%       
    dplyr::ungroup()
  
  # [Step 1] 클러스터 빈도수(Frequency) 필터링
  # 특정 Pathway가 몇 개의 클러스터에 등장하는지 계산
  pathway_counts <- deduplicated_df %>%
    dplyr::group_by(.data[[y_var]]) %>%
    dplyr::summarise(cluster_count = dplyr::n_distinct(.data[[x_var]]))
  
  # n_clusters가 지정된 경우 필터링
  if (!is.null(n_clusters)) {
    target_pathways <- pathway_counts %>%
      dplyr::filter(cluster_count %in% n_clusters) %>%
      dplyr::pull(.data[[y_var]])
    
    deduplicated_df <- deduplicated_df %>%
      dplyr::filter(.data[[y_var]] %in% target_pathways)
      
    if(nrow(deduplicated_df) == 0) {
      message("조건에 맞는 Pathway가 없습니다. n_clusters 조건을 확인하세요.")
      return(NULL)
    }
  }

  # [Step 2] 표시할 Pathway 선정 (각 클러스터별 상위 N개 합집합)
  top_pathways <- deduplicated_df %>%
    dplyr::group_by(.data[[x_var]]) %>%
    dplyr::filter(.data[[p_var]] < 0.05) %>% 
    dplyr::arrange(dplyr::desc(abs(.data[[fill_var]]))) %>%
    dplyr::slice_head(n = top_n_per_cluster) %>%
    dplyr::pull(.data[[y_var]]) %>%
    unique()
  
  # Pathway가 없으면 NULL 반환
  if (length(top_pathways) == 0) {
    message("조건에 맞는 유의한 Pathway가 없습니다.")
    return(NULL)
  }
  
  # [Step 3] 데이터 준비
  plot_data <- deduplicated_df %>%
    dplyr::filter(.data[[y_var]] %in% top_pathways) %>%
    dplyr::mutate(stars = dplyr::case_when(
      .data[[p_var]] < 0.001 ~ "***",
      .data[[p_var]] < 0.01  ~ "**",
      .data[[p_var]] < 0.05  ~ "*",
      TRUE ~ ""
    ))
  
  # 데이터가 비어있으면 NULL 반환
  if (nrow(plot_data) == 0) {
    message("플롯할 데이터가 없습니다.")
    return(NULL)
  }
  
  # 텍스트 줄바꿈
  if (!is.null(wrap_width)) {
    plot_data <- plot_data %>%
      dplyr::mutate(!!rlang::sym(y_var) := stringr::str_wrap(.data[[y_var]], width = wrap_width))
  }

  # [Step 4] Y축 정렬 (Clustering)
  mat <- plot_data %>% 
    dplyr::select(.data[[y_var]], .data[[x_var]], .data[[fill_var]]) %>%
    tidyr::pivot_wider(names_from = .data[[x_var]], values_from = .data[[fill_var]], values_fill = 0) %>%
    as.data.frame()
  
  # mat가 올바르게 생성되었는지 확인
  if (is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) {
    message("행렬 생성에 실패했습니다.")
    return(NULL)
  }
  
  # y_var 컬럼이 있는지 확인
  if (!y_var %in% colnames(mat)) {
    message("y_var 컬럼이 행렬에 없습니다.")
    return(NULL)
  }
  
  rownames(mat) <- mat[[y_var]]
  mat <- mat[, -1, drop = FALSE]
  
  # nrow, ncol가 유효한 숫자인지 확인하고 클러스터링 수행
  n_rows <- nrow(mat)
  n_cols <- ncol(mat)
  if (!is.na(n_rows) && !is.na(n_cols) && n_rows > 1 && n_cols > 1) {
      dd <- stats::dist(mat)
      hc <- stats::hclust(dd)
      pathway_order <- rownames(mat)[hc$order]
      plot_data[[y_var]] <- factor(plot_data[[y_var]], levels = pathway_order)
  } else {
    # 클러스터링 불가능한 경우 기본 순서 사용
    plot_data[[y_var]] <- factor(plot_data[[y_var]], levels = unique(plot_data[[y_var]]))
  }

  # [Step 5] 그리기
  ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[[x_var]], y = .data[[y_var]], fill = .data[[fill_var]])) +
    ggplot2::geom_tile(color = "white") + 
    ggplot2::geom_text(ggplot2::aes(label = stars), color = "black", size = 4) +
    
    ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                         name = color_label) +
    
    ggplot2::labs(x = "Cluster", y = NULL, title = "Cross-Cluster Heatmap",
         subtitle = paste0("Filtered by presence in ", 
                           ifelse(is.null(n_clusters), "any cluster", 
                                  paste(paste(n_clusters, collapse=","), "clusters"))),
         caption = "* P < 0.05, ** P < 0.01, *** P < 0.001") +
    
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 10),
      panel.grid = ggplot2::element_blank(), 
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    ggplot2::coord_fixed(ratio = 0.5)
}

#' Plot Cross-Cluster GSEA/GO Dot Plot
#'
#' Creates a dot plot visualization comparing pathway enrichment results across
#' multiple clusters. Pathways are filtered by significance and frequency across
#' clusters, then displayed as a dot plot where point size represents significance
#' and color represents effect size.
#'
#' @param merged_df Data frame containing merged pathway analysis results from
#'   multiple clusters. Must contain columns specified by \code{x_var}, \code{y_var},
#'   \code{color_var}, and \code{p_var}.
#' @param top_n_per_cluster Numeric. Number of top pathways to select per cluster
#'   (default: 10). Pathways are ranked by absolute value of \code{color_var}.
#' @param n_clusters Numeric vector or NULL. Filter pathways that appear in
#'   exactly this many clusters (e.g., \code{c(1,2,3)} for pathways in 1-3 clusters).
#'   If NULL, no filtering is applied (default: NULL).
#' @param x_var Character string. Column name for x-axis grouping (default: "cluster").
#' @param y_var Character string. Column name for y-axis pathway labels (default: "description").
#' @param color_var Character string. Column name for point color values (default: "effect_size").
#'   For GSEA results, use "ES" or "NES"; for GO/KEGG, use "effect_size".
#' @param color_label Character string or NULL. Label for the color scale legend.
#'   If NULL, uses default label (default: NULL).
#' @param p_var Character string. Column name for p-value/adjusted p-value
#'   (default: "adj_p_value"). Used for significance filtering, star annotations,
#'   and point size (via -log10 transformation).
#' @param wrap_width Numeric or NULL. Width for text wrapping of pathway names.
#'   If NULL, no wrapping is applied (default: 50).
#' @return A ggplot2 object showing the cross-cluster dot plot, or NULL if no
#'   pathways match the filtering criteria.
#' @examples
#' \dontrun{
#' # Basic usage with GSEA results
#' dot_plot <- plot_compare_gsea_dot(
#'   all_gsea_data,
#'   top_n_per_cluster = 50,
#'   color_var = "ES",
#'   color_label = "ES",
#'   wrap_width = 40
#' )
#' 
#' # Use raw p-value instead of adjusted
#' dot_plot <- plot_compare_gsea_dot(
#'   all_gsea_data,
#'   top_n_per_cluster = 50,
#'   color_var = "ES",
#'   p_var = "p_value",
#'   wrap_width = 40
#' )
#' 
#' # Filter pathways appearing in 2 clusters
#' dot_plot <- plot_compare_gsea_dot(
#'   all_go_data,
#'   top_n_per_cluster = 50,
#'   n_clusters = 2,
#'   wrap_width = 40
#' )
#' }
#' @importFrom dplyr group_by arrange slice_head ungroup filter pull mutate select
#' @importFrom dplyr n_distinct summarise case_when
#' @importFrom tidyr pivot_wider
#' @importFrom ggplot2 ggplot aes geom_tile geom_point geom_text scale_color_gradient2
#' @importFrom ggplot2 scale_size_continuous labs theme_bw theme element_text
#' @importFrom rlang .data sym
#' @importFrom stringr str_wrap
#' @importFrom stats dist hclust
#' @export
plot_compare_gsea_dot <- function(merged_df, 
                                  top_n_per_cluster = 10,
                                  n_clusters = NULL,
                                  x_var = "cluster", 
                                  y_var = "description",
                                  color_var = "effect_size", 
                                  color_label = NULL,
                                  p_var = "adj_p_value",
                                  wrap_width = 50) {
  
  # [Step 0] 중복 제거
  deduplicated_df <- merged_df %>%
    dplyr::group_by(.data[[x_var]], .data[[y_var]]) %>%
    dplyr::arrange(.data[[p_var]]) %>% 
    dplyr::slice_head(n = 1) %>%       
    dplyr::ungroup()
  
  # [Step 1] 클러스터 빈도수 필터링
  pathway_counts <- deduplicated_df %>%
    dplyr::group_by(.data[[y_var]]) %>%
    dplyr::summarise(cluster_count = dplyr::n_distinct(.data[[x_var]]))
  
  if (!is.null(n_clusters)) {
    target_pathways <- pathway_counts %>%
      dplyr::filter(cluster_count %in% n_clusters) %>%
      dplyr::pull(.data[[y_var]])
    
    deduplicated_df <- deduplicated_df %>%
      dplyr::filter(.data[[y_var]] %in% target_pathways)
    
    if(nrow(deduplicated_df) == 0) {
      message("조건에 맞는 Pathway가 없습니다.")
      return(NULL)
    }
  }
  
  # [Step 2] Pathway 선정
  top_pathways <- deduplicated_df %>%
    dplyr::group_by(.data[[x_var]]) %>%
    dplyr::filter(.data[[p_var]] < 0.05) %>% 
    dplyr::arrange(dplyr::desc(abs(.data[[color_var]]))) %>% 
    dplyr::slice_head(n = top_n_per_cluster) %>%
    dplyr::pull(.data[[y_var]]) %>%
    unique()
  
  # Pathway가 없으면 NULL 반환
  if (length(top_pathways) == 0) {
    message("조건에 맞는 유의한 Pathway가 없습니다.")
    return(NULL)
  }
  
  # [Step 3] 데이터 필터링 & 별표 & 로그 변환
  plot_data <- deduplicated_df %>%
    dplyr::filter(.data[[y_var]] %in% top_pathways) %>%
    dplyr::mutate(stars = dplyr::case_when(
      .data[[p_var]] < 0.001 ~ "***",
      .data[[p_var]] < 0.01  ~ "**",
      .data[[p_var]] < 0.05  ~ "*",
      TRUE ~ ""
    )) %>%
    dplyr::mutate(log_p = -log10(.data[[p_var]])) %>%
    dplyr::mutate(log_p = ifelse(is.infinite(log_p), 30, log_p))

  # 데이터가 비어있으면 NULL 반환
  if (nrow(plot_data) == 0) {
    message("플롯할 데이터가 없습니다.")
    return(NULL)
  }

  # 텍스트 줄바꿈
  if (!is.null(wrap_width)) {
    plot_data <- plot_data %>%
      dplyr::mutate(!!rlang::sym(y_var) := stringr::str_wrap(.data[[y_var]], width = wrap_width))
  }
  
  # [Step 4] Y축 정렬
  mat <- plot_data %>% 
    dplyr::select(.data[[y_var]], .data[[x_var]], .data[[color_var]]) %>%
    tidyr::pivot_wider(names_from = .data[[x_var]], values_from = .data[[color_var]], values_fill = 0) %>%
    as.data.frame()
  
  # mat가 올바르게 생성되었는지 확인
  if (is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) {
    message("행렬 생성에 실패했습니다.")
    return(NULL)
  }
  
  # y_var 컬럼이 있는지 확인
  if (!y_var %in% colnames(mat)) {
    message("y_var 컬럼이 행렬에 없습니다.")
    return(NULL)
  }
  
  rownames(mat) <- mat[[y_var]]
  mat <- mat[, -1, drop = FALSE]
  
  # nrow, ncol가 유효한 숫자인지 확인하고 클러스터링 수행
  n_rows <- nrow(mat)
  n_cols <- ncol(mat)
  if (!is.na(n_rows) && !is.na(n_cols) && n_rows > 1 && n_cols > 1) {
      dd <- stats::dist(mat)
      hc <- stats::hclust(dd)
      pathway_order <- rownames(mat)[hc$order]
      plot_data[[y_var]] <- factor(plot_data[[y_var]], levels = pathway_order)
  } else {
    # 클러스터링 불가능한 경우 기본 순서 사용
    plot_data[[y_var]] <- factor(plot_data[[y_var]], levels = unique(plot_data[[y_var]]))
  }
  
  # [Step 5] 그리기
  ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[[x_var]], y = .data[[y_var]])) +
    ggplot2::geom_tile(fill = NA, color = "grey95") + 
    ggplot2::geom_point(ggplot2::aes(size = log_p, color = .data[[color_var]])) +
    ggplot2::geom_text(ggplot2::aes(label = stars), color = "black", vjust = 0.75, size = 3, fontface = "bold") +
    
    ggplot2::scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, 
                          name = color_label) +
    ggplot2::scale_size_continuous(range = c(2, 8), name = "-log10(P-adj)") + 
    
    ggplot2::labs(x = "Cluster", y = NULL, title = "Cross-Cluster Dot Plot",
         subtitle = paste0("Filtered by presence in ", 
                           ifelse(is.null(n_clusters), "any cluster", 
                                  paste(paste(n_clusters, collapse=","), "clusters"))),
         caption = "* P < 0.05, ** P < 0.01, *** P < 0.001") +
    
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 10, lineheight = 0.8),
      panel.grid.major = ggplot2::element_line(color = "grey90", linetype = "dashed")
    )
}
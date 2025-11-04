#' @export
mybar <- function(data, column = NULL, bins = NULL, x_unit = NULL, y_unit = NULL, 
                  xlab = NULL, main = NULL) {
  
  # vector인 경우 처리
  if(is.numeric(data) || is.vector(data)) {
    df <- data.frame(value = data)
    column <- sym("value")
    col_name <- deparse(substitute(data))  # 변수명 추출
  } else {
    # dataframe인 경우
    col_quo <- enquo(column)
    col_name <- as_label(col_quo)
    df <- data
    column <- col_quo
  }
  
  # x, y 레이블에 단위 추가
  if(is.null(xlab)) {
    x_lab <- if(!is.null(x_unit)) sprintf("%s (%s)", col_name, x_unit) else col_name
  } else {
    x_lab <- xlab
  }
  
  y_lab <- if(!is.null(y_unit)) sprintf("Frequency (%s)", y_unit) else "Frequency"
  
  # 제목 설정
  if(is.null(main)) {
    main <- paste("Histogram of", col_name)
  }
  
  p <- ggplot(df, aes(x = !!column)) +
    # bins 인자에 따라 히스토그램 그리기
    { if(!is.null(bins)) geom_histogram(bins = bins, color = "black", fill = "skyblue") 
      else geom_histogram(color = "black", fill = "skyblue") } +
    labs(x = x_lab, y = y_lab, title = main) +
    theme_minimal()
  
  print(p)
}

#' 단일 벡터 line plot with cutoff intersection
#' @export
myline <- function(vector, main = "Sorted Values", ylab = "Value", 
                   cutoff = NULL, cutoff_col = "red", show_intersection = TRUE, ...) {
  sorted_vec <- sort(vector, decreasing = TRUE)
  x_vals <- 1:length(sorted_vec)
  
  plot(x_vals, sorted_vec, type = "l", 
       main = main, xlab = "Rank", ylab = ylab, ...)
  
  if(!is.null(cutoff)) {
    abline(h = cutoff, col = cutoff_col, lty = 2, lwd = 2)
    
    # cutoff 교점 찾기 및 표시
    if(show_intersection) {
      # 교점 찾기 (선형 보간)
      idx_above <- which(sorted_vec >= cutoff)
      idx_below <- which(sorted_vec < cutoff)
      
      if(length(idx_above) > 0 && length(idx_below) > 0) {
        idx_cross <- max(idx_above)
        if(idx_cross < length(sorted_vec)) {
          # 선형 보간으로 정확한 x 위치 계산
          y1 <- sorted_vec[idx_cross]
          y2 <- sorted_vec[idx_cross + 1]
          x_intersect <- idx_cross + (cutoff - y1) / (y2 - y1)
          
          # 교점 표시
          points(x_intersect, cutoff, col = cutoff_col, pch = 19, cex = 1.5)
          text(x_intersect, cutoff, labels = sprintf("x=%.1f", x_intersect), 
               pos = 3, col = cutoff_col, cex = 0.8)
        }
      }
    }
  }
}

#' 다중 벡터 line plot with dual y-axis support
#' @export
mylines <- function(..., vectors_right = NULL, 
                    main = "Sorted Values Comparison", 
                    ylab = "Value", ylab_right = NULL,
                    cutoff = NULL, cutoff_right = NULL,
                    cutoff_col = "red", cutoff_col_right = "blue",
                    legend_pos = "topright", show_intersection = TRUE) {
  
  # 왼쪽 y축 벡터들
  vectors_left <- list(...)
  n_left <- length(vectors_left)
  
  # 오른쪽 y축 벡터들 처리
  if(!is.null(vectors_right)) {
    if(!is.list(vectors_right)) {
      vectors_right <- list(vectors_right)
    }
  }
  n_right <- length(vectors_right)
  
  # 모든 벡터 정렬
  sorted_left <- lapply(vectors_left, function(x) sort(x, decreasing = TRUE))
  sorted_right <- if(!is.null(vectors_right)) {
    lapply(vectors_right, function(x) sort(x, decreasing = TRUE))
  } else NULL
  
  # x축 스케일 조정
  all_vectors <- c(sorted_left, sorted_right)
  max_len <- max(sapply(all_vectors, length))
  
  x_coords <- lapply(all_vectors, function(vec) {
    len <- length(vec)
    seq(1, max_len, length.out = len)
  })
  
  # 왼쪽 y축 설정 및 플롯
  par(mar = c(5, 4, 4, 4))  # 오른쪽 여백 확보
  
  y_range_left <- range(unlist(sorted_left), na.rm = TRUE)
  
  # 첫 번째 플롯 (왼쪽 y축)
  plot(x_coords[[1]], sorted_left[[1]], type = "l", 
       xlim = c(1, max_len), ylim = y_range_left,
       main = main, xlab = "Normalized Rank", ylab = ylab, 
       col = 1, lwd = 2)
  
  # 추가 왼쪽 y축 라인들
  if(n_left > 1) {
    for(i in 2:n_left) {
      lines(x_coords[[i]], sorted_left[[i]], col = i, lwd = 2)
    }
  }
  
  # 왼쪽 y축 cutoff
  if(!is.null(cutoff)) {
    abline(h = cutoff, col = cutoff_col, lty = 2, lwd = 2)
    
    # 교점 찾기 및 표시
    if(show_intersection) {
      for(i in 1:n_left) {
        find_and_mark_intersection(x_coords[[i]], sorted_left[[i]], 
                                   cutoff, cutoff_col, paste0("L", i))
      }
    }
  }
  
  # 오른쪽 y축 플롯
  if(!is.null(vectors_right)) {
    par(new = TRUE)
    
    y_range_right <- range(unlist(sorted_right), na.rm = TRUE)
    
    # 오른쪽 y축 색상은 n_left + 1부터 시작
    plot(x_coords[[n_left + 1]], sorted_right[[1]], type = "l",
         xlim = c(1, max_len), ylim = y_range_right,
         axes = FALSE, xlab = "", ylab = "",
         col = n_left + 1, lwd = 2, lty = 2)  # 점선으로 구분
    
    # 추가 오른쪽 y축 라인들
    if(n_right > 1) {
      for(i in 2:n_right) {
        lines(x_coords[[n_left + i]], sorted_right[[i]], 
              col = n_left + i, lwd = 2, lty = 2)
      }
    }
    
    # 오른쪽 y축 그리기
    axis(4)
    if(!is.null(ylab_right)) {
      mtext(ylab_right, side = 4, line = 2.5)
    }
    
    # 오른쪽 y축 cutoff
    if(!is.null(cutoff_right)) {
      abline(h = cutoff_right, col = cutoff_col_right, lty = 3, lwd = 2)
      
      if(show_intersection) {
        for(i in 1:n_right) {
          find_and_mark_intersection(x_coords[[n_left + i]], sorted_right[[i]], 
                                     cutoff_right, cutoff_col_right, paste0("R", i))
        }
      }
    }
  }
  
  # 범례
  all_labels <- c(paste("Left", 1:n_left))
  all_cols <- 1:n_left
  all_lty <- rep(1, n_left)
  
  if(!is.null(vectors_right)) {
    all_labels <- c(all_labels, paste("Right", 1:n_right))
    all_cols <- c(all_cols, (n_left + 1):(n_left + n_right))
    all_lty <- c(all_lty, rep(2, n_right))
  }
  
  legend(legend_pos, legend = all_labels, 
         col = all_cols, lty = all_lty, lwd = 2)
}

# 교점 찾기 및 표시 헬퍼 함수; myline, mylines에 사용됨.
find_and_mark_intersection <- function(x_vals, y_vals, cutoff, col, label_prefix) {
  idx_above <- which(y_vals >= cutoff)
  idx_below <- which(y_vals < cutoff)
  
  if(length(idx_above) > 0 && length(idx_below) > 0) {
    # 마지막으로 cutoff 위에 있는 인덱스
    idx_last_above <- max(which(y_vals >= cutoff))
    
    if(idx_last_above < length(y_vals)) {
      # 선형 보간
      x1 <- x_vals[idx_last_above]
      x2 <- x_vals[idx_last_above + 1]
      y1 <- y_vals[idx_last_above]
      y2 <- y_vals[idx_last_above + 1]
      
      x_intersect <- x1 + (x2 - x1) * (cutoff - y1) / (y2 - y1)
      
      points(x_intersect, cutoff, col = col, pch = 19, cex = 1.2)
      text(x_intersect, cutoff, labels = sprintf("%s:%.1f", label_prefix, x_intersect), 
           pos = 1, col = col, cex = 0.7)
    }
  }
}

#' @export
mydensity <- function(data, column = NULL, adjust = 1, x_unit = NULL, y_unit = NULL,
                      xlab = NULL, main = NULL) {
  
  # vector인 경우 처리
  if(is.numeric(data) || is.vector(data)) {
    df <- data.frame(value = data)
    column <- sym("value")
    col_name <- deparse(substitute(data))  # 변수명 추출
  } else {
    # dataframe인 경우
    col_quo <- enquo(column)
    col_name <- as_label(col_quo)
    df <- data
    column <- col_quo
  }
  
  # x, y 레이블
  if(is.null(xlab)) {
    x_lab <- if(!is.null(x_unit)) sprintf("%s (%s)", col_name, x_unit) else col_name
  } else {
    x_lab <- xlab
  }
  
  y_lab <- if(!is.null(y_unit)) sprintf("Density (%s)", y_unit) else "Density"
  
  # 제목 설정
  if(is.null(main)) {
    main <- paste("Density Plot of", col_name)
  }
  
  p <- ggplot(df, aes(x = !!column)) +
    geom_density(adjust = adjust, color = "white", fill = "tomato", alpha = 0.7) +
    labs(x = x_lab, y = y_lab, title = main) +
    theme_minimal()
  
  print(p)
}

# Boxplot ----------

.fetch_plot_data <- function(data, features, sample_col = NULL, group.by = NULL, split.by = NULL, 
                             assay = NULL, layer = "data") {
  
  `%||%` <- rlang::`%||%`
  
  # 1. Seurat 객체 처리
  if (inherits(data, "Seurat")) {
    sobj <- data
    meta_cols <- names(sobj@meta.data)
    
    meta_features <- intersect(features, meta_cols)
    assay_features <- setdiff(features, meta_cols)
    
    # 필요한 모든 컬럼 정의
    grouping_cols <- unique(c(sample_col, group.by, split.by))
    grouping_cols <- grouping_cols[!is.null(grouping_cols)]
    
    vars_to_fetch <- unique(c(grouping_cols, meta_features, assay_features))
    vars_to_fetch <- vars_to_fetch[vars_to_fetch %in% c(meta_cols, rownames(sobj[[assay %||% DefaultAssay(sobj)]]))] # 존재하는 것만
    
    # 데이터 추출
    fetched_df <- data.frame(matrix(ncol = 0, nrow = ncol(sobj)))
    rownames(fetched_df) <- colnames(sobj)
    
    if (length(vars_to_fetch) > 0) {
        if (length(assay_features) > 0) {
          assay_to_use <- assay %||% DefaultAssay(sobj)
          fetched_df <- FetchData(sobj, vars = vars_to_fetch, assay = assay_to_use, layer = layer)
        } else {
          fetched_df <- FetchData(sobj, vars = vars_to_fetch)
        }
    } else {
       warning("가져올 유효한 피처나 메타데이터가 없습니다.")
       return(data.frame())
    }

  # 2. 데이터 프레임 처리
  } else if (is.data.frame(data)) {
    fetched_df <- data
    req_cols <- unique(c(features, sample_col, group.by, split.by))
    req_cols <- req_cols[!is.null(req_cols)]
    
    missing_cols <- setdiff(req_cols, names(fetched_df))
    if (length(missing_cols) > 0) {
      stop("데이터 프레임에 다음 컬럼이 없습니다: ", paste(missing_cols, collapse=", "))
    }
  } else {
    stop("`data` 인수는 Seurat 객체 또는 data.frame이어야 합니다.")
  }
  
  return(fetched_df)
}

#' 세포 단위 피처 발현 박스플롯 (plot_gene_boxplot 스타일)
#'
#' @param data Seurat 객체 또는 data.frame
#' @param features 플로팅할 피처 (유전자 또는 메타데이터)
#' @param group.by 패널(facet)로 분리할 그룹 변수 (예: "cell_type")
#' @param split.by X축에서 나눌 그룹 변수 (예: "timepoint")
#' @param idents 'group.by'에서 포함할 특정 값들
#' @param assay Seurat 객체용 assay
#' @param layer Seurat 객체용 layer
#' @param ncol 플롯 컬럼 수
#' @param pt.size 점 크기 (0이면 숨김)
#' @param violin 바이올린 플롯 오버레이
#' @param add_stats 통계 검정 추가 (ggpubr::stat_compare_means)
#' @param ... stat_compare_means에 전달할 추가 인수 (예: paired = TRUE, method = "t.test")
#'
#' @return ggplot 객체
#' @export
mybox_cell <- function(data, 
                                      features, 
                                      group.by, 
                                      split.by, 
                                      idents = NULL, 
                                      assay = NULL, 
                                      layer = "data",
                                      ncol = 3, 
                                      pt.size = 0.2, 
                                      violin = FALSE,
                                      add_stats = TRUE,
                                      ...) {

  plot_df_base <- .fetch_plot_data(data, features, NULL, group.by, split.by, assay, layer)
  
  all_plots <- list()
  
  for (feature_name in features) {
    current_df <- plot_df_base
    
    # 1. group.by가 NULL일 경우 "Overall"로 처리
    group.by.internal <- group.by
    if (is.null(group.by)) {
      group.by.internal <- ".internal_placeholder_group"
      current_df[[group.by.internal]] <- "Overall"
    }
    
    # 2. idents 필터링 (mybox 로직과 동일)
    if (!is.null(group.by) && !is.null(idents)) {
      original_col_type <- class(current_df[[group.by.internal]])
      current_df[[group.by.internal]] <- as.character(current_df[[group.by.internal]])
      current_df <- current_df[current_df[[group.by.internal]] %in% idents, ]
      
      if (nrow(current_df) == 0) {
         warning(feature_name, "에서 idents 필터링 후 데이터가 없습니다. 건너뜁니다.")
         next
      }
      
      if ("factor" %in% original_col_type) {
        original_levels <- levels(if(inherits(data, "Seurat")) data@meta.data[[group.by]] else data[[group.by]])
        valid_subset_levels <- intersect(original_levels, idents)
        current_df[[group.by.internal]] <- factor(current_df[[group.by.internal]], levels = valid_subset_levels)
      } else {
        current_df[[group.by.internal]] <- as.factor(current_df[[group.by.internal]])
      }
    }
    
    current_df[[group.by.internal]] <- as.factor(current_df[[group.by.internal]])
    current_df[[split.by]] <- as.factor(current_df[[split.by]])
    
    # 3. ggplot 생성 (Req 3 레이아웃)
    p <- ggplot(current_df, aes(x = .data[[split.by]], y = .data[[feature_name]], fill = .data[[split.by]]))
    
    if (violin) {
      p <- p + geom_violin(alpha = 0.5, scale = "width", trim = FALSE, na.rm = TRUE)
    }
    
    # outlier.shape=NA로 설정 (Req 2)
    p <- p + geom_boxplot(na.rm = TRUE, outlier.shape = NA, alpha = 0.7)
    
    # pt.size > 0 일 때만 jitter 추가 (Req 2: 색상 "black" 고정)
    if (pt.size > 0) {
      p <- p + geom_jitter(color = "black", size = pt.size, alpha = 0.3, 
                           width = 0.2, height = 0, na.rm = TRUE)
    }
    
    # 통계 추가 (ggpubr)
    if (add_stats) {
      # `...` 인수로 paired=TRUE 등을 받을 수 있음
      p <- p + stat_compare_means(label = "p.format", ...) 
    }
    
    # 패널(facet) 적용 및 테마 설정 (Req 3)
    p <- p + facet_wrap(as.formula(paste("~", group.by.internal))) +
      labs(title = feature_name,
           x = split.by,
           y = feature_name) +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5),
            legend.position = "none") # 범례 제거
    
    all_plots[[feature_name]] <- p
  }
  
  if (length(all_plots) == 0) {
    message("생성된 플롯이 없습니다.")
    return(invisible(NULL))
  }
  
  return(patchwork::wrap_plots(all_plots, ncol = ncol))
}
#' 샘플 단위(Pseudobulk) 피처 평균 박스플롯 (mybox 로직 + plot_gene_boxplot 스타일)
#'
#' @param data Seurat 객체 또는 data.frame
#' @param features 플로팅할 피처 (유전자 또는 메타데이터)
#' @param sample_col 샘플/환자 식별 컬럼 (필수)
#' @param group.by 패널(facet)로 분리할 그룹 변수 (예: "cell_type")
#' @param split.by X축에서 나눌 그룹 변수 (예: "timepoint")
#' @param idents 'group.by'에서 포함할 특정 값들
#' @param assay Seurat 객체용 assay
#' @param layer Seurat 객체용 layer
#' @param ncol 플롯 컬럼 수
#' @param pt.size 샘플 평균 점 크기 (0이면 숨김)
#' @param violin 바이올린 플롯 오버레이
#' @param add_stats 통계 검정 추가 (요약된 데이터 대상)
#' @param ... stat_compare_means에 전달할 추가 인수 (예: paired = TRUE)
#'
#' @return ggplot 객체
#' @export
mybox_pb <- function(data, 
                                    features, 
                                    sample_col, 
                                    group.by, 
                                    split.by, 
                                    idents = NULL, 
                                    assay = NULL, 
                                    layer = "data",
                                    ncol = 3, 
                                    pt.size = 1.0, 
                                    violin = FALSE,
                                    add_stats = FALSE,
                                    ...) {

  if (is.null(sample_col)) {
    stop("`sample_col` 인수는 Pseudobulk 플롯에 필수입니다.")
  }

  plot_df_base <- .fetch_plot_data(data, features, sample_col, group.by, split.by, assay, layer)
  
  all_plots <- list()
  
  for (feature_name in features) {
    current_df <- plot_df_base
    
    # 1. group.by가 NULL일 경우 "Overall"로 처리
    group.by.internal <- group.by
    if (is.null(group.by)) {
      group.by.internal <- ".internal_placeholder_group"
      current_df[[group.by.internal]] <- "Overall"
    }
    
    # 2. idents 필터링
    if (!is.null(group.by) && !is.null(idents)) {
      original_col_type <- class(current_df[[group.by.internal]])
      current_df[[group.by.internal]] <- as.character(current_df[[group.by.internal]])
      current_df <- current_df[current_df[[group.by.internal]] %in% idents, ]
      
      if (nrow(current_df) == 0) {
         warning(feature_name, "에서 idents 필터링 후 데이터가 없습니다. 건너뜁니다.")
         next
      }
      
      if ("factor" %in% original_col_type) {
        original_levels <- levels(if(inherits(data, "Seurat")) data@meta.data[[group.by]] else data[[group.by]])
        valid_subset_levels <- intersect(original_levels, idents)
        current_df[[group.by.internal]] <- factor(current_df[[group.by.internal]], levels = valid_subset_levels)
      } else {
        current_df[[group.by.internal]] <- as.factor(current_df[[group.by.internal]])
      }
    }
    
    current_df[[group.by.internal]] <- as.factor(current_df[[group.by.internal]])
    current_df[[split.by]] <- as.factor(current_df[[split.by]])
    current_df[[sample_col]] <- as.factor(current_df[[sample_col]])
    
    # 3. *** 데이터 요약 (Aggregation) ***
    grouping_vars <- unique(c(sample_col, group.by.internal, split.by))
    
    aggregated_df <- current_df %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(grouping_vars))) %>%
      dplyr::summarise(mean_value = mean(.data[[feature_name]], na.rm = TRUE), .groups = 'drop')
    
    if (nrow(aggregated_df) == 0) {
      warning(feature_name, "에서 요약 후 데이터가 없습니다. 건너뜁니다.")
      next
    }
    
    # 4. ggplot 생성 (Req 3 레이아웃)
    p <- ggplot(aggregated_df, aes(x = .data[[split.by]], y = mean_value, fill = .data[[split.by]]))
    
    if (violin) {
      p <- p + geom_violin(alpha = 0.5, scale = "width", trim = FALSE, na.rm = TRUE)
    }
    
    # outlier.shape=NA로 설정 (Req 2)
    p <- p + geom_boxplot(na.rm = TRUE, outlier.shape = NA, alpha = 0.7)
    
    # pt.size > 0 일 때만 jitter 추가 (Req 2: 색상 "black" 고정)
    if (pt.size > 0) {
      # PB 플롯은 점이 샘플을 의미하므로 더 진하게 (alpha=0.7)
      p <- p + geom_jitter(color = "black", size = pt.size, alpha = 0.7, 
                           width = 0.2, height = 0, na.rm = TRUE)
    }

    # 통계 추가 (요약된 데이터를 대상으로 함)
    if (add_stats) {
      # paired=TRUE를 `...`로 전달받아 샘플 페어 테스트 가능
      p <- p + stat_compare_means(label = "p.format", ...) 
    }

    # 패널(facet) 적용 및 테마 설정 (Req 3)
    p <- p + facet_wrap(as.formula(paste("~", group.by.internal))) +
      labs(title = feature_name,
           x = split.by,
           y = paste("Average", feature_name)) +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5),
            legend.position = "none") # 범례 제거
            
    all_plots[[feature_name]] <- p
  }
  
  if (length(all_plots) == 0) {
    message("생성된 플롯이 없습니다.")
    return(invisible(NULL))
  }
  
  return(patchwork::wrap_plots(all_plots, ncol = ncol))
}

mybox_geomx <- function(
    data_norm,
    metadata,
    features,
    group.by,
    split.by = NULL,
    test_method = "wilcox.test",
    p_adjust_method = "holm",
    hide_ns = FALSE,
    show_points = TRUE,
    comparisons = "all"  # "all" 또는 특정 비교 리스트
) {
  
  # 1. 라이브러리 로드
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggpubr)
  library(rstatix)
  
  # 2. 입력 데이터 유효성 검사
  if (!group.by %in% colnames(metadata)) {
    stop(paste("Error: `group.by` 변수인 '", group.by, "'가 metadata에 존재하지 않습니다.", sep=""))
  }
  if (!is.null(split.by) && !split.by %in% colnames(metadata)) {
    stop(paste("Error: `split.by` 변수인 '", split.by, "'가 metadata에 존재하지 않습니다.", sep=""))
  }
  if (!test_method %in% c("t.test", "wilcox.test")) {
    stop("Error: `test_method`는 't.test' 또는 'wilcox.test' 중 하나여야 합니다.")
  }
  
  # 3. 데이터 준비
  common_aois <- intersect(colnames(data_norm), metadata$SegmentDisplayName)
  if(length(common_aois) == 0){
    stop("Error: 'data_norm'의 컬럼명과 'metadata'의 SegmentDisplayName이 일치하는 AOI가 없습니다.")
  }
  
  metadata_sub <- metadata %>% filter(SegmentDisplayName %in% common_aois)
  data_norm_sub <- data_norm[features, common_aois, drop = FALSE]
  
  plot_data <- data_norm_sub %>%
    t() %>%
    as.data.frame() %>%
    mutate(SegmentDisplayName = rownames(.)) %>%
    pivot_longer(
      cols = all_of(features),
      names_to = "feature",
      values_to = "expression"
    ) %>%
    left_join(metadata_sub, by = "SegmentDisplayName")
  
  # 결과 저장용 리스트
  stat_results <- list()
  
  # 4. 각 Feature에 대해 Boxplot 생성
  for (feat in features) {
    
    feature_data <- plot_data %>% filter(feature == feat)
    
    # Factor 변환
    if(!is.factor(feature_data[[group.by]])) {
      feature_data[[group.by]] <- as.factor(feature_data[[group.by]])
    }
    if(!is.null(split.by) && !is.factor(feature_data[[split.by]])) {
      feature_data[[split.by]] <- as.factor(feature_data[[split.by]])
    }
    
    # 기본 플롯 설정
    p <- ggplot(feature_data, aes(x = .data[[group.by]], y = expression))
    
    # Case 1: split.by가 지정된 경우
    if (!is.null(split.by)) {
      dodge_width <- 0.85
      
      p <- p +
        geom_boxplot(
          aes(fill = .data[[split.by]]), 
          position = position_dodge(width = dodge_width), 
          outlier.shape = NA, 
          alpha = 0.8
        ) +
        {if(show_points) geom_jitter(
          aes(color = .data[[split.by]]), 
          position = position_jitterdodge(dodge.width = dodge_width, jitter.width = 0.2), 
          size = 1, 
          alpha = 0.5
        )} +
        scale_color_discrete(name = split.by) +
        labs(fill = split.by)
      
      # 각 group.by 수준 내에서 split.by 수준들 간 비교
      stat_formula <- as.formula(paste("expression ~", paste0("`", split.by, "`")))
      
      base_stat_test <- if (test_method == "t.test") {
        feature_data %>% 
          group_by(.data[[group.by]]) %>% 
          t_test(stat_formula)
      } else {
        feature_data %>% 
          group_by(.data[[group.by]]) %>% 
          wilcox_test(stat_formula)
      }
      
      if (nrow(base_stat_test) > 0) {
        # p-value 조정 및 위치 계산
        stat_test <- base_stat_test %>%
          adjust_pvalue(method = p_adjust_method) %>%
          add_xy_position(x = group.by, dodge = dodge_width, fun = "max", scales = "free_y")
        
        # 통계 결과 저장
        stat_results[[feat]] <- stat_test
        
        # 플롯에 p-value 추가
        p <- p + stat_pvalue_manual(
          stat_test, 
          label = "p.adj", 
          hide.ns = hide_ns, 
          tip.length = 0.01
        )
      }
      
      # Case 2: split.by가 없는 경우
    } else {
      p <- p +
        geom_boxplot(
          aes(fill = .data[[group.by]]), 
          outlier.shape = NA, 
          alpha = 0.8
        ) +
        {if(show_points) geom_jitter(width = 0.2, size = 1, alpha = 0.6)} +
        theme(legend.position = "none")
      
      # group.by의 수준이 2개 이상일 때만 통계 검정 수행
      n_groups <- length(unique(feature_data[[group.by]]))
      
      if (n_groups >= 2) {
        stat_formula <- as.formula(paste("expression ~", paste0("`", group.by, "`")))
        
        # 그룹이 2개면 단순 비교, 3개 이상이면 pairwise
        if (n_groups == 2) {
          base_stat_test <- if (test_method == "t.test") {
            feature_data %>% t_test(stat_formula)
          } else {
            feature_data %>% wilcox_test(stat_formula)
          }
        } else {
          # 3개 이상 그룹: pairwise comparison
          if (comparisons == "all") {
            base_stat_test <- if (test_method == "t.test") {
              feature_data %>% 
                t_test(stat_formula, p.adjust.method = "none")
            } else {
              feature_data %>% 
                wilcox_test(stat_formula, p.adjust.method = "none")
            }
          } else if (is.list(comparisons)) {
            base_stat_test <- if (test_method == "t.test") {
              feature_data %>% 
                t_test(stat_formula, comparisons = comparisons, p.adjust.method = "none")
            } else {
              feature_data %>% 
                wilcox_test(stat_formula, comparisons = comparisons, p.adjust.method = "none")
            }
          }
        }
        
        if (exists("base_stat_test") && nrow(base_stat_test) > 0) {
          # p-value 조정 및 위치 계산
          stat_test <- base_stat_test %>%
            adjust_pvalue(method = p_adjust_method) %>%
            add_y_position(scales = "free_y", step.increase = 0.1)
          
          # 통계 결과 저장
          stat_results[[feat]] <- stat_test
          
          # 플롯에 p-value 추가
          p <- p + stat_pvalue_manual(
            stat_test, 
            label = "p.adj", 
            hide.ns = hide_ns, 
            tip.length = 0.01
          )
        }
      }
    }
    
    # 공통 테마 및 레이블 적용
    p <- p +
      theme_classic() +
      labs(
        title = paste("Expression of", feat),
        x = group.by,
        y = "Normalized Expression"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11)
      )
    
    # 플롯 출력
    print(p)
  }
  
  # 통계 결과 반환 (invisible로 콘솔에 자동 출력되지 않도록)
  invisible(stat_results)
}


# boxplot end

#' UpSet Plot for Multiple Gene Lists
#'
#' @description
#' 주어진 여러 유전자 벡터(list)를 이용해 UpSet 플롯을 그립니다.
#'
#' @param gene_lists Named list. 각 요소가 character vector인 gene list.
#' @param sets Optional character vector. 표시할 세트 순서를 names(gene_lists)가 아닌 다른 순서로 지정할 때 사용.
#'             기본은 names(gene_lists) 순서.
#' @param min_size Integer. 최소 교집합 크기(유전자 수)가 이 값 이상인 경우에만 표시 (default: 1).
#' @param width_ratio Numeric. 메인 플롯과 바 플롯 너비 비율 (default: 0.3).
#' @param keep_sets Character vector. 특정 세트만 플롯에 포함시키고 싶을 때 사용. 기본 NULL(모두 포함).
#' @return ggplot 객체.  
#' @examples
#' \dontrun{
#' L <- list(
#'   A = c("TP53","EGFR","MYC","BRCA1"),
#'   B = c("EGFR","MYC","PTEN","BRCA2"),
#'   C = c("MYC","PTEN","ALK")
#' )
#' p <- upset_gene_lists(L)
#' print(p)
#' }
#' @export
upset_gene_lists <- function(gene_lists,
                             sets       = NULL,
                             min_size   = 1,
                             width_ratio= 0.3,
                             keep_sets  = NULL) {
  # 준비: 데이터프레임으로 변환
  # 모든 유전자의 유니크 목록
  all_genes <- unique(unlist(gene_lists))
  df <- data.frame(Gene = all_genes, stringsAsFactors = FALSE)
  
  # 각 세트별 presence/absence 컬럼 생성
  for (nm in names(gene_lists)) {
    df[[nm]] <- all_genes %in% gene_lists[[nm]]
  }
  
  # 필요한 세트만 선택
  if (!is.null(keep_sets)) {
    stopifnot(all(keep_sets %in% names(gene_lists)))
    df <- df[, c("Gene", keep_sets), drop = FALSE]
  }
  
  # sets 순서 지정
  if (is.null(sets)) sets <- setdiff(names(df), "Gene")
  stopifnot(all(sets %in% names(df)))
  
  # UpSet 플롯
  if (requireNamespace("ComplexUpset", quietly = TRUE)) {
    library(ComplexUpset)
    p <- ComplexUpset::upset(
      df,
      intersect = sets,
      name      = "Gene",
      min_size  = min_size,
      width_ratio = width_ratio
    )
  } else {
    warning("ComplexUpset가 없습니다. 설치해주세요.")
    return(NULL)
  }
  
  return(p)
}

vln_p=function(sobj, feature, group.by, split.by, pt.size=0,ncol=4, ...){
  plist <- VlnPlot(
    object    = sobj,
    features  = feature,
    group.by  = group.by,
    split.by  = split.by,
    combine   = FALSE,
    pt.size=pt.size,
    ...
  )
  
  plist2 <- lapply(plist, function(p){
    df <- p$data
    
    p +
      ggpubr::stat_compare_means(
        data       = df,
        mapping    = aes_string(group = "split"),
        method     = "wilcox.test",
        label      = "p.signif",
        hide.ns    = TRUE,
        tip.length = 0.02
      ) +
      theme(legend.position = "none")
  })
  
  return(patchwork::wrap_plots(plist2, ncol = ncol))
}


#' Create a Proportional Bar Graph of Clusters
#'
#' This function creates a stacked bar plot showing the proportional distribution of clusters
#' across different samples or groups.
#'
#' @param sobj A Seurat object
#' @param identity Character string specifying the identity to use for clustering.
#'                Default is "seurat_clusters"
#' @param group.by Character string specifying the metadata column to group by.
#'                 Default is "sample"
#' @param idents Character vector specifying which identities to include.
#'              If NULL, includes all identities
#' @param df Logical. If TRUE, returns the data frame used for plotting instead of the plot.
#'          Default is FALSE
#' @param vlines Numeric vector specifying x-axis positions for vertical lines.
#'              Default is NULL
#' @param vline_color Character string specifying the color for vertical lines.
#'                   Default is "red"
#'
#' @return A ggplot object showing the proportional distribution of clusters, or a data frame
#'         if df = TRUE
#'
#' @examples
#' \dontrun{
#' # Create proportional bar plot for all clusters
#' p <- cmb(sobj, identity = "seurat_clusters", group.by = "sample")
#' 
#' # Create plot for specific clusters with vertical lines
#' p <- cmb(sobj, idents = c("0", "1", "2"), vlines = c(3, 6))
#' 
#' # Get the underlying data frame
#' df <- cmb(sobj, df = TRUE)
#' }
#' @export
cmb <- function(sobj, identity = "seurat_clusters", group.by = "sample", idents = NULL, df=F, vlines=NULL, vline_color = "red") {
  Idents(sobj) <- identity
  cluster_ids <- Idents(sobj)
  sample_ids <- sobj@meta.data[[group.by]]
  
  data <- data.frame(cluster = cluster_ids, sample = sample_ids)
  
  # Filter identities if idents is provided
  if(!is.null(idents)) {
    data <- data[data$cluster %in% idents,]
  }
  
  summary_data <- data %>%
    group_by(sample, cluster) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(sample) %>%
    mutate(proportion = count / sum(count)) %>%
    ungroup()
  
  # Sort the samples using the new sort_samples function
  sorted_samples <- sort_samples(unique(as.character(summary_data$sample)))
  summary_data$sample <- factor(summary_data$sample, levels = sorted_samples)
  
  output=ggplot(summary_data, aes(x = sample, y = proportion, fill = cluster)) +
    geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.2) +
    labs(title = paste0("Proportional Bar Graph of Clusters for Each ", group.by),
         x = group.by,
         y = "Proportion of Cells") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right") +
    scale_fill_viridis_d() +
    scale_y_continuous(labels = scales::percent_format())
  
  if(!is.null(vlines)){
    for (line_pos in vlines) {
      output <- output + geom_vline(xintercept = line_pos + 0.5, color = vline_color, linetype = "dashed", linewidth = 0.5)
    }
  }
  
  if(df){
    return(summary_data)
  }
  return(output)
}


#' Create an Absolute Count Bar Graph of Clusters
#'
#' This function creates a stacked bar plot showing the absolute count distribution of clusters
#' across different samples or groups.
#'
#' @param sobj A Seurat object
#' @param identity Character string specifying the identity to use for clustering.
#'                Default is "seurat_clusters"
#' @param group.by Character string specifying the metadata column to group by.
#'                 Default is "sample"
#' @param idents Character vector specifying which identities to include.
#'              If NULL, includes all identities
#' @param df Logical. If TRUE, returns the data frame used for plotting instead of the plot.
#'          Default is FALSE
#' @param vlines Numeric vector specifying x-axis positions for vertical lines.
#'              Default is NULL
#' @param vline_color Character string specifying the color for vertical lines.
#'                   Default is "red"
#'
#' @return A ggplot object showing the absolute count distribution of clusters, or a data frame
#'         if df = TRUE
#'
#' @examples
#' \dontrun{
#' # Create absolute count bar plot for all clusters
#' p <- acmb(sobj, identity = "seurat_clusters", group.by = "sample")
#' 
#' # Create plot for specific clusters with vertical lines
#' p <- acmb(sobj, idents = c("0", "1", "2"), vlines = c(3, 6))
#' 
#' # Get the underlying data frame
#' df <- acmb(sobj, df = TRUE)
#' }
#' @export
acmb <- function(sobj, identity="seurat_clusters", group.by="sample", idents = NULL, df=F, vlines=NULL, vline_color = "red") {
  Idents(sobj) <- identity
  cluster_ids <- Idents(sobj)
  sample_ids <- sobj@meta.data[[group.by]]
  
  data <- data.frame(cluster = cluster_ids, sample = sample_ids)
  
  # Filter identities if idents is provided
  if(!is.null(idents)) {
    data <- data[data$cluster %in% idents,]
  }
  
  summary_data <- data %>%
    group_by(sample, cluster) %>%
    summarise(count = n(), .groups = "drop")
  
  # Sort the samples using the new sort_samples function
  sorted_samples <- sort_samples(unique(summary_data$sample))
  summary_data$sample <- factor(summary_data$sample, levels = sorted_samples)
  
  output=ggplot(summary_data, aes(x = sample, y = count, fill = cluster)) +
    geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.2) +
    labs(title = paste0("Cumulative Bar Graph of Clusters for Each ", group.by),
         x = group.by,
         y = "Number of Cells") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right") +
    scale_fill_viridis_d()
  
  if(!is.null(vlines)){
    for (line_pos in vlines) {
      output <- output + geom_vline(xintercept = line_pos + 0.5, color = vline_color, linetype = "dashed", linewidth = 0.5)
    }
  }
  
  if(df){
    return(summary_data)
  }
  return(output)
}






#' Cumulative Line Graph of Cluster Proportions
#'
#' Produces a cumulative line plot of per‑cluster cell proportions for each
#' sample (or other grouping variable). Clusters can be ordered alphabetically
#'/numerically by name, or by their overall frequency in the dataset.
#'
#' @param sobj            A **Seurat** object containing the study metadata.
#' @param cluster_col     Metadata column with cluster IDs (default
#'                        \code{"seurat_clusters"}).
#' @param group.by        Metadata column identifying samples or other groups
#'                        to compare (default \code{"gem"}).
#' @param sort.by         Cluster ordering strategy: \code{"name"} (default)
#'                        keeps the natural/alphabetical order, whereas
#'                        \code{"frequency"} orders clusters by their global
#'                        frequency (descending).
#' @param df              Logical; return underlying summary data instead of
#'                        the plot (default \code{FALSE}).
#' @param color_palette   Optional vector of colours for groups. If \code{NULL}
#'                        a Set2 palette is generated.
#' @param n_patterns,n_shapes  How many distinct linetypes / shapes to cycle.
#'
#' @return A \code{ggplot} object (or \code{data.frame} when \code{df = TRUE}).
#' @export
#' @import dplyr ggplot2 RColorBrewer
cml <- function(sobj,
                cluster_col = "seurat_clusters",
                group.by    = "gem",
                sort.by     = c("name", "frequency"),
                df          = FALSE,
                color_palette = NULL,
                n_patterns  = 5,
                n_shapes    = 7) {
  sort.by <- match.arg(sort.by)
  stopifnot("Seurat" %in% class(sobj))
  

  # 1. Tabulate counts & proportions

  data_tbl <- sobj@meta.data %>%
    dplyr::select(cluster = !!cluster_col, group = !!group.by) %>%
    dplyr::group_by(group, cluster) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
    mutate(group=as.factor(group))
  

  # 2. Determine cluster order

  cluster_levels <- if (sort.by == "frequency") {
    data_tbl %>%
      dplyr::group_by(cluster) %>%
      dplyr::summarise(n = sum(count)) %>%
      dplyr::arrange(dplyr::desc(n)) %>%
      dplyr::pull(cluster)
  } else {
    sort(unique(as.character(data_tbl$cluster)), method = "auto")
  }
  data_tbl$cluster <- factor(data_tbl$cluster, levels = cluster_levels)
  

  # 3. Compute within‑group proportions and cumulative sums

  cum_tbl <- data_tbl %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(prop = count / sum(count)) %>%
    dplyr::arrange(cluster) %>%
    dplyr::mutate(cum_prop = cumsum(prop)) %>%
    dplyr::ungroup()
  

  # 4. Styling helpers

  if (is.null(color_palette)) {
    n_col <- length(unique(cum_tbl$group))
    color_palette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(n_col)
  }
  line_patterns <- rep(c("solid", "dashed", "dotted", "dotdash", "longdash"),
                       length.out = n_patterns)
  point_shapes  <- rep(c(16, 17, 18, 15, 8, 3, 4), length.out = n_shapes)
  
  style_df <- data.frame(
    group    = unique(cum_tbl$group),
    colour   = rep(color_palette, length.out = length(unique(cum_tbl$group))),
    linetype = rep(line_patterns , length.out = length(unique(cum_tbl$group))),
    shape    = rep(point_shapes  , length.out = length(unique(cum_tbl$group)))
  )
  cum_tbl <- dplyr::left_join(cum_tbl, style_df, by = "group")
  

  # 5. Plot

  p <- ggplot2::ggplot(cum_tbl, ggplot2::aes(x = cluster, y = cum_prop,
                                             colour = group, linetype = group,
                                             shape = group, group  = group)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    ggplot2::scale_colour_manual(values  = setNames(style_df$colour  , style_df$group)) +
    ggplot2::scale_linetype_manual(values = setNames(style_df$linetype, style_df$group)) +
    ggplot2::scale_shape_manual(values    = setNames(style_df$shape   , style_df$group)) +
    ggplot2::labs(title = "Cumulative Cluster Composition per Group",
                  x = cluster_col, y = "Cumulative Proportion") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   legend.position = "right")
  
  if (df) return(cum_tbl)
  p
}

#' Compute and Plot Cumulative Distribution Function (CDF)
#'
#' 이 함수는 주어진 데이터 프레임에서 확률, 로그잇(logit) 변환된 확률, 
#' 또는 비율에 대해 누적분포함수(CDF)를 계산하고, ggplot2를 이용해 
#' CDF 그래프를 그립니다. 필요시 파일로도 저장할 수 있으며, 요약 통계량도 반환합니다.
#'
#' @param data `data.frame` 형태의 입력 데이터.
#' @param probability_col 확률 값을 가진 컬럼명. 문자열(`"Best_Probability"`) 또는  
#'   bare name(`Best_Probability`) 형태로 지정 가능.
#' @param ratio_col 확률 비율 값을 가진 컬럼명. 문자열 또는 bare name 형태 지정.
#' @param plot_type CDF를 계산할 변수 유형.  
#'   - `"probability"`: 원본 확률 값  
#'   - `"logit"`: 확률의 로그잇 변환 (`log(p/(1-p))`)  
#'   - `"ratio"`: 확률 비율 (로그 스케일)  
#'   기본값은 `"probability"`.
#' @param output_file (선택) 결과 플롯을 저장할 파일 경로(예: `"cdf_plot.png"`).  
#'   지정하지 않으면 화면에만 출력.
#'
#' @details
#' - `probability_col`, `ratio_col` 에는 `data` 내 존재하는 컬럼명을 지정해야 합니다.  
#' - `plot_type="ratio"` 인 경우 x축이 로그 스케일로 표현됩니다.  
#' - 반환값으로는 `list(plot, summary)` 형태를 가지며,  
#'   - `plot`: ggplot2 객체  
#'   - `summary`: 평균, 중앙값, 표준편차, 최소값, 최대값을 포함한 요약 통계 데이터프레임
#'
#' @return `list`  
#' \describe{  
#'   \item{plot}{ggplot2 객체, CDF 그래프}  
#'   \item{summary}{`data.frame`, 변수별 요약 통계(Mean, Median, SD, Min, Max)}  
#' }
#'
#' @examples
#' \dontrun{
#' # bare name 사용 예
#' result <- cdf(df, Best_Probability, Probability_Ratio, plot_type = "logit")
#' print(result$plot)
#' 
#' # 문자열 사용 예
#' result2 <- cdf(df, "Best_Probability", "Probability_Ratio", plot_type = "ratio",
#'                output_file = "my_cdf.png")
#' head(result2$summary)
#' }
#'
#' @import dplyr tidyr ggplot2 scales
#' @export
cdf <- function(data, probability_col, ratio_col,
                plot_type = "probability", output_file = NULL) {
  # Validate plot_type
  valid_types <- c("probability", "logit", "ratio")
  if (!plot_type %in% valid_types) {
    stop("Invalid plot_type. Choose 'probability', 'logit', or 'ratio'.")
  }
  
  plot_data <- data %>%
    mutate(
      Probability = as.numeric(.data[[probability_col]]),
      Logit       = log(.data[[probability_col]] / (1 - .data[[probability_col]])),
      Ratio       = as.numeric(.data[[ratio_col]])
    ) %>%
    select(Probability, Logit, Ratio) %>%
    pivot_longer(everything(),
                 names_to  = "Variable",
                 values_to = "Value") %>%
    filter(Variable == case_when(
      plot_type == "probability" ~ "Probability",
      plot_type == "logit"       ~ "Logit",
      TRUE                        ~ "Ratio"
    )) %>%
    arrange(Value) %>%
    mutate(CDF = row_number() / n())
  
  # Set up x-axis scale and labels
  x_scale <- if (plot_type %in% c("probability", "logit")) {
    scale_x_continuous(labels = number_format(accuracy = 0.01))
  } else {
    scale_x_log10(labels = scientific)
  }
  x_label <- switch(plot_type,
                    probability = "Probability",
                    logit       = "Logit of Probability",
                    ratio       = "Probability Ratio (log scale)")
  
  p <- ggplot(plot_data, aes(x = Value, y = CDF)) +
    geom_step() +
    x_scale +
    labs(
      title = paste("CDF of", x_label),
      x     = x_label,
      y     = "Cumulative Probability"
    ) +
    theme_minimal()
  
  if (!is.null(output_file)) {
    ggsave(output_file, plot = p, width = 10, height = 6, dpi = 300)
  }
  
  summary_stats <- plot_data %>%
    summarise(
      Mean   = mean(Value, na.rm = TRUE),
      Median = median(Value, na.rm = TRUE),
      SD     = sd(Value, na.rm = TRUE),
      Min    = min(Value, na.rm = TRUE),
      Max    = max(Value, na.rm = TRUE)
    )
  
  list(plot = p, summary = summary_stats)
}


#' Compute and Plot CDFs Across Multiple or Grouped Datasets
#'
#' 이 함수는 여러 데이터프레임 리스트 또는 단일 통합 데이터프레임에서
#' 그룹별로 CDF를 계산 및 비교합니다. `group_by_col`로 데이터를 분할하고,
#' `sample_col`을 통해 각 데이터셋의 이름을 지정할 수 있습니다.
#'
#' @param data_list `data.frame` 또는 `list` of `data.frame`.
#'   - 단일 `data.frame`에 `group_by_col`을 지정하면 해당 컬럼 값별로 분할합니다.
#'   - `list`로 제공 시 각 요소를 별도 데이터셋으로 처리합니다.
#' @param probability_col 확률 값을 담고 있는 컬럼명 (문자열). 예: "Best_Probability"
#' @param ratio_col 확률 비율 값을 담고 있는 컬럼명 (문자열). 예: "Probability_Ratio"
#' @param plot_type CDF 변수를 선택합니다: "probability", "logit", "ratio".
#'   기본값 "probability".
#' @param group_by_col (선택) 단일 데이터프레임 입력 시 그룹핑할 컬럼명 (문자열).
#' @param sample_col (선택) 리스트 입력 시 또는 분할된 각 데이터프레임에서
#'   데이터셋 이름으로 사용할 컬럼명 (문자열). 지정하지 않으면
#'   리스트 이름 또는 그룹 값이 사용됩니다.
#' @param output_file (선택) 플롯을 저장할 파일 경로. 예: "cdf_plot.png"
#'
#' @return `list`:
#'   - `plot`: `ggplot2` 객체 (비교 CDF 그래프)
#'   - `summary`: `data.frame` (각 데이터셋 요약 통계: Mean, Median, SD, Min, Max)
#'
#' @examples
#' \dontrun{
#' # 리스트 입력 예시
#' res_list <- cdf_multi(
#'   data_list       = list(df1, df2),
#'   probability_col = "Best_Probability",
#'   ratio_col       = "Probability_Ratio",
#'   plot_type       = "logit",
#'   sample_col      = "sample"
#' )
#'
#' # 통합 데이터프레임 + 그룹 컬럼 예시: group_by_col을 명시해주어야 함.
#' res_grouped <- cdf_multi(
#'   data_list       = combined_df,
#'   probability_col = "Best_Probability",
#'   ratio_col       = "Probability_Ratio",
#'   plot_type       = "ratio",
#'   group_by_col    = "GEM",
#'   sample_col      = "sample"
#' )
#' }
#'
#' @import dplyr tidyr ggplot2 scales purrr
#' @export
cdf_multi <- function(
    data_list,
    probability_col,
    ratio_col,
    plot_type     = "probability",
    group_by_col  = NULL,
    sample_col    = NULL,
    output_file   = NULL
) {
  # plot_type 유효성 검사
  types <- c("probability", "logit", "ratio")
  if (!plot_type %in% types) {
    stop("plot_type must be one of: ", paste(types, collapse = ", "))
  }
  
  # 데이터 분할 및 이름 할당
  if (is.data.frame(data_list) && !is.null(group_by_col)) {
    if (!group_by_col %in% colnames(data_list)) {
      stop("group_by_col must be a column in data_list")
    }
    # 그룹 값별로 분할
    groups <- unique(data_list[[group_by_col]])
    data_list <- purrr::map(groups, ~ dplyr::filter(data_list, .data[[group_by_col]] == .x))
    dataset_names <- as.character(groups)
  } else if (is.list(data_list) && all(purrr::map_lgl(data_list, is.data.frame))) {
    # 리스트 입력 시 sample_col이 있으면 해당 값 사용
    if (!is.null(sample_col)) {
      dataset_names <- purrr::map_chr(data_list, function(df) {
        if (!sample_col %in% colnames(df)) {
          stop("sample_col must be a column in each data.frame of data_list")
        }
        # 첫 번째 고유값 사용
        vals <- unique(df[[sample_col]])
        as.character(vals[1])
      })
    } else {
      dataset_names <- names(data_list)
      if (is.null(dataset_names)) {
        dataset_names <- paste0("Dataset_", seq_along(data_list))
      }
    }
  } else {
    # 단일 데이터프레임 입력 시
    data_list <- list(data_list)
    if (!is.null(sample_col) && sample_col %in% colnames(data_list[[1]])) {
      dataset_names <- as.character(unique(data_list[[1]][[sample_col]])[1])
    } else {
      dataset_names <- "Dataset_1"
    }
  }
  
  # CDF 계산을 위한 데이터 준비
  plot_data <- purrr::map2_dfr(
    data_list, dataset_names,
    function(df, name) {
      df %>%
        mutate(
          Probability = as.numeric(.data[[probability_col]]),
          Logit       = log(Probability / (1 - Probability)),
          Ratio       = as.numeric(.data[[ratio_col]]),
          Dataset     = name
        ) %>%
        select(Dataset, Probability, Logit, Ratio)
    }
  ) %>%
    pivot_longer(
      cols      = c(Probability, Logit, Ratio),
      names_to  = "Variable",
      values_to = "Value"
    ) %>%
    filter(
      Variable == dplyr::case_when(
        plot_type == "probability" ~ "Probability",
        plot_type == "logit"       ~ "Logit",
        TRUE                          ~ "Ratio"
      )
    ) %>%
    group_by(Dataset) %>%
    arrange(Value) %>%
    mutate(CDF = row_number() / n()) %>%
    ungroup()
  
  # x축 스케일 및 레이블 설정
  x_scale <- if (plot_type %in% c("probability", "logit")) {
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01))
  } else {
    scale_x_log10(labels = scales::scientific)
  }
  x_label <- switch(
    plot_type,
    probability = "Probability",
    logit       = "Logit of Probability",
    ratio       = "Probability Ratio (log scale)"
  )
  
  # CDF 그래프 생성
  p <- ggplot(plot_data, aes(x = Value, y = CDF, color = Dataset)) +
    geom_step() +
    x_scale +
    labs(
      title = paste("CDF of", x_label),
      x     = x_label,
      y     = "Cumulative Probability"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # 파일 저장 옵션
  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, plot = p, width = 10, height = 6, dpi = 300)
  }
  
  # 요약 통계 계산
  summary_stats <- plot_data %>%
    group_by(Dataset) %>%
    summarise(
      Mean   = mean(Value, na.rm = TRUE),
      Median = median(Value, na.rm = TRUE),
      SD     = sd(Value, na.rm = TRUE),
      Min    = min(Value, na.rm = TRUE),
      Max    = max(Value, na.rm = TRUE)
    ) %>%
    ungroup()
  
  return(list(plot = p, summary = summary_stats))
}




#  Scatter‑smooth with colour & line options


#' Averaged Expression vs Numeric Covariate (with optional colour)
#'
#' For each sample (\code{group.by}), the function first averages gene
#' expression (\code{feature}) across cells, then scatters these averages
#' against a numeric covariate (\code{x_var}). Points can be coloured by a
#' categorical or numeric column, transparency scaled, and different smoothing
#' methods added.
#'
#' @param object      A **Seurat** object *or* a data.frame. When data.frame is
#'                    supplied it must already contain the columns listed below.
#' @param feature     Character; gene symbol (Seurat) or expression column name
#'                    (data.frame).
#' @param group.by    Column used to aggregate cells into samples (default
#'                    "sample_no").
#' @param x_var       Column providing the numeric predictor *per cell* (will be
#'                    averaged just like expression). Typical clinical scores.
#' @param transpose   Logical; if TRUE swap X and Y axes (default FALSE).
#' @param color_by    NULL (no colour) or a column name to colour points.
#' @param palette     Optional palette vector / name passed to scale_*; if NULL
#'                    defaults are used.
#' @param transparency Logical; map \code{color_by} (numeric only) to point
#'                    alpha (default TRUE).
#' @param transparency_desc Logical; when TRUE, higher values become more
#'                    transparent (inverse mapping).
#' @param fitted_line One of "linear" (default), "loess", "lasso", or NULL.
#' @return A \code{ggplot} scatter‑smooth.
#' @export
#' @import dplyr ggplot2 viridisLite
scatter_smooth_colored <- function(object,
                                     feature,
                                     group.by   = "sample_no",
                                     x_var      = "nih_change",
                                     transpose  = FALSE,
                                     color_by   = NULL,
                                     palette    = NULL,
                                     transparency    = TRUE,
                                     transparency_desc = FALSE,
                                     fitted_line = c("linear", "loess", "lasso", NULL)) {
  fitted_line <- match.arg(fitted_line)
  stopifnot(is.character(feature), length(feature) == 1)
  
  # 1. Build per-cell tibble ------------------------------------------------
  
  if (inherits(object, "Seurat")) {
    expr_vec <- Seurat::FetchData(object, vars = feature)[, 1]
    meta_df  <- tibble::as_tibble(object@meta.data)
    cell_df  <- dplyr::mutate(meta_df, !!feature := expr_vec)
  } else {
    cell_df <- tibble::as_tibble(object)
    if (!feature %in% names(cell_df))
      stop("In data.frame mode the column '", feature, "' must exist.")
  }
  
  for (col in c(group.by, x_var, color_by)) {
    if (!is.null(col) && !col %in% names(cell_df))
      stop("Column '", col, "' not found in data.")
  }
  
  
  # 2. Aggregate by group ---------------------------------------------------
  
  # *** 수정된 부분 시작 ***
  #
  # color_by 열이 숫자형인지 미리 확인합니다.
  is_color_numeric <- if (!is.null(color_by)) {
    is.numeric(cell_df[[color_by]])
  } else {
    FALSE # color_by가 NULL이면 어차피 사용되지 않음
  }
  
  agg_df <- cell_df %>%
    dplyr::group_by(.data[[group.by]]) %>%
    dplyr::summarise(
      avg_expr = mean(.data[[feature]], na.rm = TRUE),
      x_val    = mean(.data[[x_var]], na.rm = TRUE),
      # color_by의 타입에 따라 집계 방식을 다르게 합니다.
      colour   = if (!is.null(color_by)) {
        if (is_color_numeric) {
          # 숫자형이면 평균 계산
          mean(.data[[color_by]], na.rm = TRUE)
        } else {
          # 범주형(문자열/팩터)이면 그룹의 첫 번째 값 사용
          dplyr::first(.data[[color_by]], na_rm = TRUE)
        }
      } else {
        NA
      },
      .groups  = "drop"
    )
  #
  # *** 수정된 부분 끝 ***
  
  
  # 3. Aesthetics -----------------------------------------------------------
  
  x_col <- if (transpose) "avg_expr" else "x_val"
  y_col <- if (transpose) "x_val"   else "avg_expr"
  
  p <- ggplot2::ggplot(agg_df, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]]))
  
  # Point layer with colour / alpha mapping
  if (!is.null(color_by)) {
    # (원본 코드)
    # 여기서는 cell_df의 원본 열 타입을 확인하는 것이 올바릅니다.
    # 위에서 agg_df$colour의 타입을 올바르게 생성했기 때문에
    # 이 로직이 이제 정상적으로 작동합니다.
    if (is.numeric(cell_df[[color_by]])) {
      # numeric colour gradient
      p <- p + ggplot2::geom_point(ggplot2::aes(colour = colour,
                                                 alpha  = colour), size = 3)
      alpha_range <- if (transparency_desc) c(1, 0.2) else c(0.2, 1)
      if (transparency) {
        p <- p + ggplot2::scale_alpha(range = alpha_range, guide = "none")
      } else {
        p <- p + ggplot2::guides(alpha = "none")
      }
      pal <- if (is.null(palette)) viridisLite::viridis(256) else palette
      p <- p + ggplot2::scale_colour_gradientn(colours = pal, name = color_by)
    } else {
      # categorical palette
      # agg_df$colour가 이제 문자열이므로 ggplot이 자동으로 이산형으로 처리합니다.
      p <- p + ggplot2::geom_point(ggplot2::aes(colour = colour), size = 3)
      pal <- palette %||% RColorBrewer::brewer.pal(max(3, length(unique(agg_df$colour))), "Set1")
      p <- p + ggplot2::scale_colour_manual(values = pal, name = color_by)
    }
  } else {
    p <- p + ggplot2::geom_point(size = 3)
  }
  
  
  # 4. Smoothing line -------------------------------------------------------
  
  if (!is.null(fitted_line)) {
    if (fitted_line == "linear") {
      p <- p + ggplot2::geom_smooth(method = "lm", se = TRUE, colour = "black")
      fit <- stats::lm(agg_df[[y_col]] ~ agg_df[[x_col]])
      coef <- round(stats::coef(fit), 3)
      pval <- signif(summary(fit)$coefficients[2, 4], 3)
      annot <- paste0("y = ", coef[1], " + ", coef[2], " * x\np = ", pval)
      p <- p + ggplot2::annotate("text", x = min(agg_df[[x_col]], na.rm = TRUE),
                                 y = max(agg_df[[y_col]], na.rm = TRUE),
                                 label = annot, hjust = 0, vjust = 1, size = 4)
    } else if (fitted_line == "loess") {
      p <- p + ggplot2::geom_smooth(method = "loess", se = TRUE, colour = "black")
    } else if (fitted_line == "lasso") {
      if (!requireNamespace("glmnet", quietly = TRUE)) {
        warning("glmnet not installed; falling back to linear fit.")
        p <- p + ggplot2::geom_smooth(method = "lm", se = TRUE)
      } else {
        xmat <- as.matrix(agg_df[[x_col]])
        fit  <- glmnet::cv.glmnet(xmat, agg_df[[y_col]], alpha = 1)
        preds <- as.numeric(glmnet::predict.glmnet(fit$glmnet.fit, newx = xmat,
                                                   s = fit$lambda.min))
        pred_df <- agg_df %>% dplyr::mutate(pred = preds)
        p <- p + ggplot2::geom_line(data = pred_df[order(pred_df[[x_col]]), ],
                                    ggplot2::aes(x = .data[[x_col]], y = pred),
                                    colour = "red", linewidth = 1)
      }
    }
  }
  
  
  # 5. Labels & theme -------------------------------------------------------
  
  p <- p + ggplot2::theme_bw() +
    ggplot2::labs(x = if (transpose) paste("Average", feature, "expression") else x_var,
                  y = if (transpose) x_var else paste("Average", feature, "expression"),
                  colour = color_by)
  
  p
  return(p)
}






#' Create a Heatmap of Gene Set Expression per Cluster
#'
#' This function creates a heatmap showing the normalized expression of gene sets across clusters.
#' The expression values are z-score normalized for better visualization.
#'
#' @param sobj A Seurat object
#' @param group Character string specifying the identity to use for clustering.
#'             Default is "seurat_clusters"
#' @param value Character string specifying how to aggregate expression values.
#'             Options: "average" (default) or "sum"
#' @param assay Character string specifying which assay to use.
#'             Default is "SCT"
#' @param gene_sets A list of character vectors, where each vector contains gene names.
#'                 The list can be named, and names will be used as gene set labels.
#'                 If a single vector is provided, it will be converted to a list.
#' @param title Character string for the plot title.
#'             Default is "Normalized Gene Set Expression per Cluster"
#' @param x_label Character string for the x-axis label.
#'               Default is "Cluster"
#' @param y_label Character string for the y-axis label.
#'               Default is "Gene Set"
#'
#' @return A data frame containing the normalized expression values for each gene set
#'         across clusters, with an additional column indicating the assigned cell type
#'         based on highest expression
#'
#' @examples
#' \dontrun{
#' # Create heatmap with named gene sets
#' gene_sets <- list(
#'   T_cells = c("CD3D", "CD3E", "CD4", "CD8A"),
#'   B_cells = c("MS4A1", "CD19", "CD79A"),
#'   Myeloid = c("CD14", "FCGR3A", "LYZ")
#' )
#' result <- myhm_genesets2(sobj, gene_sets = gene_sets)
#' 
#' # Create heatmap with a single gene set
#' result <- myhm_genesets2(sobj, gene_sets = c("CD3D", "CD3E", "CD4"))
#' }
#' @export
myhm_genesets4 <- function(
    sobj,
    group = "seurat_clusters",
    value = "average",
    assay = "SCT",
    gene_sets = NULL,
    title="Normalized Gene Set Expression per Cluster",
    x_label="Cluster",
    y_label="Gene Set"
){
  library(Seurat)
  library(dplyr)
  library(reshape2)
  library(ggplot2)
  
  #-------------------------------'
  # (A) 유효성 체크
  #-------------------------------'
  if(is.null(gene_sets)){
    stop("gene_sets를 지정해 주세요. (ex: list(Immune=c('CD3D','CD3E'), Bcell=c('MS4A1','CD79A'))) ")
  }
  
  # 만약 gene_sets가 리스트가 아니라 벡터만 들어왔다면 리스트로 변환
  if(!is.list(gene_sets)){
    gene_sets <- list(GeneSet1 = gene_sets)
  }
  
  # 이름이 없는 리스트 원소가 있다면 자동으로 이름 부여
  if(is.null(names(gene_sets)) || any(names(gene_sets) == "")){
    for(i in seq_along(gene_sets)){
      if(is.null(names(gene_sets)[i]) || names(gene_sets)[i] == ""){
        names(gene_sets)[i] <- paste0("GeneSet", i)
      }
    }
  }
  
  #-------------------------------'
  # (B) Seurat 객체에 grouping 적용
  #-------------------------------'
  Idents(sobj) <- group
  
  #-------------------------------'
  # (C) 평균 발현량(또는 합계 등) 계산
  #-------------------------------'
  if(value == "average"){
    cluster_avg <- AverageExpression(sobj, assays = assay, slot = "data", group.by = group)[[assay]]
  } else {
    cluster_avg <- AggregateExpression(sobj, assays = assay, slot = "data", group.by = group)[[assay]]
  }
  
  # Handle case when result is a vector (single group)
  if(is.null(dim(cluster_avg)) || length(dim(cluster_avg)) == 1){
    # Convert vector to single-column matrix
    cluster_avg <- as.matrix(cluster_avg)
    if(ncol(cluster_avg) == 1 && is.null(colnames(cluster_avg))){
      # Get the unique group value
      unique_groups <- unique(as.character(Idents(sobj)))
      if(length(unique_groups) == 1){
        colnames(cluster_avg) <- unique_groups[1]
      } else {
        colnames(cluster_avg) <- "Group1"
      }
    }
  }
  
  #-------------------------------'
  # (D) Gene Set 별 발현량 계산
  #-------------------------------'
  cluster_names <- colnames(cluster_avg)
  
  # 결과를 담을 data.frame 생성
  gene_set_expression <- data.frame(Cluster = cluster_names, stringsAsFactors = FALSE)
  
  # gene_sets 각각에 대해 평균 발현량을 구함
  for(gset_name in names(gene_sets)){
    genes <- gene_sets[[gset_name]]
    genes_present <- genes[genes %in% rownames(cluster_avg)]
    
    if(length(genes_present) == 0){
      warning(paste("No genes from", gset_name, "found in the dataset."))
      gene_set_expression[[gset_name]] <- NA
      next
    }
    
    # Calculate means
    subset_data <- cluster_avg[genes_present, , drop = FALSE]
    
    # Ensure subset_data is a matrix
    if(!is.matrix(subset_data)){
      subset_data <- as.matrix(subset_data)
    }
    
    # Calculate column means
    if(nrow(subset_data) == 1){
      # If only one gene, the values are already the means
      gene_set_expression[[gset_name]] <- as.numeric(subset_data[1,])
    } else if(ncol(subset_data) == 1){
      # If only one cluster/group
      gene_set_expression[[gset_name]] <- mean(subset_data[,1])
    } else {
      # Multiple genes and multiple clusters
      gene_set_expression[[gset_name]] <- colMeans(subset_data)
    }
  }
  
  #-------------------------------'
  # (E) Z-score 정규화
  #-------------------------------'
  gene_set_expression_normalized <- gene_set_expression
  
  # Check if we have more than one cluster
  if(nrow(gene_set_expression_normalized) > 1){
    gene_set_expression_normalized[,-1] <- scale(gene_set_expression_normalized[,-1])
  } else {
    # For single cluster, z-score normalization doesn't make sense
    # Set all values to 0 (neutral)
    gene_set_expression_normalized[,-1] <- 0
  }
  
  # 각 Cluster에서 가장 높은 값을 가지는 gene set을 배정
  gene_set_expression_normalized$Assigned_CellType <- apply(
    gene_set_expression_normalized[,-c(1, ncol(gene_set_expression_normalized))], 1, 
    function(x){
      if(all(is.na(x))) return(NA)
      names(x)[which.max(x)]
    }
  )
  
  #-------------------------------'
  # (F) 클러스터 순서 정렬
  #-------------------------------'
  numeric_test <- suppressWarnings(as.numeric(gene_set_expression_normalized$Cluster))
  
  if(!all(is.na(numeric_test))){
    if(sum(is.na(numeric_test)) == 0){
      sorted_levels <- sort(unique(numeric_test))
      gene_set_expression_normalized$Cluster <- factor(
        gene_set_expression_normalized$Cluster,
        levels = as.character(sorted_levels)
      )
    } else {
      sorted_levels <- sort(unique(gene_set_expression_normalized$Cluster))
      gene_set_expression_normalized$Cluster <- factor(
        gene_set_expression_normalized$Cluster,
        levels = sorted_levels
      )
    }
  } else {
    sorted_levels <- sort(unique(gene_set_expression_normalized$Cluster))
    gene_set_expression_normalized$Cluster <- factor(
      gene_set_expression_normalized$Cluster,
      levels = sorted_levels
    )
  }
  
  #-------------------------------'
  # (G) Heatmap용 long format 만들기
  #-------------------------------'
  melted_data <- melt(
    gene_set_expression_normalized,
    id.vars = c("Cluster","Assigned_CellType"),
    variable.name = "GeneSet",
    value.name = "Zscore"
  )
  
  #-------------------------------'
  # (H) Heatmap 그리기
  #-------------------------------'
  p <- ggplot(melted_data, aes(x = Cluster, y = GeneSet, fill = Zscore)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    theme_minimal() +
    labs(
      title = title,
      x = x_label,
      y = y_label
    ) +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=12),
          axis.text.y = element_text(hjust=1, size=12),
          axis.title.x = element_text(face="bold", size=14),
          axis.title.y = element_text(face="bold", size=14),
          plot.title = element_text(size=16, face="bold", hjust=0.5))
  
  print(p)
  
  #-------------------------------'
  # (I) 결과 반환
  #-------------------------------'
  return(gene_set_expression_normalized)
}

#' Create a Heatmap of Individual Gene Expression per Cluster
#'
#' This function creates a heatmap showing the normalized expression of individual genes
#' across clusters. The expression values are z-score normalized for better visualization.
#'
#' @param sobj A Seurat object
#' @param group Character string specifying the identity to use for clustering.
#'             Default is "seurat_clusters"
#' @param value Character string specifying how to aggregate expression values.
#'             Options: "average" (default) or "sum"
#' @param assay Character string specifying which assay to use.
#'             Default is "SCT"
#' @param genes Character vector containing gene names to plot
#' @param title Character string for the plot title.
#'             Default is "Normalized Gene Expression per Cluster"
#' @param x_label Character string for the x-axis label.
#'               Default is "Cluster"
#' @param y_label Character string for the y-axis label.
#'               Default is "Genes"
#'
#' @return A data frame containing the normalized expression values for each gene
#'         across clusters
#'
#' @examples
#' \dontrun{
#' # Create heatmap for a set of genes
#' genes <- c("CD3D", "CD3E", "CD4", "CD8A", "MS4A1", "CD19")
#' result <- myhm_genes2(sobj, genes = genes)
#' 
#' # Create heatmap with custom title and labels
#' result <- myhm_genes2(sobj, 
#'                      genes = genes,
#'                      title = "T and B Cell Markers",
#'                      x_label = "Cell Clusters",
#'                      y_label = "Marker Genes")
#' }
#' @export
myhm_genes4 <- function(
    sobj,
    group = "seurat_clusters",
    value = "average",
    assay = "SCT",
    genes = NULL,
    title = "Normalized Gene Expression per Cluster",
    x_label = "Cluster",
    y_label = "Genes"
){
  library(Seurat)
  library(dplyr)
  library(reshape2)
  library(ggplot2)
  
  #-------------------------------'
  # (A) 유효성 체크
  #-------------------------------'
  if(is.null(genes)){
    stop("genes를 지정해 주세요. (예: genes=c('CD3D','CD3E','MS4A1'))")
  }
  
  #-------------------------------'
  # (B) Seurat 객체에 grouping 적용
  #-------------------------------'
  Seurat::Idents(sobj) <- group
  
  #-------------------------------'
  # (C) 평균 발현량(또는 합계 등) 계산
  #-------------------------------'
  if(value == "average"){
    cluster_avg <- Seurat::AverageExpression(sobj, assays = assay, slot = "data", group.by = group)[[assay]]
  } else {
    cluster_avg <- Seurat::AggregateExpression(sobj, assays = assay, slot = "data", group.by = group)[[assay]]
  }
  
  # Handle case when cluster_avg is a vector (single group)
  if(is.null(dim(cluster_avg)) || length(dim(cluster_avg)) == 1){
    cluster_avg <- as.matrix(cluster_avg)
    if(ncol(cluster_avg) == 1 && is.null(colnames(cluster_avg))){
      # Get the unique group value
      unique_groups <- unique(as.character(Seurat::Idents(sobj)))
      if(length(unique_groups) == 1){
        colnames(cluster_avg) <- unique_groups[1]
      } else {
        colnames(cluster_avg) <- "Group1"
      }
    }
  }
  
  #-------------------------------'
  # (D) 원하는 유전자만 필터
  #-------------------------------'
  genes_present <- genes[genes %in% rownames(cluster_avg)]
  if(length(genes_present) == 0){
    stop("지정하신 유전자 중 데이터셋에 존재하는 유전자가 없습니다.")
  }
  
  if(length(genes_present) < length(genes)){
    missing_genes <- genes[!genes %in% rownames(cluster_avg)]
    warning(paste("다음 유전자들이 데이터셋에 없습니다:", paste(missing_genes, collapse=", ")))
  }
  
  # Subset 
  gene_expression <- cluster_avg[genes_present, , drop=FALSE]
  
  # Ensure it's a matrix
  if(!is.matrix(gene_expression)){
    gene_expression <- as.matrix(gene_expression)
  }
  
  #-------------------------------'
  # (E) Z-score 정규화
  #-------------------------------'
  # Check number of clusters
  n_clusters <- ncol(gene_expression)
  
  if(n_clusters == 1){
    # For single cluster, z-score normalization across genes doesn't make sense
    # We can normalize across genes instead or set to 0
    gene_expression_scaled <- matrix(0, 
                                     nrow = length(genes_present), 
                                     ncol = 1,
                                     dimnames = list(genes_present, colnames(gene_expression)))
    
    # Alternative: scale across genes for the single cluster
    # This shows relative expression levels among genes
    if(length(genes_present) > 1){
      scaled_values <- scale(gene_expression[,1])
      gene_expression_scaled[,1] <- scaled_values
    }
  } else {
    # Multiple clusters - scale across clusters for each gene
    # This is the original behavior
    gene_expression_scaled <- t(scale(t(gene_expression)))
  }
  
  # Convert to data.frame
  gene_expression_df <- as.data.frame(t(gene_expression_scaled))
  gene_expression_df$Cluster <- rownames(gene_expression_df)
  
  #-------------------------------'
  # (F) 클러스터 순서 정렬
  #-------------------------------'
  numeric_test <- suppressWarnings(as.numeric(gene_expression_df$Cluster))
  
  if(!all(is.na(numeric_test))){
    if(sum(is.na(numeric_test)) == 0){
      sorted_levels <- sort(unique(numeric_test))
      gene_expression_df$Cluster <- factor(
        gene_expression_df$Cluster,
        levels = as.character(sorted_levels)
      )
    } else {
      sorted_levels <- sort(unique(gene_expression_df$Cluster))
      gene_expression_df$Cluster <- factor(
        gene_expression_df$Cluster,
        levels = sorted_levels
      )
    }
  } else {
    sorted_levels <- sort(unique(gene_expression_df$Cluster))
    gene_expression_df$Cluster <- factor(
      gene_expression_df$Cluster,
      levels = sorted_levels
    )
  }
  
  #-------------------------------'
  # (G) long format으로 melt
  #-------------------------------'
  melted_data <- reshape2::melt(
    gene_expression_df,
    id.vars = "Cluster",
    variable.name = "Gene",
    value.name = "Zscore"
  )
  
  # 유전자 순서를 입력 순서대로 유지
  melted_data$Gene <- factor(melted_data$Gene, levels = genes_present)
  
  #-------------------------------'
  # (H) Heatmap 그리기
  #-------------------------------'
  # Adjust plot for single cluster case
  if(n_clusters == 1){
    # For single cluster, make the heatmap narrower
    p <- ggplot2::ggplot(melted_data, ggplot2::aes(x = Cluster, y = Gene, fill = Zscore)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = title,
        x = x_label,
        y = y_label
      ) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle=45, hjust=1, size=12),
        axis.text.y = ggplot2::element_text(size=12),
        axis.title.x = ggplot2::element_text(face="bold", size=14),
        axis.title.y = ggplot2::element_text(face="bold", size=14),
        plot.title = ggplot2::element_text(size=16, face="bold", hjust=0.5),
        aspect.ratio = length(genes_present)/2  # Make plot taller for single cluster
      )
  } else {
    # Regular plot for multiple clusters
    p <- ggplot2::ggplot(melted_data, ggplot2::aes(x = Cluster, y = Gene, fill = Zscore)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = title,
        x = x_label,
        y = y_label
      ) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle=45, hjust=1, size=12),
        axis.text.y = ggplot2::element_text(size=12),
        axis.title.x = ggplot2::element_text(face="bold", size=14),
        axis.title.y = ggplot2::element_text(face="bold", size=14),
        plot.title = ggplot2::element_text(size=16, face="bold", hjust=0.5)
      )
  }
  
  print(p)
  
  #-------------------------------'
  # (I) 결과 반환
  #-------------------------------'
  return(gene_expression_df)
}


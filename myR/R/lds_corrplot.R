# ============================================================================
# LDS SVA 상관관계 플롯 Helper Functions
# ============================================================================

#' @title SVA와 메타데이터 간 상관관계 플롯 생성
#' @description SVA 객체(또는 SV 매트릭스)와 메타데이터를 받아 상관관계 플롯 생성
#' @param sva SVA 객체 또는 SV 매트릭스 (행=셀, 열=SV)
#' @param meta.data 메타데이터 (행=셀, 열=변수)
#' @param method 상관관계 계산 방법 ("spearman" 또는 "pearson"). 기본값 "spearman"
#' @param display.value 상관관계 값을 플롯에 표시할지 여부. 기본값 FALSE
#' @param clustering 계층적 클러스터링 적용 여부. 기본값 TRUE
#' @param top_n 상관관계가 높은 상위 n개 변수만 표시 (NULL이면 모두 표시). 기본값 NULL
#' @param min_corr 최소 상관관계 절댓값 (이보다 작은 값은 필터링). 기본값 0.1
#' @param output_file 출력 파일 경로 (NULL이면 플롯만 반환). 기본값 NULL
#' @param width 플롯 너비 (인치). 기본값 12
#' @param height 플롯 높이 (인치). 기본값 8
#' @param ... 추가 파라미터 (pheatmap 또는 corrplot에 전달)
#' @return 상관관계 행렬 (invisible)
#' @export
lds_corrplot <- function(sva,
                         meta.data,
                         method = c("spearman", "pearson"),
                         display.value = FALSE,
                         clustering = TRUE,
                         top_n = NULL,
                         min_corr = 0.1,
                         output_file = NULL,
                         width = 12,
                         height = 8,
                         ...) {
  
  # 1. 입력 검증
  method <- match.arg(method)
  
  # sva 객체에서 SV 추출
  if (inherits(sva, "list") && "sv" %in% names(sva)) {
    # SVA 객체인 경우
    sv_matrix <- sva$sv
    # 열 이름이 없으면 SV1, SV2, ... 부여
    if (is.null(colnames(sv_matrix))) {
      colnames(sv_matrix) <- paste0("SV", 1:ncol(sv_matrix))
    }
  } else if (is.matrix(sva) || is.data.frame(sva)) {
    # SV 매트릭스인 경우
    sv_matrix <- as.matrix(sva)
    if (is.null(colnames(sv_matrix))) {
      colnames(sv_matrix) <- paste0("SV", 1:ncol(sv_matrix))
    }
  } else {
    stop("'sva'는 SVA 객체 또는 SV 매트릭스여야 합니다.")
  }
  
  # 행 수 일치 확인
  if (nrow(sv_matrix) != nrow(meta.data)) {
    stop("SV 매트릭스와 메타데이터의 행 수가 일치하지 않습니다: ",
         nrow(sv_matrix), " vs ", nrow(meta.data))
  }
  
  # 2. 숫자형/팩터형/문자형 변수만 선택
  numeric_vars <- sapply(colnames(meta.data), function(v) {
    var_data <- meta.data[[v]]
    is.numeric(var_data) || is.factor(var_data) || is.character(var_data)
  })
  meta_numeric <- meta.data[, numeric_vars, drop = FALSE]
  
  # 문자형을 factor로 변환 후, 팩터를 숫자로 변환
  for (i in 1:ncol(meta_numeric)) {
    if (is.character(meta_numeric[[i]])) {
      meta_numeric[[i]] <- as.factor(meta_numeric[[i]])
    }
    if (is.factor(meta_numeric[[i]])) {
      meta_numeric[[i]] <- as.numeric(meta_numeric[[i]])
    }
  }
  
  if (ncol(meta_numeric) == 0) {
    stop("숫자형/팩터형 메타데이터 변수가 없습니다.")
  }
  
  # 3. 상관관계 계산 (SV가 열이 되도록)
  sv_names <- colnames(sv_matrix)
  meta_var_names <- colnames(meta_numeric)
  
  # SV가 열, 메타데이터가 행이 되도록 전치
  cor_matrix <- matrix(NA, nrow = length(meta_var_names), ncol = length(sv_names))
  rownames(cor_matrix) <- meta_var_names
  colnames(cor_matrix) <- sv_names
  
  for (i in seq_along(meta_var_names)) {
    for (j in seq_along(sv_names)) {
      sv_vec <- sv_matrix[, j]
      meta_vec <- meta_numeric[[i]]
      
      # NA 제거 후 상관관계 계산
      valid_idx <- !is.na(sv_vec) & !is.na(meta_vec)
      if (sum(valid_idx) > 1) {
        cor_matrix[i, j] <- cor(meta_vec[valid_idx], sv_vec[valid_idx], method = method)
      }
    }
  }
  
  # 4. 필터링 (min_corr, top_n) - 행(메타데이터) 기준
  if (!is.null(min_corr) && min_corr > 0) {
    abs_cor <- abs(cor_matrix)
    abs_cor[is.na(abs_cor)] <- 0
    keep_vars <- rowSums(abs_cor >= min_corr, na.rm = TRUE) > 0
    if (sum(keep_vars) > 0) {
      cor_matrix <- cor_matrix[keep_vars, , drop = FALSE]
      meta_var_names <- rownames(cor_matrix)
    } else {
      warning("min_corr 조건을 만족하는 변수가 없습니다. 모든 변수를 표시합니다.")
    }
  }
  
  if (!is.null(top_n) && top_n > 0 && top_n < nrow(cor_matrix)) {
    # 각 메타데이터 변수에 대해 가장 높은 상관관계를 가진 SV 찾기
    abs_cor <- abs(cor_matrix)
    abs_cor[is.na(abs_cor)] <- 0
    
    # 각 행(메타데이터)의 최대 절댓값 상관관계로 정렬
    max_cor_per_var <- apply(abs_cor, 1, max, na.rm = TRUE)
    top_vars <- names(sort(max_cor_per_var, decreasing = TRUE))[1:min(top_n, length(max_cor_per_var))]
    
    if (length(top_vars) > 0) {
      cor_matrix <- cor_matrix[top_vars, , drop = FALSE]
    }
  }
  
  # 5. 플롯 생성
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("'pheatmap' 패키지가 필요합니다: install.packages('pheatmap')")
  }
  
  # NA 처리 (클러스터링을 위해)
  cor_matrix_for_clust <- cor_matrix
  cor_matrix_for_clust[is.na(cor_matrix_for_clust)] <- 0
  
  # 플롯 파라미터 설정 (SV가 열이 되도록)
  plot_params <- list(
    mat = cor_matrix,  # 행=메타데이터, 열=SV
    color = colorRampPalette(c("blue", "white", "red"))(100),
    breaks = seq(-1, 1, length.out = 101),
    cluster_rows = clustering,
    cluster_cols = clustering,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize = 8,
    fontsize_row = 10,
    fontsize_col = 10,
    na_col = "gray90",
    ...
  )
  
  # display.value가 TRUE면 숫자 표시
  if (display.value) {
    plot_params$display_numbers <- TRUE
    plot_params$number_format <- "%.2f"
    plot_params$fontsize_number <- 7
  }
  
  # 출력 파일이 지정된 경우
  if (!is.null(output_file)) {
    if (grepl("\\.png$", output_file, ignore.case = TRUE)) {
      png(output_file, width = width * 100, height = height * 100, res = 100)
    } else if (grepl("\\.pdf$", output_file, ignore.case = TRUE)) {
      pdf(output_file, width = width, height = height)
    } else {
      stop("지원하는 파일 형식: .png, .pdf")
    }
  }
  
  # 플롯 생성
  tryCatch({
    do.call(pheatmap::pheatmap, plot_params)
  }, error = function(e) {
    # 클러스터링 실패 시 NA를 0으로 대체하여 재시도
    if (clustering) {
      warning("클러스터링 실패, NA를 0으로 대체하여 재시도합니다.")
      plot_params$mat <- cor_matrix_for_clust
      do.call(pheatmap::pheatmap, plot_params)
    } else {
      stop("플롯 생성 실패: ", conditionMessage(e))
    }
  })
  
  if (!is.null(output_file)) {
    dev.off()
    message("플롯 저장: ", output_file)
  }
  
  # 6. 결과 반환 (invisible)
  return(invisible(cor_matrix))
}

#' @title SVA 상관관계 요약 테이블 생성
#' @description SVA와 메타데이터 간 상관관계를 요약한 테이블 반환
#' @param sva SVA 객체 또는 SV 매트릭스
#' @param meta.data 메타데이터
#' @param method 상관관계 계산 방법 ("spearman" 또는 "pearson")
#' @param top_n 상위 n개 결과만 반환 (NULL이면 모두 반환)
#' @return data.frame (SV, 변수, 상관관계, 절댓값)
#' @export
lds_corrplot_summary <- function(sva,
                                 meta.data,
                                 method = c("spearman", "pearson"),
                                 top_n = 20) {
  
  method <- match.arg(method)
  
  # 상관관계 계산 (lds_corrplot과 동일한 로직)
  if (inherits(sva, "list") && "sv" %in% names(sva)) {
    sv_matrix <- sva$sv
    if (is.null(colnames(sv_matrix))) {
      colnames(sv_matrix) <- paste0("SV", 1:ncol(sv_matrix))
    }
  } else if (is.matrix(sva) || is.data.frame(sva)) {
    sv_matrix <- as.matrix(sva)
    if (is.null(colnames(sv_matrix))) {
      colnames(sv_matrix) <- paste0("SV", 1:ncol(sv_matrix))
    }
  } else {
    stop("'sva'는 SVA 객체 또는 SV 매트릭스여야 합니다.")
  }
  
  numeric_vars <- sapply(colnames(meta.data), function(v) {
    var_data <- meta.data[[v]]
    is.numeric(var_data) || is.factor(var_data) || is.character(var_data)
  })
  meta_numeric <- meta.data[, numeric_vars, drop = FALSE]
  
  for (i in 1:ncol(meta_numeric)) {
    if (is.character(meta_numeric[[i]])) {
      meta_numeric[[i]] <- as.factor(meta_numeric[[i]])
    }
    if (is.factor(meta_numeric[[i]])) {
      meta_numeric[[i]] <- as.numeric(meta_numeric[[i]])
    }
  }
  
  # 상관관계 계산
  sv_names <- colnames(sv_matrix)
  meta_var_names <- colnames(meta_numeric)
  
  cor_list <- list()
  for (i in seq_along(sv_names)) {
    for (j in seq_along(meta_var_names)) {
      sv_vec <- sv_matrix[, i]
      meta_vec <- meta_numeric[[j]]
      
      valid_idx <- !is.na(sv_vec) & !is.na(meta_vec)
      if (sum(valid_idx) > 1) {
        cor_val <- cor(sv_vec[valid_idx], meta_vec[valid_idx], method = method)
        if (!is.na(cor_val)) {
          cor_list[[length(cor_list) + 1]] <- data.frame(
            SV = sv_names[i],
            Variable = meta_var_names[j],
            Correlation = cor_val,
            AbsCorrelation = abs(cor_val),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  
  if (length(cor_list) == 0) {
    return(data.frame(SV = character(), Variable = character(), 
                      Correlation = numeric(), AbsCorrelation = numeric()))
  }
  
  result_df <- do.call(rbind, cor_list)
  result_df <- result_df[order(result_df$AbsCorrelation, decreasing = TRUE), ]
  
  if (!is.null(top_n) && top_n > 0) {
    result_df <- result_df[1:min(top_n, nrow(result_df)), ]
  }
  
  rownames(result_df) <- NULL
  return(result_df)
}


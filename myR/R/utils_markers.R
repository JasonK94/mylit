#' Filter and process marker genes from Seurat's FindMarkers or FindAllMarkers results
#'
#' This function processes marker gene results by filtering based on log fold change direction,
#' adjusted p-value cutoff, and adding percentage difference information.
#'
#' @param markers A data frame from Seurat's FindMarkers or FindAllMarkers results
#' @param sign Character string indicating the direction of log fold change to keep.
#'             Options: NULL (keep all), "+" (positive only), "-" (negative only)
#' @param p_cutoff Numeric value for adjusted p-value cutoff. Default is NULL (no filtering)
#' @param filter Additional filtering criteria (not used in this function)
#'
#' @return A processed data frame with filtered markers and added pct.diff column
#'
#' @examples
#' \dontrun{
#' # Filter markers to keep only positive log fold changes
#' markers_positive <- marker_trim(markers, sign = "+")
#'
#' # Filter markers with adjusted p-value < 0.05
#' markers_sig <- marker_trim(markers, p_cutoff = 0.05)
#' }
#' @export
marker_trim <- function(markers, sign = NULL, p_cutoff = NULL, filter = NULL) {
  if ("cluster" %in% names(markers)) {} else { # FindAllMarkers object has "cluster" column, FindMarkers not.
    markers$gene <- rownames(markers)
  }
  markers$pct.diff <- markers$pct.1 - markers$pct.2
  if (is.null(sign)) {} else if (sign == "+") {
    markers <- markers[markers$avg_log2FC > 0, ]
  } else if (sign == "-") {
    markers <- markers[markers$avg_log2FC < 0, ]
  } else {
    print("sign is wrong. put NULL or " + " or " - " ")
  }

  if (is.null(p_cutoff)) {} else {
    markers <- markers[markers$p_val_adj < p_cutoff, ]
  }

  return(markers)
}

#' Filter out unwanted genes from marker results
#'
#' This function removes genes matching specific patterns (ribosomal, mitochondrial,
#' hemoglobin, etc.) from marker gene results.
#'
#' @param markers A data frame from Seurat's FindMarkers or FindAllMarkers results
#' @param filter Character vector specifying which gene types to filter out.
#'               Options: "rb" (ribosomal), "mt" (mitochondrial), "hb" (hemoglobin),
#'               "AC" (AC/AL genes), "ENSG" (ENSG genes), "LINC" (LINC genes)
#'
#' @return A filtered data frame with unwanted genes removed
#'
#' @examples
#' \dontrun{
#' # Remove ribosomal and mitochondrial genes
#' markers_filtered <- marker_filter(markers, filter = c("rb", "mt"))
#'
#' # Remove all unwanted gene types
#' markers_clean <- marker_filter(markers, filter = c("rb", "mt", "hb", "AC", "ENSG", "LINC"))
#' }
#' @export
marker_filter <- function(markers, filter = c("rb", "mt", "hb", "AC", "ENSG", "LINC")) {
  rb <- mt <- AC <- ENSG <- LINC <- ""
  if ("rb" %in% filter) {
    rb <- "^RPL|^RPS"
  }
  if ("mt" %in% filter) {
    mt <- "^MT-"
  }
  if ("AC" %in% filter) {
    AC <- "^(AC\\d+|AL\\d+)"
  }
  if ("ENSG" %in% filter) {
    ENSG <- "^ENSG"
  }
  if ("LINC" %in% filter) {
    LINC <- "^LINC"
  }
  filter_pattern <- paste(rb, mt, AC, ENSG, LINC, sep = "|")
  if ("gene" %in% names(markers)) {} else { # FindAllMarkers object generates "gene" column automatically
    markers$gene <- rownames(markers)
  }
  markers <- markers[!grepl(filter_pattern, markers$gene), ]

  if ("hb" %in% filter) {
    globin_genes <- c(
      "HBA1", "HBA2", "HBM", "HBQ1", "HBQ2", "HBZ", "HBZP1", "HBZP2",
      "HBB", "HBD", "HBG1", "HBG2", "HBE1", "HBBS", "HBBP1", "HBBP2"
    )
    markers <- markers[!markers$gene %in% globin_genes, ]
  }

  return(markers)
}

#' filterout (Intercept) terms
#' @export
lrf <- function(lmm_result) {
  lmm_result$summary <- lmm_result$summary[lmm_result$summary$term != "(Intercept)", ]
  return(lmm_result)
}


#' Convert FindAllMarkers results to a list organized by cluster
#'
#' This function takes the results from Seurat's FindAllMarkers and organizes them
#' into a list where each element contains markers for a specific cluster.
#'
#' @param markers A data frame from Seurat's FindAllMarkers results
#'
#' @return A list where each element is a data frame containing markers for one cluster.
#'         Returns NULL if the input doesn't have a 'cluster' column.
#'
#' @examples
#' \dontrun{
#' # Convert markers to a list by cluster
#' marker_list <- all_markers_to_list(markers)
#'
#' # Access markers for a specific cluster
#' cluster_0_markers <- marker_list[["cluster_0"]]
#' }
#' @export
all_markers_to_list <- function(markers) {
  if (!"cluster" %in% names(markers)) {
    print("column <cluster> should be in the dataframe.")
    return(NULL)
  } else {
    marker_list <- list()
    for (i in unique(markers$cluster)) {
      name <- paste0("cluster_", i)
      marker_list[[name]] <- markers[markers$cluster == i, ]
    }
  }
  return(marker_list)
}

#' Print top marker genes for each cluster
#'
#' This function prints the top N marker genes for each cluster, either from a
#' marker list or a FindAllMarkers data frame.
#'
#' @param markers Either a list of marker data frames (from all_markers_to_list)
#'               or a FindAllMarkers data frame
#' @param n Integer specifying the number of top markers to print per cluster
#' @param cluster_to_print Character vector specifying which clusters to print.
#'                        If NULL, prints all clusters
#'
#' @return NULL (prints results to console)
#'
#' @examples
#' \dontrun{
#' # Print top 50 markers for all clusters
#' marker_print(markers, n = 50)
#'
#' # Print top 20 markers for specific clusters
#' marker_print(markers, n = 20, cluster_to_print = c("0", "1"))
#' }
#' @export
marker_print_all <- function(markers, n = 100, cluster_to_print = NULL) {
  number_to_print <- n
  if (is.null(cluster_to_print)) {} else {
    if (class(markers) == "list") {
      markers <- markers[cluster_to_print]
    } else {
      markers <- markers[markers$cluster == cluster_to_print, ]
    }
  }
  if (class(markers) == "list") {
    for (i in names(markers)) {
      print(i)
      print(paste(markers[[i]][markers[[i]]$avg_log2FC > 0, ][1:number_to_print, ]$gene, collapse = ", "))
    }
  } else {
    for (i in unique(markers$cluster)) {
      print(i)
      print(paste(markers[markers$cluster == i, ][markers$avg_log2FC > 0, ][1:number_to_print, ]$gene, collapse = ", "))
    }
  }
}

#' @export
marker_print <- function(marker, n = 100, sign = "+") {
  if (sign == "+") {
    print(paste(marker[marker$avg_log2FC > 0, ][1:n, ]$gene, collapse = ", "))
  } else {
    print(paste(marker[marker$avg_log2FC < 0, ][1:n, ]$gene, collapse = ", "))
  }
}

#' 여러 FindMarkers 결과를 종합하여 메타분석 수행
#'
#' @param marker_list FindMarkers 결과 데이터프레임의 리스트.
#'        각 데이터프레임은 'gene', 'p_val', 'avg_log2FC', 'p_val_adj' 컬럼을 포함해야 함.
#' @return 메타분석 결과 데이터프레임.
synthesize_markers <- function(marker_list, p_sig = "p_val_adj") {
  # 모든 마커 결과를 하나의 데이터프레임으로 통합
  all_genes_df <- marker_list %>%
    # 리스트의 각 데이터프레임에 그룹 이름을 부여
    imap(~ .x %>% mutate(group = .y)) %>%
    # 리스트를 하나의 데이터프레임으로 결합
    bind_rows() %>%
    # gene 컬럼이 없는 경우 rownames를 gene 컬럼으로 변환
    rownames_to_column(var = "gene_fallback") %>%
    mutate(gene = if ("gene" %in% names(.)) gene else gene_fallback) %>%
    select(-gene_fallback)

  # 유전자별로 메타분석 수행
  meta_results <- all_genes_df %>%
    group_by(gene) %>%
    summarise(
      # P-value synthesis (Fisher's method)
      # P-value가 0인 경우 작은 값으로 대체하여 log 변환 오류 방지
      meta_p = {
        p_values <- p_val[p_val > 0]
        if (length(p_values) == 0) {
          NA_real_
        } else {
          chisq_stat <- -2 * sum(log(p_values))
          pchisq(chisq_stat, df = 2 * length(p_values), lower.tail = FALSE)
        }
      },

      # Meta-analysis for Log2FC (Averaging)
      meta_log2FC = mean(avg_log2FC, na.rm = TRUE),

      # Standard deviation as a measure of variability
      sd_p = sd(p_val, na.rm = TRUE),
      sd_adj_p = sd(p_val_adj, na.rm = TRUE),
      sd_meta_log2FC = sd(avg_log2FC, na.rm = TRUE),

      # 몇 개의 그룹에서 유의했는지 카운트
      n_groups_significant = sum(!!sym(p_sig) < 0.05, na.rm = TRUE),
      n_groups_total = n(),
      .groups = "drop"
    ) %>%
    # Meta p-value에 대한 FDR 계산
    mutate(meta_adj_p = p.adjust(meta_p, method = "BH")) %>%
    # 컬럼 순서 정리
    select(
      gene, meta_p, meta_adj_p, meta_log2FC,
      sd_p, sd_adj_p, sd_meta_log2FC,
      n_groups_significant, n_groups_total
    ) %>%
    arrange(meta_adj_p, meta_p)

  return(meta_results)
}

#' 여러 FindMarkers 결과의 순위를 종합
#'
#' @param marker_list FindMarkers 결과 데이터프레임의 리스트.
#' @param rank_by 어떤 컬럼을 기준으로 순위를 매길지 지정. ('p_val' 또는 'avg_log2FC')
#' @return 유전자별 종합 순위 정보 데이터프레임.
synthesize_ranks <- function(marker_list, rank_by = "p_val") {
  # 순위 계산 함수 정의
  calculate_ranks <- function(df, col) {
    # p-value는 오름차순, log2FC는 절대값 기준 내림차순으로 순위 부여
    if (col == "p_val") {
      df %>% mutate(rank = rank(!!sym(col), ties.method = "min"))
    } else { # avg_log2FC
      df %>% mutate(rank = rank(-abs(!!sym(col)), ties.method = "min"))
    }
  }

  # 각 그룹별로 순위 계산
  ranked_list <- marker_list %>%
    map(~ calculate_ranks(.x, rank_by))

  # 순위를 종합하여 기하 평균 계산
  rank_synthesis <- ranked_list %>%
    map(~ .x %>% select(gene, rank)) %>%
    reduce(function(df1, df2) full_join(df1, df2, by = "gene", suffix = c(".1", ".2"))) %>%
    rowwise() %>%
    mutate(
      ranks = list(c_across(starts_with("rank"))),
      # 기하 평균 계산 (NA는 제외)
      geometric_mean_rank = exp(mean(log(na.omit(ranks)), na.rm = TRUE))
    ) %>%
    ungroup() %>%
    select(gene, geometric_mean_rank) %>%
    arrange(geometric_mean_rank)

  return(rank_synthesis)
}


# -------- PTML: utility. --------
#' FindAllMarkers 결과에서 Top N 유전자를 LLM용으로 출력
#'
#' @param markers_df FindAllMarkers()의 결과인 데이터프레임
#' @param top_n 각 클러스터별로 추출할 상위 유전자 개수
#' @param cluster_col 'cluster' 정보가 있는 컬럼명 (기본값: "cluster")
#' @param gene_col 'gene' 정보가 있는 컬럼명 (기본값: "gene")
#' @param p_val_col p-value 컬럼명 (기본값: "p_val")
#' @param p_val_adj_col p-value (adjusted) 컬럼명 (기본값: "p_val_adj")
#' @param avg_log2FC_col log-fold change 컬럼명 (기본값: "avg_log2FC")
#'
#' @return 콘솔에 클러스터별 Top N 유전자 목록을 출력합니다.
#'
#' @export
PrintTopMarkersForLLM <- function(markers_df,
                                  top_n = 50,
                                  cluster_col = "cluster",
                                  gene_col = "gene",
                                  p_val_col = "p_val",
                                  p_val_adj_col = "p_val_adj",
                                  avg_log2FC_col = "avg_log2FC") {
  # 1. 필수 컬럼 존재 여부 확인
  required_cols <- c(cluster_col, gene_col, avg_log2FC_col)

  if (!all(required_cols %in% colnames(markers_df))) {
    missing_cols <- required_cols[!required_cols %in% colnames(markers_df)]
    stop(
      "필수 컬럼이 markers_df에 없습니다: ", paste(missing_cols, collapse = ", "),
      "\n'cluster_col', 'gene_col', 'avg_log2FC_col' 인자를 확인하세요."
    )
  }

  # 2. 정렬을 위한 데이터 전처리 (p-value 0 처리)
  markers_processed <- markers_df

  use_pval <- p_val_col %in% colnames(markers_processed)
  use_padj <- p_val_adj_col %in% colnames(markers_processed)

  # p_val_mod 컬럼을 안전하게 추가
  if (use_pval) {
    markers_processed$p_val_mod <- ifelse(markers_processed[[p_val_col]] == 0,
      1e-300, # 0 대신 매우 작은 값
      markers_processed[[p_val_col]]
    )
  }

  # 3. 정렬 기준에 따라 Top N 마커 필터링
  markers_grouped <- markers_processed %>%
    group_by(!!sym(cluster_col))

  # 정렬 로직
  if (use_padj && use_pval) {
    message(paste("Sorting by", p_val_adj_col, ",", p_val_col, ", and", avg_log2FC_col))
    sorted_markers <- markers_grouped %>%
      arrange(!!sym(p_val_adj_col), p_val_mod, desc(!!sym(avg_log2FC_col)), .by_group = TRUE)
  } else if (use_pval) {
    message(paste("Sorting by", p_val_col, "and", avg_log2FC_col))
    sorted_markers <- markers_grouped %>%
      arrange(p_val_mod, desc(!!sym(avg_log2FC_col)), .by_group = TRUE)
  } else {
    message(paste("Sorting by", avg_log2FC_col, "only."))
    sorted_markers <- markers_grouped %>%
      arrange(desc(!!sym(avg_log2FC_col)), .by_group = TRUE)
  }

  final_markers <- sorted_markers %>%
    slice_head(n = top_n) %>%
    ungroup()

  # 4. LLM에 복사/붙여넣기 편한 형태로 콘솔에 출력

  # 클러스터 목록 (숫자형/문자형 모두 정렬되도록)
  clusters <- unique(final_markers[[cluster_col]])
  tryCatch(
    {
      clusters_sorted <- sort(as.numeric(as.character(clusters)))
    },
    warning = function(w) {
      clusters_sorted <<- sort(as.character(clusters))
    }
  )


  cat("--- Top", top_n, "Markers per Cluster ---\n\n")

  for (cl in clusters_sorted) {
    cat("--- Cluster:", cl, "---\n")

    genes <- final_markers %>%
      filter(!!sym(cluster_col) == cl) %>%
      pull(!!sym(gene_col))

    # 쉼표로 구분된 한 줄의 문자열로 출력
    cat(paste(genes, collapse = ", "))
    cat("\n\n") # 클러스터 간 구분을 위한 공백
  }

  message("--- End of Marker List ---")
}

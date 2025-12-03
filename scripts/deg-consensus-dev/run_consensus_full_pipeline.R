# ============================================================================
# Full Pipeline: Cluster-wise DEG Consensus with Covariates
# ============================================================================
# 실행 방법:
#   cd /home/user3/GJC_KDW_250721 && Rscript /home/user3/data_user3/git_repo/_wt/deg-consensus/scripts/deg-consensus-dev/run_consensus_full_pipeline.R
# ============================================================================
# 
# 전체 파이프라인:
# 1. 각 cluster별로 쪼개서 DEG consensus 실행
# 2. cluster별 결과를 aggregate_cluster_deg_consensus로 통합
# 3. Fixed effects: sex (그리고 필요시 GEM, set)
# 4. Random effects: hos_no (patient-level, nebula 방법에서만)
#
# 데이터 확장 전략:
# - level="downsampled": 빠른 테스트
# - level="full": 전체 데이터 분석
# ============================================================================

# ============================================================================
# 0. 설정 및 환경 준비
# ============================================================================

# 환경 로드 (start.R이 이미 로드되어 있다면 생략 가능)
if (!exists(".renv")) {
  if (file.exists("/home/user3/GJC_KDW_250721/start.R")) {
    source("/home/user3/GJC_KDW_250721/start.R")
  }
}

# 함수 로드
devtools::load_all("/home/user3/data_user3/git_repo/mylit/myR")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_limma.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_edger.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_deseq2.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_base.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_dream.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_standardize.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_consensus_analysis.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/run_deg_consensus.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/aggregate_cluster_deg_consensus.R")

# 설정 로드
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/docs-main/config/vars_config.R")

# ============================================================================
# 1. 파라미터 설정
# ============================================================================

# 데이터 레벨 선택: "downsampled" (빠른 테스트) 또는 "full" (전체 분석)
data_level <- "downsampled"  # TODO: "full"로 변경하여 전체 데이터 분석

# 출력 디렉토리
out_dir <- "/data/user3/sobj"
cons_dir <- file.path(out_dir, "consensus")
if (!dir.exists(cons_dir)) {
  dir.create(cons_dir, recursive = TRUE, showWarnings = FALSE)
}
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
file_prefix <- file.path(cons_dir, paste0("deg_consensus_clusterwise_", data_level, "_", timestamp))

# 중간 데이터 저장 옵션
save_intermediate_data <- TRUE  # default=TRUE: pseudobulk 등 중간 데이터 저장

# 메타데이터 설정
config <- get_dataset_config("stroke", level = data_level)
cluster_key <- config$cluster_var
sample_key <- config$patient_var
group_key <- config$target_var
contrast_str <- "2 - 1"

# Fixed effects 설정
# 주의: GEM과 hos_no는 완전포함 관계이므로 동시에 포함 시 공선성 문제 발생
# guide_KR.md에 따르면 hos_no < GEM < set 완전포함 관계
# 따라서 sex만 사용하거나, 필요시 GEM만 사용 (hos_no는 random effect로 처리)
covar_effects <- "sex"  # Fixed effect로 sex 추가
# covar_effects <- c("sex", "GEM")  # 필요시 GEM도 추가 가능 (주의: 공선성 확인 필요)

# Batch 설정 (pseudobulk 방법에서 사용)
# batch_key <- "GEM"  # 일부 방법에서 batch로 사용 가능
batch_key <- NULL  # 공선성 문제 때문에 일단 NULL로 설정

# ============================================================================
# 2. 방법론 선택
# ============================================================================

# 기본 방법론 (pseudobulk 기반)
base_methods <- c(
  "muscat-edgeR",
  "muscat-DESeq2",
  "muscat-limma-voom",
  "muscat-limma-trend",
  "limma-voom",
  "limma-trend",
  "edgeR-LRT",
  "edgeR-QLF",
  "DESeq2-Wald",
  "DESeq2-LRT"
)

# 고급 방법론 (random effects 지원)
# 주의: dream, nebula는 원래 잘 작동하지 않았음
include_dream <- FALSE  # 원래 잘 작동하지 않았음
include_nebula <- TRUE  # 별도로 전체 데이터에서 실행 (클러스터별 쪼개지 않음)
include_nebula_pb <- FALSE

optional_methods <- c(
  if (isTRUE(include_dream)) "dream",
  if (isTRUE(include_nebula)) "nebula",
  if (isTRUE(include_nebula_pb)) "nebula-pb"
)

methods_to_run <- unique(c(base_methods, optional_methods))

# ============================================================================
# 3. 데이터 로드
# ============================================================================

message("=== 데이터 로드 ===")
message(sprintf("데이터 레벨: %s", data_level))
message(sprintf("데이터 경로: %s", config$path))

is5 <- qs::qread(config$path)
message(sprintf("로드 완료: %d cells, %d genes", ncol(is5), nrow(is5)))

# 메타데이터 확인
message("\n=== 메타데이터 확인 ===")
message(sprintf("  - cluster_var (%s): %d levels", cluster_key, length(unique(is5@meta.data[[cluster_key]]))))
message(sprintf("  - target_var (%s): %s", group_key, paste(unique(is5@meta.data[[group_key]]), collapse=", ")))
message(sprintf("  - patient_var (%s): %d unique values", sample_key, length(unique(is5@meta.data[[sample_key]]))))
if (!is.null(covar_effects)) {
  for (covar in covar_effects) {
    if (covar %in% colnames(is5@meta.data)) {
      message(sprintf("  - %s: %s", covar, paste(unique(is5@meta.data[[covar]]), collapse=", ")))
    }
  }
}

# 클러스터 목록 확인
clusters <- unique(is5@meta.data[[cluster_key]])
clusters <- clusters[!is.na(clusters)]
message(sprintf("\n총 %d 개의 클러스터 발견", length(clusters)))
message(sprintf("클러스터 목록: %s", paste(head(clusters, 10), collapse=", ")))
if (length(clusters) > 10) {
  message(sprintf("  ... 외 %d 개", length(clusters) - 10))
}

# ============================================================================
# 4. 클러스터별 DEG Consensus 실행
# ============================================================================

message("\n=== 클러스터별 DEG Consensus 분석 시작 ===")
message(sprintf("총 %d 개의 클러스터에 대해 분석 실행", length(clusters)))
message(sprintf("사용 방법론: %s", paste(methods_to_run, collapse=", ")))
message(sprintf("Fixed effects: %s", paste(covar_effects, collapse=", ")))

deg_consensus_list <- list()
cluster_results_summary <- data.frame(
  cluster = character(),
  n_cells = integer(),
  n_methods_success = integer(),
  n_methods_failed = integer(),
  n_genes = integer(),
  stringsAsFactors = FALSE
)

# 건너뛴 클러스터 추적
skipped_clusters <- list()

for (i in seq_along(clusters)) {
  clust <- clusters[i]
  clust_clean <- make.names(clust)
  
  message(sprintf("\n--- [%d/%d] 클러스터: %s ---", i, length(clusters), clust))
  
  # 클러스터별 데이터 서브셋
  # Seurat subset은 동적 변수명을 직접 지원하지 않으므로, WhichCells 사용
  cell_idx <- which(is5@meta.data[[cluster_key]] == clust)
  if (length(cell_idx) == 0) {
    message(sprintf("  ⚠ 세포가 없어 건너뜁니다."))
    next
  }
  is5_clust <- is5[, cell_idx]
  n_cells_clust <- ncol(is5_clust)
  
  message(sprintf("  세포 수: %d", n_cells_clust))
  
  # 최소 세포 수 체크 (매우 작은 클러스터는 건너뛰기)
  min_cells_absolute <- 10  # 절대 최소값
  if (n_cells_clust < min_cells_absolute) {
    message(sprintf("  ⚠ 세포 수가 너무 적어 건너뜁니다 (최소 %d개 필요)", min_cells_absolute))
    skipped_clusters[[clust_clean]] <- list(
      cluster = clust,
      n_cells = n_cells_clust,
      reason = sprintf("세포 수 부족 (최소 %d개 필요)", min_cells_absolute)
    )
    next
  }
  
  # 클러스터별 DEG Consensus 실행
  result_clust <- tryCatch({
    run_deg_consensus(
      sobj = is5_clust,
      contrast = contrast_str,
      methods = methods_to_run,
      cluster_id = cluster_key,
      sample_id = sample_key,
      group_id = group_key,
      batch_id = batch_key,
      covar_effects = covar_effects,
      remove_na_groups = TRUE,
      verbose = TRUE  # 상세 로그 출력 (muscat 방법론 확인용)
    )
  }, error = function(e) {
    message(sprintf("  ✗ 오류 발생: %s", conditionMessage(e)))
    return(NULL)
  })
  
  if (is.null(result_clust)) {
    message(sprintf("  ✗ 분석 실패"))
    skipped_clusters[[clust_clean]] <- list(
      cluster = clust,
      n_cells = n_cells_clust,
      reason = "run_deg_consensus 실행 실패"
    )
    next
  }
  
  # muscat 방법론 실패 체크 및 fallback
  muscat_methods <- c("muscat-edgeR", "muscat-DESeq2", "muscat-limma-voom", "muscat-limma-trend")
  muscat_failed <- intersect(muscat_methods, result_clust$methods_failed)
  muscat_success <- intersect(muscat_methods, result_clust$methods_run)
  
  if (length(muscat_failed) > 0 && length(muscat_success) == 0) {
    # 모든 muscat 방법론이 실패한 경우
    muscat_error_msg <- if (length(muscat_failed) > 0 && !is.null(result_clust$errors)) {
      result_clust$errors[[muscat_failed[1]]]
    } else {
      "Unknown error"
    }
    
    # muscat 실패 원인 확인
    if (grepl("Specified filtering options result in no genes|no genes in any clusters", muscat_error_msg)) {
      message(sprintf("  ⚠ 모든 muscat 방법론 실패: %s", muscat_error_msg))
      message(sprintf("  → muscat은 그룹별로 충분한 pseudobulk 샘플이 필요합니다."))
      message(sprintf("  → 다른 방법론(%d개)으로 계속 진행합니다.", length(result_clust$methods_run)))
      
      # 건너뛴 것으로 기록하지 않고, muscat만 실패한 것으로 처리
      # (다른 방법론이 성공하면 계속 진행)
    }
  }
  
  # 결과 저장
  deg_consensus_list[[clust_clean]] <- result_clust
  
  # 중간 데이터 저장 (pseudobulk 등)
  if (isTRUE(save_intermediate_data)) {
    # pseudobulk 데이터 추출 시도 (가능한 경우)
    # 일단 메타데이터만 저장
    intermediate_data <- list(
      cluster = clust,
      n_cells = n_cells_clust,
      metadata = is5_clust@meta.data,
      timestamp = Sys.time()
    )
    deg_consensus_list[[clust_clean]]$intermediate_data <- intermediate_data
  }
  
  # 표준화 및 consensus 계산
  if (length(result_clust$methods_run) > 0) {
    message(sprintf("  ✓ 성공한 방법론: %d 개", length(result_clust$methods_run)))
    if (length(result_clust$methods_failed) > 0) {
      message(sprintf("  ⚠ 실패한 방법론: %s", paste(result_clust$methods_failed, collapse=", ")))
    }
    
    # 표준화
    standardized_results <- lapply(result_clust$methods_run, function(m) {
      tryCatch({
        standardize_deg_results(result_clust$results[[m]], m)
      }, error = function(e) {
        return(NULL)
      })
    })
    names(standardized_results) <- result_clust$methods_run
    standardized_results <- standardized_results[!sapply(standardized_results, is.null)]
    
    if (length(standardized_results) > 0) {
      # 행렬 구성
      deg_matrices <- build_deg_matrices(
        standardized_results,
        genes = NULL,
        fdr_threshold = 0.1,
        significance_mode = "pvalue",
        pvalue_threshold = 0.01
      )
      
      # Consensus 계산
      agreement_scores <- compute_agreement_scores(deg_matrices$significance)
      consensus_scores <- compute_consensus_scores(deg_matrices, agreement_scores)
      
      # 결과 저장
      deg_consensus_list[[clust_clean]]$standardized_results <- standardized_results
      deg_consensus_list[[clust_clean]]$deg_matrices <- deg_matrices
      deg_consensus_list[[clust_clean]]$agreement_scores <- agreement_scores
      deg_consensus_list[[clust_clean]]$consensus_scores <- consensus_scores
      
      message(sprintf("  ✓ Consensus 계산 완료: %d 유전자", nrow(consensus_scores)))
    }
    
    # 요약 정보 저장
    cluster_results_summary <- rbind(
      cluster_results_summary,
      data.frame(
        cluster = clust,
        n_cells = n_cells_clust,
        n_methods_success = length(result_clust$methods_run),
        n_methods_failed = length(result_clust$methods_failed),
        n_genes = if (exists("consensus_scores")) nrow(consensus_scores) else 0,
        stringsAsFactors = FALSE
      )
    )
  } else {
    message(sprintf("  ✗ 성공한 방법론이 없습니다"))
    skipped_clusters[[clust_clean]] <- list(
      cluster = clust,
      n_cells = n_cells_clust,
      reason = "모든 방법론 실패"
    )
  }
}

# 중간 결과 저장
message("\n=== 중간 결과 저장 ===")
qs::qsave(deg_consensus_list, paste0(file_prefix, "_clusterwise_results.qs"))
qs::qsave(cluster_results_summary, paste0(file_prefix, "_cluster_summary.qs"))

# 건너뛴 클러스터 정보 저장
if (length(skipped_clusters) > 0) {
  skipped_df <- do.call(rbind, lapply(skipped_clusters, function(x) {
    data.frame(
      cluster = x$cluster,
      n_cells = x$n_cells,
      reason = x$reason,
      stringsAsFactors = FALSE
    )
  }))
  qs::qsave(skipped_df, paste0(file_prefix, "_skipped_clusters.qs"))
  message(sprintf("\n⚠ 건너뛴 클러스터: %d 개", length(skipped_clusters)))
  print(skipped_df)
} else {
  message("\n✓ 모든 클러스터 분석 완료")
}

message(sprintf("클러스터별 분석 완료: %d 개 클러스터", length(deg_consensus_list)))
print(cluster_results_summary)

# ============================================================================
# 5. 클러스터별 결과 통합
# ============================================================================

message("\n=== 클러스터별 결과 통합 ===")

# 각 클러스터의 consensus_scores 추출
deg_list_for_aggregation <- lapply(names(deg_consensus_list), function(clust_name) {
  result_clust <- deg_consensus_list[[clust_name]]
  
  if (!is.null(result_clust$consensus_scores) && nrow(result_clust$consensus_scores) > 0) {
    df <- result_clust$consensus_scores %>%
      mutate(
        cluster = clust_name,
        p_val = meta_p,
        p_val_adj = meta_p_adj,
        avg_log2FC = mean_beta
      )
    return(df)
  } else {
    return(NULL)
  }
})

names(deg_list_for_aggregation) <- names(deg_consensus_list)
deg_list_for_aggregation <- deg_list_for_aggregation[!sapply(deg_list_for_aggregation, is.null)]

if (length(deg_list_for_aggregation) == 0) {
  stop("통합할 클러스터별 결과가 없습니다.")
}

message(sprintf("통합 대상 클러스터: %d 개", length(deg_list_for_aggregation)))

# 클러스터 가중치 계산 (선택적)
# inverse_variance 기반 가중치 계산
compute_cluster_weights <- function(deg_list, cluster_names) {
  cluster_variances <- lapply(cluster_names, function(clust) {
    df <- deg_list[[clust]]
    if (is.null(df) || nrow(df) == 0) return(Inf)
    
    meta_p_valid <- df$meta_p[!is.na(df$meta_p) & df$meta_p > 0 & df$meta_p < 1]
    var_meta_p <- if (length(meta_p_valid) > 1) {
      stats::var(meta_p_valid)
    } else {
      NA_real_
    }
    
    mean_beta_valid <- df$mean_beta[!is.na(df$mean_beta)]
    var_mean_beta <- if (length(mean_beta_valid) > 1) {
      stats::var(mean_beta_valid)
    } else {
      NA_real_
    }
    
    combined_var <- 0
    if (!is.na(var_meta_p)) combined_var <- combined_var + var_meta_p
    if (!is.na(var_mean_beta)) combined_var <- combined_var + var_mean_beta
    if (is.na(var_meta_p) && is.na(var_mean_beta)) return(Inf)
    if (combined_var == 0) combined_var <- 1e-10
    
    return(combined_var)
  })
  
  names(cluster_variances) <- cluster_names
  variances <- unlist(cluster_variances)
  epsilon <- 1e-10
  inv_variances <- 1 / (variances + epsilon)
  inv_variances[is.infinite(inv_variances) | is.na(inv_variances)] <- epsilon
  
  if (sum(inv_variances) > 0) {
    weights <- inv_variances / sum(inv_variances)
  } else {
    weights <- rep(1 / length(cluster_names), length(cluster_names))
  }
  
  names(weights) <- cluster_names
  return(weights)
}

# 클러스터 가중치 계산
cluster_weights <- compute_cluster_weights(deg_list_for_aggregation, names(deg_list_for_aggregation))
message(sprintf("클러스터 가중치 계산 완료: %d 개 클러스터", length(cluster_weights)))
message("가중치 샘플:")
print(head(cluster_weights, 5))

# 클러스터 통합
meta_deg <- aggregate_cluster_deg_consensus(
  deg_list = deg_list_for_aggregation,
  meta_p_method = "stouffer",
  mean_beta_method = "simple_mean",
  cluster_weights = cluster_weights
)

message(sprintf("통합 완료: %d 유전자", nrow(meta_deg)))

# 통합 결과 저장
qs::qsave(meta_deg, paste0(file_prefix, "_meta_deg.qs"))

# ============================================================================
# 6. 최종 결과 요약 및 저장
# ============================================================================

message("\n=== 최종 결과 요약 ===")

# 유의한 유전자 필터링
fdr_threshold <- 0.1
meta_deg_sig <- meta_deg %>%
  filter(!is.na(meta_meta_p_adj), meta_meta_p_adj < fdr_threshold) %>%
  arrange(meta_meta_p_adj)

message(sprintf("전체 유전자: %d", nrow(meta_deg)))
message(sprintf("FDR < %.2f 유전자: %d", fdr_threshold, nrow(meta_deg_sig)))

if (nrow(meta_deg_sig) > 0) {
  message("\n상위 20개 유의 유전자:")
  print(head(meta_deg_sig[, c("gene", "mean_mean_beta", "meta_meta_p", "meta_meta_p_adj", "n_clusters", "concordance")], 20))
}

# ============================================================================
# 7. NEBULA 별도 실행 (전체 데이터, 클러스터별 쪼개지 않음)
# ============================================================================

nebula_result <- NULL
if (isTRUE(include_nebula)) {
  message("\n=== NEBULA 분석 (전체 데이터) ===")
  message("주의: NEBULA는 클러스터별로 쪼개지 않고 전체 데이터에서 실행됩니다.")
  message("      공선성 문제로 인해 GEM은 covar_effects에서 제외됩니다.")
  
  nebula_result <- tryCatch({
    runNEBULA2_v1(
      sobj = is5,  # 전체 데이터
      fixed_effects = c(group_key),
      covar_effects = covar_effects,  # sex만 포함 (GEM 제외)
      patient_col = sample_key,
      offset = "nCount_RNA",
      min_count = 20,  # 더 높은 min_count로 필터링 강화
      max_genes = 5000,
      remove_na_cells = TRUE
    )
  }, error = function(e) {
    message(sprintf("NEBULA 실패: %s", conditionMessage(e)))
    message("제안: min_count를 높이거나, covar_effects를 조정하세요.")
    return(NULL)
  })
  
  if (!is.null(nebula_result)) {
    message("✓ NEBULA 분석 완료")
    # NEBULA 결과 저장
    qs::qsave(nebula_result, paste0(file_prefix, "_nebula_result.qs"))
  } else {
    message("✗ NEBULA 분석 실패")
  }
}

# 최종 결과 저장
final_result <- list(
  config = list(
    data_level = data_level,
    cluster_key = cluster_key,
    sample_key = sample_key,
    group_key = group_key,
    contrast = contrast_str,
    covar_effects = covar_effects,
    methods = methods_to_run,
    nebula_included = isTRUE(include_nebula)
  ),
  cluster_results = deg_consensus_list,
  cluster_summary = cluster_results_summary,
  skipped_clusters = if (length(skipped_clusters) > 0) {
    do.call(rbind, lapply(skipped_clusters, function(x) {
      data.frame(cluster = x$cluster, n_cells = x$n_cells, reason = x$reason, 
                 stringsAsFactors = FALSE)
    }))
  } else {
    NULL
  },
  meta_deg = meta_deg,
  meta_deg_sig = meta_deg_sig,
  nebula_result = nebula_result
)

qs::qsave(final_result, paste0(file_prefix, "_final_result.qs"))

message("\n=== 분석 완료 ===")
message(sprintf("결과 저장 위치: %s", file_prefix))
if (length(skipped_clusters) > 0) {
  message(sprintf("건너뛴 클러스터: %d 개 (자세한 내용은 skipped_clusters 파일 참조)", length(skipped_clusters)))
}


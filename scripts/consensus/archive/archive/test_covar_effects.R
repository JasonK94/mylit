# ============================================================================
# Test Script for covar_effects (sex) functionality
# ============================================================================
# 빠른 테스트를 위해 downsampled 데이터와 일부 방법론만 사용
# ============================================================================
# 실행 방법: cd /home/user3/GJC_KDW_250721 && Rscript /home/user3/data_user3/git_repo/_wt/deg-consensus/scripts/deg-consensus-dev/test_covar_effects.R
# ============================================================================

# 0. 환경 설정 (start.R이 이미 로드되어 있다면 생략 가능)
if (!exists(".renv")) {
  if (file.exists("/home/user3/GJC_KDW_250721/start.R")) {
    source("/home/user3/GJC_KDW_250721/start.R")
  }
}

# 1. 함수 로드
devtools::load_all("/home/user3/data_user3/git_repo/mylit/myR")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_limma.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_edger.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_deseq2.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_base.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_dream.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_standardize.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_consensus_analysis.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/run_deg_consensus.R")

# 2. 설정 로드
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/docs-main/config/vars_config.R")
config <- get_dataset_config("stroke", level = "downsampled")

# 3. 데이터 로드
message("데이터 로드 중...")
is5 <- qs::qread(config$path)
message(sprintf("로드 완료: %d cells, %d genes", ncol(is5), nrow(is5)))

# 4. 메타데이터 확인
message("\n메타데이터 확인:")
message(sprintf("  - cluster_var (%s): %d levels", config$cluster_var, length(unique(is5@meta.data[[config$cluster_var]]))))
message(sprintf("  - target_var (%s): %s", config$target_var, paste(unique(is5@meta.data[[config$target_var]]), collapse=", ")))
message(sprintf("  - patient_var (%s): %d unique values", config$patient_var, length(unique(is5@meta.data[[config$patient_var]]))))
message(sprintf("  - sex: %s", paste(unique(is5@meta.data$sex), collapse=", ")))

# 5. 테스트: 하나의 클러스터만 선택 (빠른 테스트)
test_cluster <- unique(is5@meta.data[[config$cluster_var]])[1]
message(sprintf("\n테스트 클러스터: %s", test_cluster))
is5_test <- subset(is5, anno3.scvi == test_cluster)
message(sprintf("테스트 데이터: %d cells", ncol(is5_test)))

# 6. 빠른 테스트를 위한 방법론 선택 (일부만)
test_methods <- c(
  "limma-voom",
  "limma-trend",
  "edgeR-LRT",
  "edgeR-QLF",
  "DESeq2-Wald"
  # muscat 방법론들은 시간이 오래 걸릴 수 있으므로 일단 제외
  # "muscat-edgeR",
  # "muscat-DESeq2",
  # "muscat-limma-voom",
  # "muscat-limma-trend"
)

# 7. DEG Consensus 실행 (covar_effects = "sex" 포함)
message("\n=== DEG Consensus 분석 시작 (covar_effects = 'sex') ===")
result_consensus <- run_deg_consensus(
  sobj = is5_test,
  contrast = "2 - 1",
  methods = test_methods,
  cluster_id = config$cluster_var,
  sample_id = config$patient_var,
  group_id = config$target_var,
  batch_id = NULL,  # 일단 batch는 제외하고 sex만 테스트
  covar_effects = "sex",  # sex를 fixed effect로 추가
  remove_na_groups = TRUE,
  verbose = TRUE
)

# 8. 결과 확인
message("\n=== 결과 요약 ===")
message(sprintf("성공한 방법론: %d 개", length(result_consensus$methods_run)))
message(sprintf("실패한 방법론: %d 개", length(result_consensus$methods_failed)))

if (length(result_consensus$methods_run) > 0) {
  message("\n성공한 방법론:")
  for (m in result_consensus$methods_run) {
    message(sprintf("  - %s", m))
  }
}

if (length(result_consensus$methods_failed) > 0) {
  message("\n실패한 방법론:")
  for (m in result_consensus$methods_failed) {
    message(sprintf("  - %s: %s", m, result_consensus$errors[[m]]))
  }
}

# 9. 결과 표준화 및 간단한 확인
if (length(result_consensus$methods_run) > 0) {
  message("\n=== 결과 표준화 ===")
  standardized_results <- lapply(result_consensus$methods_run, function(m) {
    tryCatch({
      standardize_deg_results(result_consensus$results[[m]], m)
    }, error = function(e) {
      message(sprintf("표준화 실패 (%s): %s", m, conditionMessage(e)))
      return(NULL)
    })
  })
  names(standardized_results) <- result_consensus$methods_run
  standardized_results <- standardized_results[!sapply(standardized_results, is.null)]
  
  if (length(standardized_results) > 0) {
    message(sprintf("표준화 성공: %d 개 방법론", length(standardized_results)))
    
    # 각 방법론별 결과 요약
    for (m in names(standardized_results)) {
      res <- standardized_results[[m]]
      message(sprintf("\n%s:", m))
      message(sprintf("  - 총 유전자: %d", nrow(res)))
      if ("pvalue" %in% colnames(res)) {
        sig_count <- sum(res$pvalue < 0.05, na.rm = TRUE)
        message(sprintf("  - p < 0.05: %d", sig_count))
      }
      if ("FDR" %in% colnames(res) || "pvalue_adj" %in% colnames(res)) {
        fdr_col <- if ("FDR" %in% colnames(res)) "FDR" else "pvalue_adj"
        sig_fdr <- sum(res[[fdr_col]] < 0.05, na.rm = TRUE)
        message(sprintf("  - FDR < 0.05: %d", sig_fdr))
      }
    }
  }
}

message("\n=== 테스트 완료 ===")


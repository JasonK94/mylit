# LDS 디버깅 및 테스트 스크립트
# p-value 문제와 필터링 유연성 확인

cat("=== LDS 디버깅 테스트 시작 ===\n\n")

# 1. 패키지 로드
cat("1. 패키지 로드 중...\n")
suppressPackageStartupMessages({
  library(qs)
  library(limma)
  library(edgeR)
  library(lme4)
  library(BiocParallel)
  library(sva)
  library(dplyr)
  library(Seurat)
  library(variancePartition)
})
cat("✓ 패키지 로드 완료\n\n")

# 2. myR 패키지 로드
cat("2. myR 패키지 로드 중...\n")
devtools::load_all("/home/user3/data_user3/git_repo/_wt/lds/myR", quiet = FALSE)
cat("✓ LDS 함수 확인 완료\n\n")

# 3. 데이터 로드
cat("3. 데이터 로드 중...\n")
data_path <- "/data/user3/sobj/IS_scvi_251107_ds2500.qs"
is5s <- qs::qread(data_path)
is5s_test <- is5s
is5s_test$g3_clean <- as.numeric(as.character(is5s_test$g3))
is5s_test <- is5s_test[, !is.na(is5s_test$g3_clean)]
cat("✓ 데이터 로드 완료: ", ncol(is5s_test), " 셀, ", nrow(is5s_test), " 유전자\n\n")

# 4. 설정
output_dir <- "/data/user3/sobj/lds_intermediate"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
prefix <- "lds_debug"
formula <- ~ g3_clean + (1|hos_no) + (1|GEM)

cat("4. 설정:\n")
cat("  - 포뮬러:", as.character(formula), "\n")
cat("  - 필터링: min.count=5, min.prop=0.05 (완화)\n")
cat("  - 중간 결과 저장:", output_dir, "\n\n")

# ============================================================================
# 단계별 실행 및 검증
# ============================================================================

cat("=== 단계별 실행 시작 ===\n\n")

# --- 단계 1: 데이터 추출 ---
cat(">>> 단계 1: 데이터 추출\n")
step1 <- lds_01_extract_data(
  sobj = is5s_test,
  save_intermediate = TRUE,
  output_dir = output_dir,
  prefix = prefix
)
cat("✓ 완료\n\n")

# --- 단계 1b: NA 필터링 ---
cat(">>> 단계 1b: NA 필터링\n")
step1b <- lds_01b_filter_na(
  counts_matrix = step1$counts_matrix,
  meta.data = step1$meta.data,
  formula = formula,
  remove_na = TRUE,
  save_intermediate = TRUE,
  output_dir = output_dir,
  prefix = prefix
)
cat("✓ 완료: ", step1b$n_removed, " 개 셀 제거\n\n")

# --- 단계 2: 포뮬러 파싱 ---
cat(">>> 단계 2: 포뮬러 파싱\n")
step2 <- lds_02_parse_formula(
  formula = formula,
  save_intermediate = TRUE,
  output_dir = output_dir,
  prefix = prefix
)
cat("✓ 완료\n\n")

# --- 단계 3: DGEList 전처리 (필터링 유연성 테스트) ---
cat(">>> 단계 3: DGEList 전처리 (필터링 옵션 테스트)\n")
step3 <- lds_03_preprocess_dge(
  counts_matrix = step1b$counts_matrix,
  meta.data = step1b$meta.data,
  fixed_effects_formula = step2$fixed_effects_formula,
  min.count = 5,           # 완화된 필터링
  min.total.count = 10,
  min.prop = 0.05,         # 5%로 완화
  large.n = 10,
  save_intermediate = TRUE,
  output_dir = output_dir,
  prefix = prefix
)
cat("✓ 완료: ", step3$n_genes_after, " / ", step3$n_genes_before, " 유전자 통과\n")
cat("  필터링 비율: ", round(100 * step3$n_genes_after / step3$n_genes_before, 1), "%\n\n")

# --- 단계 4: SVA 실행 ---
cat(">>> 단계 4: SVA 실행\n")
step4 <- lds_04_run_sva(
  dge = step3$dge,
  meta.data = step1b$meta.data,
  fixed_effects_formula = step2$fixed_effects_formula,
  n_sv = NULL,
  sv_var_cutoff = 0.5,
  save_intermediate = TRUE,
  output_dir = output_dir,
  prefix = prefix
)
cat("✓ 완료: ", step4$n_sv_final, " 개 SV 사용\n\n")

# --- 단계 5: 최종 포뮬러 생성 ---
cat(">>> 단계 5: 최종 포뮬러 생성\n")
step5 <- lds_05_build_final_formula(
  original_formula = step2$original_formula,
  svs_final = step4$svs_final,
  save_intermediate = TRUE,
  output_dir = output_dir,
  prefix = prefix
)
cat("✓ 완료\n\n")

# --- 단계 6: Dream 실행 (p-value 문제 확인) ---
cat(">>> 단계 6: Dream 실행\n")
step6 <- lds_06_run_dream(
  dge = step3$dge,
  final_formula = step5$final_formula,
  meta.data = step4$meta.data_with_sv,
  n_cores = 4,
  save_intermediate = TRUE,
  output_dir = output_dir,
  prefix = prefix
)
cat("✓ 완료\n\n")

# ============================================================================
# p-value 검증
# ============================================================================

cat("=== p-value 검증 ===\n\n")

fit <- step6$fit_ebayes

# fit 객체 구조 확인
cat("1. fit 객체 구조:\n")
cat("  - 클래스:", class(fit), "\n")
cat("  - coefficients 컬럼:", paste(colnames(fit$coefficients), collapse = ", "), "\n")
if (!is.null(fit$p.value)) {
  cat("  - p.value 컬럼:", paste(colnames(fit$p.value), collapse = ", "), "\n")
  # p.value가 NA인지 확인
  for (col in colnames(fit$p.value)) {
    n_na <- sum(is.na(fit$p.value[, col]))
    n_total <- nrow(fit$p.value)
    cat("  - ", col, ": NA=", n_na, "/", n_total, " (", round(100*n_na/n_total, 1), "%)\n", sep = "")
  }
} else {
  cat("  - p.value: NULL\n")
}

# topTable 테스트
cat("\n2. topTable 테스트:\n")
if (length(step6$fixed_vars) > 0) {
  coef_name <- step6$fixed_vars[1]
  cat("  coef 사용:", coef_name, "\n")
  
  # coef 없이 시도
  cat("  a) coef 없이 topTable:\n")
  top_no_coef <- tryCatch({
    limma::topTable(fit, number = 5)
  }, error = function(e) {
    cat("    오류:", conditionMessage(e), "\n")
    return(NULL)
  })
  
  if (!is.null(top_no_coef)) {
    if ("P.Value" %in% colnames(top_no_coef)) {
      n_valid <- sum(!is.na(top_no_coef$P.Value))
      cat("    P.Value 컬럼 존재, 유효한 값:", n_valid, "/", nrow(top_no_coef), "\n")
      if (n_valid > 0) {
        cat("    첫 번째 P.Value:", top_no_coef$P.Value[1], "\n")
      }
    } else {
      cat("    P.Value 컬럼 없음\n")
    }
  }
  
  # coef 지정하여 시도
  cat("  b) coef 지정하여 topTable:\n")
  top_with_coef <- tryCatch({
    limma::topTable(fit, coef = coef_name, number = 5)
  }, error = function(e) {
    cat("    오류:", conditionMessage(e), "\n")
    return(NULL)
  })
  
  if (!is.null(top_with_coef)) {
    if ("P.Value" %in% colnames(top_with_coef)) {
      n_valid <- sum(!is.na(top_with_coef$P.Value))
      cat("    P.Value 컬럼 존재, 유효한 값:", n_valid, "/", nrow(top_with_coef), "\n")
      if (n_valid > 0) {
        cat("    첫 번째 P.Value:", top_with_coef$P.Value[1], "\n")
        cat("    첫 번째 adj.P.Val:", top_with_coef$adj.P.Val[1], "\n")
        cat("\n    상위 3개 유전자:\n")
        print(head(top_with_coef, 3))
      }
    } else {
      cat("    P.Value 컬럼 없음\n")
    }
  }
  
  # variancePartition::topTable 시도
  cat("  c) variancePartition::topTable:\n")
  top_vp <- tryCatch({
    variancePartition::topTable(fit, coef = coef_name, number = 5)
  }, error = function(e) {
    cat("    오류:", conditionMessage(e), "\n")
    return(NULL)
  })
  
  if (!is.null(top_vp)) {
    if ("P.Value" %in% colnames(top_vp)) {
      n_valid <- sum(!is.na(top_vp$P.Value))
      cat("    P.Value 컬럼 존재, 유효한 값:", n_valid, "/", nrow(top_vp), "\n")
    } else {
      cat("    P.Value 컬럼 없음\n")
    }
  }
  
} else {
  cat("  고정 효과 변수가 없어 topTable 테스트를 스킵합니다.\n")
}

# ============================================================================
# 결과 저장
# ============================================================================

cat("\n=== 결과 저장 ===\n\n")

result <- list(
  fit = step6$fit_ebayes,
  voom = step6$v_dream,
  sva_obj = step4$sva_obj,
  svs_used = step4$svs_final,
  final_formula = step5$final_formula,
  dge = step3$dge,
  fixed_vars = step6$fixed_vars,
  contrast_applied = step6$contrast_applied,
  n_genes_filtered = step3$n_genes_after,
  n_genes_total = step3$n_genes_before
)

output_path <- "/data/user3/sobj/test_lds_debug_result.qs"
qs::qsave(result, output_path)
cat("✓ 최종 결과 저장:", output_path, "\n")

cat("\n=== 테스트 완료 ===\n")


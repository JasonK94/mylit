# LDS 모듈화 함수 테스트 및 디버깅 스크립트
# 각 단계별로 실행하고 중간 결과를 저장하여 디버깅 가능
# R 세션에서 직접 실행: source("scripts/lds/test_lds_modular.R")

cat("=== LDS 모듈화 함수 테스트 시작 ===\n\n")

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
})
cat("✓ 패키지 로드 완료\n\n")

# 2. myR 패키지 로드
cat("2. myR 패키지 로드 중...\n")
# 워크트리에서 직접 로드
devtools::load_all("/home/user3/data_user3/git_repo/_wt/lds/myR", quiet = FALSE)

# 모듈화된 함수 확인
required_functions <- c("LDS", "lds_01_extract_data", "lds_01b_filter_na", 
                        "lds_02_parse_formula", "lds_03_preprocess_dge", 
                        "lds_04_run_sva", "lds_05_build_final_formula", 
                        "lds_06_run_dream", "lds_07_analyze_sva_correlation")

missing_functions <- required_functions[!sapply(required_functions, exists)]
if (length(missing_functions) > 0) {
  stop("필수 함수를 찾을 수 없습니다: ", paste(missing_functions, collapse = ", "))
}
cat("✓ 모든 LDS 함수 확인 완료\n\n")

# 3. 설정
cat("3. 설정:\n")
output_dir <- "/data/user3/sobj/lds_intermediate"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("  중간 결과 디렉터리 생성:", output_dir, "\n")
} else {
  cat("  중간 결과 디렉터리:", output_dir, "\n")
}
prefix <- "lds_test"
save_intermediate <- TRUE  # 중간 결과 저장 여부
cat("  접두사:", prefix, "\n")
cat("  중간 결과 저장:", ifelse(save_intermediate, "예", "아니오"), "\n\n")

# 4. 데이터 로드
cat("4. 데이터 로드 중...\n")
data_path <- "/data/user3/sobj/IS_scvi_251107_ds2500.qs"
if (!file.exists(data_path)) {
  stop("데이터 파일을 찾을 수 없습니다: ", data_path)
}
is5s <- qs::qread(data_path)
cat("✓ 데이터 로드 완료\n")
cat("  - 셀 수:", ncol(is5s), "\n")
cat("  - 유전자 수:", nrow(is5s), "\n\n")

# 5. 데이터 전처리
cat("5. 데이터 전처리 중...\n")
is5s_test <- is5s
is5s_test$g3_clean <- as.numeric(as.character(is5s_test$g3))
n_before <- ncol(is5s_test)
is5s_test <- is5s_test[, !is.na(is5s_test$g3_clean)]
n_after <- ncol(is5s_test)
cat("✓ 전처리 완료\n")
cat("  - NA 제거: ", n_before, " -> ", n_after, " 셀\n", sep = "")
cat("  - g3_clean 값:", paste(sort(unique(is5s_test$g3_clean)), collapse = ", "), "\n")
cat("  - hos_no 샘플 수:", length(unique(is5s_test@meta.data$hos_no)), "\n")
cat("  - GEM 배치 수:", length(unique(is5s_test@meta.data$GEM)), "\n\n")

# 6. 포뮬러 설정
formula <- ~ g3_clean + (1|hos_no) + (1|GEM)
cat("6. 포뮬러:", as.character(formula), "\n\n")

# ============================================================================
# 단계별 실행 (각 단계를 개별적으로 테스트 가능)
# ============================================================================

cat("=== 단계별 실행 시작 ===\n\n")

# --- 단계 1: 데이터 추출 ---
cat(">>> 단계 1: 데이터 추출\n")
step1 <- tryCatch({
  lds_01_extract_data(
    sobj = is5s_test,
    layer = "counts",
    save_intermediate = save_intermediate,
    output_dir = output_dir,
    prefix = prefix
  )
}, error = function(e) {
  cat("❌ 오류:", conditionMessage(e), "\n")
  traceback()
  stop("단계 1 실패")
})
cat("✓ 완료: ", step1$n_cells, " 셀, ", step1$n_genes, " 유전자\n\n")

# --- 단계 1b: NA 필터링 ---
cat(">>> 단계 1b: NA 필터링\n")
step1b <- tryCatch({
  lds_01b_filter_na(
    counts_matrix = step1$counts_matrix,
    meta.data = step1$meta.data,
    formula = formula,
    remove_na = TRUE,
    save_intermediate = save_intermediate,
    output_dir = output_dir,
    prefix = prefix
  )
}, error = function(e) {
  cat("❌ 오류:", conditionMessage(e), "\n")
  traceback()
  stop("단계 1b 실패")
})
cat("✓ 완료: ", step1b$n_removed, " 개 셀 제거\n\n")

# --- 단계 2: 포뮬러 파싱 ---
cat(">>> 단계 2: 포뮬러 파싱\n")
step2 <- tryCatch({
  lds_02_parse_formula(
    formula = formula,
    save_intermediate = save_intermediate,
    output_dir = output_dir,
    prefix = prefix
  )
}, error = function(e) {
  cat("❌ 오류:", conditionMessage(e), "\n")
  traceback()
  stop("단계 2 실패")
})
cat("✓ 완료: 고정 효과 포뮬러:", as.character(step2$fixed_effects_formula), "\n\n")

# --- 단계 3: DGEList 전처리 ---
cat(">>> 단계 3: DGEList 전처리\n")
step3 <- tryCatch({
  lds_03_preprocess_dge(
    counts_matrix = step1b$counts_matrix,
    meta.data = step1b$meta.data,
    fixed_effects_formula = step2$fixed_effects_formula,
    min.count = 5,
    min.total.count = 10,
    min.prop = 0.05,
    large.n = 10,
    save_intermediate = save_intermediate,
    output_dir = output_dir,
    prefix = prefix
  )
}, error = function(e) {
  cat("❌ 오류:", conditionMessage(e), "\n")
  traceback()
  stop("단계 3 실패")
})
cat("✓ 완료: ", step3$n_genes_after, " / ", step3$n_genes_before, " 유전자 통과\n\n")

# --- 단계 4: SVA 실행 ---
cat(">>> 단계 4: SVA 실행\n")
step4 <- tryCatch({
  lds_04_run_sva(
    dge = step3$dge,
    meta.data = step1b$meta.data,
    fixed_effects_formula = step2$fixed_effects_formula,
    n_sv = NULL,
    sv_var_cutoff = 0.5,
    save_intermediate = save_intermediate,
    output_dir = output_dir,
    prefix = prefix
  )
}, error = function(e) {
  cat("❌ 오류:", conditionMessage(e), "\n")
  traceback()
  stop("단계 4 실패")
})
cat("✓ 완료: ", step4$n_sv_final, " 개 SV 사용\n\n")

# --- 단계 5: 최종 포뮬러 생성 ---
cat(">>> 단계 5: 최종 포뮬러 생성\n")
step5 <- tryCatch({
  lds_05_build_final_formula(
    original_formula = step2$original_formula,
    svs_final = step4$svs_final,
    save_intermediate = save_intermediate,
    output_dir = output_dir,
    prefix = prefix
  )
}, error = function(e) {
  cat("❌ 오류:", conditionMessage(e), "\n")
  traceback()
  stop("단계 5 실패")
})
cat("✓ 완료: ", step5$final_formula_str, "\n\n")

# --- 단계 6: Dream 실행 ---
cat(">>> 단계 6: Dream 실행\n")
step6 <- tryCatch({
  lds_06_run_dream(
    dge = step3$dge,
    final_formula = step5$final_formula,
    meta.data = step4$meta.data_with_sv,
    n_cores = 4,
    save_intermediate = save_intermediate,
    output_dir = output_dir,
    prefix = prefix
  )
}, error = function(e) {
  cat("❌ 오류:", conditionMessage(e), "\n")
  traceback()
  stop("단계 6 실패")
})
cat("✓ 완료: Contrast 적용", ifelse(step6$contrast_applied, "예", "아니오"), "\n\n")

# --- 단계 7: SVA 상관관계 분석 ---
cat(">>> 단계 7: SVA 상관관계 분석\n")
step7 <- tryCatch({
  lds_07_analyze_sva_correlation(
    svs_final = step4$svs_final,
    meta.data = step4$meta.data_with_sv,
    plot_sva_correlation = TRUE,
    save_intermediate = save_intermediate,
    output_dir = output_dir,
    prefix = prefix
  )
}, error = function(e) {
  cat("❌ 오류:", conditionMessage(e), "\n")
  traceback()
  stop("단계 7 실패")
})
if (!is.null(step7)) {
  cat("✓ 완료: ", length(step7$metadata_vars), " 개 메타데이터 변수와 상관관계 계산\n\n")
} else {
  cat("✓ 완료: 상관관계 분석 스킵\n\n")
}

# ============================================================================
# 결과 확인 및 저장
# ============================================================================

cat("=== 결과 확인 ===\n\n")

# 최종 결과 객체 구성
result <- list(
  fit = step6$fit_ebayes,
  voom = step6$v_dream,
  sva_obj = step4$sva_obj,
  svs_used = step4$svs_final,
  final_formula = step5$final_formula,
  dge = step3$dge,
  fixed_vars = step6$fixed_vars,
  contrast_applied = step6$contrast_applied
)

if (!is.null(step7)) {
  result$sv_correlation <- step7
}

cat("7. TopTable 테스트:\n")

# dream의 경우 coef를 지정해야 p-value가 나옴
coef_name <- NULL
if (length(step6$fixed_vars) > 0) {
  coef_name <- step6$fixed_vars[1]
  cat("  coef 사용:", coef_name, "\n")
  top_genes <- tryCatch({
    limma::topTable(result$fit, coef = coef_name, number = 10)
  }, error = function(e) {
    cat("  coef 지정 실패, coef 없이 시도\n")
    limma::topTable(result$fit, number = 10)
  })
} else {
  cat("  coef 없이 topTable 실행\n")
  top_genes <- limma::topTable(result$fit, number = 10)
}

cat("  상위 5개 유전자:\n")
print(head(top_genes, 5))

# p-value 확인
if ("P.Value" %in% colnames(top_genes)) {
  n_valid_p <- sum(!is.na(top_genes$P.Value))
  cat("  ✓ P.Value 컬럼 존재\n")
  cat("  NA가 아닌 P.Value 개수:", n_valid_p, "/", nrow(top_genes), "\n")
} else {
  cat("  ⚠ P.Value 컬럼 없음\n")
}

# 결과 저장
cat("\n8. 최종 결과 저장 중...\n")
output_path <- "/data/user3/sobj/test_lds_result_modular.qs"
qs::qsave(result, output_path)
cat("  ✓ 결과 저장:", output_path, "\n")

# 전체 TopTable 저장
if (!is.null(coef_name)) {
  top_all <- limma::topTable(result$fit, coef = coef_name, number = Inf)
} else {
  top_all <- limma::topTable(result$fit, number = Inf)
}
output_csv <- "/data/user3/sobj/test_lds_topTable_modular.csv"
write.csv(top_all, output_csv, row.names = TRUE)
cat("  ✓ TopTable 저장:", output_csv, "\n")

# SVA 상관관계 확인
if (!is.null(result$sv_correlation)) {
  cat("\n9. SVA 상관관계 확인:\n")
  cor_mat <- result$sv_correlation$correlation_matrix
  cat("  상관관계 행렬 크기:", nrow(cor_mat), "x", ncol(cor_mat), "\n")
  cat("  메타데이터 변수:", length(result$sv_correlation$metadata_vars), "개\n")
  if (length(cor_mat) > 0 && !all(is.na(cor_mat))) {
    abs_cor <- abs(cor_mat)
    max_idx <- which(abs_cor == max(abs_cor, na.rm = TRUE), arr.ind = TRUE)
    if (nrow(max_idx) > 0) {
      cat("  상위 상관관계 (절댓값 기준):\n")
      for (i in 1:min(5, nrow(max_idx))) {
        sv_name <- rownames(cor_mat)[max_idx[i, 1]]
        meta_name <- colnames(cor_mat)[max_idx[i, 2]]
        cor_val <- cor_mat[max_idx[i, 1], max_idx[i, 2]]
        cat(sprintf("    %s <-> %s: %.3f\n", sv_name, meta_name, cor_val))
      }
    }
  }
}

cat("\n=== 테스트 완료 ===\n")
cat("중간 결과 저장 위치:", output_dir, "\n")
cat("각 단계별 결과 파일:\n")
cat("  - ", prefix, "_01_extract_data.qs\n", sep = "")
cat("  - ", prefix, "_01b_filter_na.qs\n", sep = "")
cat("  - ", prefix, "_02_parse_formula.qs\n", sep = "")
cat("  - ", prefix, "_03_preprocess_dge.qs\n", sep = "")
cat("  - ", prefix, "_04_run_sva.qs\n", sep = "")
cat("  - ", prefix, "_05_build_final_formula.qs\n", sep = "")
cat("  - ", prefix, "_06_run_dream.qs\n", sep = "")
cat("  - ", prefix, "_07_sva_correlation.qs\n", sep = "")


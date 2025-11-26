# LDS 함수 테스트 스크립트
# 
# 사용법:
#   source("/home/user3/data_user3/git_repo/_wt/lds/scripts/lds/test_lds.R")
#
# 또는 R 세션에서:
#   devtools::load_all("/home/user3/data_user3/git_repo/mylit/myR")
#   source("/home/user3/data_user3/git_repo/_wt/lds/scripts/lds/test_lds.R")

# ============================================================================
# 환경 설정
# ============================================================================

# 패키지 로드
cat("=== 패키지 로드 ===\n")

# 필수 패키지 확인
required_packages <- c("qs", "limma", "edgeR", "lme4", "BiocParallel", "sva", "dplyr", "Seurat")
missing_packages <- c()

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    missing_packages <- c(missing_packages, pkg)
    cat("⚠ 패키지 누락:", pkg, "\n")
  } else {
    cat("✓ 패키지 확인:", pkg, "\n")
  }
}

if (length(missing_packages) > 0) {
  cat("\n누락된 패키지 설치 필요:\n")
  cat("BiocManager::install(c('", paste(missing_packages, collapse = "', '"), "'))\n", sep = "")
  stop("필수 패키지가 설치되지 않았습니다.")
}

# 패키지 로드
library(qs)
library(limma)
library(edgeR)
library(lme4)
library(BiocParallel)
library(sva)
library(dplyr)
library(Seurat)

cat("\n=== myR 패키지 로드 ===\n")
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::load_all("/home/user3/data_user3/git_repo/mylit/myR")

# LDS 함수 확인
if (!exists("LDS")) {
  cat("\n⚠ LDS 함수를 찾을 수 없습니다. test_working.R을 확인하세요.\n")
  cat("직접 소스 로드 시도...\n")
  source("/home/user3/data_user3/git_repo/mylit/myR/R/test_working.R")
  
  if (!exists("LDS")) {
    stop("LDS 함수를 로드할 수 없습니다.")
  }
}
cat("✓ LDS 함수 확인 완료\n")

# ============================================================================
# 데이터 로드
# ============================================================================

cat("=== 데이터 로드 ===\n")

# 다운샘플링된 데이터 (테스트용)
data_path_ds <- "/data/user3/sobj/IS_scvi_251107_ds2500.qs"
data_path_full <- "/data/user3/sobj/IS_scvi_251107.qs"

if (file.exists(data_path_ds)) {
  cat("다운샘플링 데이터 로드 중...\n")
  is5s <- qs::qread(data_path_ds)
  cat("셀 수:", ncol(is5s), "\n")
  cat("유전자 수:", nrow(is5s), "\n")
} else {
  warning("다운샘플링 데이터를 찾을 수 없습니다: ", data_path_ds)
  is5s <- NULL
}

if (file.exists(data_path_full)) {
  cat("전체 데이터 로드 중...\n")
  is5 <- qs::qread(data_path_full)
  cat("셀 수:", ncol(is5), "\n")
  cat("유전자 수:", nrow(is5), "\n")
} else {
  warning("전체 데이터를 찾을 수 없습니다: ", data_path_full)
  is5 <- NULL
}

# ============================================================================
# 메타데이터 확인
# ============================================================================

cat("\n=== 메타데이터 확인 ===\n")

if (!is.null(is5s)) {
  cat("\n다운샘플링 데이터 메타데이터:\n")
  cat("컬럼:", paste(head(colnames(is5s@meta.data), 10), collapse = ", "), "...\n")
  
  if ("hos_no" %in% colnames(is5s@meta.data)) {
    cat("샘플 수 (hos_no):", length(unique(is5s@meta.data$hos_no)), "\n")
  }
  
  if ("g3" %in% colnames(is5s@meta.data)) {
    cat("g3 분포:\n")
    print(table(is5s@meta.data$g3, useNA = "ifany"))
  }
  
  if ("GEM" %in% colnames(is5s@meta.data)) {
    cat("GEM 분포:\n")
    print(table(is5s@meta.data$GEM, useNA = "ifany"))
  }
}

# ============================================================================
# 데이터 전처리
# ============================================================================

cat("\n=== 데이터 전처리 ===\n")

if (!is.null(is5s)) {
  # g3에서 "NA" 제거
  is5s_test <- is5s
  
  # g3를 숫자로 변환
  is5s_test$g3_clean <- as.numeric(as.character(is5s_test$g3))
  
  # NA 제거
  n_before <- ncol(is5s_test)
  is5s_test <- is5s_test[, !is.na(is5s_test$g3_clean)]
  n_after <- ncol(is5s_test)
  
  cat("NA 제거: ", n_before, " -> ", n_after, " 셀\n", sep = "")
  
  # g3_clean 분포 확인
  cat("g3_clean 분포:\n")
  print(table(is5s_test$g3_clean, useNA = "ifany"))
  
  # hos_no와 GEM 확인
  if ("hos_no" %in% colnames(is5s_test@meta.data)) {
    cat("샘플 수 (hos_no):", length(unique(is5s_test@meta.data$hos_no)), "\n")
  }
  if ("GEM" %in% colnames(is5s_test@meta.data)) {
    cat("배치 수 (GEM):", length(unique(is5s_test@meta.data$GEM)), "\n")
  }
} else {
  stop("테스트할 데이터가 없습니다.")
}

# ============================================================================
# LDS 테스트 1: 기본 테스트 (자동 SV 결정)
# ============================================================================

cat("\n=== LDS 테스트 1: 기본 테스트 (자동 SV 결정) ===\n")

test_lds_basic <- function() {
  cat("포뮬러: ~ g3_clean + (1|hos_no) + (1|GEM)\n")
  cat("SV: 자동 결정 (sv_var_cutoff = 0.5)\n")
  
  # 데이터 검증
  cat("\n데이터 검증:\n")
  cat("- 셀 수:", ncol(is5s_test), "\n")
  cat("- 유전자 수:", nrow(is5s_test), "\n")
  cat("- g3_clean 값:", paste(unique(is5s_test$g3_clean), collapse = ", "), "\n")
  cat("- hos_no 샘플 수:", length(unique(is5s_test@meta.data$hos_no)), "\n")
  cat("- GEM 배치 수:", length(unique(is5s_test@meta.data$GEM)), "\n")
  
  # 필수 변수 확인
  required_vars <- c("g3_clean", "hos_no", "GEM")
  missing_vars <- required_vars[!required_vars %in% colnames(is5s_test@meta.data)]
  if (length(missing_vars) > 0) {
    stop("필수 변수가 없습니다: ", paste(missing_vars, collapse = ", "))
  }
  
  result <- tryCatch({
    cat("\nLDS 실행 시작...\n")
    LDS(
      sobj = is5s_test,
      formula = ~ g3_clean + (1|hos_no) + (1|GEM),
      n_sv = NULL,  # 자동 결정
      sv_var_cutoff = 0.5,
      n_cores = 4
    )
  }, error = function(e) {
    cat("\n❌ 오류 발생:\n")
    cat("  메시지:", conditionMessage(e), "\n")
    cat("  클래스:", class(e), "\n")
    if (inherits(e, "simpleError")) {
      cat("  호출:\n")
      print(conditionCall(e))
    }
    traceback()
    return(NULL)
  })
  
  if (!is.null(result)) {
    cat("\n✓ LDS 실행 성공\n")
    
    # 결과 구조 확인
    cat("\n결과 구조:\n")
    cat("- fit:", class(result$fit), "\n")
    cat("- voom:", class(result$voom), "\n")
    cat("- sva_obj:", class(result$sva_obj), "\n")
    
    # SV 정보
    if (!is.null(result$svs_used)) {
      cat("- 사용된 SV 개수:", ncol(result$svs_used), "\n")
      cat("- SV 컬럼:", paste(colnames(result$svs_used), collapse = ", "), "\n")
    } else {
      cat("- 사용된 SV: 없음\n")
    }
    
    # 최종 포뮬러
    cat("- 최종 포뮬러:", as.character(result$final_formula), "\n")
    
    # topTable 테스트
    cat("\nTopTable 테스트:\n")
    top_genes <- topTable(result$fit, number = 10)
    cat("상위 10개 유전자:\n")
    print(head(top_genes, 5))
    
    return(result)
  } else {
    return(NULL)
  }
}

result_lds_basic <- test_lds_basic()

# ============================================================================
# LDS 테스트 2: SV 개수 지정
# ============================================================================

if (!is.null(result_lds_basic) && !is.null(result_lds_basic$svs_used)) {
  cat("\n=== LDS 테스트 2: SV 개수 지정 ===\n")
  
  n_sv_to_use <- min(3, ncol(result_lds_basic$svs_used))
  cat("SV 개수:", n_sv_to_use, "\n")
  
  result_lds_manual <- tryCatch({
    LDS(
      sobj = is5s_test,
      formula = ~ g3_clean + (1|hos_no) + (1|GEM),
      n_sv = n_sv_to_use,
      n_cores = 4
    )
  }, error = function(e) {
    cat("오류 발생:", conditionMessage(e), "\n")
    return(NULL)
  })
  
  if (!is.null(result_lds_manual)) {
    cat("✓ SV 개수 지정 테스트 성공\n")
    if (!is.null(result_lds_manual$svs_used)) {
      cat("사용된 SV 개수:", ncol(result_lds_manual$svs_used), "\n")
    }
  }
}

# ============================================================================
# 결과 저장
# ============================================================================

if (!is.null(result_lds_basic)) {
  cat("\n=== 결과 저장 ===\n")
  
  output_path <- "/data/user3/sobj/test_lds_result.qs"
  qs::qsave(result_lds_basic, output_path)
  cat("결과 저장:", output_path, "\n")
  
  # topTable 결과도 저장
  if (exists("result_lds_basic")) {
    top_all <- topTable(result_lds_basic$fit, number = Inf)
    output_csv <- "/data/user3/sobj/test_lds_topTable.csv"
    write.csv(top_all, output_csv, row.names = TRUE)
    cat("TopTable 저장:", output_csv, "\n")
  }
}

# ============================================================================
# 요약
# ============================================================================

cat("\n=== 테스트 완료 ===\n")
cat("기본 테스트:", ifelse(!is.null(result_lds_basic), "성공", "실패"), "\n")

if (!is.null(result_lds_basic)) {
  cat("\n다음 단계:\n")
  cat("1. 결과 확인: qs::qread('/data/user3/sobj/test_lds_result.qs')\n")
  cat("2. TopTable 확인: limma::topTable(result$fit, number = 100)\n")
  cat("3. 전체 데이터로 확장 테스트\n")
}


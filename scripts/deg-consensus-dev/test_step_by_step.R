# ============================================================================
# 단계별 테스트 및 디버깅 스크립트
# ============================================================================
# 각 단계를 개별적으로 테스트하여 문제를 빠르게 발견
# ============================================================================

cat("================================================================================\n")
cat("Multi-Model DEG Consensus Engine - 단계별 테스트\n")
cat("================================================================================\n\n")

# 환경 설정
if (!requireNamespace("qs", quietly = TRUE)) {
  stop("qs 패키지가 필요합니다.")
}

# ============================================================================
# Step 1: 함수 로드
# ============================================================================
cat("Step 1: 함수 로드\n")
cat("--------------------------------------------------------------------------------\n")

load_ok <- TRUE

tryCatch({
  source("/home/user3/data_user3/git_repo/mylit/myR/R/test_analysis.R")
  cat("✓ test_analysis.R (runMUSCAT2_v1, runNEBULA2_v1)\n")
}, error = function(e) {
  cat("✗ test_analysis.R:", conditionMessage(e), "\n")
  load_ok <<- FALSE
})

tryCatch({
  source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_limma.R")
  cat("✓ deg_methods_limma.R\n")
}, error = function(e) {
  cat("✗ deg_methods_limma.R:", conditionMessage(e), "\n")
  load_ok <<- FALSE
})

tryCatch({
  source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_edger.R")
  cat("✓ deg_methods_edger.R\n")
}, error = function(e) {
  cat("✗ deg_methods_edger.R:", conditionMessage(e), "\n")
  load_ok <<- FALSE
})

tryCatch({
  source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_deseq2.R")
  cat("✓ deg_methods_deseq2.R\n")
}, error = function(e) {
  cat("✗ deg_methods_deseq2.R:", conditionMessage(e), "\n")
  load_ok <<- FALSE
})

tryCatch({
  source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_standardize.R")
  cat("✓ deg_standardize.R\n")
}, error = function(e) {
  cat("✗ deg_standardize.R:", conditionMessage(e), "\n")
  load_ok <<- FALSE
})

tryCatch({
  source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_consensus_analysis.R")
  cat("✓ deg_consensus_analysis.R\n")
}, error = function(e) {
  cat("✗ deg_consensus_analysis.R:", conditionMessage(e), "\n")
  load_ok <<- FALSE
})

tryCatch({
  source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/run_deg_consensus.R")
  cat("✓ run_deg_consensus.R\n")
}, error = function(e) {
  cat("✗ run_deg_consensus.R:", conditionMessage(e), "\n")
  load_ok <<- FALSE
})

if (!load_ok) {
  stop("함수 로드 실패")
}

cat("\n")

# ============================================================================
# Step 2: 데이터 로드
# ============================================================================
cat("Step 2: 데이터 로드\n")
cat("--------------------------------------------------------------------------------\n")

is5s <- tryCatch({
  qs::qread("/data/user3/sobj/IS6_sex_added_251110_ds2500.qs")
}, error = function(e) {
  cat("✗ 데이터 로드 실패:", conditionMessage(e), "\n")
  NULL
})

if (is.null(is5s)) {
  stop("데이터 로드 실패")
}

cat(sprintf("✓ 데이터 로드: %d cells, %d genes\n", ncol(is5s), nrow(is5s)))

# 메타데이터 확인
if ("anno3.scvi" %in% colnames(is5s@meta.data)) {
  cluster_key <- "anno3.scvi"
} else if ("seurat_clusters" %in% colnames(is5s@meta.data)) {
  cluster_key <- "seurat_clusters"
} else {
  stop("cluster_id 컬럼을 찾을 수 없습니다.")
}

sample_key <- "hos_no"
group_key <- "g3"
batch_key <- if ("GEM" %in% colnames(is5s@meta.data)) "GEM" else NULL
contrast_str <- "2 - 1"

cat(sprintf("  cluster_id: %s (%d levels)\n", cluster_key, length(unique(is5s@meta.data[[cluster_key]]))))
cat(sprintf("  sample_id: %s (%d samples)\n", sample_key, length(unique(is5s@meta.data[[sample_key]]))))
cat(sprintf("  group_id: %s (%s)\n", group_key, paste(unique(is5s@meta.data[[group_key]]), collapse=", ")))
cat(sprintf("  batch_id: %s\n", ifelse(is.null(batch_key), "NULL", batch_key)))
cat(sprintf("  contrast: %s\n", contrast_str))
cat("\n")

# ============================================================================
# Step 3: 개별 방법론 테스트 (하나씩)
# ============================================================================
cat("Step 3: 개별 방법론 테스트\n")
cat("--------------------------------------------------------------------------------\n")

test_methods <- c("muscat-edgeR", "limma-voom", "edgeR-LRT")
individual_results <- list()

for (method in test_methods) {
  cat(sprintf("\n[%s] 테스트 중...\n", method))
  
  if (method == "muscat-edgeR") {
    result <- tryCatch({
      runMUSCAT2_v1(
        sobj = is5s,
        cluster_id = cluster_key,
        sample_id = sample_key,
        group_id = group_key,
        batch_id = batch_key,
        contrast = contrast_str,
        method = "edgeR",
        remove_na_groups = TRUE
      )
    }, error = function(e) {
      cat(sprintf("  ✗ 실패: %s\n", conditionMessage(e)))
      traceback()
      NULL
    })
  } else if (method == "limma-voom") {
    result <- tryCatch({
      runLIMMA_voom_v1(
        sobj = is5s,
        cluster_id = cluster_key,
        sample_id = sample_key,
        group_id = group_key,
        batch_id = batch_key,
        contrast = contrast_str,
        remove_na_groups = TRUE
      )
    }, error = function(e) {
      cat(sprintf("  ✗ 실패: %s\n", conditionMessage(e)))
      traceback()
      NULL
    })
  } else if (method == "edgeR-LRT") {
    result <- tryCatch({
      runEDGER_LRT_v1(
        sobj = is5s,
        cluster_id = cluster_key,
        sample_id = sample_key,
        group_id = group_key,
        batch_id = batch_key,
        contrast = contrast_str,
        remove_na_groups = TRUE
      )
    }, error = function(e) {
      cat(sprintf("  ✗ 실패: %s\n", conditionMessage(e)))
      traceback()
      NULL
    })
  }
  
  if (!is.null(result)) {
    cat(sprintf("  ✓ 완료: %d rows\n", nrow(result)))
    cat(sprintf("  컬럼: %s\n", paste(head(colnames(result), 10), collapse=", ")))
    individual_results[[method]] <- result
  }
}

if (length(individual_results) == 0) {
  stop("모든 방법론이 실패했습니다.")
}

cat(sprintf("\n✓ %d 개 방법론 성공\n", length(individual_results)))
cat("\n")

# ============================================================================
# Step 4: 표준화 테스트
# ============================================================================
cat("Step 4: 결과 표준화 테스트\n")
cat("--------------------------------------------------------------------------------\n")

standardized_results <- list()
for (method in names(individual_results)) {
  cat(sprintf("[%s] 표준화 중...\n", method))
  standardized <- tryCatch({
    standardize_deg_results(individual_results[[method]], method)
  }, error = function(e) {
    cat(sprintf("  ✗ 실패: %s\n", conditionMessage(e)))
    if (interactive()) traceback()
    NULL
  })
  
  if (!is.null(standardized)) {
    cat(sprintf("  ✓ 완료: %d rows\n", nrow(standardized)))
    cat(sprintf("  컬럼: %s\n", paste(colnames(standardized), collapse=", ")))
    standardized_results[[method]] <- standardized
  }
}

if (length(standardized_results) == 0) {
  stop("표준화 실패")
}

cat(sprintf("\n✓ %d 개 방법론 표준화 완료\n", length(standardized_results)))
cat("\n")

# ============================================================================
# Step 5: 행렬 구성 테스트
# ============================================================================
cat("Step 5: 행렬 구성 테스트\n")
cat("--------------------------------------------------------------------------------\n")

deg_matrices <- tryCatch({
  build_deg_matrices(standardized_results, genes = NULL, fdr_threshold = 0.05)
}, error = function(e) {
  cat(sprintf("✗ 행렬 구성 실패: %s\n", conditionMessage(e)))
  if (interactive()) traceback()
  NULL
})

if (is.null(deg_matrices)) {
  stop("행렬 구성 실패")
}

cat(sprintf("✓ 행렬 구성 완료\n"))
cat(sprintf("  beta: %d genes × %d methods\n", nrow(deg_matrices$beta), ncol(deg_matrices$beta)))
cat(sprintf("  pvalue: %d genes × %d methods\n", nrow(deg_matrices$pvalue), ncol(deg_matrices$pvalue)))
cat(sprintf("  significance: %d genes × %d methods\n", nrow(deg_matrices$significance), ncol(deg_matrices$significance)))
cat("\n")

# ============================================================================
# Step 6: Consensus 분석 테스트
# ============================================================================
cat("Step 6: Consensus 분석 테스트\n")
cat("--------------------------------------------------------------------------------\n")

# Agreement scores
agreement_scores <- tryCatch({
  compute_agreement_scores(deg_matrices$significance)
}, error = function(e) {
  cat(sprintf("✗ Agreement scores 계산 실패: %s\n", conditionMessage(e)))
  NULL
})

if (!is.null(agreement_scores)) {
  cat(sprintf("✓ Agreement scores: %d genes\n", length(agreement_scores)))
  cat(sprintf("  범위: %.3f ~ %.3f\n", min(agreement_scores, na.rm=TRUE), max(agreement_scores, na.rm=TRUE)))
  cat(sprintf("  평균: %.3f\n", mean(agreement_scores, na.rm=TRUE)))
}

# Consensus scores
consensus_scores <- tryCatch({
  compute_consensus_scores(deg_matrices, agreement_scores)
}, error = function(e) {
  cat(sprintf("✗ Consensus scores 계산 실패: %s\n", conditionMessage(e)))
  NULL
})

if (!is.null(consensus_scores)) {
  cat(sprintf("✓ Consensus scores: %d genes\n", nrow(consensus_scores)))
  cat(sprintf("  상위 5개 유전자:\n"))
  print(head(consensus_scores[order(-consensus_scores$consensus_score), 
                              c("gene", "agreement", "consensus_score", "n_significant")], 5))
}

# Consensus DEG list
consensus_deg_list <- tryCatch({
  generate_consensus_deg_list(
    consensus_scores,
    fdr_threshold = 0.05,
    agreement_threshold = 0.5,
    min_methods = ceiling(length(standardized_results) * 0.5)
  )
}, error = function(e) {
  cat(sprintf("✗ Consensus DEG list 생성 실패: %s\n", conditionMessage(e)))
  NULL
})

if (!is.null(consensus_deg_list)) {
  cat(sprintf("\n✓ Consensus DEG list: %d genes\n", nrow(consensus_deg_list)))
}

cat("\n")

# ============================================================================
# Step 7: 결과 저장
# ============================================================================
cat("Step 7: 결과 저장\n")
cat("--------------------------------------------------------------------------------\n")

final_result <- list(
  individual_results = individual_results,
  standardized_results = standardized_results,
  deg_matrices = deg_matrices,
  agreement_scores = agreement_scores,
  consensus_scores = consensus_scores,
  consensus_deg_list = consensus_deg_list
)

tryCatch({
  qs::qsave(final_result, "/data/user3/sobj/test_deg_consensus_step_by_step.qs")
  cat("✓ 결과 저장: /data/user3/sobj/test_deg_consensus_step_by_step.qs\n")
}, error = function(e) {
  cat(sprintf("✗ 저장 실패: %s\n", conditionMessage(e)))
})

cat("\n")
cat("================================================================================\n")
cat("단계별 테스트 완료\n")
cat("================================================================================\n")


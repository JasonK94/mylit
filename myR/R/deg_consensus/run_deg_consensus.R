# ============================================================================
# Multi-Model DEG Consensus Engine - Main Integration Function
# ============================================================================
# Phase 1: 통합 함수 구현
# ============================================================================
# 
# Note: This function requires:
# - runMUSCAT2_v1 and runNEBULA2_v1 (from test_analysis.R)
# - runLIMMA_voom_v1, runLIMMA_trend_v1 (from deg_methods_limma.R)
# - runEDGER_LRT_v1, runEDGER_QLF_v1 (from deg_methods_edger.R)
# - runDESEQ2_Wald_v1, runDESEQ2_LRT_v1 (from deg_methods_deseq2.R)
# These functions should be loaded before calling run_deg_consensus
# ============================================================================

#' Run Multiple DEG Methods and Collect Results
#'
#' @description
#' Runs multiple differential expression analysis methods on the same dataset
#' and collects results for consensus analysis. This is the main integration
#' function that orchestrates multiple DEG methods.
#'
#' @param sobj Seurat object
#' @param contrast Contrast string (e.g., "2 - 1" for g3==2 vs g3==1)
#' @param methods Character vector of methods to run. Supported methods:
#'   - "muscat-edgeR", "muscat-DESeq2", "muscat-limma-voom", "muscat-limma-trend"
#'   - "nebula"
#'   - "limma-voom", "limma-trend", "limma-wt", "dream"
#'   - "edgeR-LRT", "edgeR-QLF", "edgeR-robust"
#'   - "DESeq2-Wald", "DESeq2-LRT"
#' @param cluster_id Column name for cell type/cluster (default: "anno3.scvi")
#' @param sample_id Column name for sample ID (default: "hos_no")
#' @param group_id Column name for group/condition (default: "g3")
#' @param batch_id Optional column name for batch variable (default: "GEM")
#' @param pb_min_cells Minimum cells per pseudobulk sample (default: 3)
#' @param filter_genes Filtering method for muscat (default: "edgeR")
#' @param keep_clusters Optional vector of cluster IDs to keep
#' @param remove_na_groups Remove cells with NA in group_id (default: TRUE)
#' @param parallel Logical. If TRUE, run methods in parallel (default: FALSE)
#' @param n_cores Number of cores for parallel execution (default: 4)
#' @param verbose Logical. Print progress messages (default: TRUE)
#'
#' @return A list with components:
#'   \describe{
#'     \item{results}{Named list of DEG results from each method}
#'     \item{methods_run}{Character vector of successfully run methods}
#'     \item{methods_failed}{Character vector of failed methods}
#'     \item{errors}{Named list of error messages for failed methods}
#'   }
#'
#' @export
run_deg_consensus <- function(
  sobj,
  contrast = "2 - 1",
  methods = c("muscat-edgeR", "muscat-DESeq2", "muscat-limma-voom", 
              "muscat-limma-trend", "nebula"),
  cluster_id = "anno3.scvi",
  sample_id = "hos_no",
  group_id = "g3",
  batch_id = "GEM",
  pb_min_cells = 3,
  filter_genes = "edgeR",
  keep_clusters = NULL,
  remove_na_groups = TRUE,
  parallel = FALSE,
  n_cores = 4,
  verbose = TRUE
) {
  
  # --- 0. Validation ---
  if (is.null(contrast)) {
    stop("'contrast'를 지정하세요. 예: '2 - 1'")
  }
  
  if (length(methods) == 0) {
    stop("최소 하나의 방법론을 지정해야 합니다.")
  }
  
  # Check required packages
  req_pkgs <- c("Seurat")
  miss_pkgs <- req_pkgs[!vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss_pkgs) > 0) {
    stop("필요 패키지 설치: ", paste(miss_pkgs, collapse = ", "))
  }
  
  # --- 1. Method mapping ---
  # Map method names to their handler functions
  method_handlers <- list(
    # muscat methods
    "muscat-edgeR" = function(...) {
      runMUSCAT2_v1(..., method = "edgeR")
    },
    "muscat-DESeq2" = function(...) {
      runMUSCAT2_v1(..., method = "DESeq2")
    },
    "muscat-limma-voom" = function(...) {
      runMUSCAT2_v1(..., method = "limma-voom")
    },
    "muscat-limma-trend" = function(...) {
      runMUSCAT2_v1(..., method = "limma-trend")
    },
    # nebula
    "nebula" = function(...) {
      runNEBULA2_v1(...)
    },
    # limma methods
    "limma-voom" = function(...) {
      if (!exists("runLIMMA_voom_v1")) {
        stop("runLIMMA_voom_v1 함수가 로드되지 않았습니다. deg_methods_limma.R를 source하세요.")
      }
      runLIMMA_voom_v1(...)
    },
    "limma-trend" = function(...) {
      if (!exists("runLIMMA_trend_v1")) {
        stop("runLIMMA_trend_v1 함수가 로드되지 않았습니다. deg_methods_limma.R를 source하세요.")
      }
      runLIMMA_trend_v1(...)
    },
    "limma-wt" = function(...) {
      if (!exists("runLIMMA_wt_v1")) {
        stop("runLIMMA_wt_v1 함수가 아직 구현되지 않았습니다.")
      }
      runLIMMA_wt_v1(...)
    },
    "dream" = function(...) {
      if (!exists("runDREAM_v1")) {
        stop("runDREAM_v1 함수가 아직 구현되지 않았습니다.")
      }
      runDREAM_v1(...)
    },
    # edgeR methods
    "edgeR-LRT" = function(...) {
      if (!exists("runEDGER_LRT_v1")) {
        stop("runEDGER_LRT_v1 함수가 로드되지 않았습니다. deg_methods_edger.R를 source하세요.")
      }
      runEDGER_LRT_v1(...)
    },
    "edgeR-QLF" = function(...) {
      if (!exists("runEDGER_QLF_v1")) {
        stop("runEDGER_QLF_v1 함수가 로드되지 않았습니다. deg_methods_edger.R를 source하세요.")
      }
      runEDGER_QLF_v1(...)
    },
    "edgeR-robust" = function(...) {
      if (!exists("runEDGER_robust_v1")) {
        stop("runEDGER_robust_v1 함수가 아직 구현되지 않았습니다.")
      }
      runEDGER_robust_v1(...)
    },
    # DESeq2 methods
    "DESeq2-Wald" = function(...) {
      if (!exists("runDESEQ2_Wald_v1")) {
        stop("runDESEQ2_Wald_v1 함수가 로드되지 않았습니다. deg_methods_deseq2.R를 source하세요.")
      }
      runDESEQ2_Wald_v1(...)
    },
    "DESeq2-LRT" = function(...) {
      if (!exists("runDESEQ2_LRT_v1")) {
        stop("runDESEQ2_LRT_v1 함수가 로드되지 않았습니다. deg_methods_deseq2.R를 source하세요.")
      }
      runDESEQ2_LRT_v1(...)
    }
  )
  
  # Check if all requested methods are available
  missing_methods <- setdiff(methods, names(method_handlers))
  if (length(missing_methods) > 0) {
    warning("다음 방법론은 지원되지 않습니다: ", paste(missing_methods, collapse = ", "))
    methods <- setdiff(methods, missing_methods)
  }
  
  if (length(methods) == 0) {
    stop("실행 가능한 방법론이 없습니다.")
  }
  
  # --- 2. Prepare common arguments ---
  # Arguments that are common across methods
  common_args <- list(
    sobj = sobj,
    contrast = contrast,
    cluster_id = cluster_id,
    sample_id = sample_id,
    group_id = group_id,
    batch_id = batch_id,
    pb_min_cells = pb_min_cells,
    filter_genes = filter_genes,
    keep_clusters = keep_clusters,
    remove_na_groups = remove_na_groups
  )
  
  # Method-specific arguments
  muscat_args <- c("pb_min_cells", "filter_genes", "keep_clusters", "cluster_label_map")
  nebula_args <- c("layer", "fixed_effects", "covar_effects", "patient_col", 
                   "offset", "min_count", "remove_na_cells")
  
  # --- 3. Run methods ---
  results <- list()
  methods_run <- character(0)
  methods_failed <- character(0)
  errors <- list()
  
  if (verbose) {
    message(sprintf("=== DEG Consensus Analysis 시작 ==="))
    message(sprintf("총 %d 개의 방법론 실행 예정", length(methods)))
    message(sprintf("Contrast: %s", contrast))
  }
  
  for (method in methods) {
    if (verbose) {
      message(sprintf("\n--- [%s] 실행 중 ---", method))
    }
    
    # Check if handler function exists
    handler <- method_handlers[[method]]
    if (is.null(handler)) {
      warning(sprintf("방법론 '%s'에 대한 핸들러가 없습니다. 건너뜁니다.", method))
      methods_failed <- c(methods_failed, method)
      errors[[method]] <- "Handler function not found"
      next
    }
    
    # Prepare method-specific arguments
    method_args <- common_args
    
    # For nebula, adjust arguments
    if (method == "nebula") {
      method_args$fixed_effects <- c(group_id)
      method_args$patient_col <- sample_id
      method_args$remove_na_cells <- remove_na_groups
      # Remove muscat-specific args that nebula doesn't use
      method_args$pb_min_cells <- NULL
      method_args$filter_genes <- NULL
      method_args$keep_clusters <- NULL
      method_args$cluster_id <- NULL
      method_args$batch_id <- NULL
      method_args$contrast <- NULL  # nebula doesn't use contrast string
    }
    
    # For muscat methods, ensure cluster_id is provided and fix filter_genes
    if (startsWith(method, "muscat-")) {
      if (is.null(method_args$cluster_id)) {
        stop(sprintf("muscat 방법론은 'cluster_id'가 필요합니다."))
      }
      # filter_genes should be one of "both", "genes", "samples", "none"
      if (!is.null(method_args$filter_genes) && method_args$filter_genes == "edgeR") {
        method_args$filter_genes <- "both"  # edgeR은 "both"로 변경
      }
    }
    
    # For non-muscat methods, remove filter_genes if present
    if (!startsWith(method, "muscat-") && !startsWith(method, "nebula")) {
      method_args$filter_genes <- NULL
    }
    
    # Try to run the method
    tryCatch({
      result <- do.call(handler, method_args)
      
      # Store result with method name
      results[[method]] <- result
      methods_run <- c(methods_run, method)
      
      if (verbose) {
        message(sprintf("✓ [%s] 완료", method))
      }
    }, error = function(e) {
      error_msg <- conditionMessage(e)
      methods_failed <- c(methods_failed, method)
      errors[[method]] <- error_msg
      
      if (verbose) {
        message(sprintf("✗ [%s] 실패: %s", method, error_msg))
      }
    })
  }
  
  # --- 4. Summary ---
  if (verbose) {
    message(sprintf("\n=== 실행 완료 ==="))
    message(sprintf("성공: %d 개", length(methods_run)))
    message(sprintf("실패: %d 개", length(methods_failed)))
    if (length(methods_failed) > 0) {
      message(sprintf("실패한 방법론: %s", paste(methods_failed, collapse = ", ")))
    }
  }
  
  # --- 5. Return results ---
  return(list(
    results = results,
    methods_run = methods_run,
    methods_failed = methods_failed,
    errors = errors,
    contrast = contrast,
    metadata = list(
      cluster_id = cluster_id,
      sample_id = sample_id,
      group_id = group_id,
      batch_id = batch_id
    )
  ))
}


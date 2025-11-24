# ==============================================================================
# TML Core Implementation
# ==============================================================================

TML7 <- function(
  l1_signatures,
  holdout_data,
  target_var,
  l2_methods = c(
    "glm", "ranger", "xgbTree",
    "glmnet", "svmRadial", "mlp",
    "mlpKerasDropout", "nnet", "earth"
  ),
  cv_folds = 5,
  cv_method = "cv",
  repeats = 1,
  metric = c("AUC", "ROC", "Accuracy", "Kappa"),
  fgs_seed = 42,
  layer = "data",
  allow_parallel = FALSE,
  parallel_workers = NULL,
  cv_group_var = "emrid"
) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  # ---- NEW: 표준화 유틸 ----
  as_signature <- function(sig) {
    # NULL 또는 빈 signature 처리
    if (is.null(sig) || (is.list(sig) && length(sig) == 0)) {
      stop("Signature is NULL or empty.")
    }

    # 문자형 유전자 벡터
    if (is.character(sig)) {
      if (length(sig) == 0) stop("Character signature is empty.")
      return(structure(rep(1, length(sig)), names = sig))
    }

    # 이름달린 numeric (가중치)
    if (is.numeric(sig) && !is.null(names(sig))) {
      if (length(sig) == 0) stop("Numeric signature is empty.")
      return(sig)
    }

    # list(up=..., down=...)
    if (is.list(sig) && all(c("up", "down") %in% names(sig))) {
      up <- sig$up
      down <- sig$down
      if (!(is.character(up) || is.null(up)) || !(is.character(down) || is.null(down))) {
        stop("Signature up/down must be character vectors or NULL.")
      }
      up <- up %||% character()
      down <- down %||% character()
      if (length(up) == 0 && length(down) == 0) {
        stop("Signature up and down are both empty.")
      }
      nm <- c(up, down)
      wt <- c(rep(1, length(up)), rep(-1, length(down)))
      return(structure(wt, names = nm))
    }

    # list(genes=..., weights=...)
    if (is.list(sig) && all(c("genes", "weights") %in% names(sig))) {
      g <- sig$genes
      w <- sig$weights
      if (!is.character(g)) {
        stop(sprintf("Signature genes must be character, got: %s", class(g)[1]))
      }
      if (!is.numeric(w)) {
        stop(sprintf("Signature weights must be numeric, got: %s", class(w)[1]))
      }
      if (length(g) == 0 || length(w) == 0) {
        stop("Signature genes or weights is empty.")
      }
      if (length(g) != length(w)) {
        stop(sprintf(
          "Signature genes/weights length mismatch: genes=%d, weights=%d",
          length(g), length(w)
        ))
      }
      names(w) <- g
      return(w)
    }

    # data.frame/tibble: gene/feature + weight/w/score
    if (is.data.frame(sig)) {
      if (nrow(sig) == 0) stop("Data.frame signature is empty.")
      cn <- tolower(colnames(sig))
      gene_col <- which(cn %in% c("gene", "genes", "feature", "features", "symbol", "id"))[1]
      weight_col <- which(cn %in% c("weight", "weights", "w", "score", "scores", "coef", "coefs"))[1]
      if (!is.na(gene_col) && !is.na(weight_col)) {
        g <- as.character(sig[[gene_col]])
        w <- as.numeric(sig[[weight_col]])
        ok <- !is.na(g) & !is.na(w)
        w <- w[ok]
        g <- g[ok]
        if (length(g) == 0) stop("No valid gene-weight pairs in data.frame signature.")
        names(w) <- g
        return(w)
      }
      stop("Data.frame signature missing required columns (gene/genes/feature and weight/weights/w/score).")
    }

    stop(sprintf(
      "Unsupported signature format. Type: %s, Class: %s",
      typeof(sig), paste(class(sig), collapse = ", ")
    ))
  }

  .score_signature <- function(expr_data, signature, normalize = TRUE) {
    # signature를 표준화해 항상 "이름달린 가중치 numeric"으로 맞춘다
    w <- as_signature(signature)

    genes <- base::intersect(rownames(expr_data), names(w))
    if (length(genes) == 0) stop("Signature has no overlap with expression data.")

    ww <- w[genes]
    # (표준) 가중평균
    s <- as.numeric(Matrix::t(expr_data[genes, , drop = FALSE]) %*% ww)
    den <- sum(abs(ww))
    if (!is.finite(den) || den == 0) den <- length(ww)
    s <- s / den

    if (normalize) {
      # zero-variance 체크: scale()이 NA를 반환하는 것을 방지
      s_sd <- stats::sd(s, na.rm = TRUE)
      if (!is.finite(s_sd) || s_sd == 0) {
        # 모든 값이 동일하면 normalize를 건너뛰고 0으로 설정
        s <- rep(0, length(s))
      } else {
        s_scaled <- as.numeric(scale(s))
        # scale()이 NA를 반환할 수 있으므로 체크
        if (any(!is.finite(s_scaled))) {
          warning("Some signature scores became NA/Inf after scaling. Using original scores.")
          # s는 이미 normalize 전 상태이므로 그대로 사용
        } else {
          s <- s_scaled
        }
      }
    }
    s
  }

  .is_binary <- function(y) length(levels(y)) == 2

  metric_choices <- c("AUC", "ROC", "Accuracy", "Kappa")
  if (!is.null(cv_group_var)) {
    cv_group_var <- as.character(cv_group_var)[1]
    if (is.na(cv_group_var)) cv_group_var <- NULL
  }

  .metric_map <- function(user_metric, is_binary) {
    if (is_binary) {
      # caret twoClassSummary → columns: ROC, Sens, Spec
      if (user_metric %in% c("AUC", "ROC")) {
        return(list(train_metric = "ROC", summary = "twoClassSummary"))
      }
      if (user_metric %in% c("Accuracy", "Kappa")) {
        return(list(train_metric = user_metric, summary = "twoClassSummary"))
      }
    } else {
      # 멀티클래스 → defaultSummary (Accuracy, Kappa)
      if (user_metric == "AUC") {
        message("WARNING: 'AUC/ROC' not defined for multi-class defaultSummary. Falling back to 'Accuracy'.")
        return(list(train_metric = "Accuracy", summary = "defaultSummary"))
      }
      return(list(train_metric = user_metric, summary = "defaultSummary"))
    }
    list(
      train_metric = if (is_binary) "ROC" else "Accuracy",
      summary = if (is_binary) "twoClassSummary" else "defaultSummary"
    )
  }

  .package_ok <- function(pkg) {
    suppressWarnings(suppressMessages(requireNamespace(pkg, quietly = TRUE)))
  }

  supported_l2_methods <- c(
    "glm", "ranger", "xgbTree",
    "glmnet", "svmRadial", "mlp",
    "mlpKerasDropout", "nnet", "earth"
  )
  method_dependencies <- c(
    ranger = "ranger",
    xgbTree = "xgboost",
    glmnet = "glmnet",
    svmRadial = "kernlab",
    mlp = "RSNNS",
    mlpKerasDropout = "keras",
    nnet = "nnet",
    earth = "earth"
  )
  default_l2_methods <- c("glm", "ranger", "xgbTree")
  l2_methods <- l2_methods %||% default_l2_methods
  l2_methods <- base::unique(l2_methods)

  unsupported <- base::setdiff(l2_methods, supported_l2_methods)
  if (length(unsupported) > 0) {
    warning(sprintf(
      "Unsupported L2 method(s) removed: %s. Supported methods: %s",
      paste(unsupported, collapse = ", "),
      paste(supported_l2_methods, collapse = ", ")
    ))
    l2_methods <- base::intersect(l2_methods, supported_l2_methods)
  }
  if (length(l2_methods) == 0) {
    stop("No supported L2 methods specified.")
  }

  if (!.package_ok("caret")) stop("caret package required.")
  if (!.package_ok("Matrix")) stop("Matrix package required.")

  register_sequential <- function() {
    if (.package_ok("foreach")) {
      try(foreach::registerDoSEQ(), silent = TRUE)
    }
  }

  # CRITICAL: Set global xgboost config to single-threaded BEFORE any operations
  # This must be done early to prevent xgboost from spawning threads
  if (requireNamespace("xgboost", quietly = TRUE)) {
    tryCatch(
      {
        if (exists("xgb.set.config", envir = asNamespace("xgboost"))) {
          xgboost::xgb.set.config(nthread = 1, verbosity = 0)
        }
        options(xgboost.nthread = 1)
      },
      error = function(e) {
        # If xgb.set.config fails, at least set the option
        options(xgboost.nthread = 1)
      }
    )
    # Also set environment variable that xgboost checks at C level
    Sys.setenv(XGBOOST_NTHREAD = "1")
  }

  # Strict CPU core limiting to prevent cascade parallelization
  # Especially important when child processes load start.R which may spawn more workers
  max_cores_limit <- 16L
  if (requireNamespace("parallel", quietly = TRUE)) {
    detected_cores <- tryCatch(parallel::detectCores(logical = FALSE), error = function(e) NA_integer_)
    if (!is.na(detected_cores)) {
      max_cores_limit <- min(max_cores_limit, as.integer(detected_cores))
    }
  }

  # Set environment variables BEFORE any parallel processing
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1",
    # Prevent future workers from loading start.R
    KDW_START_LOAD_ALL_PACKAGES = "FALSE",
    KDW_START_AUTOLOAD_QS = "FALSE"
  )
  options(mc.cores = 1L)

  allow_parallel <- isTRUE(allow_parallel)
  worker_count <- parallel_workers
  if (is.null(worker_count)) {
    fallback <- getOption("mylit.meta_learner.workers", 4L)
    cores <- tryCatch(parallel::detectCores(logical = FALSE), error = function(e) NA_integer_)
    cores <- cores %||% fallback
    # Limit to max_cores_limit to prevent excessive CPU usage
    worker_count <- max(1L, min(max_cores_limit, as.integer(cores)))
  } else {
    # Also cap user-specified workers
    worker_count <- min(max_cores_limit, as.integer(worker_count))
  }

  # Disable parallel processing by default to prevent cascade effects
  # Even if allow_parallel=TRUE, we'll use sequential execution to avoid
  # child processes spawning more workers via start.R
  allow_parallel <- FALSE
  register_sequential()

  # Note: We disable parallel processing to prevent cascade parallelization
  # when child processes load start.R. If parallel processing is needed,
  # it should be handled at a higher level with proper resource management.

  set.seed(fgs_seed)
  RNGkind(sample.kind = "Rejection")

  # --- 데이터 추출 ---
  if (inherits(holdout_data, "Seurat")) {
    if (!.package_ok("Seurat")) stop("Seurat package required for Seurat input.")
    expr_mat <- Seurat::GetAssayData(holdout_data, layer = layer)
    meta_data <- holdout_data@meta.data
    if (!target_var %in% colnames(meta_data)) {
      stop(sprintf("Target column '%s' not found in Seurat meta.data.", target_var))
    }
    l2_target <- meta_data[[target_var]]
  } else {
    if (is.data.frame(holdout_data)) holdout_data <- as.matrix(holdout_data)
    if (!is.matrix(holdout_data)) stop("holdout_data must be Seurat or a matrix/data.frame [features x cells].")
    expr_mat <- holdout_data
    l2_target <- target_var
  }
  if (is.null(rownames(expr_mat))) stop("Expression matrix must have rownames (feature IDs).")
  if (is.null(colnames(expr_mat))) colnames(expr_mat) <- paste0("cell_", seq_len(ncol(expr_mat)))

  # --- L2 특성 구성 (시그니처 점수 계산) ---
  if (is.null(names(l1_signatures))) names(l1_signatures) <- paste0("L1_", seq_along(l1_signatures))

  # Filter out invalid signatures (failed FGS results, empty signatures, etc.)
  original_sig_count <- length(l1_signatures)
  valid_signatures <- list()
  skipped_signatures <- character(0)

  for (sig_name in names(l1_signatures)) {
    sig <- l1_signatures[[sig_name]]

    # Check if it's a failed FGS result
    if (is.list(sig) && "error" %in% names(sig)) {
      skipped_signatures <- c(skipped_signatures, sig_name)
      next
    }

    # Check if it's an empty result
    if (is.list(sig) && "genes" %in% names(sig) && "weights" %in% names(sig)) {
      if (length(sig$genes) == 0 || length(sig$weights) == 0) {
        skipped_signatures <- c(skipped_signatures, sig_name)
        next
      }
    }

    # Try to validate signature format (will fail if unsupported)
    sig_valid <- tryCatch(
      {
        test_w <- as_signature(sig)
        if (is.null(test_w) || length(test_w) == 0) {
          FALSE
        } else {
          TRUE
        }
      },
      error = function(e) {
        FALSE
      }
    )

    if (sig_valid) {
      valid_signatures[[sig_name]] <- sig
    } else {
      skipped_signatures <- c(skipped_signatures, sig_name)
    }
  }

  if (length(valid_signatures) == 0) {
    stop(sprintf(
      "No valid signatures found after filtering. Input: %d, Skipped: %d (%s)",
      original_sig_count,
      length(skipped_signatures),
      paste(skipped_signatures, collapse = ", ")
    ))
  }

  if (length(skipped_signatures) > 0) {
    message(sprintf(
      "Filtered out %d invalid signature(s): %s",
      length(skipped_signatures),
      paste(skipped_signatures, collapse = ", ")
    ))
    message(sprintf(
      "Proceeding with %d valid signature(s): %s",
      length(valid_signatures),
      paste(names(valid_signatures), collapse = ", ")
    ))
  }

  l1_signatures <- valid_signatures

  l2_features_list <- mapply(function(sig, sig_name) {
    tryCatch(
      {
        .score_signature(expr_mat, sig, normalize = TRUE)
      },
      error = function(e) {
        warning(sprintf(
          "Failed to score signature '%s': %s",
          sig_name, e$message
        ))
        return(rep(NA_real_, ncol(expr_mat)))
      }
    )
  }, l1_signatures, names(l1_signatures), SIMPLIFY = FALSE)
  l2_train_df <- as.data.frame(do.call(cbind, l2_features_list))
  colnames(l2_train_df) <- make.names(names(l1_signatures))

  # 디버깅: 각 signature의 분산 체크
  sig_variances <- apply(l2_train_df, 2, function(x) {
    x_clean <- x[is.finite(x)]
    if (length(x_clean) == 0) {
      return(NA_real_)
    }
    stats::var(x_clean, na.rm = TRUE)
  })
  zero_var_sigs <- names(sig_variances)[is.na(sig_variances) | sig_variances == 0]
  if (length(zero_var_sigs) > 0 && length(zero_var_sigs) < ncol(l2_train_df)) {
    message(sprintf(
      "Note: %d signature(s) have zero variance before nearZeroVar: %s",
      length(zero_var_sigs), paste(zero_var_sigs, collapse = ", ")
    ))
  }

  # --- 타깃 준비/정리 ---
  if (!is.factor(l2_target)) l2_target <- factor(l2_target)
  orig_lvls <- levels(l2_target)
  safe_lvls <- make.names(orig_lvls)
  if (!all(orig_lvls == safe_lvls)) {
    message(sprintf(
      "Sanitizing target levels: '%s' -> '%s'",
      paste(orig_lvls, collapse = "'/'"),
      paste(safe_lvls, collapse = "'/'")
    ))
    levels(l2_target) <- safe_lvls
  }

  # --- 그룹 컬럼 준비 (Seurat + meta.data 컬럼 있을 때만) ---
  cv_group <- NULL
  if (inherits(holdout_data, "Seurat") && !is.null(cv_group_var) && nzchar(cv_group_var)) {
    meta_data <- holdout_data@meta.data
    if (cv_group_var %in% colnames(meta_data)) {
      cv_group <- as.vector(meta_data[[cv_group_var]])
    } else {
      message(sprintf(
        "cv_group_var '%s' not found in Seurat meta.data; using standard cell-wise CV.",
        cv_group_var
      ))
    }
  } else if (!inherits(holdout_data, "Seurat") && !is.null(cv_group_var) && nzchar(cv_group_var)) {
    message("cv_group_var requires Seurat metadata; using standard cell-wise CV.")
  }

  if (!is.null(cv_group)) {
    if (all(is.na(cv_group))) {
      message(sprintf("cv_group_var '%s' is entirely NA; using standard cell-wise CV.", cv_group_var))
      cv_group <- NULL
    } else if (anyNA(cv_group)) {
      message(sprintf(
        "cv_group_var '%s' contains %d NA value(s); using standard cell-wise CV.",
        cv_group_var, sum(is.na(cv_group))
      ))
      cv_group <- NULL
    }
  }

  keep <- !is.na(l2_target)
  if (!all(keep)) {
    message(sprintf("Removing %d rows with NA target.", sum(!keep)))
    l2_target <- l2_target[keep]
    l2_train_df <- l2_train_df[keep, , drop = FALSE]
    if (!is.null(cv_group)) cv_group <- cv_group[keep]
  }

  row_ok <- stats::complete.cases(l2_train_df) &
    apply(l2_train_df, 1, function(r) all(is.finite(r)))
  if (!all(row_ok)) {
    message(sprintf("Removing %d rows with NA/NaN/Inf features.", sum(!row_ok)))
    l2_target <- l2_target[row_ok]
    l2_train_df <- l2_train_df[row_ok, , drop = FALSE]
    if (!is.null(cv_group)) cv_group <- cv_group[row_ok]
  }
  l2_target <- base::droplevels(l2_target)

  nzv <- caret::nearZeroVar(l2_train_df, saveMetrics = FALSE)
  if (length(nzv) > 0) {
    message(sprintf(
      "Removing %d zero-variance features: %s",
      length(nzv), paste(colnames(l2_train_df)[nzv], collapse = ", ")
    ))
    l2_train_df <- l2_train_df[, -nzv, drop = FALSE]
  }

  if (nrow(l2_train_df) == 0 || ncol(l2_train_df) == 0) {
    # 더 자세한 에러 메시지 제공
    n_sigs_input <- length(l1_signatures)
    n_sigs_after_scoring <- length(l2_features_list)
    n_rows_after_na_removal <- sum(row_ok)
    n_cols_after_nzv <- ncol(l2_train_df)

    msg <- sprintf(
      "No usable data remains after cleaning.\n  Input signatures: %d\n  Signatures after scoring: %d\n  Rows after NA/Inf removal: %d\n  Columns after zero-variance removal: %d\n  Possible causes: (1) All signature scores are identical (zero variance), (2) All rows contain NA/Inf, (3) Signatures have no overlap with expression data.",
      n_sigs_input, n_sigs_after_scoring, n_rows_after_na_removal, n_cols_after_nzv
    )
    stop(msg)
  }

  # --- metric 매핑 & group-wise CV index 구성 ---
  is_bin <- .is_binary(l2_target)
  metric <- match.arg(metric, metric_choices)
  map <- .metric_map(metric, is_bin)
  caret_metric <- map$train_metric

  summary_fun <- if (map$summary == "twoClassSummary") caret::twoClassSummary else caret::defaultSummary
  positive_class <- if (length(levels(l2_target)) > 0) tail(levels(l2_target), 1) else NA_character_

  index <- indexOut <- NULL
  log_group_fold_details <- function(index_list, index_out_list, group_vector) {
    if (is.null(index_list) || is.null(index_out_list)) {
      return()
    }
    for (fold in seq_along(index_list)) {
      in_idx <- index_list[[fold]]
      out_idx <- index_out_list[[fold]]
      in_groups <- sort(base::unique(as.character(group_vector[in_idx])))
      out_groups <- sort(base::unique(as.character(group_vector[out_idx])))
      holdout_preview <- if (length(out_groups) > 5) {
        paste0(paste(head(out_groups, 5), collapse = ", "), ", ...")
      } else if (length(out_groups) == 0) {
        "<none>"
      } else {
        paste(out_groups, collapse = ", ")
      }
      message(sprintf(
        "Group CV fold %d: index=%d cells (%d groups), indexOut=%d cells (%d groups), hold-out groups: %s",
        fold,
        length(in_idx),
        length(in_groups),
        length(out_idx),
        length(out_groups),
        holdout_preview
      ))
    }
  }
  # --- CV Folds Generation ---
  set.seed(fgs_seed)

  # Sample size check for LOGO recommendation
  if (!is.null(cv_group)) {
    n_groups <- length(unique(cv_group))
    if (n_groups < 20 && cv_method != "LOGO") {
      message(sprintf("Suggestion: Sample size (groups) is small (n=%d). Consider using cv_method='LOGO' (Leave-One-Group-Out) for more robust validation.", n_groups))
    }
  }

  if (cv_method == "LOGO") {
    if (is.null(cv_group)) stop("cv_method='LOGO' requires valid cv_group_var.")
    unique_groups <- unique(cv_group)
    message(sprintf("Using Leave-One-Group-Out CV (%d groups).", length(unique_groups)))

    # Manual LOGO implementation to ensure correctness
    index <- list()
    for (i in seq_along(unique_groups)) {
      g <- unique_groups[i]
      # Train indices: rows where cv_group is NOT g
      train_idx <- which(cv_group != g)
      index[[paste0("Fold", i)]] <- train_idx
    }
    indexOut <- lapply(index, function(train_idx) base::setdiff(seq_along(l2_target), train_idx))

    log_group_fold_details(index, indexOut, cv_group)
  } else if (cv_method == "repeatedcv") {
    message(sprintf("Using Repeated CV (%d folds, %d repeats).", cv_folds, repeats))
    if (!is.null(cv_group)) {
      warning("Repeated CV with groups is not fully supported. Using standard Repeated CV (ignoring groups).")
    }
    index <- caret::createMultiFolds(l2_target, k = cv_folds, times = repeats)
    indexOut <- lapply(index, function(train_idx) base::setdiff(seq_along(l2_target), train_idx))
  } else {
    # Standard CV or Group CV
    use_group_cv <- !is.null(cv_group) && !all(is.na(cv_group))

    if (use_group_cv) {
      group_factor <- factor(cv_group)
      unique_groups <- levels(group_factor)

      if (length(unique_groups) < cv_folds) {
        message(sprintf(
          "cv_group_var '%s' has only %d unique values (< cv_folds=%d); using standard cell-wise CV.",
          cv_group_var, length(unique_groups), cv_folds
        ))
        use_group_cv <- FALSE
      } else {
        # Manual Group K-Fold to ensure balance if needed, or use caret::groupKFold
        # Existing logic used manual sampling. Let's stick to it or use caret.
        # Caret's groupKFold is simple.
        index <- caret::groupKFold(cv_group, k = cv_folds)
        indexOut <- lapply(index, function(train_idx) base::setdiff(seq_along(l2_target), train_idx))

        message(sprintf(
          "Using group-wise CV on '%s' (%d groups) for L2 (cv_folds=%d).",
          cv_group_var, length(unique_groups), cv_folds
        ))
        log_group_fold_details(index, indexOut, group_factor)
      }
    }

    if (!use_group_cv) {
      message("Using standard cell-wise CV folds (explicitly generated).")
      index <- caret::createFolds(l2_target, k = cv_folds, returnTrain = TRUE)
      indexOut <- lapply(index, function(train_idx) base::setdiff(seq_along(l2_target), train_idx))
    }
  }

  # Always disable parallel processing in caret to prevent cascade parallelization
  # Child processes from future/doParallel may load start.R and spawn more workers
  ctrl <- caret::trainControl(
    method = "cv", number = cv_folds,
    classProbs = is_bin,
    summaryFunction = summary_fun,
    savePredictions = "final",
    allowParallel = FALSE, # Always FALSE to prevent cascade parallelization
    index = index,
    indexOut = indexOut
  )

  drop_if_missing <- function(method_name, pkg) {
    if (method_name %in% l2_methods && !.package_ok(pkg)) {
      warning(sprintf("%s package required for '%s'; dropping.", pkg, method_name))
      l2_methods <<- base::setdiff(l2_methods, method_name)
    }
  }

  deps_to_check <- base::intersect(names(method_dependencies), l2_methods)
  if (length(deps_to_check) > 0) {
    for (m in deps_to_check) {
      drop_if_missing(m, method_dependencies[[m]])
    }
  }

  if ("xgbTree" %in% l2_methods) {
    if (!.package_ok("xgboost")) {
      warning("xgboost not available; dropping 'xgbTree'.")
      l2_methods <- base::setdiff(l2_methods, "xgbTree")
    } else {
      ok <- TRUE
      tryCatch(utils::packageVersion("xgboost"), error = function(e) ok <<- FALSE)
      if (!ok) {
        warning("xgboost seems broken; dropping 'xgbTree'.")
        l2_methods <- base::setdiff(l2_methods, "xgbTree")
      } else {
        # 가능하면 전역 verbosity를 0으로 낮춰 C-level 경고(특히 ntree_limit)를 숨긴다
        if (isTRUE(.package_ok("xgboost"))) {
          xgb_ns <- asNamespace("xgboost")
          if (exists("xgb.set.config", envir = xgb_ns, mode = "function")) {
            try(xgboost::xgb.set.config(verbosity = 0), silent = TRUE)
          }
        }
      }
    }
  }
  if (length(l2_methods) == 0) stop("No usable models remain.")

  message(sprintf(
    "Validated L2 method set (%d): %s",
    length(l2_methods),
    paste(l2_methods, collapse = ", ")
  ))

  # Ensure no parallel backends are active before training
  if (requireNamespace("foreach", quietly = TRUE)) {
    try(foreach::registerDoSEQ(), silent = TRUE)
  }
  if (requireNamespace("doParallel", quietly = TRUE)) {
    try(doParallel::stopImplicitCluster(), silent = TRUE)
  }
  if (requireNamespace("doMC", quietly = TRUE)) {
    try(doMC::registerDoMC(cores = 1), silent = TRUE)
  }

  # Progress tracking for TML7
  tml_timing_file <- file.path(getwd(), ".tml7_timing_cache.rds")
  tml_timing_cache <- if (file.exists(tml_timing_file)) {
    tryCatch(readRDS(tml_timing_file), error = function(e) list())
  } else {
    list()
  }

  total_l2_methods <- length(l2_methods)
  completed_l2_methods <- 0
  tml_start_time <- Sys.time()

  model_list <- list()
  for (m_idx in seq_along(l2_methods)) {
    m <- l2_methods[m_idx]

    # Estimate time based on previous runs
    method_key <- paste0("tml7_l2_", m)
    estimated_sec <- if (!is.null(tml_timing_cache[[method_key]])) {
      mean(tml_timing_cache[[method_key]], na.rm = TRUE)
    } else {
      NA_real_
    }

    if (!is.na(estimated_sec) && estimated_sec > 0) {
      estimated_str <- sprintf("%.1f분", estimated_sec / 60)
      message(sprintf(
        "Training L2 candidate: %s [%d/%d] (예상: %s, metric=%s)",
        m, m_idx, total_l2_methods, estimated_str, caret_metric
      ))
    } else {
      message(sprintf(
        "Training L2 candidate: %s [%d/%d] (metric=%s)",
        m, m_idx, total_l2_methods, caret_metric
      ))
    }

    l2_method_start_time <- Sys.time()

    # === STRATEGY 1: Force cleanup of all parallel backends before caret::train ===
    # This prevents caret from using any registered parallel backends
    if (requireNamespace("doParallel", quietly = TRUE)) {
      tryCatch(
        {
          # Stop any implicit cluster
          doParallel::stopImplicitCluster()
        },
        error = function(e) {
          # Ignore if no cluster exists
        }
      )
      tryCatch(
        {
          # Force sequential backend
          doParallel::registerDoSEQ()
        },
        error = function(e) {
          # Ignore errors
        }
      )
    }

    if (requireNamespace("foreach", quietly = TRUE)) {
      tryCatch(
        {
          foreach::registerDoSEQ()
        },
        error = function(e) {}
      )
    }

    if (requireNamespace("doMC", quietly = TRUE)) {
      tryCatch(
        {
          doMC::registerDoMC(cores = 1)
        },
        error = function(e) {}
      )
    }

    # === STRATEGY 2: Re-set environment variables before each model ===
    # Ensure BLAS/LAPACK threads are still limited
    Sys.setenv(
      OMP_NUM_THREADS = "1",
      OPENBLAS_NUM_THREADS = "1",
      MKL_NUM_THREADS = "1",
      VECLIB_MAXIMUM_THREADS = "1",
      NUMEXPR_NUM_THREADS = "1"
    )

    # === STRATEGY 3: Model-specific thread limits (ENHANCED) ===
    # CRITICAL: Model packages (xgboost, ranger) use C/C++ level parallel processing
    # Even with nthread=1, they may spawn threads. We need to set global configs.
    if (m == "xgbTree" && requireNamespace("xgboost", quietly = TRUE)) {
      # xgboost: Set global config BEFORE any xgboost operations
      # This must be done at the C level, not just R level
      tryCatch(
        {
          # Set global xgboost config to single-threaded
          if (exists("xgb.set.config", envir = asNamespace("xgboost"))) {
            xgboost::xgb.set.config(nthread = 1, verbosity = 0)
          }
          # Also set R-level option as backup
          options(xgboost.nthread = 1)
        },
        error = function(e) {
          # If xgb.set.config fails, at least set the option
          options(xgboost.nthread = 1)
        }
      )
      # CRITICAL: Also set environment variable that xgboost checks at C level
      Sys.setenv(OMP_NUM_THREADS = "1")
      Sys.setenv(XGBOOST_NTHREAD = "1")
    } else if (m == "ranger" && requireNamespace("ranger", quietly = TRUE)) {
      # ranger: Set global option (ranger checks this at C level)
      options(ranger.num.threads = 1)
      # CRITICAL: ranger also respects OMP_NUM_THREADS
      Sys.setenv(OMP_NUM_THREADS = "1")
    } else if (m %in% c("glmnet", "nnet", "earth", "mlp", "mlpKerasDropout")) {
      # These use OpenMP, ensure it's limited
      Sys.setenv(OMP_NUM_THREADS = "1")
    }

    # === STRATEGY 4: Garbage collection before training ===
    # Clean up memory to reduce overhead
    gc(verbose = FALSE)

    # === Train model with all safeguards ===
    # CRITICAL: Wrap caret::train in a function that enforces single-threading
    # This ensures that even if caret spawns child processes, they inherit our settings
    train_with_strict_limits <- function(...) {
      # Re-assert environment variables inside the function
      # This ensures child processes (if any) inherit these settings
      old_env <- Sys.getenv(c(
        "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"
      ))
      on.exit(
        {
          # Restore old environment (though we want to keep limits)
          do.call(Sys.setenv, as.list(old_env))
        },
        add = TRUE
      )

      # Force single-threading environment
      Sys.setenv(
        OMP_NUM_THREADS = "1",
        OPENBLAS_NUM_THREADS = "1",
        MKL_NUM_THREADS = "1",
        VECLIB_MAXIMUM_THREADS = "1",
        NUMEXPR_NUM_THREADS = "1"
      )

      # Call caret::train
      caret::train(...)
    }

    if (m == "xgbTree") {
      # R-level warning은 suppressWarnings로, C-level 로그는 위에서 xgb.set.config(verbosity=0)로 최대한 억제
      fit <- try(
        suppressWarnings(
          train_with_strict_limits(
            x = l2_train_df,
            y = l2_target,
            method = m,
            trControl = ctrl,
            metric = caret_metric,
            tuneLength = 5,
            nthread = 1 # Force single-threaded xgboost (also set globally above)
          )
        ),
        silent = TRUE
      )
    } else if (m == "ranger") {
      fit <- try(
        train_with_strict_limits(
          x = l2_train_df,
          y = l2_target,
          method = m,
          trControl = ctrl,
          metric = caret_metric,
          tuneLength = 5,
          num.threads = 1, # Force single-threaded ranger (also set globally above)
          importance = "permutation" # Enable importance calculation for compute_meta_gene_importance
        ),
        silent = TRUE
      )
    } else {
      fit <- try(
        train_with_strict_limits(
          x = l2_train_df,
          y = l2_target,
          method = m,
          trControl = ctrl,
          metric = caret_metric,
          tuneLength = 5
        ),
        silent = TRUE
      )
    }

    # === STRATEGY 5: Cleanup after training ===
    # Force cleanup again after training to prevent accumulation
    if (requireNamespace("doParallel", quietly = TRUE)) {
      tryCatch(
        {
          doParallel::stopImplicitCluster()
          doParallel::registerDoSEQ()
        },
        error = function(e) {}
      )
    }
    if (requireNamespace("foreach", quietly = TRUE)) {
      tryCatch(
        {
          foreach::registerDoSEQ()
        },
        error = function(e) {}
      )
    }
    gc(verbose = FALSE)
    l2_method_end_time <- Sys.time()
    elapsed_sec <- as.numeric(difftime(l2_method_end_time, l2_method_start_time, units = "secs"))

    if (inherits(fit, "try-error")) {
      warning(sprintf("Failed to train '%s': %s", m, as.character(fit)))
      completed_l2_methods <- completed_l2_methods + 1
      message(sprintf(
        "✗ %s 실패: %.1f초 (%.1f분) | 진행: %d/%d",
        m, elapsed_sec, elapsed_sec / 60, completed_l2_methods, total_l2_methods
      ))
    } else {
      model_list[[m]] <- fit
      completed_l2_methods <- completed_l2_methods + 1

      # Update timing cache
      method_key <- paste0("tml7_l2_", m)
      if (is.null(tml_timing_cache[[method_key]])) {
        tml_timing_cache[[method_key]] <- numeric()
      }
      tml_timing_cache[[method_key]] <- c(tml_timing_cache[[method_key]], elapsed_sec)
      if (length(tml_timing_cache[[method_key]]) > 5) {
        tml_timing_cache[[method_key]] <- tail(tml_timing_cache[[method_key]], 5)
      }

      # Save timing cache
      tryCatch(
        {
          saveRDS(tml_timing_cache, tml_timing_file)
        },
        error = function(e) {
          # Ignore
        }
      )

      # Improved progress estimation: use actual elapsed time from first completed method
      # After first method completes, we can estimate remaining time more accurately
      elapsed_total <- as.numeric(difftime(l2_method_end_time, tml_start_time, units = "secs"))
      remaining_methods <- total_l2_methods - completed_l2_methods

      if (completed_l2_methods == 1 && remaining_methods > 0) {
        # First method completed: estimate based on this method's time
        # For CV with cv_folds, we can estimate: this method took elapsed_sec for all folds
        # So remaining methods will take approximately: elapsed_sec * remaining_methods
        estimated_remaining <- elapsed_sec * remaining_methods
        message(sprintf(
          "✓ %s 완료: %.1f초 (%.1f분) | 진행: %d/%d | 예상 남은 시간: %.1f분 (첫 모델 기준)",
          m, elapsed_sec, elapsed_sec / 60, completed_l2_methods, total_l2_methods, estimated_remaining / 60
        ))
      } else if (completed_l2_methods > 1 && remaining_methods > 0) {
        # Multiple methods completed: use average time
        avg_time_per_method <- elapsed_total / completed_l2_methods
        estimated_remaining <- avg_time_per_method * remaining_methods

        # Also show per-fold estimate if we have CV fold info
        # Note: caret doesn't expose fold-level timing, but we can estimate based on cv_folds
        if (!is.null(cv_folds) && cv_folds > 1) {
          # Estimate: if this method took elapsed_sec for cv_folds, each fold takes ~elapsed_sec/cv_folds
          # For remaining methods: (elapsed_sec/cv_folds) * cv_folds * remaining_methods = elapsed_sec * remaining_methods
          # But we use the average across completed methods for better accuracy
          avg_per_fold <- avg_time_per_method / cv_folds
          message(sprintf(
            "✓ %s 완료: %.1f초 (%.1f분) | 진행: %d/%d | 예상 남은 시간: %.1f분 (평균 %.1f초/모델, %.2f초/fold)",
            m, elapsed_sec, elapsed_sec / 60, completed_l2_methods, total_l2_methods,
            estimated_remaining / 60, avg_time_per_method, avg_per_fold
          ))
        } else {
          message(sprintf(
            "✓ %s 완료: %.1f초 (%.1f분) | 진행: %d/%d | 예상 남은 시간: %.1f분 (평균 %.1f초/모델)",
            m, elapsed_sec, elapsed_sec / 60, completed_l2_methods, total_l2_methods,
            estimated_remaining / 60, avg_time_per_method
          ))
        }
      } else {
        # Last method or no remaining methods
        message(sprintf(
          "✓ %s 완료: %.1f초 (%.1f분) | 진행: %d/%d | 완료!",
          m, elapsed_sec, elapsed_sec / 60, completed_l2_methods, total_l2_methods
        ))
      }
    }
  }

  # Final TML7 summary
  tml_end_time <- Sys.time()
  tml_total_elapsed <- as.numeric(difftime(tml_end_time, tml_start_time, units = "secs"))
  message(sprintf(
    "\n=== TML7 L2 학습 완료: 총 %d개 모델, %.1f분 소요 ===",
    completed_l2_methods, tml_total_elapsed / 60
  ))
  if (length(model_list) == 0) stop("No L2 models were successfully trained.")

  if (length(model_list) == 1) {
    best_name <- names(model_list)[1]
    best_fit <- model_list[[1]]
    best_val <- suppressWarnings(max(best_fit$results[[caret_metric]], na.rm = TRUE))
    message(sprintf(
      "Only one model trained. Selected '%s' (CV %s=%.4f).",
      best_name, caret_metric, best_val
    ))
    resamp <- NULL
  } else {
    resamp <- caret::resamples(model_list)
    vals <- sapply(model_list, function(f) suppressWarnings(max(f$results[[caret_metric]], na.rm = TRUE)))
    best_name <- names(which.max(vals))
    best_fit <- model_list[[best_name]]
    message(sprintf("Best model: %s (CV %s=%.4f).", best_name, caret_metric, max(vals, na.rm = TRUE)))
  }

  list(
    best_model = best_fit,
    best_model_name = best_name,
    best_metric_name = caret_metric,
    model_comparison = resamp,
    trained_models = model_list,
    positive_class = positive_class,
    l2_train = data.frame(l2_train_df, .target = l2_target),
    l1_signatures = l1_signatures,
    cv_folds = list(index = index, indexOut = indexOut)
  )
}


#' Derive gene-level importance from a trained meta learner
#'
#' Combines L2 signature importances with L1 signature weights to obtain
#' per-gene contribution scores. When the selected meta learner is a logistic
#' regression, signed contributions indicate whether higher gene expression
#' increases (`> 0`) or decreases (`< 0`) the odds of the positive class.
#' For non-linear models the magnitude is still meaningful, but signs should be
#' interpreted cautiously because the L2 model can capture non-linear effects.
#'
#' @param meta_result A result list returned by [TML7()].
#' @param normalize If `TRUE`, scale contributions based on `normalization_method`.
#' @param normalization_method Method for normalizing signature importances.
#'   Options: "max_abs" (default), "min_max", "softmax", "rank", "z_score".
#' @param target_model Optional. Name of the model to use (e.g. "ranger", "glm").
#'   If NULL, uses the best model from TML7.
#' @return A list with elements:
#'   \describe{
#'     \item{signature_importance}{Named numeric vector of L1 signature importances.}
#'     \item{gene_importance}{Data frame with columns `signature`, `gene`,
#'           `contribution` (signed), and `abs_contribution`.}
#'     \item{gene_summary}{Aggregated contributions per gene across signatures.}
#'     \item{positive_class}{Name of the positive class (usually `X2`).}
#'     \item{model_type}{Name of the selected L2 model.}
#'     \item{normalization_method}{Method used for normalization.}
#'   }
#' @export

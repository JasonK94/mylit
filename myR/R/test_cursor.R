# ============================================================================
# DEPRECATED: 이 파일의 모든 함수들은 test_analysis.R로 통합되었습니다.
# 최신 버전은 test_analysis.R의 다음 함수들을 사용하세요:
# - runNEBULA_v1 (또는 runNEBULA2_v1)
# - runMAST_v1
# - runMUSCAT_v5 (또는 runMUSCAT2_v1)
# ============================================================================

#' @title NEBULA 파이프라인 함수 (DEPRECATED)
#' @description 이 함수는 deprecated되었습니다. test_analysis.R의 runNEBULA_v1을 사용하세요.
#' @export 
runNEBULA_cursor_legacy <- function(...) {
  .Deprecated("runNEBULA_v1", package = "myR", msg = "runNEBULA is deprecated. Use runNEBULA_v1 from test_analysis.R instead.")
  if (exists("runNEBULA_v1", envir = asNamespace("myR"), inherits = FALSE)) {
    fun <- get("runNEBULA_v1", envir = asNamespace("myR"))
    return(fun(...))
  } else {
    stop("runNEBULA_v1 from test_analysis.R not found. Please ensure test_analysis.R is loaded.")
  }
}

#' @title MAST 파이프라인 함수 (DEPRECATED)
#' @description 이 함수는 deprecated되었습니다. test_analysis.R의 runMAST_v1을 사용하세요.
#' @export
runMAST_cursor_legacy <- function(...) {
  .Deprecated("runMAST_v1", package = "myR", msg = "runMAST_v1 in test_cursor.R is deprecated. Use runMAST_v1 from test_analysis.R instead.")
  if (exists("runMAST_v1", envir = asNamespace("myR"), inherits = FALSE)) {
    fun <- get("runMAST_v1", envir = asNamespace("myR"))
    return(fun(...))
  } else {
    stop("runMAST_v1 from test_analysis.R not found. Please ensure test_analysis.R is loaded.")
  }
}

#' @title Muscat 파이프라인 함수 v1 (DEPRECATED)
#' @description 이 함수는 deprecated되었습니다. test_analysis.R의 runMUSCAT_v5 또는 runMUSCAT2_v1을 사용하세요.
#' @export
runMUSCAT_cursor_legacy <- function(...) {
  .Deprecated("runMUSCAT_v5", package = "myR", msg = "runMUSCAT_v1 is deprecated. Use runMUSCAT_v5 or runMUSCAT2_v1 from test_analysis.R instead.")
  if (exists("runMUSCAT_v5", envir = asNamespace("myR"), inherits = FALSE)) {
    fun <- get("runMUSCAT_v5", envir = asNamespace("myR"))
    return(fun(...))
  } else {
    stop("runMUSCAT_v5 from test_analysis.R not found. Please ensure test_analysis.R is loaded.")
  }
}

#' @title MAST 파이프라인 함수 (DEPRECATED)
#' @description 이 함수는 deprecated되었습니다. test_analysis.R의 runMAST_v1을 사용하세요.
#' @export
runMAST <- function(...) {
  .Deprecated("runMAST_v1", package = "myR", msg = "runMAST is deprecated. Use runMAST_v1 from test_analysis.R instead.")
  if (exists("runMAST_v1", envir = asNamespace("myR"), inherits = FALSE)) {
    fun <- get("runMAST_v1", envir = asNamespace("myR"))
    return(fun(...))
  } else {
    stop("runMAST_v1 from test_analysis.R not found. Please ensure test_analysis.R is loaded.")
  }
}

#' @title Muscat 파이프라인 함수 (DEPRECATED)
#' @description 이 함수는 deprecated되었습니다. test_analysis.R의 runMUSCAT_v5 또는 runMUSCAT2_v1을 사용하세요.
#' @export
runMUSCAT <- function(...) {
  .Deprecated("runMUSCAT_v5", package = "myR", msg = "runMUSCAT is deprecated. Use runMUSCAT_v5 or runMUSCAT2_v1 from test_analysis.R instead.")
  if (exists("runMUSCAT_v5", envir = asNamespace("myR"), inherits = FALSE)) {
    fun <- get("runMUSCAT_v5", envir = asNamespace("myR"))
    return(fun(...))
  } else {
    stop("runMUSCAT_v5 from test_analysis.R not found. Please ensure test_analysis.R is loaded.")
  }
}

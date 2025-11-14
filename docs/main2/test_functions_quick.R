# ============================================================================
# 빠른 테스트 스크립트: 함수 시그니처 및 기본 검증
# ============================================================================
# 이 스크립트는 새로 개발된 함수들의 기본적인 검증을 수행합니다.
# 실제 데이터 없이 함수들이 제대로 로드되고 정의되어 있는지 확인합니다.
# ============================================================================

# 함수 소스 로드
source("myR/R/test_analysis.R")

# ============================================================================
# 1. 함수 존재 확인
# ============================================================================
message("========================================")
message("함수 존재 확인...")
message("========================================")

functions_to_check <- c("runMUSCAT2_v1", "runNEBULA2_v1", "runNEBULA2_v1_with_pseudobulk")

for (func_name in functions_to_check) {
  if (exists(func_name)) {
    message(sprintf("✓ %s: 존재함", func_name))
    
    # 함수 시그니처 확인
    func <- get(func_name)
    if (is.function(func)) {
      sig <- formals(func)
      message(sprintf("  - 파라미터 수: %d", length(sig)))
      message(sprintf("  - 파라미터: %s", paste(names(sig), collapse=", ")))
    }
  } else {
    message(sprintf("✗ %s: 존재하지 않음", func_name))
  }
}

# ============================================================================
# 2. 함수 타입 확인
# ============================================================================
message("\n========================================")
message("함수 타입 확인...")
message("========================================")

for (func_name in functions_to_check) {
  if (exists(func_name)) {
    func <- get(func_name)
    message(sprintf("%s:", func_name))
    message(sprintf("  - 타입: %s", class(func)[1]))
    message(sprintf("  - 함수인가: %s", is.function(func)))
    
    # 함수 본문 확인 (너무 길면 생략)
    body_text <- deparse(body(func))
    message(sprintf("  - 본문 라인 수: %d", length(body_text)))
  }
}

# ============================================================================
# 3. 기본 파라미터 검증
# ============================================================================
message("\n========================================")
message("기본 파라미터 검증...")
message("========================================")

# runMUSCAT2_v1 파라미터 확인
if (exists("runMUSCAT2_v1")) {
  muscat2_params <- formals(runMUSCAT2_v1)
  message("runMUSCAT2_v1 파라미터:")
  for (param_name in names(muscat2_params)) {
    param_value <- muscat2_params[[param_name]]
    if (is.null(param_value)) {
      message(sprintf("  - %s: NULL", param_name))
    } else if (is.name(param_value) && as.character(param_value) == "") {
      message(sprintf("  - %s: (필수)", param_name))
    } else {
      tryCatch({
        param_str <- paste(deparse(param_value), collapse = " ")
        message(sprintf("  - %s: %s", param_name, param_str))
      }, error = function(e) {
        message(sprintf("  - %s: (기본값)", param_name))
      })
    }
  }
}

# runNEBULA2_v1 파라미터 확인
if (exists("runNEBULA2_v1")) {
  nebula2_params <- formals(runNEBULA2_v1)
  message("\nrunNEBULA2_v1 파라미터:")
  for (param_name in names(nebula2_params)) {
    param_value <- nebula2_params[[param_name]]
    if (is.null(param_value)) {
      message(sprintf("  - %s: NULL", param_name))
    } else if (is.name(param_value) && as.character(param_value) == "") {
      message(sprintf("  - %s: (필수)", param_name))
    } else {
      tryCatch({
        param_str <- paste(deparse(param_value), collapse = " ")
        message(sprintf("  - %s: %s", param_name, param_str))
      }, error = function(e) {
        message(sprintf("  - %s: (기본값)", param_name))
      })
    }
  }
}

# runNEBULA2_v1_with_pseudobulk 파라미터 확인
if (exists("runNEBULA2_v1_with_pseudobulk")) {
  nebula2_pb_params <- formals(runNEBULA2_v1_with_pseudobulk)
  message("\nrunNEBULA2_v1_with_pseudobulk 파라미터:")
  for (param_name in names(nebula2_pb_params)) {
    param_value <- nebula2_pb_params[[param_name]]
    if (is.null(param_value)) {
      message(sprintf("  - %s: NULL", param_name))
    } else if (is.name(param_value) && as.character(param_value) == "") {
      message(sprintf("  - %s: (필수)", param_name))
    } else {
      tryCatch({
        param_str <- paste(deparse(param_value), collapse = " ")
        message(sprintf("  - %s: %s", param_name, param_str))
      }, error = function(e) {
        message(sprintf("  - %s: (기본값)", param_name))
      })
    }
  }
}

# ============================================================================
# 4. 함수 도움말 확인
# ============================================================================
message("\n========================================")
message("함수 도움말 확인...")
message("========================================")

for (func_name in functions_to_check) {
  if (exists(func_name)) {
    message(sprintf("\n%s 도움말:", func_name))
    tryCatch({
      help_text <- help(func_name, package = "myR", try.all.packages = FALSE)
      if (!is.null(help_text)) {
        message("  - 도움말 페이지 존재")
      }
    }, error = function(e) {
      message(sprintf("  - 도움말 페이지 없음 (로컬 함수이므로 정상)")
    })
  }
}

# ============================================================================
# 완료
# ============================================================================
message("\n========================================")
message("빠른 테스트 완료")
message("========================================")
message("\n다음 단계:")
message("1. 실제 데이터로 테스트: test_functions.R 실행")
message("2. 각 테스트 플래그를 TRUE로 변경")
message("3. 다운샘플링 데이터로 먼저 테스트")


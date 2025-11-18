# 데이터 전처리 테스트 스크립트
# 
# 이 스크립트는 데이터 전처리 함수를 테스트하고 전처리된 데이터를 저장합니다.

# 라이브러리 로드
if (!requireNamespace("qs", quietly = TRUE)) {
  stop("qs package is required. Install with: install.packages('qs')")
}

# 전처리 함수 로드
# Assume script is run from project root or scripts/pseudotime directory
if (file.exists("preprocess_data.R")) {
  source("preprocess_data.R")
} else if (file.exists("scripts/pseudotime/preprocess_data.R")) {
  source("scripts/pseudotime/preprocess_data.R")
} else {
  stop("Cannot find preprocess_data.R. Please run from project root or scripts/pseudotime directory.")
}

# --- 1. 다운샘플 데이터 전처리 ---
message("=== Testing preprocessing on downsampled data ===")
input_file_ds <- "/data/user3/sobj/IS_scvi_251107_ds2500.qs"
output_file_ds <- "/data/user3/sobj/IS_scvi_251107_ds2500_preprocessed.qs"

if (file.exists(input_file_ds)) {
  message("Loading downsampled data...")
  sobj_ds_preprocessed <- preprocess_pseudotime_data(
    input_file = input_file_ds,
    output_file = output_file_ds
  )
  message("Downsampled data preprocessing completed!")
  message("Saved to: ", output_file_ds)
} else {
  warning("Downsampled data file not found: ", input_file_ds)
}

# --- 2. 원본 데이터 전처리 ---
message("\n=== Testing preprocessing on full data ===")
input_file_full <- "/data/user3/sobj/IS_scvi_251107.qs"
output_file_full <- "/data/user3/sobj/IS_scvi_251107_preprocessed.qs"

if (file.exists(input_file_full)) {
  message("Loading full data...")
  sobj_full_preprocessed <- preprocess_pseudotime_data(
    input_file = input_file_full,
    output_file = output_file_full
  )
  message("Full data preprocessing completed!")
  message("Saved to: ", output_file_full)
} else {
  warning("Full data file not found: ", input_file_full)
}

message("\n=== IMPORTANT ===")
message("Please use the preprocessed files for all future analyses:")
message("  - Downsampled: ", output_file_ds)
message("  - Full: ", output_file_full)
message("==================\n")


#' 데이터 전처리 스크립트
#'
#' Seurat 객체의 metadata를 전처리하여 분석에 적합한 형태로 변환합니다.
#' 특히 sex 변수와 time 변수들을 올바르게 변환합니다.
#'
#' @param seurat_obj Seurat 객체 또는 qs 파일 경로
#' @param input_file Optional. qs 파일 경로 (seurat_obj가 경로인 경우)
#' @param output_file Optional. 전처리된 데이터를 저장할 경로
#'   기본값: `/data/user3/sobj/IS_scvi_251107_preprocessed.qs` (원본 파일명에 _preprocessed 추가)
#'
#' @return 전처리된 Seurat 객체
#'
#' @import Seurat
#' @importFrom lubridate parse_date_time ymd_hms as.duration
#' @export
#'
#' @examples
#' \dontrun{
#' # 방법 1: Seurat 객체 직접 전달
#' sobj_preprocessed <- preprocess_pseudotime_data(seurat_obj = sobj)
#'
#' # 방법 2: 파일 경로로 로드 및 전처리
#' sobj_preprocessed <- preprocess_pseudotime_data(
#'   input_file = "/data/user3/sobj/IS_scvi_251107_ds2500.qs",
#'   output_file = "/data/user3/sobj/IS_scvi_251107_ds2500_preprocessed.qs"
#' )
#' }
preprocess_pseudotime_data <- function(seurat_obj = NULL,
                                      input_file = NULL,
                                      output_file = NULL) {
  
  # --- 1. 데이터 로드 ---
  if (is.null(seurat_obj) && is.null(input_file)) {
    stop("Either seurat_obj or input_file must be provided.", call. = FALSE)
  }
  
  if (is.character(seurat_obj) && is.null(input_file)) {
    # seurat_obj가 경로인 경우
    input_file <- seurat_obj
    seurat_obj <- NULL
  }
  
  if (!is.null(input_file)) {
    message("Loading data from: ", input_file)
    if (!file.exists(input_file)) {
      stop("Input file not found: ", input_file, call. = FALSE)
    }
    
    if (requireNamespace("qs", quietly = TRUE)) {
      seurat_obj <- qs::qread(input_file)
    } else {
      seurat_obj <- readRDS(input_file)
    }
    
    # output_file이 지정되지 않았으면 자동 생성
    if (is.null(output_file)) {
      base_name <- tools::file_path_sans_ext(basename(input_file))
      dir_path <- dirname(input_file)
      output_file <- file.path(dir_path, paste0(base_name, "_preprocessed.qs"))
    }
  }
  
  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object.", call. = FALSE)
  }
  
  message("Starting data preprocessing...")
  message("Original cells: ", ncol(seurat_obj))
  message("Original genes: ", nrow(seurat_obj))
  
  # --- 2. sex 변수 전처리 (매우 중요) ---
  message("\n--- Step 1: Preprocessing sex variable ---")
  if ("sex" %in% colnames(seurat_obj@meta.data)) {
    sex_original <- seurat_obj@meta.data$sex
    message("Original sex values: ", paste(unique(sex_original), collapse = ", "))
    
    # 다양한 encoding을 M, F로 변환
    sex_clean <- as.character(sex_original)
    
    # 대소문자 통일 및 일반적인 변환
    sex_clean <- toupper(sex_clean)
    sex_clean[grepl("^M", sex_clean, ignore.case = TRUE)] <- "M"
    sex_clean[grepl("^F", sex_clean, ignore.case = TRUE)] <- "F"
    sex_clean[grepl("^1|MALE", sex_clean, ignore.case = TRUE)] <- "M"
    sex_clean[grepl("^2|FEMALE", sex_clean, ignore.case = TRUE)] <- "F"
    
    # NA 처리
    sex_clean[is.na(sex_clean) | sex_clean == "" | sex_clean == "NA"] <- NA_character_
    
    # 유효하지 않은 값 확인
    valid_sex <- sex_clean %in% c("M", "F", NA_character_)
    if (any(!valid_sex, na.rm = TRUE)) {
      invalid_values <- unique(sex_clean[!valid_sex])
      warning("Found invalid sex values: ", paste(invalid_values, collapse = ", "),
              ". These will be set to NA.", call. = FALSE)
      sex_clean[!valid_sex] <- NA_character_
    }
    
    seurat_obj@meta.data$sex <- factor(sex_clean, levels = c("M", "F"))
    message("Converted sex values: ", paste(levels(seurat_obj@meta.data$sex), collapse = ", "))
    message("Sex distribution: M=", sum(seurat_obj@meta.data$sex == "M", na.rm = TRUE),
            ", F=", sum(seurat_obj@meta.data$sex == "F", na.rm = TRUE),
            ", NA=", sum(is.na(seurat_obj@meta.data$sex)))
  } else {
    warning("'sex' column not found in metadata.", call. = FALSE)
  }
  
  # --- 3. g3 변수 전처리 (매우 중요) ---
  message("\n--- Step 2: Preprocessing g3 variable ---")
  if ("g3" %in% colnames(seurat_obj@meta.data)) {
    g3_original <- seurat_obj@meta.data$g3
    message("Original g3 class: ", class(g3_original))
    message("Original g3 values: ", paste(unique(g3_original[!is.na(g3_original)]), collapse = ", "))
    
    # factor로 변환 (numeric으로 인식되는 것을 방지)
    g3_clean <- factor(g3_original, levels = c("1", "2"))
    seurat_obj@meta.data$g3 <- g3_clean
    
    message("Converted g3 to factor with levels: ", paste(levels(g3_clean), collapse = ", "))
    message("G3 distribution: ", paste(table(g3_clean, useNA = "ifany"), collapse = ", "))
  } else {
    warning("'g3' column not found in metadata.", call. = FALSE)
  }
  
  # --- 4. 시간 변수 전처리 (매우 중요) ---
  message("\n--- Step 3: Preprocessing time variables ---")
  
  time_vars <- c("icu_adm_dt", "ia_start", "ia_end")
  time_parsed <- list()
  
  for (tv in time_vars) {
    if (tv %in% colnames(seurat_obj@meta.data)) {
      message("Processing ", tv, "...")
      time_raw <- seurat_obj@meta.data[[tv]]
      
      # 숫자 형식인지 확인
      if (is.numeric(time_raw)) {
        message("  ", tv, " is numeric. Attempting to parse as datetime...")
        # 숫자를 datetime으로 변환 시도 (여러 형식 시도)
        if (requireNamespace("lubridate", quietly = TRUE)) {
          # Unix timestamp인지 확인 (일반적으로 10자리 또는 13자리)
          if (all(time_raw > 1e9 & time_raw < 1e13, na.rm = TRUE)) {
            # Unix timestamp로 해석
            time_parsed[[tv]] <- lubridate::as_datetime(time_raw)
            message("  Parsed as Unix timestamp")
          } else {
            # Excel serial date로 해석 시도
            time_parsed[[tv]] <- lubridate::as_datetime((time_raw - 25569) * 86400, tz = "UTC")
            message("  Parsed as Excel serial date")
          }
        } else {
          # lubridate가 없으면 as.POSIXct 시도
          time_parsed[[tv]] <- as.POSIXct(time_raw, origin = "1970-01-01")
          message("  Parsed using as.POSIXct")
        }
      } else if (is.character(time_raw)) {
        message("  ", tv, " is character. Attempting to parse...")
        if (requireNamespace("lubridate", quietly = TRUE)) {
          # 여러 형식 시도
          time_parsed[[tv]] <- lubridate::parse_date_time(
            time_raw,
            orders = c("ymd HMS", "ymd HM", "ymd", "Ymd HMS", "Ymd HM", "Ymd",
                     "dmy HMS", "dmy HM", "dmy", "mdy HMS", "mdy HM", "mdy"),
            quiet = TRUE
          )
        } else {
          time_parsed[[tv]] <- as.POSIXct(time_raw)
        }
      } else if (inherits(time_raw, "POSIXct") || inherits(time_raw, "Date")) {
        message("  ", tv, " is already a datetime object.")
        time_parsed[[tv]] <- time_raw
      } else {
        warning("  Cannot parse ", tv, " (class: ", class(time_raw), "). Skipping.", call. = FALSE)
        next
      }
      
      # 변환된 시간 변수 저장
      seurat_obj@meta.data[[tv]] <- time_parsed[[tv]]
      message("  Successfully parsed ", tv, ". Valid values: ", 
              sum(!is.na(time_parsed[[tv]])), "/", length(time_parsed[[tv]]))
    }
  }
  
  # --- 5. 시간 차이 계산 ---
  message("\n--- Step 4: Calculating time differences ---")
  
  if (all(c("icu_adm_dt", "ia_start") %in% names(time_parsed))) {
    if (requireNamespace("lubridate", quietly = TRUE)) {
      diff_icu_ia_start <- lubridate::as.duration(time_parsed$ia_start - time_parsed$icu_adm_dt)
      seurat_obj@meta.data$icu_to_ia_start_hours <- as.numeric(diff_icu_ia_start, "hours")
      message("Calculated icu_to_ia_start_hours")
    } else {
      seurat_obj@meta.data$icu_to_ia_start_hours <- as.numeric(
        difftime(time_parsed$ia_start, time_parsed$icu_adm_dt, units = "hours")
      )
      message("Calculated icu_to_ia_start_hours (using base R)")
    }
  }
  
  if (all(c("icu_adm_dt", "ia_end") %in% names(time_parsed))) {
    if (requireNamespace("lubridate", quietly = TRUE)) {
      diff_icu_ia_end <- lubridate::as.duration(time_parsed$ia_end - time_parsed$icu_adm_dt)
      seurat_obj@meta.data$icu_to_ia_end_hours <- as.numeric(diff_icu_ia_end, "hours")
      message("Calculated icu_to_ia_end_hours")
    } else {
      seurat_obj@meta.data$icu_to_ia_end_hours <- as.numeric(
        difftime(time_parsed$ia_end, time_parsed$icu_adm_dt, units = "hours")
      )
      message("Calculated icu_to_ia_end_hours (using base R)")
    }
  }
  
  if (all(c("ia_start", "ia_end") %in% names(time_parsed))) {
    if (requireNamespace("lubridate", quietly = TRUE)) {
      diff_ia <- lubridate::as.duration(time_parsed$ia_end - time_parsed$ia_start)
      seurat_obj@meta.data$ia_duration_hours <- as.numeric(diff_ia, "hours")
      message("Calculated ia_duration_hours")
    } else {
      seurat_obj@meta.data$ia_duration_hours <- as.numeric(
        difftime(time_parsed$ia_end, time_parsed$ia_start, units = "hours")
      )
      message("Calculated ia_duration_hours (using base R)")
    }
  }
  
  # --- 6. 숫자 변수 확인 및 변환 ---
  message("\n--- Step 5: Validating numeric variables ---")
  numeric_vars <- c("nih_change", "arrival_gcs_score", "age", "ini_nih", "nih_end1")
  
  for (nv in numeric_vars) {
    if (nv %in% colnames(seurat_obj@meta.data)) {
      if (!is.numeric(seurat_obj@meta.data[[nv]])) {
        message("Converting ", nv, " to numeric...")
        seurat_obj@meta.data[[nv]] <- as.numeric(as.character(seurat_obj@meta.data[[nv]]))
      }
      message(nv, ": ", sum(!is.na(seurat_obj@meta.data[[nv]])), " valid values")
    }
  }
  
  # --- 7. 결과 저장 ---
  if (!is.null(output_file)) {
    message("\n--- Step 6: Saving preprocessed data ---")
    message("Saving to: ", output_file)
    
    # 디렉토리 생성
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
    
    if (requireNamespace("qs", quietly = TRUE)) {
      qs::qsave(seurat_obj, output_file)
      message("Data saved successfully using qs::qsave()")
    } else {
      saveRDS(seurat_obj, sub("\\.qs$", ".rds", output_file))
      message("Data saved successfully using saveRDS()")
    }
    
    message("\n=== IMPORTANT ===")
    message("Preprocessed data has been saved to: ", output_file)
    message("Please use this file for all future analyses.")
    message("The following variables have been preprocessed:")
    message("  - sex: Converted to factor (M, F)")
    message("  - g3: Converted to factor (1, 2)")
    message("  - Time variables: Parsed as datetime objects")
    message("  - Time differences: Calculated and added")
    message("==================\n")
  }
  
  message("Preprocessing completed successfully!")
  return(seurat_obj)
}


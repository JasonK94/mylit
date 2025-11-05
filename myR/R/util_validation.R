#' Input Validation Utilities
#'
#' This module provides common validation functions used across the package.
#' It standardizes error messages and validation logic.
#'
#' @name validation
NULL

#' @param sobj Seurat object
#' @param patient Column name for patient identifier
#' @param treatment Column name for treatment/drug identifier
#' @param timepoint Column name for timepoint identifier
#' @param tissue Column name for tissue identifier
#' @param check_key Column name for the key to check
#' @param check_values Values to check for parity
#' @return A list containing:
#'   \itemize{
#'     \item `detailed`: The original metadata with added diagnostic columns.
#'     \item `diag_table`: A one-row-per-group summary table.
#'     \item `message`: A summary message.
#'     \item `paired_ids`: A vector of `SegmentDisplayName`s for samples that passed the checks.
#'   }
#' @examples
#' \dontrun{
#' diagnosis_parity(sobj, patient="EMRID", treatment="drug", timepoint="pre_post", tissue="ck_str", check_key="pre_post", check_values=c("pre", "post"))
#' }
#' @export
diagnosis_parity <- function(sobj,
                                      patient=NULL, treatment=NULL, timepoint=NULL, tissue=NULL,   # 문자열 컬럼명 허용
                                      check_key, check_values,
                                      id_col="SegmentDisplayName") {
  # 동적 진단 열 이름: "pre","post" -> "n_pre","n_post"
  safe_counts_names <- paste0("n_", make.names(check_values, unique = TRUE))
  
  make_unique_id <- function(meta, patient=NULL, treatment=NULL, timepoint=NULL, tissue=NULL, check_key=NULL) {
    # unique_id를 만들 때 쓸 후보
    cols <- c(patient, treatment, timepoint, tissue)
  
    # check_key가 들어있으면 제외
    cols <- cols[!cols %in% check_key]
  
    # NULL 제거
    cols <- cols[!sapply(cols, is.null)]
  
    # 조합
    if (length(cols) == 0) {
      meta$unique_id <- seq_len(nrow(meta))  # 전부 NULL이면 그냥 row별 고유 번호
    } else {
      meta <- meta %>%
        mutate(across(all_of(cols), ~ tidyr::replace_na(as.character(.x), ""))) %>%
        mutate(unique_id = do.call(paste0, .[cols]))
    }
    meta
  }
  
  meta <- sobj@meta.data %>%
    make_unique_id(patient=patient,
                 treatment=treatment,
                 timepoint=timepoint,
                 tissue=tissue,
                 check_key=check_key) %>%
    mutate(
      .val_chr  = as.character(.data[[check_key]])
    )

  # unique_id별 전체 개수
  totals <- meta %>% count(unique_id, name = "total_n")

  # 선택값 카운트 (없어도 0으로 생성되도록 complete + pivot_wider)
  counts_sel <- meta %>%
    mutate(.val_chr = factor(.val_chr, levels = check_values)) %>%
    count(unique_id, .val_chr) %>%
    complete(unique_id, .val_chr, fill = list(n = 0)) %>%
    pivot_wider(names_from = .val_chr, values_from = n, values_fill = 0)

  # pivot에서 특정 값 컬럼 자체가 안 생긴 경우 대비해서 강제 추가
  for (v in check_values) {
    if (!hasName(counts_sel, v)) counts_sel[[v]] <- 0L
  }

  # 읽기 쉬운 이름으로 통일: pre -> n_pre, post -> n_post ...
  counts_sel <- counts_sel %>%
    rename_with(~ paste0("n_", make.names(.x, unique = TRUE)),
                .cols = all_of(check_values))

  # 진단 테이블 계산 (rowwise 없이)
  diag_table <- totals %>%
    left_join(counts_sel, by = "unique_id") %>%
    mutate(across(c(total_n, all_of(safe_counts_names)),
                  ~ replace_na(.x, 0))) %>%
    mutate(
      # 선택값 합계
      n_selected = purrr::reduce(across(all_of(safe_counts_names)), `+`),
      # 모든 선택값이 최소 1개 이상인가?
      has_all    = if_all(all_of(safe_counts_names), ~ .x > 0),
      # 선택값들의 개수가 모두 동일한가? (pre==post==... 형태)
      parity_eq  = pmap_lgl(across(all_of(safe_counts_names)),
                            ~ { v <- c(...); length(unique(v)) == 1 }),
      # 기타 값 개수
      n_other    = total_n - n_selected,
      # 최종 통과 조건(원래 의도 유지)
      ok         = has_all & parity_eq & n_other == 0
    ) %>%
    dplyr::select(unique_id, total_n, all_of(safe_counts_names),
           n_selected, n_other, has_all, parity_eq, ok)

  # 원본 롱형 테이블에 진단 열 조인
  detailed <- meta %>% left_join(diag_table, by = "unique_id")

  # 실패 id 메시지
  no_pass_ids <- diag_table %>% filter(!ok) %>% pull(unique_id)
  msg <- if (length(no_pass_ids) == 0) {
    "✅ All groups passed parity & exclusivity checks."
  } else {
    sprintf("checkout %s", paste(no_pass_ids, collapse = ", "))
  }
  message(msg)
  message("View($detailed)")
  message("enter View(sobj$diag_table) to check which row has the problem")
  ok_id=detailed%>%
    filter(ok)%>%
    pull(.data[[id_col]])
  
  list(
    detailed    = detailed,    # 원본 행 + n_pre/n_post/... + ok 등
    diag_table  = diag_table,  # id별 한 줄 요약
    message     = msg,
    paired_ids = ok_id
  )
}


#' Validate Seurat Object
#'
#' Checks if the input is a valid Seurat object with expected components.
#'
#' @param obj Object to validate
#' @param assay Optional assay name to check for
#' @param reduction Optional reduction name to check for
#' @param min_cells Minimum number of cells required
#' @param min_features Minimum number of features required
#'
#' @return TRUE if valid, stops with error message otherwise
#' @keywords internal
#' @export
validate_seurat <- function(obj, 
                            assay = NULL, 
                            reduction = NULL,
                            min_cells = 0,
                            min_features = 0) {
  
  # Check if it's a Seurat object
  if (!inherits(obj, "Seurat")) {
    stop("Input must be a Seurat object. Current class: ", 
         paste(class(obj), collapse = ", "))
  }
  
  # Check cell count
  n_cells <- ncol(obj)
  if (n_cells < min_cells) {
    stop("Seurat object has ", n_cells, " cells, but ", min_cells, " required")
  }
  
  # Check feature count
  n_features <- nrow(obj)
  if (n_features < min_features) {
    stop("Seurat object has ", n_features, " features, but ", min_features, " required")
  }
  
  # Check assay if specified
  if (!is.null(assay)) {
    available_assays <- Seurat::Assays(obj)
    if (!assay %in% available_assays) {
      stop("Assay '", assay, "' not found. Available assays: ", 
           paste(available_assays, collapse = ", "))
    }
  }
  
  # Check reduction if specified
  if (!is.null(reduction)) {
    available_reductions <- names(obj@reductions)
    if (!reduction %in% available_reductions) {
      # Try case-insensitive match
      matched_reduction <- available_reductions[tolower(available_reductions) == tolower(reduction)]
      if (length(matched_reduction) == 1) {
        message("Note: Using case-insensitive match: '", matched_reduction, "' for '", reduction, "'")
        return(TRUE)
      }
      stop("Reduction '", reduction, "' not found. Available reductions: ", 
           paste(available_reductions, collapse = ", "))
    }
  }
  
  return(TRUE)
}

#' Validate Metadata Column
#'
#' Checks if a metadata column exists and optionally validates its type.
#'
#' @param obj Seurat object or data frame
#' @param column_name Column name to validate
#' @param required_type Optional required type ("numeric", "factor", "character")
#' @param allow_na Whether NA values are allowed
#'
#' @return TRUE if valid, stops with error message otherwise
#' @keywords internal
#' @export
validate_metadata_column <- function(obj, 
                                     column_name, 
                                     required_type = NULL,
                                     allow_na = TRUE) {
  
  # Get metadata
  if (inherits(obj, "Seurat")) {
    metadata <- obj@meta.data
  } else if (is.data.frame(obj)) {
    metadata <- obj
  } else {
    stop("Input must be a Seurat object or data frame")
  }
  
  # Check if column exists
  if (!column_name %in% colnames(metadata)) {
    stop("Column '", column_name, "' not found in metadata. Available columns: ",
         paste(head(colnames(metadata), 10), collapse = ", "),
         if (ncol(metadata) > 10) "..." else "")
  }
  
  column_data <- metadata[[column_name]]
  
  # Check for NA values if not allowed
  if (!allow_na && any(is.na(column_data))) {
    n_na <- sum(is.na(column_data))
    stop("Column '", column_name, "' contains ", n_na, " NA values, but NA not allowed")
  }
  
  # Check type if specified
  if (!is.null(required_type)) {
    is_correct_type <- switch(
      required_type,
      "numeric" = is.numeric(column_data),
      "factor" = is.factor(column_data),
      "character" = is.character(column_data),
      stop("Unknown required_type: ", required_type)
    )
    
    if (!is_correct_type) {
      stop("Column '", column_name, "' must be of type '", required_type, 
           "', but is: ", paste(class(column_data), collapse = ", "))
    }
  }
  
  return(TRUE)
}

#' Validate Gene List
#'
#' Checks if genes exist in the object and optionally filters to valid genes.
#'
#' @param obj Seurat object or character vector of available genes
#' @param genes Character vector of genes to validate
#' @param min_present Minimum number of genes that must be present (default: all)
#' @param assay Assay to check genes in (for Seurat objects)
#' @param warn_missing Whether to warn about missing genes
#'
#' @return Character vector of valid genes present in the object
#' @keywords internal
#' @export
validate_genes <- function(obj, 
                          genes, 
                          min_present = NULL,
                          assay = NULL,
                          warn_missing = TRUE) {
  
  if (!is.character(genes) || length(genes) == 0) {
    stop("genes must be a non-empty character vector")
  }
  
  # Get available genes
  if (inherits(obj, "Seurat")) {
    if (is.null(assay)) {
      assay <- Seurat::DefaultAssay(obj)
    }
    available_genes <- rownames(obj[[assay]])
  } else if (is.character(obj)) {
    available_genes <- obj
  } else {
    stop("obj must be a Seurat object or character vector of gene names")
  }
  
  # Find present and missing genes
  genes_present <- intersect(genes, available_genes)
  genes_missing <- setdiff(genes, available_genes)
  
  # Set default minimum
  if (is.null(min_present)) {
    min_present <- length(genes)
  }
  
  # Check if enough genes are present
  if (length(genes_present) < min_present) {
    stop("Only ", length(genes_present), " of ", length(genes), 
         " genes found, but ", min_present, " required. ",
         "First missing genes: ", 
         paste(head(genes_missing, 5), collapse = ", "),
         if (length(genes_missing) > 5) "..." else "")
  }
  
  # Warn about missing genes if requested
  if (warn_missing && length(genes_missing) > 0) {
    warning("Genes not found (", length(genes_missing), "/", length(genes), "): ",
            paste(head(genes_missing, 10), collapse = ", "),
            if (length(genes_missing) > 10) "..." else "")
  }
  
  return(genes_present)
}

#' Validate Numeric Range
#'
#' Checks if a numeric parameter is within acceptable range.
#'
#' @param value Numeric value to validate
#' @param param_name Name of parameter (for error messages)
#' @param min Minimum allowed value (inclusive)
#' @param max Maximum allowed value (inclusive)
#' @param allow_na Whether NA is allowed
#'
#' @return TRUE if valid, stops with error message otherwise
#' @keywords internal
#' @export
validate_numeric_range <- function(value, 
                                   param_name, 
                                   min = -Inf, 
                                   max = Inf,
                                   allow_na = FALSE) {
  
  # Check for NA
  if (is.na(value)) {
    if (allow_na) {
      return(TRUE)
    } else {
      stop(param_name, " cannot be NA")
    }
  }
  
  # Check if numeric
  if (!is.numeric(value) || length(value) != 1) {
    stop(param_name, " must be a single numeric value")
  }
  
  # Check range
  if (value < min || value > max) {
    stop(param_name, " must be between ", min, " and ", max, 
         ", but is: ", value)
  }
  
  return(TRUE)
}

#' Validate Choice
#'
#' Validates that a value is one of allowed choices (like match.arg but more informative).
#'
#' @param value Value to validate
#' @param param_name Name of parameter (for error messages)
#' @param choices Vector of allowed values
#' @param multiple Whether multiple choices are allowed
#'
#' @return The validated value (or values if multiple=TRUE)
#' @keywords internal
#' @export
validate_choice <- function(value, param_name, choices, multiple = FALSE) {
  
  if (is.null(value)) {
    stop(param_name, " cannot be NULL")
  }
  
  if (multiple) {
    if (!all(value %in% choices)) {
      invalid <- setdiff(value, choices)
      stop(param_name, " contains invalid choices: ", 
           paste(invalid, collapse = ", "),
           ". Allowed choices: ", 
           paste(choices, collapse = ", "))
    }
  } else {
    if (length(value) != 1) {
      stop(param_name, " must be a single value, not ", length(value))
    }
    if (!value %in% choices) {
      stop(param_name, " must be one of: ", 
           paste(choices, collapse = ", "),
           ". Got: ", value)
    }
  }
  
  return(value)
}

#' Validate File Path
#'
#' Checks if a file path exists and optionally validates extension.
#'
#' @param path File path to validate
#' @param must_exist Whether file must already exist
#' @param extensions Optional vector of allowed extensions (e.g., c("csv", "txt"))
#' @param type Type of path ("file" or "directory")
#'
#' @return Normalized path if valid, stops with error otherwise
#' @keywords internal
#' @export
validate_path <- function(path, 
                         must_exist = TRUE, 
                         extensions = NULL,
                         type = c("file", "directory")) {
  
  type <- match.arg(type)
  
  if (!is.character(path) || length(path) != 1) {
    stop("path must be a single character string")
  }
  
  # Check existence if required
  if (must_exist) {
    if (type == "file") {
      if (!file.exists(path)) {
        stop("File not found: ", path)
      }
      if (dir.exists(path)) {
        stop("Expected a file but got a directory: ", path)
      }
    } else if (type == "directory") {
      if (!dir.exists(path)) {
        stop("Directory not found: ", path)
      }
    }
  }
  
  # Check extension if specified
  if (!is.null(extensions) && type == "file") {
    file_ext <- tools::file_ext(path)
    if (!file_ext %in% extensions) {
      stop("File extension must be one of: ", 
           paste(extensions, collapse = ", "),
           ". Got: ", file_ext)
    }
  }
  
  # Return normalized path
  return(normalizePath(path, mustWork = must_exist))
}

#' Create Informative Error Message
#'
#' Helper to create consistent, informative error messages.
#'
#' @param context Context where error occurred (e.g., function name)
#' @param message Main error message
#' @param suggestion Optional suggestion for fixing the error
#'
#' @return Formatted error message
#' @keywords internal
#' @export
create_error_message <- function(context, message, suggestion = NULL) {
  msg <- paste0("[", context, "] ", message)
  if (!is.null(suggestion)) {
    msg <- paste0(msg, "\nSuggestion: ", suggestion)
  }
  return(msg)
}

#' Check Package Dependencies
#'
#' Checks if required packages are installed and optionally loads them.
#'
#' @param packages Character vector of package names
#' @param load Whether to load the packages (default: FALSE)
#'
#' @return TRUE if all packages available, stops with error otherwise
#' @keywords internal
#' @export
check_packages <- function(packages, load = FALSE) {
  
  missing <- character(0)
  
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing <- c(missing, pkg)
    } else if (load) {
      library(pkg, character.only = TRUE)
    }
  }
  
  if (length(missing) > 0) {
    stop("Required packages not installed: ", 
         paste(missing, collapse = ", "),
         "\nInstall with: install.packages(c('", 
         paste(missing, collapse = "', '"), "'))")
  }
  
  return(TRUE)
}



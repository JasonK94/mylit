#' Pseudo-bulk 데이터 생성 및 edgeR 분석 준비 함수
#'
#' Seurat 객체로부터 pseudo-bulk 데이터를 생성하고, edgeR 분석을 위한
#' DGEList 객체, 디자인 매트릭스 등을 준비합니다.
#'
#' @param seurat_obj Seurat 객체
#' @param assay 사용할 assay 이름 (기본값: "SCT")
#' @param slot 사용할 slot 이름 (기본값: "counts"). AggregateExpression은 raw count 사용 권장.
#' @param sample_col 메타데이터에서 샘플을 구분하는 컬럼 이름
#' @param cluster_col 메타데이터에서 클러스터를 구분하는 컬럼 이름
#' @param group_col 메타데이터에서 비교할 그룹을 구분하는 컬럼 이름
#' @param min_count filterByExpr에서 사용할 최소 카운트 (기본값: 10)
#' @param norm_method calcNormFactors에서 사용할 정규화 방법 (기본값: "TMM")
#' @param design_formula 모델 디자인 ფორმულა. 예: ~ group_col 또는 ~ 0 + group_col.
#'                       NULL이면 기본값(~ group_col) 사용. 그룹 컬럼 이름이 group_col 변수값으로 자동 대체됩니다.
#' @param verbose 진행 상황 출력 여부 (기본값: TRUE)
#'
#' @return 리스트:
#'         - pb: Pseudo-bulk 발현 매트릭스 (genes x sample_cluster)
#'         - meta: Pseudo-bulk 샘플에 대한 메타데이터 (patient, ctype, group 포함)
#'         - dge: 필터링 및 정규화된 DGEList 객체
#'         - design: 생성된 디자인 매트릭스
#'         - contrast_levels: 그룹 컬럼의 레벨 (contrast 설정 시 참고용)
#'
#' @importFrom Seurat AggregateExpression VariableFeatures GetAssayData
#' @importFrom edgeR DGEList filterByExpr calcNormFactors estimateDisp glmQLFit glmQLFTest topTags
#' @importFrom dplyr %>% mutate select filter distinct left_join bind_rows pull
#' @importFrom tidyr separate
#' @importFrom tibble tibble rownames_to_column column_to_rownames
#' @importFrom stats model.matrix relevel
#'
#' @examples
#' \dontrun{
#' # 예제 데이터 생성 (실제 사용 시 실제 Seurat 객체로 대체)
#' counts <- matrix(rpois(10000, lambda = 10), ncol=100, nrow=100)
#' rownames(counts) <- paste0("Gene", 1:100)
#' colnames(counts) <- paste0("Cell", 1:100)
#' metadata <- data.frame(
#'   exp_sample_no = sample(paste0("Sample", 1:4), 100, replace = TRUE),
#'   annotation2_big = sample(paste0("Cluster", 1:3), 100, replace = TRUE),
#'   group3 = NA,
#'   row.names = paste0("Cell", 1:100)
#' )
#' metadata$group3[metadata$exp_sample_no %in% c("Sample1", "Sample2")] <- "Control"
#' metadata$group3[metadata$exp_sample_no %in% c("Sample3", "Sample4")] <- "Treatment"
#' seurat_obj_example <- CreateSeuratObject(counts = counts, meta.data = metadata)
#' # SCT assay가 없으므로 RNA assay 사용 예시
#' seurat_obj_example[["SCT"]] <- seurat_obj_example[["RNA"]] # 예시를 위해 복사
#'
#' prep_data <- prepare_pseudobulk_edgeR(
#'   seurat_obj = seurat_obj_example,
#'   assay = "SCT", # 또는 "RNA" 등
#'   slot = "counts",
#'   sample_col = "exp_sample_no",
#'   cluster_col = "annotation2_big",
#'   group_col = "group3"
#' )
#'
#' print(prep_data$contrast_levels) # 그룹 레벨 확인
#'
#' # edgeR 분석 수행 (예: 그룹 간 비교)
#' # contrast 설정 (예: Treatment vs Control)
#' # prep_data$contrast_levels 결과 보고 순서 맞게 설정
#' # 예를 들어 contrast_levels가 c("Control", "Treatment") 라면
#' my_contrast <- c(-1, 1)
#'
#' dge <- estimateDisp(prep_data$dge, prep_data$design)
#' fit <- glmQLFit(dge, prep_data$design)
#' qlf <- glmQLFTest(fit, contrast = my_contrast)
#' results <- topTags(qlf, n = Inf)$table
#' head(results)
#' }
prepare_pseudobulk_edgeR <- function(seurat_obj,
                                     assay = "SCT",
                                     slot = "counts",
                                     sample_col,
                                     cluster_col,
                                     group_col,
                                     min_count = 10,
                                     norm_method = "TMM",
                                     design_formula = NULL,
                                     verbose = TRUE) {
  
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat 패키지가 필요합니다.")
  if (!requireNamespace("edgeR", quietly = TRUE)) stop("edgeR 패키지가 필요합니다.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr 패키지가 필요합니다.")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("tidyr 패키지가 필요합니다.")
  if (!requireNamespace("tibble", quietly = TRUE)) stop("tibble 패키지가 필요합니다.")
  
  # 1. Check input columns exist in metadata
  meta_cols <- colnames(seurat_obj@meta.data)
  required_cols <- c(sample_col, cluster_col, group_col)
  if (!all(required_cols %in% meta_cols)) {
    stop("다음 컬럼이 Seurat 객체 메타데이터에 없습니다: ",
         paste(required_cols[!required_cols %in% meta_cols], collapse = ", "))
  }
  
  # Ensure group_col is a factor for reliable level ordering
  seurat_obj@meta.data[[group_col]] <- factor(seurat_obj@meta.data[[group_col]])
  
  # 2. Aggregate Expression
  if (verbose) message("1. Aggregating expression to pseudo-bulk...")
  pb <- AggregateExpression(
    seurat_obj,
    assays = assay,
    slot = slot,
    group.by = c(sample_col, cluster_col),
    return.seurat = FALSE
  )[[assay]] # Get the specific assay matrix
  
  if (nrow(pb) == 0 || ncol(pb) == 0) {
    stop("AggregateExpression 결과 pseudo-bulk 매트릭스가 비어있습니다. group.by 컬럼들을 확인하세요.")
  }
  
  # 3. Create Metadata for pseudo-bulk samples
  if (verbose) message("2. Creating metadata for pseudo-bulk samples...")
  # Determine separator based on column names (Seurat sometimes changes _)
  col_example <- colnames(pb)[1]
  sep <- ifelse(grepl("_", col_example), "_", "-") # Check if underscore exists
  
  meta_pb <- tryCatch({
    tibble(pb_sample_id = colnames(pb)) %>%
      tidyr::separate(pb_sample_id, into = c("patient", "ctype"), sep = sep, extra = "merge", remove = FALSE)
  }, error = function(e) {
    stop("Pseudo-bulk 컬럼 이름 ('", col_example, "')을 'patient'", sep, "'ctype' 형태로 분리할 수 없습니다. ",
         "AggregateExpression의 group.by 순서나 컬럼 내용을 확인하세요. 에러: ", e$message)
  })
  
  # Map group information from original Seurat metadata
  # sample_group_map 생성 로직 수정 시작
  
  # 원본 메타데이터에서 sample_col과 group_col을 가져옵니다.
  original_meta_subset <- seurat_obj@meta.data[, c(sample_col, group_col), drop = FALSE]
  
  # sample_col을 문자형으로 변환하고 앞뒤 공백을 제거합니다.
  original_meta_subset[[sample_col]] <- trimws(as.character(original_meta_subset[[sample_col]]))
  
  if (verbose) {
    message("Debug: Raw '", sample_col, "' values from Seurat metadata before 'g' prefixing (first 10): ", 
            paste(head(unique(original_meta_subset[[sample_col]]), 10), collapse=", "))
  }
  
  # AggregateExpression의 동작을 모방하여, 숫자처럼 보이는 sample_col 값에 'g'를 붙입니다.
  # 이 로직은 현재 디버그 결과(숫자형 ID에 'g'가 붙음)에 기반합니다.
  # 이미 'g'로 시작하는 경우 중복으로 붙지 않도록 합니다.
  needs_g_prefix_flags <- grepl("^[0-9]+$", original_meta_subset[[sample_col]]) & 
    !startsWith(original_meta_subset[[sample_col]], "g")
  
  if (any(needs_g_prefix_flags) && verbose) {
    message("Debug: Identified numeric-like values in '", sample_col, "' that need 'g' prefix. Example: '",
            original_meta_subset[[sample_col]][Position(function(x) needs_g_prefix_flags[which(original_meta_subset[[sample_col]] == x)[1]], original_meta_subset[[sample_col]])],
            "' will become 'g",
            original_meta_subset[[sample_col]][Position(function(x) needs_g_prefix_flags[which(original_meta_subset[[sample_col]] == x)[1]], original_meta_subset[[sample_col]])],
            "'.")
  }
  original_meta_subset[[sample_col]][needs_g_prefix_flags] <- paste0("g", original_meta_subset[[sample_col]][needs_g_prefix_flags])
  
  sample_group_map <- original_meta_subset %>% distinct()
  
  if (verbose) {
    message("Debug: Unique values in '", sample_col, "' for sample_group_map after potential 'g' prefixing (first 10): ", 
            paste(head(unique(sample_group_map[[sample_col]]), 10), collapse=", "))
  }
  # sample_group_map 생성 로직 수정 끝
  
  
  # join 키 설정: meta_pb의 'patient' 컬럼과 sample_group_map의 'sample_col' (이제 'g'가 붙은 ID)을 매칭합니다.
  join_key_map <- stats::setNames(sample_col, "patient")
  
  meta_pb <- meta_pb %>%
    left_join(sample_group_map, by = join_key_map)
  
  # ---- 디버깅을 위한 추가 출력 (기존과 동일하거나 유사하게 유지) ----
  if (sum(is.na(meta_pb[[group_col]]))) {
    message("Debug: WARNING - Group mapping failed for some samples AFTER 'g' prefix adjustment.")
    failed_patients <- meta_pb %>% filter(is.na(!!sym(group_col))) %>% pull(patient) %>% unique()
    message("Debug: Patient IDs from aggregated names that failed to map: ", paste(head(failed_patients, 10), collapse=", "))
  } else {
    message("Debug: Group mapping appears successful for all pseudo-bulk samples after 'g' prefix adjustment.")
  }
  
  # Ensure factor levels are consistent and set reference level if needed (optional, depends on formula)
  meta_pb[[group_col]] <- factor(meta_pb[[group_col]], levels = levels(seurat_obj@meta.data[[group_col]]))
  
  # Make sure row order of meta_pb matches column order of pb
  meta_pb <- meta_pb[match(colnames(pb), meta_pb$pb_sample_id),]
  rownames(meta_pb) <- meta_pb$pb_sample_id # Set rownames for edgeR convenience
  
  # 4. Prepare for edgeR
  if (verbose) message("3. Preparing DGEList and design matrix...")
  dge <- DGEList(counts = pb, group = meta_pb[[group_col]], samples = meta_pb)
  
  # Filtering
  keep <- filterByExpr(dge, group = meta_pb[[group_col]], min.count = min_count)
  if (verbose) message("   - Filtering: Kept ", sum(keep), " out of ", nrow(dge), " genes.")
  if (sum(keep) == 0) {
    warning("filterByExpr 결과 남은 유전자가 없습니다. min.count 값을 낮추거나 데이터를 확인하세요.")
  }
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  
  # Normalization
  dge <- calcNormFactors(dge, method = norm_method)
  
  # Design Matrix
  contrast_levels <- levels(meta_pb[[group_col]])
  if (is.null(design_formula)) {
    # Default design: ~ group_col (intercept + group effects)
    formula_str <- paste("~", group_col)
  } else {
    # User-provided formula string (replace 'group_col' placeholder if present)
    formula_str <- gsub("group_col", group_col, deparse(design_formula))
  }
  
  design <- tryCatch({
    model.matrix(as.formula(formula_str), data = meta_pb)
  }, error = function(e) {
    stop("디자인 매트릭스 생성 실패: ", e$message,
         "\n   Formula: ", formula_str,
         "\n   Metadata head:\n", paste(utils::capture.output(head(meta_pb)), collapse="\n"))
  })
  
  # Clean up design matrix column names (remove backticks or invalid chars)
  colnames(design) <- make.names(colnames(design))
  
  if (verbose) message("4. Preparation complete.")
  
  return(list(
    pb = pb,
    meta = meta_pb,
    dge = dge,
    design = design,
    contrast_levels = contrast_levels
  ))
}

prepare_pseudobulk_edgeR <- function(seurat_obj,
                                     assay = "SCT",
                                     slot = "counts",
                                     sample_col,
                                     cluster_col,
                                     group_col,
                                     min_count = 10,
                                     norm_method = "TMM",
                                     design_formula = NULL,
                                     verbose = TRUE) {
  
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat 패키지가 필요합니다.")
  if (!requireNamespace("edgeR", quietly = TRUE)) stop("edgeR 패키지가 필요합니다.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr 패키지가 필요합니다.")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("tidyr 패키지가 필요합니다.")
  if (!requireNamespace("tibble", quietly = TRUE)) stop("tibble 패키지가 필요합니다.")
  
  # 1. Check input columns exist in metadata
  meta_cols <- colnames(seurat_obj@meta.data)
  required_cols <- c(sample_col, cluster_col, group_col)
  if (!all(required_cols %in% meta_cols)) {
    stop("다음 컬럼이 Seurat 객체 메타데이터에 없습니다: ",
         paste(required_cols[!required_cols %in% meta_cols], collapse = ", "))
  }
  
  # Ensure group_col is a factor for reliable level ordering
  seurat_obj@meta.data[[group_col]] <- factor(seurat_obj@meta.data[[group_col]])
  
  # 2. Aggregate Expression
  if (verbose) message("1. Aggregating expression to pseudo-bulk...")
  pb <- AggregateExpression(
    seurat_obj,
    assays = assay,
    slot = slot,
    group.by = c(sample_col, cluster_col),
    return.seurat = FALSE
  )[[assay]] # Get the specific assay matrix
  
  if (nrow(pb) == 0 || ncol(pb) == 0) {
    stop("AggregateExpression 결과 pseudo-bulk 매트릭스가 비어있습니다. group.by 컬럼들을 확인하세요.")
  }
  
  # 3. Create Metadata for pseudo-bulk samples
  if (verbose) message("2. Creating metadata for pseudo-bulk samples...")
  # Determine separator based on column names (Seurat sometimes changes _)
  col_example <- colnames(pb)[1]
  sep <- ifelse(grepl("_", col_example), "_", "-") # Check if underscore exists
  
  meta_pb <- tryCatch({
    tibble(pb_sample_id = colnames(pb)) %>%
      tidyr::separate(pb_sample_id, into = c("patient", "ctype"), sep = sep, extra = "merge", remove = FALSE)
  }, error = function(e) {
    stop("Pseudo-bulk 컬럼 이름 ('", col_example, "')을 'patient'", sep, "'ctype' 형태로 분리할 수 없습니다. ",
         "AggregateExpression의 group.by 순서나 컬럼 내용을 확인하세요. 에러: ", e$message)
  })
  
  # Modify this section:
  # Map group information from original Seurat metadata
  # sample_group_map 생성 로직 수정 시작
  
  # 원본 메타데이터에서 sample_col과 group_col을 가져옵니다.
  # Ensure sample_col in original_meta_subset is character and trimmed
  original_meta_subset <- seurat_obj@meta.data %>%
    dplyr::select(dplyr::all_of(c(sample_col, group_col))) %>% # Removed cluster_col unless needed for debugging here
    dplyr::mutate(!!rlang::sym(sample_col) := trimws(as.character(.data[[sample_col]])))
  
  # Apply g-prefixing logic consistently to sample_col in original_meta_subset
  needs_g_prefix_flags <- grepl("^[0-9]+$", original_meta_subset[[sample_col]]) &
    !startsWith(original_meta_subset[[sample_col]], "g")
  
  if (any(needs_g_prefix_flags) && verbose) {
    message("Debug (prepare_pseudobulk_edgeR): Applying 'g' prefix to numeric-like values in '", sample_col, "' for sample_group_map creation.")
  }
  original_meta_subset[[sample_col]][needs_g_prefix_flags] <- paste0("g", original_meta_subset[[sample_col]][needs_g_prefix_flags])
  
  # Check for samples mapped to multiple distinct groups, which is problematic
  sample_to_multiple_groups_check <- original_meta_subset %>%
    dplyr::group_by(!!rlang::sym(sample_col)) %>%
    dplyr::summarise(n_distinct_groups = dplyr::n_distinct(!!rlang::sym(group_col)), .groups = "drop") %>%
    dplyr::filter(n_distinct_groups > 1)
  
  if(nrow(sample_to_multiple_groups_check) > 0 && verbose) {
    warning("prepare_pseudobulk_edgeR: Multiple distinct groups (in '", group_col,
            "') found for the same sample ID (in '", sample_col, "') in input Seurat object. Affected sample(s): ",
            paste(utils::head(sample_to_multiple_groups_check[[sample_col]], 5), collapse=", "),
            ". This can lead to ambiguous group assignments. `sample_group_map` will use the first encountered group for each such sample.")
  }
  
  # Create sample_group_map, ensuring one group per sample_col value
  sample_group_map <- original_meta_subset %>%
    dplyr::select(!!rlang::sym(sample_col), !!rlang::sym(group_col)) %>%
    dplyr::distinct(!!rlang::sym(sample_col), .keep_all = TRUE) # Critical: ensures one row per sample_col
  
  if (verbose) {
    message("Debug (prepare_pseudobulk_edgeR): `sample_group_map` created with ", nrow(sample_group_map), " unique '", sample_col, "' entries.")
    # message("Debug (prepare_pseudobulk_edgeR): Head of sample_group_map:\n", paste(utils::capture.output(utils::head(sample_group_map)), collapse="\n"))
    # message("Debug (prepare_pseudobulk_edgeR): Head of meta_pb before join (patient column):\n", paste(utils::capture.output(utils::head(meta_pb[, "patient", drop=FALSE])), collapse="\n"))
  }
  # sample_group_map 생성 로직 수정 끝
  
  # join 키 설정: meta_pb의 'patient' 컬럼과 sample_group_map의 'sample_col'
  join_key_map <- stats::setNames(sample_col, "patient") # meta_pb$patient = sample_group_map[[sample_col]]
  
  # Perform the join. If `sample_group_map` is correctly unique by `sample_col`,
  # this should be a many-meta_pb-rows-to-one-sample_group_map-row relationship
  meta_pb_joined <- meta_pb %>%
    dplyr::left_join(sample_group_map, by = join_key_map #, relationship = "many-to-one" # Can be explicit if dplyr version supports
    )
  
  meta_pb <- meta_pb_joined # Replace original meta_pb
  
  # ---- 디버깅을 위한 추가 출력 (기존과 동일하거나 유사하게 유지) ----
  if (sum(is.na(meta_pb[[group_col]]))) {
    message("Debug: WARNING - Group mapping resulted in NAs for some pseudo-bulk samples AFTER join modifications.")
    failed_patients <- meta_pb %>% dplyr::filter(is.na(!!rlang::sym(group_col))) %>% dplyr::pull(patient) %>% unique()
    message("Debug: Patient IDs from aggregated names (meta_pb$patient) that failed to map to a group: ", paste(utils::head(failed_patients, 10), collapse=", "))
    # Further debugging: check if these failed_patients exist in sample_group_map[[sample_col]]
    if (verbose && length(failed_patients) > 0) {
      unmapped_in_sgm <- setdiff(failed_patients, sample_group_map[[sample_col]])
      if (length(unmapped_in_sgm) > 0) {
        message("Debug: Of these failed patients, the following were NOT found in the prepared sample_group_map[['",sample_col,"']]: ", paste(utils::head(unmapped_in_sgm, 5), collapse=", "))
      }
      mapped_but_NA_group <- intersect(failed_patients, sample_group_map[[sample_col]][is.na(sample_group_map[[group_col]])])
      if (length(mapped_but_NA_group) > 0) {
        message("Debug: Of these failed patients, the following WERE found in sample_group_map but had an NA '",group_col,"' value: ", paste(utils::head(mapped_but_NA_group, 5), collapse=", "))
      }
    }
    
  } else {
    if (verbose) message("Debug (prepare_pseudobulk_edgeR): Group mapping appears successful for all pseudo-bulk samples after join modifications (no NAs introduced in group column).")
  }
  
  # Ensure factor levels are consistent and set reference level if needed (optional, depends on formula)
  # This re-factors based on levels from original Seurat object's group_col,
  # which might not match the actual levels present in meta_pb[[group_col]] if some groups are absent.
  # It's generally safer to factor based on actual values in meta_pb[[group_col]]
  # and ensure the reference level is set correctly if the design formula implies one.
  # The `cluster_pseudobulk_deg` already factors `temp_group_col` with GROUP2_REFERENCE as the ref.
  # Here, we just need to ensure it's a factor.
  if (!is.factor(meta_pb[[group_col]])) {
    meta_pb[[group_col]] <- factor(meta_pb[[group_col]])
  }
  # If a specific reference level from the original data is needed for other reasons:
  # meta_pb[[group_col]] <- factor(meta_pb[[group_col]], levels = levels(seurat_obj@meta.data[[group_col]]))
  
  
  # Make sure row order of meta_pb matches column order of pb
  meta_pb <- meta_pb[match(colnames(pb), meta_pb$pb_sample_id),]
  
  # Replace `rownames(meta_pb) <- meta_pb$pb_sample_id` to avoid deprecation warning
  if ("pb_sample_id" %in% colnames(meta_pb)) {
    meta_pb <- tibble::column_to_rownames(meta_pb, var = "pb_sample_id")
  } else if (verbose) {
    warning("prepare_pseudobulk_edgeR: 'pb_sample_id' column not found in meta_pb after join. Cannot set rownames for DGEList samples.")
  }
  
  # 4. Prepare for edgeR
  if (verbose) message("3. Preparing DGEList and design matrix...")
  dge <- DGEList(counts = pb, group = meta_pb[[group_col]], samples = meta_pb)
  
  # Filtering
  keep <- filterByExpr(dge, group = meta_pb[[group_col]], min.count = min_count)
  if (verbose) message("   - Filtering: Kept ", sum(keep), " out of ", nrow(dge), " genes.")
  if (sum(keep) == 0) {
    warning("filterByExpr 결과 남은 유전자가 없습니다. min.count 값을 낮추거나 데이터를 확인하세요.")
  }
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  
  # Normalization
  dge <- calcNormFactors(dge, method = norm_method)
  
  # Design Matrix
  contrast_levels <- levels(meta_pb[[group_col]])
  if (is.null(design_formula)) {
    # Default design: ~ group_col (intercept + group effects)
    formula_str <- paste("~", group_col)
  } else {
    # User-provided formula string (replace 'group_col' placeholder if present)
    formula_str <- gsub("group_col", group_col, deparse(design_formula))
  }
  
  design <- tryCatch({
    model.matrix(as.formula(formula_str), data = meta_pb)
  }, error = function(e) {
    stop("디자인 매트릭스 생성 실패: ", e$message,
         "\n   Formula: ", formula_str,
         "\n   Metadata head:\n", paste(utils::capture.output(head(meta_pb)), collapse="\n"))
  })
  
  # Clean up design matrix column names (remove backticks or invalid chars)
  colnames(design) <- make.names(colnames(design))
  
  if (verbose) message("4. Preparation complete.")
  
  return(list(
    pb = pb,
    meta = meta_pb,
    dge = dge,
    design = design,
    contrast_levels = contrast_levels
  ))
}

#' Pseudo-bulk Differential Gene Expression Analysis
#'
#' 주어진 Seurat 객체 또는 준비된 pseudo-bulk 데이터에 대해
#' 다양한 수준의 edgeR 기반 DEG 분석을 수행합니다.
#'
#' @param analysis_level 분석 수준:
#'                       - "overall": 모든 클러스터를 통합하여 전체 그룹 간 비교 (코드 1 방식).
#'                       - "per_cluster": 각 클러스터 내에서 그룹 간 비교 (코드 2 방식).
#'                       - "specific_cluster": 지정된 단일 클러스터 내에서 그룹 간 비교.
#' @param contrast edgeR의 glmQLFTest에 전달할 contrast 벡터. 그룹 레벨 순서에 맞게 지정해야 함.
#'                 예: 그룹 레벨이 c("Control", "Treatment")일 때 Treatment vs Control 비교는 c(-1, 1).
#'                 prepare_pseudobulk_edgeR 결과의 contrast_levels 참고.
#' @param target_cluster analysis_level이 "specific_cluster"일 때 분석할 클러스터 이름.
#' @param min_samples_per_group_cluster 클러스터별 분석 시, 각 그룹이 가져야 하는 최소 샘플 수.
#'                                      미달 시 해당 클러스터 분석 건너뜀 (기본값: 2).
#' @param ... prepare_pseudobulk_edgeR 함수에 전달할 인자들 (seurat_obj, assay, slot, sample_col, cluster_col, group_col 등)
#'            또는 prepare_pseudobulk_edgeR 함수의 반환값인 리스트 (pb, meta, dge, design 포함).
#' @param verbose 진행 상황 출력 여부 (기본값: TRUE)
#'
#' @return DEG 결과 테이블 (tibble).
#'         - "overall" 분석: 단일 tibble.
#'         - "per_cluster" 분석: 각 클러스터 결과가 합쳐진 tibble (cluster 컬럼 포함).
#'         - "specific_cluster" 분석: 해당 클러스터 결과만 담긴 tibble (cluster 컬럼 포함).
#'         결과 테이블에는 logFC, logCPM, F, PValue, FDR 등의 edgeR 결과와 gene 컬럼이 포함됩니다.
#'
#' @importFrom edgeR estimateDisp glmQLFit glmQLFTest topTags
#' @importFrom dplyr %>% filter bind_rows mutate rename select group_by summarise n case_when pull
#' @importFrom tibble tibble rownames_to_column
#' @importFrom stats setNames
#'
#' @examples
#' \dontrun{
#' # 예제 데이터 준비 (prepare_pseudobulk_edgeR 예제 참고)
#' counts <- matrix(rpois(10000, lambda = 10), ncol=100, nrow=100)
#' rownames(counts) <- paste0("Gene", 1:100)
#' colnames(counts) <- paste0("Cell", 1:100)
#' metadata <- data.frame(
#'   sample_id = sample(paste0("Sample", 1:4), 100, replace = TRUE),
#'   cell_type = sample(paste0("Cluster", 1:3), 100, replace = TRUE),
#'   condition = NA,
#'   row.names = paste0("Cell", 1:100)
#' )
#' metadata$condition[metadata$sample_id %in% c("Sample1", "Sample2")] <- "Control"
#' metadata$condition[metadata$sample_id %in% c("Sample3", "Sample4")] <- "Treatment"
#' seurat_obj_example <- CreateSeuratObject(counts = counts, meta.data = metadata)
#' seurat_obj_example[["SCT"]] <- seurat_obj_example[["RNA"]] # 예시용
#'
#' # 사용례 1: 전체 그룹 간 비교 (Overall)
#' overall_results <- run_pseudobulk_deg(
#'   analysis_level = "overall",
#'   contrast = c(-1, 1), # Treatment vs Control (Control이 첫번째 레벨이라고 가정)
#'   seurat_obj = seurat_obj_example,
#'   assay = "SCT", slot = "counts",
#'   sample_col = "sample_id", cluster_col = "cell_type", group_col = "condition"
#' )
#' head(overall_results)
#'
#' # 사용례 2: 클러스터별 그룹 간 비교 (Per Cluster)
#' per_cluster_results <- run_pseudobulk_deg(
#'   analysis_level = "per_cluster",
#'   contrast = c(-1, 1),
#'   seurat_obj = seurat_obj_example,
#'   assay = "SCT", slot = "counts",
#'   sample_col = "sample_id", cluster_col = "cell_type", group_col = "condition"
#' )
#' head(per_cluster_results)
#' table(per_cluster_results$cluster) # 각 클러스터별 결과 확인
#'
#' # 사용례 3: 특정 클러스터("Cluster2") 내 그룹 간 비교 (Specific Cluster)
#' specific_cluster_results <- run_pseudobulk_deg(
#'   analysis_level = "specific_cluster",
#'   target_cluster = "Cluster2",
#'   contrast = c(-1, 1),
#'   seurat_obj = seurat_obj_example,
#'   assay = "SCT", slot = "counts",
#'   sample_col = "sample_id", cluster_col = "cell_type", group_col = "condition"
#' )
#' head(specific_cluster_results)
#' }
run_pseudobulk_deg <- function(analysis_level = c("overall", "per_cluster", "specific_cluster"),
                               contrast,
                               target_cluster = NULL, # 기본값을 NULL로 유지
                               min_samples_per_group_cluster = 2,
                               min_count = 10,
                               ...,
                               verbose = TRUE) {
  
  analysis_level <- match.arg(analysis_level)
  dot_args <- list(...)
  
  # --- 1. 데이터 준비 (prepare_pseudobulk_edgeR 호출 또는 기존 데이터 사용) ---
  # <<< START MODIFIED SECTION >>>
  prep_data <- NULL
  group_col <- NULL
  
  if ("prepared_data" %in% names(dot_args) &&
      is.list(dot_args$prepared_data) &&
      all(c("pb", "meta", "dge", "design", "contrast_levels") %in% names(dot_args$prepared_data))) {
    
    if (verbose) message("Using explicitly provided 'prepared_data' argument.")
    prep_data <- dot_args$prepared_data
    # Extract necessary components from prep_data
    pb   <- prep_data$pb
    meta <- prep_data$meta # This meta is meta_pb from prepare_pseudobulk_edgeR
    dge  <- prep_data$dge
    design <- prep_data$design
    contrast_levels <- prep_data$contrast_levels
    
    # Determine group_col from the design formula or meta structure
    # This logic attempts to find the name of the grouping column used to create the design matrix.
    # It's often the primary variable in the design formula (excluding intercept).
    formula_terms <- try(terms(design), silent = TRUE) # Get terms from design matrix attributes
    if (inherits(formula_terms, "try-error") || is.null(attr(formula_terms, "term.labels"))) {
      # Fallback: try to infer from meta columns (excluding known ones)
      group_col_inferred <- setdiff(colnames(meta), c("pb_sample_id", "patient", "ctype"))[1]
      if (length(group_col_inferred) == 1 && nzchar(group_col_inferred) && group_col_inferred %in% colnames(meta)) {
        group_col <- group_col_inferred
        if (verbose) message("Inferred 'group_col' from meta structure: ", group_col)
      } else {
        stop("Cannot reliably determine grouping column from 'prepared_data$design' attributes or 'prepared_data$meta' structure.")
      }
    } else {
      term_labels <- attr(formula_terms, "term.labels")
      # Remove interactions or function calls to get base variable names
      term_labels_cleaned <- gsub("factor\\((.*?)\\)", "\\1", term_labels) # Remove factor()
      term_labels_cleaned <- gsub(":.*", "", term_labels_cleaned) # Remove interaction terms part
      term_labels_cleaned <- unique(term_labels_cleaned)
      
      # Try to find a match in meta columns
      possible_group_cols <- intersect(term_labels_cleaned, colnames(meta))
      if (length(possible_group_cols) == 1) {
        group_col <- possible_group_cols[1]
        if (verbose) message("Determined 'group_col' from design formula: ", group_col)
      } else if (length(possible_group_cols) > 1) {
        group_col <- possible_group_cols[1] # Default to the first one if multiple
        if (verbose) message("Multiple possible 'group_col's from design formula (", paste(possible_group_cols, collapse=", "), "), using first: ", group_col)
      } else {
        # If still not found, use the previous fallback
        group_col_inferred <- setdiff(colnames(meta), c("pb_sample_id", "patient", "ctype"))[1]
        if (length(group_col_inferred) == 1 && nzchar(group_col_inferred) && group_col_inferred %in% colnames(meta)) {
          group_col <- group_col_inferred
          if (verbose) message("Inferred 'group_col' from meta structure (fallback): ", group_col)
        } else {
          stop("Cannot determine grouping column from 'prepared_data'. Ensure design matrix attributes are set or meta structure is as expected.")
        }
      }
    }
    if (is.null(group_col) || !group_col %in% colnames(meta)) {
      stop(paste0("Failed to determine a valid 'group_col' from 'prepared_data' that exists in 'prepared_data$meta'. Determined: '",
                  ifelse(is.null(group_col), "NULL", group_col), "'. Meta columns: ", paste(colnames(meta), collapse=", ")))
    }
    
    
  } else if (all(c("pb", "meta", "dge", "design", "contrast_levels") %in% names(dot_args))) {
    # This case is when pb, meta, etc., are passed as top-level arguments in ...
    if (verbose) message("Using pre-prepared pseudo-bulk data components passed directly via ...")
    prep_data <- dot_args # The whole dot_args is the list of components
    pb   <- prep_data$pb
    meta <- prep_data$meta
    dge  <- prep_data$dge
    design <- prep_data$design
    contrast_levels <- prep_data$contrast_levels
    group_col <- setdiff(colnames(meta), c("pb_sample_id", "patient", "ctype"))[1] # Best guess
    if(length(group_col) != 1 || !group_col %in% colnames(meta)) {
      group_col <- names(attr(design, "contrasts"))[1] # Another guess
      if(is.null(group_col)) stop("Cannot determine grouping column from direct components.")
    }
  } else if ("seurat_obj" %in% names(dot_args)) {
    if (verbose) message("Preparing pseudo-bulk data using prepare_pseudobulk_edgeR...")
    req_args_prep <- c("seurat_obj", "sample_col", "cluster_col", "group_col")
    if (!all(req_args_prep %in% names(dot_args))) {
      stop("Seurat 객체를 사용할 경우, 다음 인수가 필요합니다: ", paste(req_args_prep, collapse=", "))
    }
    group_col <- dot_args$group_col # Crucial: get group_col for prepare_pseudobulk_edgeR
    
    all_args_for_prep <- c(dot_args, list(verbose = verbose))
    # Pass min_count from run_pseudobulk_deg to prepare_pseudobulk_edgeR if not already in dot_args
    if (!"min_count" %in% names(all_args_for_prep)) {
      all_args_for_prep$min_count <- min_count # Use min_count from run_pseudobulk_deg's signature
    }
    
    prep_data_list_internal <- do.call(prepare_pseudobulk_edgeR, all_args_for_prep)
    pb   <- prep_data_list_internal$pb
    meta <- prep_data_list_internal$meta
    dge  <- prep_data_list_internal$dge
    design <- prep_data_list_internal$design
    contrast_levels <- prep_data_list_internal$contrast_levels
    # group_col is already defined from dot_args$group_col for this path
  } else {
    stop("Seurat 객체(seurat_obj) 또는 준비된 pseudo-bulk 데이터 리스트('prepared_data' argument or direct components)를 제공해야 합니다.")
  }
  
  if (missing(contrast)) stop("'contrast' 인수는 필수입니다.")
  
  # --- 2. 분석 수준에 따른 DEG 수행 ---
  results <- NULL
  loop_over_these_clusters <- NULL # 반복할 클러스터 목록을 저장할 변수
  
  if (analysis_level == "overall") {
    if (verbose) message("Performing 'overall' DEG analysis...")
    # contrast 유효성 검사는 이전에 있었음 (필요시 여기에 더 명시적으로 추가 가능)
    dge_overall <- estimateDisp(dge, design) # dge, design은 전체 데이터에 대한 것
    fit_overall <- glmQLFit(dge_overall, design)
    qlf_overall <- glmQLFTest(fit_overall, contrast = contrast)
    results <- topTags(qlf_overall, n = Inf)$table %>%
      rownames_to_column("gene") %>%
      tibble()
    if (verbose) message("  Overall analysis complete. Found ", nrow(results), " genes.")
    
  } else if (analysis_level == "per_cluster") {
    # 'per_cluster'는 meta 데이터에 있는 모든 고유한 클러스터에 대해 수행
    # target_cluster 인수는 이 경우 무시됨
    all_meta_ctypes_cleaned <- trimws(as.character(meta$ctype)) # meta$ctype은 하이픈 사용 (AggregateExpression 결과)
    loop_over_these_clusters <- unique(all_meta_ctypes_cleaned)
    
    if (verbose) {
      message("Performing 'per_cluster' analysis. Iterating over clusters found in data: '",
              paste(sort(loop_over_these_clusters), collapse="', '"), "'")
    }
    
  } else if (analysis_level == "specific_cluster") {
    # 'specific_cluster'는 제공된 단일 target_cluster에 대해서만 수행
    if (is.null(target_cluster) || length(target_cluster) != 1 || !is.character(target_cluster)) {
      stop("'specific_cluster' 분석 시 'target_cluster' 인수를 단일 클러스터 이름(문자열)으로 정확히 제공해야 합니다.")
    }
    
    target_cluster_cleaned <- trimws(as.character(target_cluster))
    # AggregateExpression 출력 형식(하이픈)에 맞추기 위해 밑줄(_)을 하이픈(-)으로 변경
    target_cluster_normalized <- gsub("_", "-", target_cluster_cleaned)
    
    all_meta_ctypes_cleaned <- trimws(as.character(meta$ctype))
    unique_meta_ctypes_cleaned <- unique(all_meta_ctypes_cleaned)
    
    if (verbose) {
      message("Debug specific_cluster: Original target_cluster supplied (cleaned): '", target_cluster_cleaned, "'")
      message("Debug specific_cluster: Normalized target_cluster for matching (hyphenated): '", target_cluster_normalized, "'")
      message("Debug specific_cluster: Available unique ctype values in meta (cleaned): '", paste(sort(unique_meta_ctypes_cleaned), collapse="', '"), "'")
    }
    
    if (!target_cluster_normalized %in% unique_meta_ctypes_cleaned) {
      stop("지정된 'target_cluster' (정규화 후: '", target_cluster_normalized, 
           "')가 메타데이터의 ctype 컬럼에 없습니다. 확인된 ctype 값 (정리 후): '", 
           paste(sort(unique_meta_ctypes_cleaned), collapse="', '"), "'")
    }
    loop_over_these_clusters <- target_cluster_normalized # 단일 클러스터 이름
    if (verbose) {
      message("Performing 'specific_cluster' analysis for target cluster: '", loop_over_these_clusters, "'")
    }
  }
  
  # --- 3. "per_cluster" 또는 "specific_cluster" 경우 실제 루프 실행 ---
  if (analysis_level == "per_cluster" || analysis_level == "specific_cluster") {
    if (is.null(loop_over_these_clusters) || length(loop_over_these_clusters) == 0) {
      warning("분석할 클러스터가 결정되지 않았거나 없습니다.")
      return(tibble()) # 빈 tibble 반환
    }
    
    res_list <- list()
    for (cl_name in loop_over_these_clusters) { # 루프 변수 이름 변경 (예: cl_name)
      if (verbose) message("  Processing cluster: ", cl_name)
      
      # meta$ctype은 이미 하이픈으로 정리된 형태, cl_name도 하이픈으로 정리된 형태임
      ix <- meta$ctype == cl_name 
      if (sum(ix) == 0) {
        warning("클러스터 '", cl_name, "'에 해당하는 pseudo-bulk 샘플이 없습니다. 건너<0xEB>니다.")
        next
      }
      pb_sub <- pb[, ix, drop = FALSE]
      md_sub <- meta[ix, , drop = FALSE]
      
      sample_counts <- md_sub %>% count(!!sym(group_col))
      if(any(sample_counts$n < min_samples_per_group_cluster)) {
        warning("클러스터 '", cl_name, "'는 그룹별 최소 샘플 수(", min_samples_per_group_cluster, ")를 만족하지 못합니다: ",
                paste(sample_counts[[group_col]], "=", sample_counts$n, collapse=", "), ". 건너<0xEB>니다.")
        next
      }
      md_sub[[group_col]] <- factor(md_sub[[group_col]], levels = contrast_levels)
      
      dge_sub <- DGEList(counts = pb_sub, group = md_sub[[group_col]], samples = md_sub)
      
      # 클러스터 부분집합에 대한 디자인 매트릭스 생성
      # group_col 인수가 올바르게 함수 범위 내에 있는지 확인 (함수 인수로 직접 받거나 dot_args에서 가져옴)
      formula_sub_str <- paste("~", group_col) 
      design_sub <- tryCatch({
        model.matrix(as.formula(formula_sub_str), data = md_sub)
      }, error = function(e) {
        warning("클러스터 '", cl_name, "'에 대한 디자인 매트릭스 생성 실패: ", e$message, ". 건너<0xEB>니다.")
        return(NULL) 
      })
      if(is.null(design_sub)) next
      colnames(design_sub) <- make.names(colnames(design_sub))
      
      # 이 부분에서 min_count는 run_pseudobulk_deg 함수의 인수로 받은 min_count를 사용합니다.
      keep_sub <- filterByExpr(dge_sub, design = design_sub, group = md_sub[[group_col]], min.count = min_count) 
      if (verbose) message("    - Filtering for '", cl_name, "': Kept ", sum(keep_sub), " out of ", nrow(dge_sub), " genes.")
      if(sum(keep_sub) == 0) {
        warning("클러스터 '", cl_name, "'에서 filterByExpr 후 남은 유전자가 없습니다. 건너<0xEB>니다.")
        next
      }
      dge_sub <- dge_sub[keep_sub, , keep.lib.sizes = FALSE]
      dge_sub <- calcNormFactors(dge_sub)
      
      # contrast 유효성 검사 (클러스터 부분집합의 디자인 매트릭스에 대해)
      expected_contrast_len_sub <- ncol(design_sub) - (!grepl("~ 0 +|~0+", deparse(attributes(design_sub)$terms)))
      if (length(contrast) != expected_contrast_len_sub) {
        warning("클러스터 '", cl_name, "' 분석에서 contrast 길이(", length(contrast), 
                ")가 디자인 매트릭스 컬럼 수(", ncol(design_sub), 
                ", 인터셉트 고려 시 ", expected_contrast_len_sub,")와 다를 수 있습니다. ",
                "Design columns: ", paste(colnames(design_sub), collapse=", "))
      }
      
      dge_sub_disp <- tryCatch({ estimateDisp(dge_sub, design_sub) }, 
                               warning = function(w){
                                 warning("Dispersion estimation 경고 (클러스터 ", cl_name,"): ", w$message); estimateDisp(dge_sub, design_sub, robust=TRUE) }, 
                               error = function(e){ warning("Dispersion estimation 에러 (클러스터 ", cl_name,"): ", e$message, ". 건너<0xEB>니다."); return(NULL) })
      if(is.null(dge_sub_disp)) next
      
      fit_sub <- tryCatch({ glmQLFit(dge_sub_disp, design_sub) }, 
                          error = function(e){ warning("glmQLFit 에러 (클러스터 ", cl_name, "): ", e$message, ". 건너<0xEB>니다."); return(NULL) })
      if(is.null(fit_sub)) next
      
      qlf_sub <- tryCatch({ glmQLFTest(fit_sub, contrast = contrast) }, 
                          error = function(e){ warning("glmQLFTest 에러 (클러스터 ", cl_name, "): ", e$message, ". 건너<0xEB>니다."); return(NULL) })
      if(is.null(qlf_sub)) next
      
      res_sub <- topTags(qlf_sub, n = Inf)$table %>%
        rownames_to_column("gene") %>%
        # mutate의 cluster 값에 현재 루프의 클러스터 이름(cl_name) 사용
        mutate(cluster = cl_name, .before = 1) %>% 
        tibble()
      
      res_list[[cl_name]] <- res_sub
      if (verbose) message("    - Cluster '", cl_name, "' analysis complete. Found ", nrow(res_sub), " genes.")
    } # end for loop
    
    if (length(res_list) > 0) {
      results <- bind_rows(res_list)
    } else {
      warning("분석 수준 '", analysis_level, "'에 대한 결과를 생성하지 못했습니다 (처리된 클러스터 없음).")
      results <- tibble() 
    }
  } # end if per_cluster or specific_cluster
  
  return(results)
}





#' @title 유전자 발현 데이터에 대한 그룹별 또는 전체 선형 회귀 분석 수행 (최적화됨)
#' @description
#' 이 함수는 Seurat 객체, 유전자 목록, (선택적) 그룹 지정 메타데이터 열,
#' 샘플 식별 메타데이터 열, 그리고 수치형 예측 변수 메타데이터 열을 입력으로 받습니다.
#' 먼저 전체 샘플에 대한 유사벌크 평균 발현량을 계산하고 메타데이터와 병합합니다.
#' 그 후, 각 그룹별(또는 전체)로 각 유전자에 대해 수치형 예측 변수를 사용한 선형 회귀를 수행합니다.
#'
#' @param sobj Seurat 객체.
#' @param genes 분석할 유전자 이름의 문자형 벡터.
#' @param sample_col 샘플 또는 유사벌크 단위를 식별하는 메타데이터 열 이름. 기본값은 "sample".
#' @param numeric_predictor 선형 모델에서 예측 변수(독립 변수, X)로 사용될 수치형
#'   메타데이터 열 이름. 기본값은 "severity_score".
#' @param group_col (선택 사항) 그룹별 분석을 위한 메타데이터 열 이름.
#'   `NULL` (기본값)이면 전체 샘플에 대해 분석을 수행합니다.
#' @param p_adjust_method `stats::p.adjust`에 사용될 p-값 보정 방법. 기본값은 "BH".
#' @param min_samples_per_group 그룹별 또는 전체 분석 시 필요한 최소 고유 샘플 수. 기본값은 3.
#' @param min_distinct_predictor_values 회귀 분석을 위해 필요한 `numeric_predictor`의 최소 고유 값 수. 기본값은 2.
#'
#' @return 각 유전자(및 각 그룹)에 대한 선형 모델 결과 (절편, 기울기, 기울기 표준오차,
#'   모델에 사용된 샘플 수, p-값, 보정된 p-값, R-제곱)를 포함하는 데이터 프레임.
#'   `group_col`이 지정된 경우, 결과에 'group' 열이 추가됩니다.
#'
#' @importFrom dplyr group_by summarise across all_of select distinct left_join arrange bind_rows filter n_distinct ungroup
#' @importFrom rlang sym !! :=
#' @importFrom stats lm coef summary.lm p.adjust as.formula na.omit
#' @import Seurat
#' @importFrom tibble rownames_to_column
#'
#' @export
pseudobulk_linear_fit <- function(sobj,
                                  genes,
                                  sample_col = "sample",
                                  numeric_predictor = "severity_score",
                                  group_col = NULL,
                                  p_adjust_method = "BH",
                                  min_samples_per_group = 3,
                                  min_distinct_predictor_values = 2) {
  
  # --- 1. 입력값 유효성 검사 ---
  if (!inherits(sobj, "Seurat")) stop("`sobj`는 Seurat 객체여야 합니다.")
  if (!is.character(genes) || length(genes) == 0) stop("`genes`는 하나 이상의 유전자 이름을 포함하는 문자형 벡터여야 합니다.")
  if (!all(genes %in% rownames(sobj))) {
    missing_genes <- genes[!genes %in% rownames(sobj)]
    stop(paste0("일부 유전자가 Seurat 객체에 없습니다: ", paste(missing_genes, collapse=", ")))
  }
  meta_cols <- colnames(sobj@meta.data)
  if (!sample_col %in% meta_cols) stop(paste0("`sample_col` '", sample_col, "'이 메타데이터에 없습니다."))
  if (!numeric_predictor %in% meta_cols) stop(paste0("`numeric_predictor` '", numeric_predictor, "'이 메타데이터에 없습니다."))
  if (!is.null(group_col) && !group_col %in% meta_cols) stop(paste0("`group_col` '", group_col, "'이 메타데이터에 없습니다."))
  
  predictor_data_orig <- sobj@meta.data[[numeric_predictor]]
  if (!is.numeric(predictor_data_orig)) {
    warning_msg <- paste0("`numeric_predictor` '", numeric_predictor, "'이(가) ", class(predictor_data_orig)[1], "형입니다. 수치형으로 변환을 시도합니다.")
    original_NAs <- is.na(predictor_data_orig)
    if (is.factor(predictor_data_orig)) {
      converted_numeric <- suppressWarnings(as.numeric(as.character(predictor_data_orig)))
    } else if (is.character(predictor_data_orig)) {
      converted_numeric <- suppressWarnings(as.numeric(predictor_data_orig))
    } else {
      stop(paste0("`numeric_predictor` '", numeric_predictor, "'의 타입(", class(predictor_data_orig)[1], ")을 수치형으로 자동 변환할 수 없습니다."))
    }
    new_NAs_introduced <- any(is.na(converted_numeric) & !original_NAs)
    if (all(is.na(converted_numeric)) && !all(original_NAs)) stop(paste0("`numeric_predictor` 변환 실패: 모든 값이 NA가 되었습니다."))
    if (new_NAs_introduced) warning(paste0(warning_msg, " 일부 값이 NA로 변환되었습니다.")) else warning(warning_msg)
    sobj@meta.data[[numeric_predictor]] <- converted_numeric
    if (!is.numeric(sobj@meta.data[[numeric_predictor]])) stop(paste0("`numeric_predictor`를 수치형으로 최종 변환 실패."))
  }
  
  # --- 2. 사전 그룹별/전체 샘플 수 및 예측 변수 유효성 검사 ---
  meta_df_check <- sobj@meta.data %>%
    select(all_of(c(sample_col, numeric_predictor, group_col))) %>%
    filter(!is.na(.data[[numeric_predictor]])) 
  
  valid_groups <- c()
  if (is.null(group_col)) {
    n_unique_samples <- n_distinct(meta_df_check[[sample_col]])
    n_unique_preds <- n_distinct(meta_df_check[[numeric_predictor]])
    if (n_unique_samples >= min_samples_per_group && n_unique_preds >= min_distinct_predictor_values) {
      valid_groups <- "all_samples"
    } else {
      stop(paste0("전체 샘플에 대해 유효한 샘플 수(", n_unique_samples, "<", min_samples_per_group, ") 또는 예측 변수 고유값 수(", n_unique_preds, "<", min_distinct_predictor_values,")가 부족합니다."))
    }
  } else {
    meta_df_check_grouped <- meta_df_check %>%
      group_by(!!sym(group_col)) %>%
      summarise(
        n_unique_samples = n_distinct(.data[[sample_col]]),
        n_unique_preds = n_distinct(.data[[numeric_predictor]]),
        .groups = "drop"
      ) %>%
      filter(n_unique_samples >= min_samples_per_group, n_unique_preds >= min_distinct_predictor_values)
    
    valid_groups <- unique(as.character(meta_df_check_grouped[[group_col]]))
    all_defined_groups <- unique(as.character(sobj@meta.data[[group_col]]))
    # NA 그룹 ID는 제외하고 비교
    invalid_groups <- setdiff(all_defined_groups[!is.na(all_defined_groups)], valid_groups)
    if (length(invalid_groups) > 0) {
      warning(paste0("다음 그룹들은 샘플 수 또는 예측 변수 다양성 부족으로 분석에서 제외됩니다: ", paste(invalid_groups, collapse=", ")))
    }
    if (length(valid_groups) == 0) {
      stop("모든 그룹이 분석을 위한 최소 샘플 수 또는 예측 변수 다양성 기준을 충족하지 못합니다.")
    }
  }
  
  # --- 3. 전체 샘플에 대한 유사벌크 및 메타데이터 준비 ---
  # GetAssayData의 colnames는 sobj의 colnames와 일치
  cell_barcodes <- colnames(sobj) 
  # 분석할 유전자만, 그리고 현재 Seurat 객체에 존재하는 cell만 필터링
  expr_data <- GetAssayData(sobj, assay = DefaultAssay(sobj), slot = "data")[genes, cell_barcodes, drop = FALSE]
  
  expr_df_transposed <- as.data.frame(t(as.matrix(expr_data))) 
  # expr_df_transposed의 행 이름(cell_barcodes)을 사용하여 sample_col 정보 매칭
  expr_df_transposed[[sample_col]] <- sobj@meta.data[rownames(expr_df_transposed), sample_col]
  
  avg_expr_all_samples <- expr_df_transposed %>%
    group_by(!!sym(sample_col)) %>%
    summarise(across(all_of(genes), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  
  cols_to_select <- c(sample_col, numeric_predictor)
  if (!is.null(group_col)) cols_to_select <- c(cols_to_select, group_col)
  
  meta_data_full <- sobj@meta.data %>%
    select(all_of(unique(cols_to_select))) %>% # unique() to handle cases where sample_col could be same as group_col
    distinct()
  
  df_merged_full <- left_join(avg_expr_all_samples, meta_data_full, by = sample_col)
  
  if (!is.numeric(df_merged_full[[numeric_predictor]])) {
    df_merged_full[[numeric_predictor]] <- as.numeric(as.character(df_merged_full[[numeric_predictor]]))
    if (all(is.na(df_merged_full[[numeric_predictor]]))) {
      stop("최종 병합 데이터에서 numeric_predictor가 모두 NA가 되었습니다.")
    }
  }
  
  # --- 4. 내부 함수: 단일 그룹 (또는 전체)에 대한 분석 수행 ---
  .run_lm_on_merged_data <- function(current_df_merged, current_genes, current_numeric_predictor, current_sample_col_name, min_samples_thresh, min_preds_thresh) {
    if (nrow(current_df_merged) == 0) return(NULL)
    
    results_list <- lapply(current_genes, function(gene) {
      formula_str <- paste0("`", gene, "` ~ `", current_numeric_predictor, "`")
      model_data <- current_df_merged[, c(gene, current_numeric_predictor, current_sample_col_name), drop = FALSE] # sample_col도 포함하여 고유 샘플 수 확인
      model_data_complete <- stats::na.omit(model_data) # NA 제거
      
      # NA 제거 후 실제 사용될 샘플 수 및 예측 변수 다양성 확인
      n_samples_for_model <- n_distinct(model_data_complete[[current_sample_col_name]])
      n_distinct_preds_in_model_data <- n_distinct(model_data_complete[[current_numeric_predictor]])
      
      if (n_samples_for_model < min_samples_thresh || n_distinct_preds_in_model_data < min_preds_thresh) {
        # 이 경고는 이미 사전 체크에서 대부분 걸러졌겠지만, 유전자별 NA 때문에 추가 발생 가능
        warning(paste0("유전자 '", gene, "'에 대해 NA 제거 후 유효 데이터 부족 (샘플수:", n_samples_for_model, 
                       ", 예측변수 고유값:", n_distinct_preds_in_model_data, ")으로 모델을 적합할 수 없습니다."))
        return(data.frame(
          gene = gene, intercept = NA_real_, slope = NA_real_, slope_se = NA_real_,
          n_samples_in_model = n_samples_for_model, p_value = NA_real_, r_squared = NA_real_
        ))
      }
      
      model <- lm(as.formula(formula_str), data = model_data_complete) # model_data_complete에서 sample_col 제외하고 사용해도 무방
      summary_model <- summary(model)
      coefs <- coef(summary_model)
      
      intercept_val <- NA_real_; slope_val <- NA_real_; slope_se_val <- NA_real_; p_val <- NA_real_
      if ("(Intercept)" %in% rownames(coefs)) intercept_val <- coefs["(Intercept)", "Estimate"]
      
      predictor_in_model <- names(coef(model))[2]
      if (!is.na(predictor_in_model) && predictor_in_model %in% rownames(coefs)) {
        slope_val <- coefs[predictor_in_model, "Estimate"]
        slope_se_val <- coefs[predictor_in_model, "Std. Error"]
        p_val <- coefs[predictor_in_model, "Pr(>|t|)"]
      } else if (nrow(coefs) >= 2) {
        warning(paste0("유전자 '", gene, "' 모델에서 예측변수 '", current_numeric_predictor, 
                       "'의 계수 이름을 명시적으로 찾지 못했습니다. 모델의 두 번째 계수를 사용합니다. (실제 모델 계수명: '", rownames(coefs)[2], "')"))
        slope_val <- coefs[2, "Estimate"]; slope_se_val <- coefs[2, "Std. Error"]; p_val <- coefs[2, "Pr(>|t|)"]
      }
      
      data.frame(
        gene = gene, intercept = intercept_val, slope = slope_val, slope_se = slope_se_val,
        n_samples_in_model = n_samples_for_model, # 여기서는 n_distinct(model_data_complete[[current_sample_col_name]]) 사용
        p_value = p_val, r_squared = summary_model$r.squared
      )
    })
    bind_rows(results_list)
  }
  
  # --- 5. 분석 실행 ---
  final_results_list <- list()
  for (grp_val in valid_groups) {
    current_data_subset <- df_merged_full
    current_group_name_for_output <- grp_val # 기본값
    
    if (!is.null(group_col)) { # 그룹별 분석일 때
      if (grp_val == "all_samples" && !"all_samples" %in% unique(sobj@meta.data[[group_col]])) {
        # group_col이 NULL이 아니지만 valid_groups에 "all_samples"가 들어간 비정상적 상황 방지
        # (이런 경우는 위 로직상 발생 안해야함)
        # 이 부분은 group_col이 있을 때 "all_samples"라는 그룹명이 실제 메타데이터에 있을 경우를 위한 것.
        # 일반적으로 group_col이 지정되면 valid_groups는 실제 그룹명만 포함.
        current_data_subset <- df_merged_full %>% filter(!!sym(group_col) == grp_val)
      } else if (grp_val != "all_samples") {
        current_data_subset <- df_merged_full %>% filter(!!sym(group_col) == grp_val)
      }
      # "all_samples"는 group_col이 NULL일 때만 사용되는 플레이스홀더
    } else { # group_col이 NULL일 때, grp_val은 "all_samples" 임
      current_group_name_for_output <- "all_samples"
    }
    
    # 최종 데이터셋에 대한 마지막 체크 (이미 위에서 필터링 했지만, 안전장치)
    if(n_distinct(na.omit(current_data_subset[[numeric_predictor]])) < min_distinct_predictor_values ||
       n_distinct(na.omit(current_data_subset[[sample_col]])) < min_samples_per_group ) {
      warning(paste0(ifelse(is.null(group_col), "전체 데이터셋", paste0("그룹 '", grp_val, "'")),
                     "에서 최종 분석 데이터의 예측변수 다양성 또는 샘플 수가 부족하여 건너뜁니다."))
      next 
    }
    
    group_results_df <- .run_lm_on_merged_data(current_data_subset, genes, numeric_predictor, sample_col, min_samples_per_group, min_distinct_predictor_values)
    if (!is.null(group_results_df) && nrow(group_results_df) > 0) {
      group_results_df$group <- current_group_name_for_output
      final_results_list[[current_group_name_for_output]] <- group_results_df # 리스트 이름도 고유하게
    }
  }
  final_results_df <- bind_rows(final_results_list)
  
  # --- 6. 결과 정리 및 보정된 p-값 계산 ---
  if (is.null(final_results_df) || nrow(final_results_df) == 0) {
    warning("분석 결과 데이터 프레임이 비어있습니다.")
    return(data.frame(group=character(), gene=character(), intercept=numeric(), slope=numeric(), 
                      slope_se=numeric(), n_samples_in_model=integer(), p_value=numeric(), 
                      adj_p_value=numeric(), r_squared=numeric()))
  }
  
  valid_p_indices <- !is.na(final_results_df$p_value)
  if (any(valid_p_indices)) {
    final_results_df$adj_p_value <- NA_real_
    final_results_df$adj_p_value[valid_p_indices] <- stats::p.adjust(
      final_results_df$p_value[valid_p_indices], method = p_adjust_method
    )
  } else {
    final_results_df$adj_p_value <- NA_real_
  }
  
  cols_ordered <- c("group", "gene", "intercept", "slope", "slope_se", "n_samples_in_model", "p_value", "adj_p_value", "r_squared")
  # group 열이 없는 경우 (is.null(group_col)이었고, all_samples로 추가 안 했다면) 대비 -> group 열은 항상 존재하게 수정됨
  final_results_df <- final_results_df[, intersect(cols_ordered, names(final_results_df)), drop = FALSE]
  
  return(final_results_df)
}




#' @title 그룹 간 선형 회귀 기울기 사후 비교
#' @description
#' `pseudobulk_linear_fit` 함수의 결과를 사용하여 각 유전자에 대해 그룹 간 기울기를 비교합니다.
#' 두 그룹의 경우 t-검정을 수행하고, 세 그룹 이상인 경우 모든 쌍에 대해 t-검정을 수행하고 p-값을 보정합니다.
#'
#' @param results_df `pseudobulk_linear_fit`에서 반환된 데이터 프레임.
#' @param gene_col `results_df`에서 유전자 이름을 포함하는 열의 이름. 기본값 "gene".
#' @param group_col `results_df`에서 그룹 식별자를 포함하는 열의 이름. 기본값 "group".
#' @param slope_col `results_df`에서 기울기 값을 포함하는 열의 이름. 기본값 "slope".
#' @param se_col `results_df`에서 기울기의 표준 오차를 포함하는 열의 이름. 기본값 "slope_se".
#' @param n_samples_col `results_df`에서 모델 피팅에 사용된 샘플 수를 포함하는 열의 이름. 기본값 "n_samples_in_model".
#' @param p_adjust_method 다중 비교를 위한 p-값 보정 방법. 기본값 "BH".
#' @param adjustment_scope p-값 보정 범위. "global" (모든 비교에 대해 전체적으로 보정) 또는
#'                         "per_gene" (각 유전자 내의 비교들에 대해서만 보정). 기본값 "global".
#'
#' @return 각 유전자 내 그룹 쌍 간의 기울기 비교 결과를 담은 데이터 프레임.
#'
#' @importFrom dplyr group_by do filter arrange mutate ungroup slice
#' @importFrom rlang sym !! .data
#' @importFrom utils combn
#' @importFrom stats pt p.adjust
#'
#' @export
post_hoc_slope_comparison <- function(results_df,
                                      gene_col = "gene",
                                      group_col = "group",
                                      slope_col = "slope",
                                      se_col = "slope_se",
                                      n_samples_col = "n_samples_in_model",
                                      p_adjust_method = "BH",
                                      adjustment_scope = "global") { # 새 파라미터 및 기본값
  
  # --- 입력값 및 필수 열 확인 ---
  required_cols <- c(gene_col, group_col, slope_col, se_col, n_samples_col)
  if (!all(required_cols %in% names(results_df))) {
    missing_cols <- required_cols[!required_cols %in% names(results_df)]
    stop(paste0("필수 열 중 일부가 results_df에 없습니다: ", paste(missing_cols, collapse=", ")))
  }
  
  if (!adjustment_scope %in% c("global", "per_gene")) {
    stop("`adjustment_scope`는 'global' 또는 'per_gene' 이어야 합니다.")
  }
  
  # --- 중복된 (유전자-그룹) 항목 처리 ---
  results_df_deduplicated <- results_df %>%
    group_by(!!sym(gene_col), !!sym(group_col)) %>%
    slice(1) %>% 
    ungroup()
  
  if(nrow(results_df_deduplicated) < nrow(results_df)){
    warning("입력 `results_df`에 중복된 유전자-그룹 항목이 발견되었습니다. 각 조합의 첫 번째 항목을 사용합니다.")
  }
  
  # --- NA 값 및 유효하지 않은 샘플 수/SE 처리 ---
  results_df_complete <- results_df_deduplicated %>%
    filter(
      !is.na(.data[[slope_col]]), 
      !is.na(.data[[se_col]]),
      .data[[se_col]] > 0, 
      !is.na(.data[[n_samples_col]]),
      .data[[n_samples_col]] >= 3 
    )
  
  if (nrow(results_df_complete) == 0) {
    warning("NA 값, 유효하지 않은 SE 또는 샘플 수 필터링 후 비교할 데이터가 없습니다.")
    return(data.frame(gene=character(), group1=character(), group2=character(), 
                      slope1=numeric(), slope2=numeric(), se1=numeric(), se2=numeric(),
                      n1=integer(), n2=integer(), t_statistic=numeric(), df=numeric(),
                      p_value_raw=numeric(), adj_p_value=numeric()))
  }
  
  # --- 유전자별 그룹 간 비교 수행 ---
  # comparison_results_from_do 변수명 변경 (이전 답변에서 _list 로 끝났었음)
  comparison_results_from_do <- results_df_complete %>%
    group_by(!!sym(gene_col)) %>%
    do({
      gene_data_for_current_gene <- .
      
      if (nrow(gene_data_for_current_gene) < 2) {
        return(data.frame()) 
      }
      
      group_combinations_matrix <- combn(gene_data_for_current_gene[[group_col]], 2)
      
      individual_pair_results_list_for_gene <- lapply(1:ncol(group_combinations_matrix), function(col_idx) {
        g1_name <- group_combinations_matrix[1, col_idx]
        g2_name <- group_combinations_matrix[2, col_idx]
        
        d1 <- gene_data_for_current_gene[gene_data_for_current_gene[[group_col]] == g1_name, ]
        d2 <- gene_data_for_current_gene[gene_data_for_current_gene[[group_col]] == g2_name, ]
        
        b1 <- d1[[slope_col]]; se1 <- d1[[se_col]]; n1 <- d1[[n_samples_col]]
        b2 <- d2[[slope_col]]; se2 <- d2[[se_col]]; n2 <- d2[[n_samples_col]]
        
        df_reg1 <- n1 - 2 
        df_reg2 <- n2 - 2
        
        t_stat <- NA_real_; df_ws <- NA_real_; p_val <- NA_real_
        
        t_stat <- (b1 - b2) / sqrt(se1^2 + se2^2)
        numerator_ws <- (se1^2 + se2^2)^2
        denominator_ws <- (se1^4 / df_reg1) + (se2^4 / df_reg2) 
        
        df_ws <- numerator_ws / denominator_ws
        p_val <- 2 * stats::pt(abs(t_stat), df = df_ws, lower.tail = FALSE) 
        
        data.frame(
          group1 = g1_name, group2 = g2_name,
          slope1 = b1, slope2 = b2,
          se1 = se1, se2 = se2,
          n1 = n1, n2 = n2,
          t_statistic = t_stat, df = df_ws, p_value_raw = p_val,
          stringsAsFactors = FALSE
        )
      }) 
      
      if (length(individual_pair_results_list_for_gene) > 0) {
        do.call(rbind, individual_pair_results_list_for_gene)
      } else {
        data.frame() 
      }
    }) %>% 
    ungroup() %>% 
    filter(nrow(.) > 0) 
  
  if (nrow(comparison_results_from_do) == 0) { # 변수명 수정
    warning("최종적으로 비교 결과가 생성되지 않았습니다.")
    return(data.frame(gene=character(), group1=character(), group2=character(), 
                      slope1=numeric(), slope2=numeric(), se1=numeric(), se2=numeric(),
                      n1=integer(), n2=integer(), t_statistic=numeric(), df=numeric(),
                      p_value_raw=numeric(), adj_p_value=numeric()))
  }
  
  # p-값 보정 (adjustment_scope에 따라 범위 결정)
  if (adjustment_scope == "per_gene") {
    final_adjusted_results <- comparison_results_from_do %>% # 변수명 수정
      group_by(!!sym(gene_col)) %>%
      mutate(adj_p_value = stats::p.adjust(.data$p_value_raw, method = p_adjust_method)) %>%
      ungroup()
  } else { # "global" (기본값)
    final_adjusted_results <- comparison_results_from_do %>% # 변수명 수정
      mutate(adj_p_value = stats::p.adjust(.data$p_value_raw, method = p_adjust_method))
  }
  
  final_adjusted_results <- final_adjusted_results %>%
    # adj_p_value가 NA가 될 수 있으므로, NA를 마지막으로 정렬
    arrange(!!sym(gene_col), desc(is.na(.data$adj_p_value)), .data$adj_p_value, desc(is.na(.data$p_value_raw)), .data$p_value_raw)
  
  
  return(final_adjusted_results)
}


#' Identify Cluster Markers using Pseudo-bulk DEG
#'
#' This function identifies marker genes for specified clusters by performing
#' pseudo-bulk differential expression analysis. It is similar to Seurat's
#' FindMarkers/FindAllMarkers but uses a pseudo-bulk approach based on edgeR
#' or t-tests.
#'
#' @param sobj A Seurat object.
#' @param sample_col Character string. The column name in `sobj@meta.data` that
#'   identifies individual biological samples/replicates. This is crucial for
#'   forming pseudo-bulks.
#' @param group.by Character string. The metadata column to use for defining clusters
#'   (e.g., "seurat_clusters"). Default is "seurat_clusters".
#' @param assay Character string. Assay to use. If NULL, uses DefaultAssay(sobj).
#' @param layer Character string. Layer to use from the assay (Seurat v5).
#'   Default is "data". If using Seurat v3/4, consider using the `slot` argument.
#' @param slot Character string. Slot to use from the assay (Seurat v3/v4).
#'   If `layer` is also provided, `layer` takes precedence. If NULL and `layer` is also NULL,
#'   defaults to "data" for `pct.1`/`pct.2` and "counts" for DEG methods.
#' @param features Character vector. Genes to test. If NULL, tests all genes.
#' @param ident.1 Character or NULL. Identity of the first group. If NULL,
#'   finds markers for all clusters against the rest (like FindAllMarkers).
#' @param ident.2 Character or NULL. Identity of the second group. If NULL and
#'   `ident.1` is specified, `ident.1` is compared against all other clusters.
#' @param logfc.threshold Numeric. Limit testing to genes which show, on average,
#'   at least X-fold difference (log-scale). Default is 0.1.
#' @param test.use Character. Statistical test to use. Options are "wilcox" (uses
#'   edgeR via `run_pseudobulk_deg`) or "t-test" (performs t-tests on pseudo-bulk averages).
#'   Default is "wilcox".
#' @param min.pct Numeric. Only test genes that are detected in a minimum fraction
#'   of cells in either of the two populations. Default is 0.1.
#' @param min.diff.pct Numeric. Only test genes that show a minimum difference in
#'   the fraction of detection between the two groups. Default is -Inf (no filter).
#' @param only.pos Logical. Only return positive markers (positive avg_log2FC). Default is FALSE.
#' @param norm_method_edgeR Character. Normalization method for edgeR if `test.use = "wilcox"`.
#'   Default is "TMM". Passed to `prepare_pseudobulk_edgeR`.
#' @param min_count_edgeR Numeric. Minimum count for `filterByExpr` if `test.use = "wilcox"`.
#'   Default is 10. Passed to `prepare_pseudobulk_edgeR`.
#' @param p_adjust_method Character. Method for p-value adjustment (e.g., "BH", "bonferroni"). Default is "BH".
#' @param verbose Logical. Print progress messages. Default is TRUE.
#' @param ... Additional arguments passed to downstream functions if any (currently not extensively used).
#'
#' @return A data frame (or a list of data frames if `ident.1` is NULL and multiple
#'   clusters are tested) containing marker genes, log2 fold changes, p-values,
#'   adjusted p-values, and percentage of cells expressing the gene in each group.
#'   The columns are named to be similar to Seurat's `FindMarkers` output.
#'
#' @importFrom Seurat DefaultAssay GetAssayData Idents AggregateExpression
#' @importFrom dplyr %>% filter select mutate rename group_by summarise arrange bind_rows left_join distinct
#' @importFrom tibble tibble rownames_to_column column_to_rownames
#' @importFrom stats p.adjust t.test shapiro.test as.formula na.omit
#' @importFrom car leveneTest
#' @importFrom rlang sym
#'
#' @export
cluster_pseudobulk_deg <- function(sobj,
                                   sample_col,
                                   group.by = "seurat_clusters",
                                   assay = NULL,
                                   layer = "data",
                                   slot = NULL,
                                   features = NULL,
                                   ident.1 = NULL,
                                   ident.2 = NULL,
                                   logfc.threshold = 0.1,
                                   test.use = "wilcox",
                                   min.pct = 0.1,
                                   min.diff.pct = -Inf,
                                   only.pos = FALSE,
                                   norm_method_edgeR = "TMM",
                                   min_count_edgeR = 10,
                                   p_adjust_method = "BH",
                                   verbose = TRUE,
                                   ...) {
  
  # --- 0. Preliminaries & Input Checks ---
  if (verbose) message("Starting cluster_pseudobulk_deg analysis...")
  
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat package is required.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr package is required.")
  if (test.use == "wilcox" && !requireNamespace("edgeR", quietly = TRUE)) {
    stop("edgeR package is required for test.use = 'wilcox'.")
  }
  if (test.use == "t-test" && !requireNamespace("car", quietly = TRUE)) {
    stop("car package is required for Levene's test when test.use = 't-test'.")
  }
  
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(sobj)
    if (verbose) message("Using default assay: ", assay)
  }
  if (!assay %in% Seurat::Assays(sobj)) {
    stop("Specified assay '", assay, "' not found in Seurat object.")
  }
  
  # Handle layer vs slot for Seurat version compatibility
  # Prioritize 'layer' if provided (Seurat v5), else use 'slot'
  data_source_param <- if (!is.null(layer)) layer else slot
  if (is.null(data_source_param)) data_source_param <- "data" # Default if both are NULL
  
  # For DEG methods (edgeR/t-test), raw counts are often preferred for pseudo-bulking
  counts_source_param <- "counts" # Typically, pseudo-bulk methods start from counts
  
  if (!sample_col %in% colnames(sobj@meta.data)) {
    stop("`sample_col` '", sample_col, "' not found in Seurat object metadata.")
  }
  if (!group.by %in% colnames(sobj@meta.data)) {
    stop("`group.by` column '", group.by, "' not found in Seurat object metadata.")
  }
  
  # Ensure the group.by column is a factor
  sobj@meta.data[[group.by]] <- factor(sobj@meta.data[[group.by]])
  all_group_by_levels <- levels(sobj@meta.data[[group.by]])
  
  all_results_list <- list()
  clusters_to_iterate <- NULL
  comparison_mode <- ""
  
  # --- 1. Determine comparison mode and clusters to iterate ---
  if (is.null(ident.1)) {
    comparison_mode <- "FindAllMarkers"
    clusters_to_iterate <- all_group_by_levels
    if (verbose) message("Mode: FindAllMarkers. Iterating through all clusters in '", group.by, "'.")
  } else {
    if (!ident.1 %in% all_group_by_levels) {
      stop("`ident.1` ('", ident.1, "') not found in `group.by` column ('", group.by, "'). Available levels: ", paste(all_group_by_levels, collapse=", "))
    }
    if (is.null(ident.2)) {
      comparison_mode <- "OneVsRest"
      clusters_to_iterate <- ident.1
      if (verbose) message("Mode: OneVsRest. Comparing cluster '", ident.1, "' vs all others.")
    } else {
      if (!ident.2 %in% all_group_by_levels) {
        stop("`ident.2` ('", ident.2, "') not found in `group.by` column ('", group.by, "'). Available levels: ", paste(all_group_by_levels, collapse=", "))
      }
      if(ident.1 == ident.2) stop("`ident.1` and `ident.2` cannot be the same.")
      comparison_mode <- "GroupVsGroup"
      clusters_to_iterate <- ident.1 # Process once for this specific comparison
      if (verbose) message("Mode: GroupVsGroup. Comparing cluster '", ident.1, "' vs '", ident.2, "'.")
    }
  }
  
  # --- 2. Loop through target clusters for DEG ---
  for (current_target_cluster_name in clusters_to_iterate) {
    if (verbose) message("\nProcessing comparison for: ", current_target_cluster_name)
    
    # Create a temporary Seurat object or a copy of metadata for modification
    # Using a copy of metadata is safer if sobj is large
    temp_meta <- sobj@meta.data
    temp_group_col <- paste0("temp_comparison_status_", make.names(Sys.time())) # Unique temp column name
    
    group1_label <- "GROUP1_TARGET"
    group2_label <- "GROUP2_REFERENCE"
    
    label_for_ident.1_in_output <- current_target_cluster_name # This will be used for the 'cluster' column
    label_for_ident.2_in_output <- "Rest" # Default for OneVsRest/FindAllMarkers
    
    if (comparison_mode == "FindAllMarkers" || comparison_mode == "OneVsRest") {
      temp_meta[[temp_group_col]] <- ifelse(
        temp_meta[[group.by]] == current_target_cluster_name,
        group1_label,
        group2_label
      )
    } else { # GroupVsGroup
      # Filter cells to only those in ident.1 or ident.2 for this specific comparison
      cells_for_comparison <- rownames(temp_meta[temp_meta[[group.by]] %in% c(ident.1, ident.2), ])
      if (length(cells_for_comparison) == 0) {
        warning("No cells found for ident.1 ('", ident.1, "') or ident.2 ('", ident.2, "'). Skipping this comparison.")
        next
      }
      temp_meta_subset <- temp_meta[cells_for_comparison, ]
      temp_meta_subset[[temp_group_col]] <- ifelse(
        temp_meta_subset[[group.by]] == ident.1,
        group1_label,
        group2_label
      )
      # For GroupVsGroup, ensure we operate on the subset of cells
      # Create a temporary sobj for this specific comparison to pass to prepare_pseudobulk_edgeR
      sobj_for_deg <- subset(sobj, cells = cells_for_comparison)
      sobj_for_deg@meta.data <- temp_meta_subset # Assign the modified metadata for this subset
      
      label_for_ident.1_in_output <- ident.1
      label_for_ident.2_in_output <- ident.2
    }
    
    # If not GroupVsGroup, use the full sobj with the added temp_group_col
    if (comparison_mode != "GroupVsGroup") {
      sobj@meta.data[[temp_group_col]] <- temp_meta[[temp_group_col]]
      sobj_for_deg <- sobj
    }
    
    # Ensure temp_group_col is a factor with REFERENCE as the first level for edgeR contrast
    # If GROUP1 is target, GROUP2 should be reference.
    sobj_for_deg@meta.data[[temp_group_col]] <- factor(sobj_for_deg@meta.data[[temp_group_col]], levels = c(group2_label, group1_label))
    
    deg_results_final_for_loop <- NULL
    
    # --- 3. Perform DEG ---
    if (test.use == "wilcox") { # Assuming "wilcox" implies using the edgeR backend from run_pseudobulk_deg
      if (verbose) message("  Using edgeR via run_pseudobulk_deg...")
      
      # Check for sufficient samples in each group within temp_group_col *after pseudo-bulking*
      # This is implicitly handled by prepare_pseudobulk_edgeR's internal checks or edgeR itself.
      # We can add an explicit check if needed on the output of prepare_pseudobulk_edgeR$meta.
      
      # Prepare data for edgeR
      # `cluster_col` in `prepare_pseudobulk_edgeR` should be the original `group.by`
      # to form pseudo-bulks per original cluster before grouping them into GROUP1/GROUP2.
      # This ensures that the variability within original clusters is part of the pseudo-bulk.
      prep_args <- list(
        seurat_obj = sobj_for_deg,
        assay = assay,
        slot = counts_source_param,
        sample_col = sample_col,       # True biological sample ID
        cluster_col = group.by,        # Original clusters/groups specified by user
        group_col = temp_group_col,    # This is our GROUP1 vs GROUP2 comparison column
        min_count = min_count_edgeR,
        norm_method = norm_method_edgeR,
        design_formula = as.formula(paste("~", temp_group_col)),
        verbose = FALSE
      )
      
      prepared_data_for_deg <- tryCatch(
        do.call(prepare_pseudobulk_edgeR, prep_args),
        error = function(e) {
          message("    Error in prepare_pseudobulk_edgeR for ", current_target_cluster_name, ": ", e$message)
          return(NULL)
        }
      )
      
      if (is.null(prepared_data_for_deg)) next # Skip if preparation failed
      
      # --- Enhanced Check for Sample Distribution Before DEG ---
      grouping_col_in_meta_pb <- temp_group_col # This is "GROUP1_TARGET" / "GROUP2_REFERENCE" column
      patient_col_in_meta_pb <- "patient"     # As defined in prepare_pseudobulk_edgeR's separate step
      
      if (!grouping_col_in_meta_pb %in% colnames(prepared_data_for_deg$meta)) {
        message("    CRITICAL WARNING for '", current_target_cluster_name,
                "': Expected comparison grouping column '", grouping_col_in_meta_pb,
                "' NOT FOUND in pseudo-bulk metadata (prepared_data_for_deg$meta). Available columns: ",
                paste(colnames(prepared_data_for_deg$meta), collapse=", "),
                ". This indicates a problem in `prepare_pseudobulk_edgeR` or how `temp_group_col` was defined/passed.",
                " Skipping DEG for this cluster.")
        if (comparison_mode != "GroupVsGroup") sobj@meta.data[[temp_group_col]] <- NULL # Cleanup
        next
      }
      if (!patient_col_in_meta_pb %in% colnames(prepared_data_for_deg$meta)) {
        message("    CRITICAL WARNING for '", current_target_cluster_name,
                "': Expected patient identifier column '", patient_col_in_meta_pb,
                "' NOT FOUND in pseudo-bulk metadata (prepared_data_for_deg$meta). Available columns: ",
                paste(colnames(prepared_data_for_deg$meta), collapse=", "),
                ". Cannot reliably check sample counts per group. Skipping DEG for this cluster to be safe.")
        if (comparison_mode != "GroupVsGroup") sobj@meta.data[[temp_group_col]] <- NULL # Cleanup
        next
      }
      
      sample_counts_per_comparison_group <- prepared_data_for_deg$meta %>%
        dplyr::group_by(!!rlang::sym(grouping_col_in_meta_pb)) %>%
        dplyr::summarise(
          n_pseudo_bulks = dplyr::n(),
          n_unique_samples = dplyr::n_distinct(!!rlang::sym(patient_col_in_meta_pb)),
          .groups = "drop"
        )
      
      min_reps_needed_per_group <- 2 # Standard minimum for edgeR comparisons
      
      if (verbose) {
        message("    Pseudo-bulk sample distribution for comparison involving '", label_for_ident.1_in_output,
                "' (vs '", label_for_ident.2_in_output, "'):")
        print(sample_counts_per_comparison_group)
      }
      
      # Check 1: Are there at least two groups to compare? (e.g., GROUP1_TARGET and GROUP2_REFERENCE)
      if (nrow(sample_counts_per_comparison_group) < 2) {
        message("    SKIPPING DEG for '", current_target_cluster_name,
                "': Less than two comparison groups found in pseudo-bulk data (expected ",
                group1_label, " and ", group2_label, "; found ",
                nrow(sample_counts_per_comparison_group), " groups: ",
                paste(sample_counts_per_comparison_group[[grouping_col_in_meta_pb]], collapse=", "),
                "). This might occur if the target cluster '", current_target_cluster_name,
                "' or the 'Rest' group has no cells/samples after aggregation and filtering.")
        if (comparison_mode != "GroupVsGroup") sobj@meta.data[[temp_group_col]] <- NULL # Cleanup
        next
      }
      
      # Check 2: Do all present groups have enough unique biological samples?
      if (any(sample_counts_per_comparison_group$n_unique_samples < min_reps_needed_per_group, na.rm = TRUE)) {
        problematic_group_details <- sample_counts_per_comparison_group %>%
          dplyr::filter(is.na(n_unique_samples) | n_unique_samples < min_reps_needed_per_group)
        message("    SKIPPING DEG for '", current_target_cluster_name,
                "': One or more comparison groups have fewer than ", min_reps_needed_per_group,
                " unique biological samples (derived from '", sample_col,
                "'). This is a common cause for 'design matrix not full rank' errors in edgeR.")
        message("    Problematic group counts for this comparison:")
        print(problematic_group_details)
        if (comparison_mode != "GroupVsGroup") sobj@meta.data[[temp_group_col]] <- NULL # Cleanup
        next
      }
      # --- End of Enhanced Check ---
      
      deg_results_raw <- tryCatch({
        run_pseudobulk_deg(
          analysis_level = "overall",
          contrast = current_contrast, 
          prepared_data = prepared_data_for_deg,
          verbose = FALSE
        )
      }, error = function(e) {
        message("    Error in run_pseudobulk_deg for ", current_target_cluster_name, ": ", e$message)
        return(NULL)
      })
      
      if (is.null(deg_results_raw) || nrow(deg_results_raw) == 0) {
        if (verbose) message("    No DEG results from run_pseudobulk_deg for ", current_target_cluster_name)
        # Clean up temp column before next iteration if it was added to original sobj
        if (comparison_mode != "GroupVsGroup") sobj@meta.data[[temp_group_col]] <- NULL
        next
      }
      
      deg_results_final_for_loop <- deg_results_raw %>%
        dplyr::rename(avg_log2FC = logFC, p_val = PValue, p_val_adj = FDR) %>%
        dplyr::select(gene, avg_log2FC, p_val, p_val_adj, dplyr::everything())
      
    } else if (test.use == "t-test") {
      if (verbose) message("  Using t-tests on pseudo-bulk averages...")
      
      # 1. Create pseudo-bulks: sum counts per sample per (GROUP1/GROUP2)
      # `AggregateExpression` sums counts by default if slot="counts"
      pb_for_ttest <- Seurat::AggregateExpression(
        sobj_for_deg, # This has the temp_group_col and is subsetted for GroupVsGroup
        assays = assay,
        slot = counts_source_param, 
        group.by = c(sample_col, temp_group_col),
        return.seurat = FALSE
      )[[assay]]
      
      if (is.null(pb_for_ttest) || ncol(pb_for_ttest) == 0) {
        message("    No pseudo-bulk data generated for t-test for ", current_target_cluster_name)
        if (comparison_mode != "GroupVsGroup") sobj@meta.data[[temp_group_col]] <- NULL
        next
      }
      
      # Create metadata for these pseudo-bulks
      pb_meta_for_ttest <- data.frame(pb_sample_id = colnames(pb_for_ttest)) %>%
        tidyr::separate(pb_sample_id, into = c("sample_id_temp", temp_group_col), 
                        sep = "_", extra = "merge", remove = FALSE) %>%
        dplyr::select(pb_sample_id, !!rlang::sym(temp_group_col))
      
      # Genes to test
      genes_to_test_ttest <- if (!is.null(features)) intersect(features, rownames(pb_for_ttest)) else rownames(pb_for_ttest)
      if (length(genes_to_test_ttest) == 0) {
        message("    No valid features to test with t-test for ", current_target_cluster_name)
        if (comparison_mode != "GroupVsGroup") sobj@meta.data[[temp_group_col]] <- NULL
        next
      }
      
      pb_for_ttest_filtered <- pb_for_ttest[genes_to_test_ttest, , drop = FALSE]
      
      ttest_results_list <- lapply(rownames(pb_for_ttest_filtered), function(g) {
        gene_expr_df <- data.frame(
          expression = as.numeric(pb_for_ttest_filtered[g, ]),
          group = pb_meta_for_ttest[[temp_group_col]]
        )
        
        group1_vals <- gene_expr_df$expression[gene_expr_df$group == group1_label]
        group2_vals <- gene_expr_df$expression[gene_expr_df$group == group2_label]
        
        # Ensure enough data points per group for t-test
        if (length(group1_vals) < 2 || length(group2_vals) < 2) {
          return(NULL)
        }
        
        # Normalize (e.g., CPM or log2(CPM+1)) - edgeR does library size normalization internally
        # For t-test on pseudo-bulk, it's common to use log-transformed normalized counts (e.g., log2(CPM+1))
        # Here, using raw pseudo-bulk sums; consider normalization if variance is an issue.
        # For simplicity, let's do log2(sum + 1) for FC calculation.
        
        logFC <- mean(log2(group1_vals + 1), na.rm = TRUE) - mean(log2(group2_vals + 1), na.rm = TRUE)
        
        # Levene's Test (homogeneity of variances)
        levene_res <- car::leveneTest(expression ~ group, data = gene_expr_df)
        levene_pval <- levene_res$"Pr(>F)"[1]
        
        # Shapiro-Wilk Test for normality (on each group)
        shapiro_g1_pval <- if(length(group1_vals) >= 3) stats::shapiro.test(group1_vals)$p.value else NA
        shapiro_g2_pval <- if(length(group2_vals) >= 3) stats::shapiro.test(group2_vals)$p.value else NA
        
        if (verbose) {
          cat("  Gene:", g, "- Levene p-val:", format(levene_pval, digits=3),
              "| Shapiro G1 p-val:", format(shapiro_g1_pval, digits=3),
              "| Shapiro G2 p-val:", format(shapiro_g2_pval, digits=3), "\n")
        }
        
        t_res <- stats::t.test(group1_vals, group2_vals, var.equal = (levene_pval > 0.05))
        
        data.frame(
          gene = g,
          avg_log2FC = logFC, # This is log2(mean(grp1)/mean(grp2)) if using geometric means, or diff of log-means
          p_val = t_res$p.value,
          t_statistic = t_res$statistic,
          levene_pval = levene_pval,
          shapiro_pval_g1 = shapiro_g1_pval,
          shapiro_pval_g2 = shapiro_g2_pval
        )
      })
      
      deg_results_final_for_loop <- dplyr::bind_rows(ttest_results_list)
      if (nrow(deg_results_final_for_loop) > 0) {
        deg_results_final_for_loop$p_val_adj <- stats::p.adjust(deg_results_final_for_loop$p_val, method = p_adjust_method)
      }
    } else {
      stop("Unsupported test.use method. Choose 'wilcox' (for edgeR) or 't-test'.")
    }
    
    if (is.null(deg_results_final_for_loop) || nrow(deg_results_final_for_loop) == 0) {
      if (verbose) message("    No DEG results for ", current_target_cluster_name)
      if (comparison_mode != "GroupVsGroup") sobj@meta.data[[temp_group_col]] <- NULL
      next
    }
    
    # --- 4. Calculate pct.1 and pct.2 ---
    # These should reflect the original cell populations, not pseudo-bulks
    genes_for_pct_calc <- deg_results_final_for_loop$gene
    
    # Cells belonging to the target group (ident.1 or current_target_cluster_name)
    cells_group1_orig <- rownames(sobj@meta.data[sobj@meta.data[[group.by]] == current_target_cluster_name, ])
    
    # Cells belonging to the reference group (ident.2 or "Rest")
    if (comparison_mode == "FindAllMarkers" || comparison_mode == "OneVsRest") {
      other_levels_for_group2 <- setdiff(all_group_by_levels, current_target_cluster_name)
      cells_group2_orig <- rownames(sobj@meta.data[sobj@meta.data[[group.by]] %in% other_levels_for_group2, ])
    } else { # GroupVsGroup
      cells_group2_orig <- rownames(sobj@meta.data[sobj@meta.data[[group.by]] == ident.2, ])
    }
    
    # Use the specified layer/slot for pct calculation (e.g., "data" for log-normalized, "counts" for raw)
    assay_data_for_pct <- Seurat::GetAssayData(sobj, assay = assay, layer = data_source_param)
    
    # Ensure genes_for_pct_calc exist in the assay
    genes_present_in_assay <- intersect(genes_for_pct_calc, rownames(assay_data_for_pct))
    if(length(genes_present_in_assay) == 0) {
      message("    No genes for pct calculation found in the assay for ", current_target_cluster_name)
      deg_results_final_for_loop$pct.1 <- NA
      deg_results_final_for_loop$pct.2 <- NA
    } else {
      data_group1_pct <- assay_data_for_pct[genes_present_in_assay, intersect(cells_group1_orig, colnames(assay_data_for_pct)), drop = FALSE]
      data_group2_pct <- assay_data_for_pct[genes_present_in_assay, intersect(cells_group2_orig, colnames(assay_data_for_pct)), drop = FALSE]
      
      pct.1_values <- if(ncol(data_group1_pct) > 0) rowSums(data_group1_pct > 0) / ncol(data_group1_pct) else rep(0, length(genes_present_in_assay))
      pct.2_values <- if(ncol(data_group2_pct) > 0) rowSums(data_group2_pct > 0) / ncol(data_group2_pct) else rep(0, length(genes_present_in_assay))
      
      names(pct.1_values) <- genes_present_in_assay
      names(pct.2_values) <- genes_present_in_assay
      
      pct_df_for_merge <- data.frame(
        gene = genes_present_in_assay,
        pct.1 = pct.1_values[genes_present_in_assay],
        pct.2 = pct.2_values[genes_present_in_assay]
      )
      
      deg_results_final_for_loop <- deg_results_final_for_loop %>%
        dplyr::left_join(pct_df_for_merge, by = "gene")
    }
    
    deg_results_final_for_loop$cluster <- label_for_ident.1_in_output
    if (comparison_mode == "GroupVsGroup") {
      deg_results_final_for_loop$ident.2 <- label_for_ident.2_in_output
    }
    
    all_results_list[[current_target_cluster_name]] <- deg_results_final_for_loop
    
    # Clean up temporary metadata column from the original sobj if it was modified
    if (comparison_mode != "GroupVsGroup") {
      sobj@meta.data[[temp_group_col]] <- NULL
    }
    if (verbose) message("  Finished comparison for: ", current_target_cluster_name)
  } # End loop over clusters_to_iterate
  
  # --- 5. Combine results and apply final filters ---
  if (length(all_results_list) == 0) {
    if (verbose) message("No results generated.")
    return(data.frame())
  }
  
  final_df <- dplyr::bind_rows(all_results_list)
  
  # Ensure p_val_adj is present (might be missing if t-test path had issues)
  if (!"p_val_adj" %in% names(final_df) && "p_val" %in% names(final_df)) {
    final_df <- final_df %>%
      dplyr::group_by(cluster) %>% # Adjust per comparison group
      dplyr::mutate(p_val_adj = stats::p.adjust(p_val, method = p_adjust_method)) %>%
      dplyr::ungroup()
  }
  
  # Apply feature filter if provided
  if (!is.null(features)) {
    final_df <- final_df %>% dplyr::filter(gene %in% features)
  }
  
  # Apply post-hoc filters
  if (!is.null(final_df$avg_log2FC)) { # Check if column exists
    final_df <- final_df %>% dplyr::filter(abs(avg_log2FC) >= logfc.threshold)
  }
  if (!is.null(final_df$pct.1) && !is.null(final_df$pct.2)) { # Check if columns exist
    final_df <- final_df %>% dplyr::filter(pct.1 >= min.pct | pct.2 >= min.pct)
    final_df <- final_df %>% dplyr::filter(abs(pct.1 - pct.2) >= min.diff.pct)
  }
  
  if (only.pos && !is.null(final_df$avg_log2FC)) {
    final_df <- final_df %>% dplyr::filter(avg_log2FC > 0)
  }
  
  # Select and order columns
  cols_ordered <- c("gene", "cluster", "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj")
  # Add other common/useful columns if they exist
  for (col_add in c("logCPM", "F", "t_statistic", "levene_pval", "shapiro_pval_g1", "shapiro_pval_g2", "ident.2")) {
    if (col_add %in% names(final_df)) {
      cols_ordered <- c(cols_ordered, col_add)
    }
  }
  # Add any remaining columns not yet selected
  cols_ordered <- unique(c(cols_ordered, names(final_df)))
  
  final_df <- final_df %>% 
    dplyr::select(dplyr::all_of(intersect(cols_ordered, names(final_df)))) %>%
    dplyr::arrange(cluster, p_val_adj, desc(abs(avg_log2FC)))
  
  if (verbose) message("Analysis complete. Returning ", nrow(final_df), " marker genes.")
  return(final_df)
}




#' Perform Pseudobulk Differential Gene Expression Analysis
#'
#' This function takes a Seurat object, aggregates counts to create pseudobulks
#' based on specified sample and bulk category keys, and then performs DEG analysis
#' between groups defined by a comparison key. It supports edgeR, DESeq2,
#' Wilcoxon rank-sum test, Student's t-test, and ROC analysis.
#'
#' @param sobj A Seurat object.
#' @param sample_col Character string. The name of the metadata column identifying individual samples or patients.
#' @param compare_col Character string. The name of the metadata column used for grouping samples for DEG comparison (e.g., prognosis factor, treatment group).
#' @param bulk_col Character string. The name of the metadata column that defines the categories within which pseudobulking and DEG analysis will be performed (e.g., cell type). DEGs are found for each unique value in this column separately.
#' @param idents.1 Character vector or NULL. Specifies the identity (or identities) in `compare_col` for the first group.
#'                 If NULL (default), performs a "FindAllMarkers" style analysis, comparing each identity in `compare_col` against all others within each `bulk_col` category.
#' @param idents.2 Character vector or NULL. Specifies the identity (or identities) in `compare_col` for the second group.
#'                 If NULL (default) and `idents.1` is specified, `idents.1` is compared against all other identities.
#'                 If `idents.1` is also NULL, this parameter is ignored.
#' @param test.use Character string. The statistical test to use. Options are:
#'                 "edgeR", "DESeq2", "wilcox" (default), "t.test", "roc".
#' @param assay Character string or NULL. The Seurat assay to use for fetching counts. If NULL, defaults to `DefaultAssay(sobj)`.
#' @param slot Character string. The slot from which to pull data for aggregation (e.g., "counts"). Defaults to "counts".
#' @param min.cells.gene Numeric. For 'wilcox', 't.test', 'roc': minimum number of pseudobulk samples where a gene must be detected to be included in the analysis. This is applied per comparison (i.e., gene must be in `min.cells.gene` samples in group1 OR group2). Seurat's `FindMarkers` uses this for single cells.
#' @param min.pct Numeric. For 'wilcox', 't.test', 'roc': only test genes that are detected in at least this fraction of pseudobulk samples in either of the two groups.
#' @param logfc.threshold Numeric. For 'wilcox', 't.test', 'roc': limit testing to genes which show at least this an absolute (log2) fold change between the two groups of pseudobulk samples. For 'roc', this is applied to the AUC value.
#' @param min.samples.group Numeric. Minimum number of pseudobulk samples required in each group for a comparison to proceed. Defaults to 2.
#' @param edgeR.filterByExpr.min.count Numeric. Passed to `min.count` in `edgeR::filterByExpr`. Default is 5.
#' @param edgeR.filterByExpr.min.total.count Numeric. Passed to `min.total.count` in `edgeR::filterByExpr`. Default is 10.
#' @param DESeq2.useT Logical. For DESeq2, passed to `useT` in `DESeq` function if initial fit fails. Defaults to FALSE.
#' @param DESeq2.fitType Character. For DESeq2, passed to `fitType` in `DESeq` function. Default is "parametric". If error, "local" may be tried.
#' @param BPPARAM BiocParallel::BiocParallelParam instance, for parallel execution in DESeq2. Default is `BiocParallel::SerialParam()`.
#' @param ... Additional arguments passed to the underlying statistical test functions (e.g., `exact` for `wilcox.test`).
#'
#' @return A data frame with DEG results, similar to Seurat's `FindMarkers` output. Columns include:
#'         `gene`, `cluster` (the `idents.1` group or current group in `FindAllMarkers` mode),
#'         `bulk_category` (from `bulk_col`), `avg_log2FC`, `p_val`, `p_val_adj`, `pct.1`, `pct.2`.
#'         For `test.use = "roc"`, `avg_log2FC` column contains the AUC.
#'
#' @import Seurat
#' @import dplyr
#' @import tidyr
#' @importFrom Matrix.utils aggregate.Matrix
#' @importFrom stats p.adjust t.test wilcox.test model.matrix na.omit sd
#' @importFrom BiocParallel SerialParam bplapply
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming 'pbmc_small' is a Seurat object with relevant metadata
#' # pbmc_small$sample <- rep(paste0("sample", 1:2), length.out = ncol(pbmc_small))
#' # pbmc_small$cell_type <- pbmc_small$RNA_snn_res.0.8
#' # pbmc_small$condition <- rep(c("A", "B"), length.out = ncol(pbmc_small))
#' # Ensure 'condition' is consistent per 'sample'
#' # sample_conditions <- pbmc_small@meta.data %>% distinct(sample, condition)
#' # pbmc_small <- pbmc_small[, !(pbmc_small$sample %in% names(which(table(sample_conditions$sample) > 1)))]
#' # temp_meta <- pbmc_small@meta.data
#' # temp_meta <- temp_meta %>% group_by(sample) %>% mutate(condition = first(condition)) %>% ungroup()
#' # rownames(temp_meta) <- colnames(pbmc_small)
#' # pbmc_small@meta.data <- temp_meta
#'
#' # deg_results_wilcox <- pseudobulk_deg(
#' #   sobj = pbmc_small,
#' #   sample_col = "sample",
#' #   compare_col = "condition", # e.g. Control vs Treatment
#' #   bulk_col = "cell_type",    # e.g. B cells, T cells
#' #   test.use = "wilcox"
#' # )
#'
#' # To use edgeR or DESeq2, ensure they are installed
#' # if (requireNamespace("edgeR", quietly = TRUE)) {
#' #   deg_results_edgeR <- pseudobulk_deg(
#' #     sobj = pbmc_small,
#' #     sample_col = "sample",
#' #     compare_col = "condition",
#' #     bulk_col = "cell_type",
#' #     test.use = "edgeR"
#' #  )
#' # }
#' # if (requireNamespace("DESeq2", quietly = TRUE) &&
#' #     requireNamespace("BiocParallel", quietly = TRUE)) {
#' #   deg_results_DESeq2 <- pseudobulk_deg(
#' #     sobj = pbmc_small,
#' #     sample_col = "sample",
#' #     compare_col = "condition",
#' #     bulk_col = "cell_type",
#' #     test.use = "DESeq2"
#' #   )
#' # }
#' }
pseudobulk_deg <- function(sobj,
                           sample_col,
                           compare_col,
                           bulk_col,
                           idents.1 = NULL,
                           idents.2 = NULL,
                           test.use = "wilcox",
                           assay = NULL,
                           slot = "counts",
                           min.cells.gene = 3,
                           min.pct = 0.1,
                           logfc.threshold = 0.25,
                           min.samples.group = 2,
                           edgeR.filterByExpr.min.count = 5,
                           edgeR.filterByExpr.min.total.count = 10,
                           DESeq2.useT = FALSE,
                           DESeq2.fitType = "parametric",
                           BPPARAM = BiocParallel::SerialParam(),
                           ...) {
  
  # --- 0. Load necessary libraries and Argument Checks ---
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Package 'Seurat' needed for this function to work. Please install it.", call. = FALSE)
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' needed for this function to work. Please install it.", call. = FALSE)
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("Package 'tidyr' needed for this function to work. Please install it.", call. = FALSE)
  
  
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(sobj)
  }
  meta_cols <- c(sample_col, compare_col, bulk_col)
  if (!all(meta_cols %in% colnames(sobj@meta.data))) {
    missing_cols <- meta_cols[!meta_cols %in% colnames(sobj@meta.data)]
    stop("Metadata column(s) not found: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  valid_tests <- c("wilcox", "t.test", "roc", "edgeR", "DESeq2")
  if (!test.use %in% valid_tests) {
    stop("Invalid test.use. Choose from: ", paste(valid_tests, collapse = ", "), call. = FALSE)
  }
  if (test.use == "edgeR" && !requireNamespace("edgeR", quietly = TRUE)) {
    stop("Package 'edgeR' needed for test.use='edgeR'. Please install it.", call. = FALSE)
  }
  if (test.use == "DESeq2") {
    if(!requireNamespace("DESeq2", quietly = TRUE)) stop("Package 'DESeq2' needed for test.use='DESeq2'. Please install it.", call. = FALSE)
    if(!requireNamespace("BiocParallel", quietly = TRUE)) stop("Package 'BiocParallel' needed for DESeq2 parallel execution. Please install it.", call. = FALSE)
  }
  if (test.use == "roc" && !requireNamespace("pROC", quietly = TRUE)) {
    stop("Package 'pROC' needed for test.use='roc'. Please install it.", call. = FALSE)
  }
  
  all_results_list <- list()
  unique_bulk_categories <- unique(as.character(sobj@meta.data[[bulk_col]]))
  unique_bulk_categories <- sort(na.omit(unique_bulk_categories))
  
  
  # --- 1. Loop over each bulk category (e.g., cell type) ---
  for (current_bulk_id in unique_bulk_categories) {
    message(paste0("\nProcessing bulk category: '", current_bulk_id, "'"))
    
    cells_in_bulk <- rownames(sobj@meta.data[which(sobj@meta.data[[bulk_col]] == current_bulk_id), ])
    if (length(cells_in_bulk) == 0) {
      message(paste0("  No cells found for bulk category: '", current_bulk_id, "'. Skipping."))
      next
    }
    sobj_bulk_subset <- subset(sobj, cells = cells_in_bulk)
    
    if (length(unique(sobj_bulk_subset@meta.data[[sample_col]])) < min.samples.group) {
      message(paste0("  Not enough unique samples (", length(unique(sobj_bulk_subset@meta.data[[sample_col]])),
                     ") in bulk category: '", current_bulk_id, "' for comparison (min required: ", min.samples.group, "). Skipping."))
      next
    }
    
    # --- 2. Create pseudobulk counts for the current bulk category ---
    current_pb_counts <- tryCatch({
      Seurat::AggregateExpression(sobj_bulk_subset,
                                  group.by = sample_col,
                                  assays = assay,
                                  slot = slot,
                                  return.seurat = FALSE)[[assay]]
    }, error = function(e) {
      message(paste0("  Error during AggregateExpression for bulk category '", current_bulk_id, "': ", e$message))
      return(NULL)
    })
    
    if (is.null(current_pb_counts) || ncol(current_pb_counts) == 0) {
      message(paste0("  Aggregation resulted in no pseudobulk samples for bulk category: '", current_bulk_id, "'. Skipping."))
      next
    }
    # Ensure counts are integers for DESeq2/edgeR
    if (test.use %in% c("edgeR", "DESeq2")) {
      current_pb_counts <- round(current_pb_counts)
    }
    
    
    # Create metadata for pseudobulks
    # Important: Ensure compare_col value is unique per sample_col. If not, this takes the first one.
    # This assumes that compare_col is a sample-level attribute.
    current_pb_meta <- sobj_bulk_subset@meta.data %>%
      dplyr::distinct(!!dplyr::sym(sample_col), .keep_all = TRUE) %>%
      dplyr::select(!!dplyr::sym(sample_col), !!dplyr::sym(compare_col)) %>%
      dplyr::filter(!!dplyr::sym(sample_col) %in% colnames(current_pb_counts)) %>%
      dplyr::mutate(!!dplyr::sym(compare_col) := as.character(!!dplyr::sym(compare_col))) # Ensure character for consistency
    
    
    shared_samples <- intersect(colnames(current_pb_counts), current_pb_meta[[sample_col]])
    if (length(shared_samples) < min.samples.group) {
      message(paste0("  Not enough shared samples (",length(shared_samples),") after aligning counts and metadata for '", current_bulk_id, "'. Skipping."))
      next
    }
    current_pb_counts <- current_pb_counts[, shared_samples, drop = FALSE]
    current_pb_meta <- current_pb_meta %>% dplyr::filter(!!dplyr::sym(sample_col) %in% shared_samples)
    rownames(current_pb_meta) <- current_pb_meta[[sample_col]]
    
    # --- 3. Define comparison groups ---
    all_compare_values <- na.omit(unique(current_pb_meta[[compare_col]]))
    all_compare_values <- sort(all_compare_values)
    
    comparisons_to_make <- list()
    if (is.null(idents.1)) { # FindAllMarkers style
      if (length(all_compare_values) < 2) {
        message(paste0("  Not enough groups (", length(all_compare_values), ") to compare within '", current_bulk_id, "'. Skipping."))
        next
      }
      for (id1_val_loop in all_compare_values) {
        # Check if id1_val_loop has enough samples and if "others" group has enough samples
        n_group1 = sum(current_pb_meta[[compare_col]] == id1_val_loop, na.rm = TRUE)
        n_group_other = sum(current_pb_meta[[compare_col]] != id1_val_loop, na.rm = TRUE)
        if (n_group1 >= min.samples.group && n_group_other >= min.samples.group) {
          comparisons_to_make[[length(comparisons_to_make) + 1]] <- list(ident1 = id1_val_loop, ident2 = NULL, name = id1_val_loop)
        } else {
          message(paste0("    Skipping FindAllMarkers-style comparison for '", id1_val_loop, "' vs rest in '", current_bulk_id,
                         "' due to insufficient samples in one or both groups (", n_group1, " vs ", n_group_other,"). Min per group: ", min.samples.group))
        }
      }
    } else { # idents.1 is specified
      ident1_val_input <- idents.1
      n_group1 = sum(current_pb_meta[[compare_col]] %in% ident1_val_input, na.rm = TRUE)
      
      if (n_group1 < min.samples.group) {
        message(paste0("  Ident1 group '", paste(ident1_val_input, collapse=", "), "' has insufficient samples (", n_group1,
                       ") in '", current_bulk_id, "'. Min per group: ", min.samples.group, ". Skipping this comparison."))
      } else {
        if (is.null(idents.2)) { # One vs Rest
          n_group_other = sum(!current_pb_meta[[compare_col]] %in% ident1_val_input, na.rm = TRUE)
          if (n_group_other >= min.samples.group) {
            comparisons_to_make[[1]] <- list(ident1 = ident1_val_input, ident2 = NULL, name = paste(ident1_val_input, collapse="_vs_"))
          } else {
            message(paste0("  'Rest' group for ident1 '", paste(ident1_val_input, collapse=", "), "' has insufficient samples (", n_group_other,
                           ") in '", current_bulk_id, "'. Min per group: ", min.samples.group, ". Skipping comparison."))
          }
        } else { # One vs One
          ident2_val_input <- idents.2
          if(any(ident1_val_input %in% ident2_val_input)){
            message(paste0("  Ident1 ('", paste(ident1_val_input, collapse=", "), "') and Ident2 ('", paste(ident2_val_input, collapse=", "),
                           "') overlap. Skipping this comparison in '", current_bulk_id, "'."))
          } else {
            n_group2 = sum(current_pb_meta[[compare_col]] %in% ident2_val_input, na.rm = TRUE)
            if (n_group2 >= min.samples.group) {
              comparisons_to_make[[1]] <- list(ident1 = ident1_val_input, ident2 = ident2_val_input, name = paste(ident1_val_input, collapse="_vs_"))
            } else {
              message(paste0("  Ident2 group '", paste(ident2_val_input, collapse=", "), "' has insufficient samples (", n_group2,
                             ") in '", current_bulk_id, "'. Min per group: ", min.samples.group, ". Skipping comparison."))
            }
          }
        }
      }
    }
    
    if (length(comparisons_to_make) == 0) {
      message(paste0("  No valid comparisons to make for bulk category: '", current_bulk_id, "' based on current parameters."))
      next
    }
    
    # --- 4. Loop for each comparison to be made ---
    for (comp_idx in seq_along(comparisons_to_make)) {
      comp_params <- comparisons_to_make[[comp_idx]]
      id1_value <- comp_params$ident1
      id2_value <- comp_params$ident2
      cluster_name_for_output <- comp_params$name
      
      
      message(paste0("    Running DEG for: cluster '", cluster_name_for_output, "' (", paste(id1_value, collapse=", "), " vs ",
                     if(is.null(id2_value)) "all others" else paste(id2_value, collapse=", "), ")"))
      
      samples1_names <- current_pb_meta %>% dplyr::filter(!!dplyr::sym(compare_col) %in% id1_value) %>% dplyr::pull(!!dplyr::sym(sample_col))
      if (is.null(id2_value)) {
        samples2_names <- current_pb_meta %>% dplyr::filter(!(!!dplyr::sym(compare_col) %in% id1_value)) %>% dplyr::pull(!!dplyr::sym(sample_col))
      } else {
        samples2_names <- current_pb_meta %>% dplyr::filter(!!dplyr::sym(compare_col) %in% id2_value) %>% dplyr::pull(!!dplyr::sym(sample_col))
      }
      
      # Already checked min.samples.group at comparison_to_make stage, but double check here after sample name extraction.
      if (length(samples1_names) < min.samples.group || length(samples2_names) < min.samples.group) {
        message(paste0("      Skipping comparison for '", cluster_name_for_output, "' due to insufficient samples in one group after final selection (",
                       length(samples1_names), " vs ", length(samples2_names), "). Min per group: ", min.samples.group))
        next
      }
      
      comp_specific_pb_counts <- current_pb_counts[, c(samples1_names, samples2_names), drop = FALSE]
      comp_specific_pb_meta <- current_pb_meta[c(samples1_names, samples2_names), , drop = FALSE]
      comp_specific_pb_meta$de_group <- factor(ifelse(comp_specific_pb_meta[[sample_col]] %in% samples1_names, "group1", "group2"), levels = c("group2", "group1")) # group1 is case, group2 is control
      
      res_df <- NULL # Initialize results data frame for this comparison
      
      # --- 5. Perform DEG analysis based on test.use ---
      if (test.use == "edgeR") {
        if (!requireNamespace("edgeR", quietly = TRUE)) { message("edgeR not found, skipping."); next }
        y <- edgeR::DGEList(counts = comp_specific_pb_counts, group = comp_specific_pb_meta$de_group)
        keep <- edgeR::filterByExpr(y, group = comp_specific_pb_meta$de_group,
                                    min.count = edgeR.filterByExpr.min.count,
                                    min.total.count = edgeR.filterByExpr.min.total.count)
        y <- y[keep, , keep.lib.sizes = FALSE]
        
        if (nrow(y$counts) == 0) { message("      No genes left after edgeR::filterByExpr. Skipping."); next; }
        y <- edgeR::calcNormFactors(y)
        design_edgeR <- stats::model.matrix(~ de_group, data = comp_specific_pb_meta)
        
        if (nrow(design_edgeR) < ncol(design_edgeR) || qr(design_edgeR)$rank < ncol(design_edgeR)) {
          message("      Design matrix for edgeR is not full rank or has too few observations. Skipping.")
          next
        }
        y <- edgeR::estimateDisp(y, design_edgeR, robust = TRUE)
        fit <- edgeR::glmQLFit(y, design_edgeR, robust = TRUE)
        qlf <- edgeR::glmQLFTest(fit, coef = ncol(design_edgeR)) # Compares de_groupgroup1 vs de_groupgroup2 (ref)
        res_edgeR_table <- edgeR::topTags(qlf, n = Inf)$table
        
        if (nrow(res_edgeR_table) > 0) {
          res_df <- data.frame(
            gene = rownames(res_edgeR_table),
            p_val = res_edgeR_table$PValue,
            avg_log2FC = res_edgeR_table$logFC,
            p_val_adj = res_edgeR_table$FDR,
            stringsAsFactors = FALSE
          )
          # Calculate pct.1 and pct.2 based on original counts for these genes
          pct_counts_subset <- comp_specific_pb_counts[res_df$gene, , drop = FALSE]
          res_df$pct.1 <- rowMeans(pct_counts_subset[, samples1_names, drop = FALSE] > 0)
          res_df$pct.2 <- rowMeans(pct_counts_subset[, samples2_names, drop = FALSE] > 0)
        }
        
      } else if (test.use == "DESeq2") {
        if (!requireNamespace("DESeq2", quietly = TRUE)) { message("DESeq2 not found, skipping."); next }
        
        dds <- DESeq2::DESeqDataSetFromMatrix(countData = comp_specific_pb_counts,
                                              colData = comp_specific_pb_meta,
                                              design = ~ de_group)
        # Optional: Pre-filter low count genes, though DESeq2 does its own independent filtering.
        # keep_dds <- rowSums(DESeq2::counts(dds) >= 5) >= min.samples.group
        # dds <- dds[keep_dds,]
        if (nrow(dds) == 0) { message("      No genes left after initial filtering for DESeq2. Skipping."); next; }
        
        dds_fit <- tryCatch({
          DESeq2::DESeq(dds, test="Wald", fitType = DESeq2.fitType, quiet = TRUE, parallel = (inherits(BPPARAM, "SerialParam") == FALSE), BPPARAM = BPPARAM)
        }, error = function(e_deseq) {
          message(paste0("      DESeq2 standard fit failed: ", e_deseq$message))
          if (grepl("every gene contains at least one zero|all samples have 0 counts|matrix contains parameters that are not estimable|model matrix not full rank", e_deseq$message)) {
            if (DESeq2.fitType != "local") {
              message("      Attempting DESeq2 with fitType='local'.")
              return(tryCatch(DESeq2::DESeq(dds, test="Wald", fitType = "local", quiet = TRUE, parallel = (inherits(BPPARAM, "SerialParam") == FALSE), BPPARAM = BPPARAM), error=function(e_local) {message(paste0("        DESeq2 fitType='local' also failed: ", e_local$message)); NULL}))
            }
          }
          return(NULL) # Return NULL if all attempts fail
        })
        
        if (is.null(dds_fit) || !"results" %in% mcols(dds_fit)$type) {
          message("      DESeq2 analysis failed to produce results. Skipping.")
          next
        }
        # Contrast group1 vs group2 (group2 is reference)
        res_deseq_table <- DESeq2::results(dds_fit, contrast = c("de_group", "group1", "group2"), parallel = (inherits(BPPARAM, "SerialParam") == FALSE), BPPARAM = BPPARAM)
        
        if (nrow(res_deseq_table) > 0) {
          res_df <- data.frame(
            gene = rownames(res_deseq_table),
            p_val = res_deseq_table$pvalue,
            avg_log2FC = res_deseq_table$log2FoldChange,
            p_val_adj = res_deseq_table$padj,
            stringsAsFactors = FALSE
          )
          res_df <- res_df[!is.na(res_df$p_val), ] # DESeq2 can produce NAs for p_val/padj
          if (nrow(res_df)>0){
            pct_counts_subset <- comp_specific_pb_counts[res_df$gene, , drop = FALSE]
            res_df$pct.1 <- rowMeans(pct_counts_subset[, samples1_names, drop = FALSE] > 0)
            res_df$pct.2 <- rowMeans(pct_counts_subset[, samples2_names, drop = FALSE] > 0)
          }
        }
        
      } else if (test.use %in% c("wilcox", "t.test", "roc")) {
        # Normalize for these tests: log2(CPM+1) or similar
        if (requireNamespace("edgeR", quietly = TRUE)) {
          pb_norm <- edgeR::cpm(comp_specific_pb_counts, log = TRUE, prior.count = 1)
        } else { # Fallback basic normalization if edgeR not available
          pb_norm <- log2(sweep(comp_specific_pb_counts, 2, colSums(comp_specific_pb_counts)/1e6, FUN="*") + 1)
        }
        
        # Filter genes based on min.cells.gene (applied to pseudobulks)
        # Gene must be detected in min.cells.gene samples in EITHER group1 OR group2
        genes_to_test_mask <- sapply(rownames(comp_specific_pb_counts), function(g_name) {
          (sum(comp_specific_pb_counts[g_name, samples1_names] > 0) >= min.cells.gene ||
             sum(comp_specific_pb_counts[g_name, samples2_names] > 0) >= min.cells.gene)
        })
        pb_norm_filtered <- pb_norm[genes_to_test_mask, , drop = FALSE]
        counts_filtered_for_pct <- comp_specific_pb_counts[genes_to_test_mask, , drop = FALSE]
        
        if (nrow(pb_norm_filtered) == 0) { message("      No genes left after min.cells.gene filtering for wilcox/t.test/roc. Skipping."); next; }
        
        results_list_genes <- BiocParallel::bplapply(rownames(pb_norm_filtered), FUN = function(g) {
          vec1 <- pb_norm_filtered[g, samples1_names]
          vec2 <- pb_norm_filtered[g, samples2_names]
          
          p_val_gene <- NA
          auc_val_gene <- NA
          
          # Initial check for t-test variance issues
          if (test.use == "t.test" && (stats::sd(vec1) == 0 && stats::sd(vec2) == 0 && mean(vec1) == mean(vec2))) { # Both groups are identical constants
            p_val_gene <- 1.0 # No difference
          } else if (test.use == "t.test" && (stats::sd(vec1) == 0 || stats::sd(vec2) == 0)) { # One group constant, other varies (or also constant but different mean)
            # t.test might still work or give an error; Wilcoxon is safer here.
            # For now, let it try, but this is a known edge case.
          }
          
          
          if (test.use == "wilcox") {
            p_val_gene <- tryCatch(stats::wilcox.test(vec1, vec2, ...)$p.value, error = function(e) NA)
          } else if (test.use == "t.test") {
            p_val_gene <- tryCatch(stats::t.test(vec1, vec2, ...)$p.value, error = function(e) NA)
          } else if (test.use == "roc") {
            if (!requireNamespace("pROC", quietly = TRUE)) { message("pROC not found for ROC test."); return(NULL) }
            response_roc <- c(rep("group1", length(samples1_names)), rep("group2", length(samples2_names)))
            predictor_roc <- c(vec1, vec2)
            
            if(any(!is.finite(predictor_roc))) {
              auc_val_gene <- NA
            } else {
              # Ensure levels are correctly specified for pROC; group1 should be "case"
              roc_obj <- tryCatch(pROC::roc(response = factor(response_roc, levels=c("group2", "group1")),
                                            predictor = predictor_roc,
                                            quiet = TRUE, direction = "<"), # "<": higher predictor means more likely group1
                                  error = function(e) NULL)
              if (!is.null(roc_obj)) auc_val_gene <- as.numeric(pROC::auc(roc_obj)) else auc_val_gene <- NA
            }
            # p-value for ROC: often use Wilcoxon as Seurat does
            p_val_gene <- tryCatch(stats::wilcox.test(vec1, vec2)$p.value, error = function(e) NA)
          }
          
          avg_lfc_gene <- if (test.use == "roc") auc_val_gene else (mean(vec1, na.rm=TRUE) - mean(vec2, na.rm=TRUE))
          
          # Use original counts for pct calculation (from counts_filtered_for_pct)
          pct1_gene <- mean(counts_filtered_for_pct[g, samples1_names] > 0, na.rm=TRUE)
          pct2_gene <- mean(counts_filtered_for_pct[g, samples2_names] > 0, na.rm=TRUE)
          
          # Apply min.pct and logfc.threshold filters
          pass_min_pct_filter <- (pct1_gene >= min.pct || pct2_gene >= min.pct)
          if (test.use == "roc") {
            pass_logfc_filter <- (!is.na(avg_lfc_gene) && avg_lfc_gene >= logfc.threshold) # AUC should be > threshold
          } else {
            pass_logfc_filter <- (!is.na(avg_lfc_gene) && abs(avg_lfc_gene) >= logfc.threshold)
          }
          
          if (pass_min_pct_filter && pass_logfc_filter) {
            return(data.frame(
              gene = g, p_val = p_val_gene, avg_log2FC = avg_lfc_gene,
              pct.1 = pct1_gene, pct.2 = pct2_gene,
              stringsAsFactors = FALSE
            ))
          } else {
            return(NULL)
          }
        }, BPPARAM = BPPARAM) # End bplapply over genes
        
        res_df <- do.call(rbind, Filter(NROW, results_list_genes))
        if (!is.null(res_df) && nrow(res_df) > 0) {
          if(any(!is.na(res_df$p_val))){
            res_df$p_val_adj <- stats::p.adjust(res_df$p_val, method = "BH")
          } else {
            res_df$p_val_adj <- NA
          }
        }
      } # End if/else for test.use
      
      # --- 6. Store results for the current comparison ---
      if (!is.null(res_df) && nrow(res_df) > 0) {
        res_df$cluster <- cluster_name_for_output
        res_df$bulk_category <- current_bulk_id
        all_results_list[[length(all_results_list) + 1]] <- res_df
      } else {
        message(paste0("      No significant DEG results found for '", cluster_name_for_output, "' in '", current_bulk_id, "' with test '", test.use, "'."))
      }
    } # End loop for comparisons_to_make
  } # End loop for unique_bulk_categories
  
  
  # --- 7. Combine all results and format ---
  if (length(all_results_list) == 0) {
    message("\nNo DEG results generated from any comparison.")
    return(data.frame(gene=character(), cluster=character(), bulk_category=character(),
                      avg_log2FC=numeric(), p_val=numeric(), p_val_adj=numeric(),
                      pct.1=numeric(), pct.2=numeric()))
  }
  
  final_results_df <- do.call(rbind, all_results_list)
  
  # Ensure correct column order and sort
  cols_ordered <- c("gene", "cluster", "bulk_category", "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj")
  # Add any missing columns with NA if they weren't generated
  for(col_name in cols_ordered){
    if(!col_name %in% colnames(final_results_df)){
      final_results_df[[col_name]] <- NA_real_ # Use type-appropriate NA
      if(col_name %in% c("gene", "cluster", "bulk_category")) final_results_df[[col_name]] <- NA_character_
    }
  }
  final_results_df <- final_results_df[, cols_ordered] %>%
    dplyr::arrange(bulk_category, cluster, p_val_adj, p_val)
  
  message("\nDone with all comparisons.")
  return(final_results_df)
}




#' Advanced Pseudobulk FindMarkers function with covariates
#'
#' @param seurat_obj Seurat object
#' @param ident.1 Identity class to define markers for
#' @param ident.2 Identity class(es) to compare against (default: all others)
#' @param group.by Column name in metadata to use for identity classes
#' @param sample.by Column name in metadata defining biological replicates/samples
#' @param covariates Character vector of additional metadata columns to include in the model
#' @param aggregate.by Additional column to create finer pseudobulk aggregations (e.g., subcluster)
#' @param min.cells Minimum number of cells per pseudobulk sample
#' @param method DE method to use: "DESeq2" or "edgeR"
#' @param design.formula Custom design formula (overrides automatic formula generation)
#' @param contrast.type Type of contrast: "simple" or "interaction"
#' @param logfc.threshold Log fold change threshold
#' @param min.pct Minimum percentage of samples expressing the gene
#' @param assay Assay to use
#' @param slot Slot to use (counts recommended for pseudobulk)
#' @return Data frame with DE results
FindMarkers_pseudobulk <- function(seurat_obj,
                                   ident.1,
                                   ident.2 = NULL,
                                   group.by = "seurat_clusters",
                                   sample.by,
                                   covariates = NULL,
                                   aggregate.by = NULL,
                                   min.cells = 10,
                                   method = "DESeq2",
                                   design.formula = NULL,
                                   contrast.type = "simple",
                                   logfc.threshold = 0.25,
                                   min.pct = 0.1,
                                   assay = DefaultAssay(seurat_obj),
                                   slot = "counts") {
  
  # Set identities
  Idents(seurat_obj) <- group.by
  
  # Get cells for each group
  cells.1 <- WhichCells(seurat_obj, idents = ident.1)
  if (is.null(ident.2)) {
    cells.2 <- WhichCells(seurat_obj, idents = setdiff(levels(Idents(seurat_obj)), ident.1))
  } else {
    cells.2 <- WhichCells(seurat_obj, idents = ident.2)
  }
  
  # Subset object
  cells.use <- c(cells.1, cells.2)
  seurat_subset <- subset(seurat_obj, cells = cells.use)
  
  # Create pseudobulk matrix with advanced aggregation
  pb_data <- create_pseudobulk_matrix_advanced(
    seurat_obj = seurat_subset,
    sample.by = sample.by,
    group.by = group.by,
    aggregate.by = aggregate.by,
    covariates = covariates,
    assay = assay,
    slot = slot,
    min.cells = min.cells
  )
  
  pb_matrix <- pb_data$matrix
  pb_metadata <- pb_data$metadata
  
  # Filter for minimum expression
  keep_genes <- rowSums(pb_matrix > 0) >= (ncol(pb_matrix) * min.pct)
  pb_matrix <- pb_matrix[keep_genes, ]
  
  # Run differential expression
  if (method == "DESeq2") {
    de_results <- run_DESeq2_pseudobulk_advanced(
      pb_matrix = pb_matrix,
      metadata = pb_metadata,
      group.by = group.by,
      ident.1 = ident.1,
      ident.2 = ident.2,
      covariates = covariates,
      design.formula = design.formula,
      logfc.threshold = logfc.threshold
    )
  } else if (method == "edgeR") {
    de_results <- run_edgeR_pseudobulk_advanced(
      pb_matrix = pb_matrix,
      metadata = pb_metadata,
      group.by = group.by,
      ident.1 = ident.1,
      ident.2 = ident.2,
      covariates = covariates,
      design.formula = design.formula,
      contrast.type = contrast.type,
      logfc.threshold = logfc.threshold
    )
  } else {
    stop("Method must be either 'DESeq2' or 'edgeR'")
  }
  
  return(de_results)
}

#' Create advanced pseudobulk expression matrix with flexible aggregation
#'
#' @param seurat_obj Seurat object
#' @param sample.by Column name defining samples
#' @param group.by Column name defining groups for comparison
#' @param aggregate.by Additional column for finer aggregation
#' @param covariates Additional metadata columns to include
#' @param assay Assay to use
#' @param slot Slot to use
#' @param min.cells Minimum cells per pseudobulk sample
#' @return List with pseudobulk matrix and metadata
create_pseudobulk_matrix_advanced <- function(seurat_obj,
                                              sample.by,
                                              group.by,
                                              aggregate.by = NULL,
                                              covariates = NULL,
                                              assay = DefaultAssay(seurat_obj),
                                              slot = "counts",
                                              min.cells = 10) {
  
  # Get expression matrix
  expr_matrix <- GetAssayData(seurat_obj, assay = assay, slot = slot)
  
  # Prepare metadata columns
  meta_cols <- c(sample.by, group.by)
  if (!is.null(aggregate.by)) meta_cols <- c(meta_cols, aggregate.by)
  if (!is.null(covariates)) meta_cols <- c(meta_cols, covariates)
  
  # Get metadata
  metadata <- seurat_obj@meta.data[, meta_cols, drop = FALSE]
  metadata$cell_id <- colnames(seurat_obj)
  
  # Create aggregation ID
  if (!is.null(aggregate.by)) {
    # Include aggregate.by in the pseudobulk sample ID
    metadata$pb_sample_id <- paste(metadata[[sample.by]], 
                                   metadata[[group.by]], 
                                   metadata[[aggregate.by]], 
                                   sep = "_")
  } else {
    metadata$pb_sample_id <- paste(metadata[[sample.by]], 
                                   metadata[[group.by]], 
                                   sep = "_")
  }
  
  # Aggregate expression
  pb_list <- list()
  sample_metadata <- list()
  
  for (pb_id in unique(metadata$pb_sample_id)) {
    cells <- metadata$cell_id[metadata$pb_sample_id == pb_id]
    
    if (length(cells) >= min.cells) {
      # Sum expression across cells
      if (length(cells) == 1) {
        pb_expr <- expr_matrix[, cells, drop = FALSE]
      } else {
        pb_expr <- Matrix::rowSums(expr_matrix[, cells, drop = FALSE])
      }
      
      pb_list[[pb_id]] <- pb_expr
      
      # Store metadata for this pseudobulk sample
      sample_info <- metadata[metadata$pb_sample_id == pb_id, ][1, ]
      sample_meta <- data.frame(
        pb_sample_id = pb_id,
        sample = sample_info[[sample.by]],
        group = sample_info[[group.by]],
        n_cells = length(cells),
        stringsAsFactors = FALSE
      )
      
      # Add aggregate.by if present
      if (!is.null(aggregate.by)) {
        sample_meta[[aggregate.by]] <- sample_info[[aggregate.by]]
      }
      
      # Add covariates if present
      if (!is.null(covariates)) {
        for (cov in covariates) {
          sample_meta[[cov]] <- sample_info[[cov]]
        }
      }
      
      sample_metadata[[pb_id]] <- sample_meta
    }
  }
  
  # Create pseudobulk matrix
  pb_matrix <- do.call(cbind, pb_list)
  pb_metadata <- do.call(rbind, sample_metadata)
  rownames(pb_metadata) <- pb_metadata$pb_sample_id
  
  return(list(matrix = pb_matrix, metadata = pb_metadata))
}

#' Run advanced edgeR analysis on pseudobulk data
run_edgeR_pseudobulk_advanced <- function(pb_matrix,
                                          metadata,
                                          group.by,
                                          ident.1,
                                          ident.2 = NULL,
                                          covariates = NULL,
                                          design.formula = NULL,
                                          contrast.type = "simple",
                                          logfc.threshold = 0.25) {
  
  # Prepare metadata
  metadata$group <- factor(metadata$group)
  
  # Set conditions
  if (is.null(ident.2)) {
    metadata$condition <- ifelse(metadata$group == ident.1, "target", "reference")
  } else {
    # Only keep ident.1 and ident.2 samples
    keep_samples <- metadata$group %in% c(ident.1, ident.2)
    pb_matrix <- pb_matrix[, keep_samples]
    metadata <- metadata[keep_samples, ]
    metadata$condition <- ifelse(metadata$group == ident.1, "target", "reference")
  }
  
  metadata$condition <- factor(metadata$condition, levels = c("reference", "target"))
  
  # Create DGEList
  y <- DGEList(counts = round(pb_matrix))
  
  # Filter low expression genes
  keep <- filterByExpr(y)
  y <- y[keep, , keep.lib.sizes = FALSE]
  
  # Normalize
  y <- calcNormFactors(y)
  
  # Create design matrix
  if (!is.null(design.formula)) {
    # Use custom formula
    design <- model.matrix(design.formula, data = metadata)
  } else if (!is.null(covariates)) {
    # Build formula with covariates
    formula_str <- paste("~ 0 + condition", paste(covariates, collapse = " + "), sep = " + ")
    design <- model.matrix(as.formula(formula_str), data = metadata)
  } else {
    # Simple design
    design <- model.matrix(~ 0 + condition, data = metadata)
  }
  
  # Ensure we have condition columns in the design matrix
  if (!any(grepl("condition", colnames(design)))) {
    stop("Design matrix must include condition factor")
  }
  
  # Estimate dispersion
  y <- estimateDisp(y, design)
  
  # Fit model
  fit <- glmQLFit(y, design)
  
  # Make contrast based on contrast type
  if (contrast.type == "simple") {
    # Find the columns for conditions
    target_col <- grep("conditiontarget", colnames(design))
    ref_col <- grep("conditionreference", colnames(design))
    
    if (length(target_col) == 0 || length(ref_col) == 0) {
      # Try alternative naming
      contrast <- makeContrasts(
        conditiontarget - conditionreference, 
        levels = design
      )
    } else {
      contrast_vector <- numeric(ncol(design))
      contrast_vector[target_col] <- 1
      contrast_vector[ref_col] <- -1
      contrast <- contrast_vector
    }
  } else if (contrast.type == "interaction") {
    # For interaction terms, user should provide custom contrast
    stop("For interaction contrasts, please provide a custom design formula")
  }
  
  # Perform test
  qlf <- glmQLFTest(fit, contrast = contrast)
  
  # Get results
  de_results <- topTags(qlf, n = Inf)$table
  
  # Format results
  de_results$avg_log2FC <- de_results$logFC
  de_results$p_val <- de_results$PValue
  de_results$p_val_adj <- de_results$FDR
  de_results$avg_expr <- de_results$logCPM
  
  # Filter by log fold change
  de_results <- de_results[abs(de_results$avg_log2FC) > logfc.threshold, ]
  
  # Select columns
  de_results <- de_results[, c("avg_log2FC", "p_val", "p_val_adj", "avg_expr")]
  
  # Sort by adjusted p-value
  de_results <- de_results[order(de_results$p_val_adj), ]
  
  return(de_results)
}

#' Run advanced DESeq2 analysis on pseudobulk data
run_DESeq2_pseudobulk_advanced <- function(pb_matrix,
                                           metadata,
                                           group.by,
                                           ident.1,
                                           ident.2 = NULL,
                                           covariates = NULL,
                                           design.formula = NULL,
                                           logfc.threshold = 0.25) {
  
  # Prepare metadata
  metadata$group <- factor(metadata$group)
  
  # Set conditions
  if (is.null(ident.2)) {
    metadata$condition <- ifelse(metadata$group == ident.1, "target", "reference")
  } else {
    keep_samples <- metadata$group %in% c(ident.1, ident.2)
    pb_matrix <- pb_matrix[, keep_samples]
    metadata <- metadata[keep_samples, ]
    metadata$condition <- ifelse(metadata$group == ident.1, "target", "reference")
  }
  
  metadata$condition <- factor(metadata$condition, levels = c("reference", "target"))
  
  # Create design
  if (!is.null(design.formula)) {
    design <- design.formula
  } else if (!is.null(covariates)) {
    formula_str <- paste("~ condition", paste(covariates, collapse = " + "), sep = " + ")
    design <- as.formula(formula_str)
  } else {
    design <- ~ condition
  }
  
  # Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(
    countData = round(pb_matrix),
    colData = metadata,
    design = design
  )
  
  # Filter low counts
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep, ]
  
  # Run DESeq2
  dds <- DESeq(dds)
  
  # Get results
  res <- results(dds, contrast = c("condition", "target", "reference"))
  
  # Format results
  de_results <- as.data.frame(res)
  de_results$avg_log2FC <- de_results$log2FoldChange
  de_results$p_val <- de_results$pvalue
  de_results$p_val_adj <- de_results$padj
  
  # Filter
  de_results <- de_results[!is.na(de_results$p_val_adj), ]
  de_results <- de_results[abs(de_results$avg_log2FC) > logfc.threshold, ]
  
  # Select columns
  de_results <- de_results[, c("avg_log2FC", "p_val", "p_val_adj", "baseMean")]
  names(de_results)[4] <- "avg_expr"
  
  # Sort
  de_results <- de_results[order(de_results$p_val_adj), ]
  
  return(de_results)
}

#' Pseudobulk FindAllMarkers function
#'
#' @param seurat_obj Seurat object
#' @param group.by Column name in metadata to use for identity classes
#' @param sample.by Column name in metadata defining biological replicates/samples
#' @param covariates Character vector of additional metadata columns to include in the model
#' @param aggregate.by Additional column to create finer pseudobulk aggregations
#' @param min.cells Minimum number of cells per sample to include
#' @param method DE method to use: "DESeq2" or "edgeR"
#' @param design.formula Custom design formula
#' @param logfc.threshold Log fold change threshold
#' @param min.pct Minimum percentage of samples expressing the gene
#' @param only.pos Only return positive markers
#' @param assay Assay to use
#' @param slot Slot to use (counts recommended for pseudobulk)
#' @return Data frame with DE results for all clusters
FindAllMarkers_pseudobulk <- function(seurat_obj,
                                      group.by = "seurat_clusters",
                                      sample.by,
                                      covariates = NULL,
                                      aggregate.by = NULL,
                                      min.cells = 10,
                                      method = "DESeq2",
                                      design.formula = NULL,
                                      logfc.threshold = 0.25,
                                      min.pct = 0.1,
                                      only.pos = FALSE,
                                      assay = DefaultAssay(seurat_obj),
                                      slot = "counts") {
  
  # Get all unique identities
  Idents(seurat_obj) <- group.by
  all_idents <- levels(Idents(seurat_obj))
  
  # Run FindMarkers for each identity
  all_markers <- list()
  
  for (ident in all_idents) {
    message(paste0("Finding markers for ", group.by, " ", ident))
    
    tryCatch({
      markers <- FindMarkers_pseudobulk(
        seurat_obj = seurat_obj,
        ident.1 = ident,
        ident.2 = NULL,
        group.by = group.by,
        sample.by = sample.by,
        covariates = covariates,
        aggregate.by = aggregate.by,
        min.cells = min.cells,
        method = method,
        design.formula = design.formula,
        logfc.threshold = logfc.threshold,
        min.pct = min.pct,
        assay = assay,
        slot = slot
      )
      
      if (nrow(markers) > 0) {
        markers$gene <- rownames(markers)
        markers$cluster <- ident
        all_markers[[as.character(ident)]] <- markers
      }
    }, error = function(e) {
      warning(paste0("Failed to find markers for ", group.by, " ", ident, ": ", e$message))
    })
  }
  
  # Combine all results
  all_markers_df <- do.call(rbind, all_markers)
  rownames(all_markers_df) <- NULL
  
  # Filter for only positive markers if requested
  if (only.pos) {
    all_markers_df <- all_markers_df[all_markers_df$avg_log2FC > 0, ]
  }
  
  # Sort by cluster and p-value
  all_markers_df <- all_markers_df[order(all_markers_df$cluster, all_markers_df$p_val_adj), ]
  
  return(all_markers_df)
}

# Example usage:
# 1. Simple comparison with sample aggregation
# markers_cluster1 <- FindMarkers_pseudobulk(
#   seurat_obj = pbmc,
#   ident.1 = "1",
#   sample.by = "patient_id",
#   group.by = "seurat_clusters",
#   method = "edgeR"
# )
#
# 2. Comparison with finer aggregation by subclusters
# markers_with_subcluster <- FindMarkers_pseudobulk(
#   seurat_obj = pbmc,
#   ident.1 = "NK_cells",
#   ident.2 = "T_cells",
#   sample.by = "patient_id",
#   group.by = "cell_type",
#   aggregate.by = "subcluster",  # Creates separate pseudobulk for each patient-celltype-subcluster
#   method = "edgeR"
# )
#
# 3. Including batch effect as covariate
# markers_batch_corrected <- FindMarkers_pseudobulk(
#   seurat_obj = pbmc,
#   ident.1 = "1",
#   sample.by = "patient_id",
#   group.by = "seurat_clusters",
#   covariates = c("batch", "sex"),  # Include batch and sex in the model
#   method = "edgeR"
# )
#
# 4. Custom design formula for complex comparisons
# markers_interaction <- FindMarkers_pseudobulk(
#   seurat_obj = pbmc,
#   ident.1 = "treated",
#   ident.2 = "control",
#   sample.by = "patient_id",
#   group.by = "treatment",
#   covariates = c("timepoint"),
#   design.formula = ~ condition + timepoint + condition:timepoint,
#   method = "edgeR"
# )




#' Advanced Pseudobulk FindMarkers function with covariates
#'
#' @param seurat_obj Seurat object
#' @param ident.1 Identity class to define markers for
#' @param ident.2 Identity class(es) to compare against (default: all others)
#' @param group.by Column name in metadata to use for identity classes
#' @param sample.by Column name in metadata defining biological replicates/samples
#' @param covariates Character vector of additional metadata columns to include in the model
#' @param aggregate.by Additional column to create finer pseudobulk aggregations (e.g., subcluster)
#' @param min.cells Minimum number of cells per pseudobulk sample
#' @param method DE method to use: "DESeq2" or "edgeR"
#' @param design.formula Custom design formula (overrides automatic formula generation)
#' @param contrast.type Type of contrast: "simple" or "interaction"
#' @param logfc.threshold Log fold change threshold
#' @param min.pct Minimum percentage of samples expressing the gene
#' @param assay Assay to use
#' @param slot Slot to use (counts recommended for pseudobulk)
#' @return Data frame with DE results
FindMarkers_pseudobulk <- function(seurat_obj,
                                   ident.1,
                                   ident.2 = NULL,
                                   group.by = "seurat_clusters",
                                   sample.by,
                                   covariates = NULL,
                                   aggregate.by = NULL,
                                   min.cells = 10,
                                   method = "DESeq2",
                                   design.formula = NULL,
                                   contrast.type = "simple",
                                   logfc.threshold = 0.25,
                                   min.pct = 0.1,
                                   assay = DefaultAssay(seurat_obj),
                                   slot = "counts") {
  
  # Set identities
  Idents(seurat_obj) <- group.by
  
  # Get cells for each group
  cells.1 <- WhichCells(seurat_obj, idents = ident.1)
  if (is.null(ident.2)) {
    cells.2 <- WhichCells(seurat_obj, idents = setdiff(levels(Idents(seurat_obj)), ident.1))
  } else {
    cells.2 <- WhichCells(seurat_obj, idents = ident.2)
  }
  
  # Subset object
  cells.use <- c(cells.1, cells.2)
  seurat_subset <- subset(seurat_obj, cells = cells.use)
  
  # Create pseudobulk matrix with advanced aggregation
  pb_data <- create_pseudobulk_matrix_advanced(
    seurat_obj = seurat_subset,
    sample.by = sample.by,
    group.by = group.by,
    aggregate.by = aggregate.by,
    covariates = covariates,
    assay = assay,
    slot = slot,
    min.cells = min.cells
  )
  
  pb_matrix <- pb_data$matrix
  pb_metadata <- pb_data$metadata
  
  # Filter for minimum expression
  keep_genes <- rowSums(pb_matrix > 0) >= (ncol(pb_matrix) * min.pct)
  pb_matrix <- pb_matrix[keep_genes, ]
  
  # Run differential expression
  if (method == "DESeq2") {
    de_results <- run_DESeq2_pseudobulk_advanced(
      pb_matrix = pb_matrix,
      metadata = pb_metadata,
      group.by = group.by,
      ident.1 = ident.1,
      ident.2 = ident.2,
      covariates = covariates,
      design.formula = design.formula,
      logfc.threshold = logfc.threshold
    )
  } else if (method == "edgeR") {
    de_results <- run_edgeR_pseudobulk_advanced(
      pb_matrix = pb_matrix,
      metadata = pb_metadata,
      group.by = group.by,
      ident.1 = ident.1,
      ident.2 = ident.2,
      covariates = covariates,
      design.formula = design.formula,
      contrast.type = contrast.type,
      logfc.threshold = logfc.threshold
    )
  } else {
    stop("Method must be either 'DESeq2' or 'edgeR'")
  }
  
  return(de_results)
}

#' Create advanced pseudobulk expression matrix with flexible aggregation
#'
#' @param seurat_obj Seurat object
#' @param sample.by Column name defining samples
#' @param group.by Column name defining groups for comparison
#' @param aggregate.by Additional column for finer aggregation
#' @param covariates Additional metadata columns to include
#' @param assay Assay to use
#' @param slot Slot to use
#' @param min.cells Minimum cells per pseudobulk sample
#' @return List with pseudobulk matrix and metadata
create_pseudobulk_matrix_advanced <- function(seurat_obj,
                                              sample.by,
                                              group.by,
                                              aggregate.by = NULL,
                                              covariates = NULL,
                                              assay = DefaultAssay(seurat_obj),
                                              slot = "counts",
                                              min.cells = 10) {
  
  # Get expression matrix
  expr_matrix <- GetAssayData(seurat_obj, assay = assay, slot = slot)
  
  # Prepare metadata columns
  meta_cols <- c(sample.by, group.by)
  if (!is.null(aggregate.by)) meta_cols <- c(meta_cols, aggregate.by)
  if (!is.null(covariates)) meta_cols <- c(meta_cols, covariates)
  
  # Get metadata
  metadata <- seurat_obj@meta.data[, meta_cols, drop = FALSE]
  metadata$cell_id <- colnames(seurat_obj)
  
  # Create aggregation ID
  if (!is.null(aggregate.by)) {
    # Include aggregate.by in the pseudobulk sample ID
    metadata$pb_sample_id <- paste(metadata[[sample.by]], 
                                   metadata[[group.by]], 
                                   metadata[[aggregate.by]], 
                                   sep = "_")
  } else {
    metadata$pb_sample_id <- paste(metadata[[sample.by]], 
                                   metadata[[group.by]], 
                                   sep = "_")
  }
  
  # Aggregate expression
  pb_list <- list()
  sample_metadata <- list()
  
  for (pb_id in unique(metadata$pb_sample_id)) {
    cells <- metadata$cell_id[metadata$pb_sample_id == pb_id]
    
    if (length(cells) >= min.cells) {
      # Sum expression across cells
      if (length(cells) == 1) {
        pb_expr <- expr_matrix[, cells, drop = FALSE]
      } else {
        pb_expr <- Matrix::rowSums(expr_matrix[, cells, drop = FALSE])
      }
      
      pb_list[[pb_id]] <- pb_expr
      
      # Store metadata for this pseudobulk sample
      sample_info <- metadata[metadata$pb_sample_id == pb_id, ][1, ]
      sample_meta <- data.frame(
        pb_sample_id = pb_id,
        sample = sample_info[[sample.by]],
        group = sample_info[[group.by]],
        n_cells = length(cells),
        stringsAsFactors = FALSE
      )
      
      # Add aggregate.by if present
      if (!is.null(aggregate.by)) {
        sample_meta[[aggregate.by]] <- sample_info[[aggregate.by]]
      }
      
      # Add covariates if present
      if (!is.null(covariates)) {
        for (cov in covariates) {
          sample_meta[[cov]] <- sample_info[[cov]]
        }
      }
      
      sample_metadata[[pb_id]] <- sample_meta
    }
  }
  
  # Create pseudobulk matrix
  pb_matrix <- do.call(cbind, pb_list)
  pb_metadata <- do.call(rbind, sample_metadata)
  rownames(pb_metadata) <- pb_metadata$pb_sample_id
  
  return(list(matrix = pb_matrix, metadata = pb_metadata))
}

#' Run advanced edgeR analysis on pseudobulk data
run_edgeR_pseudobulk_advanced <- function(pb_matrix,
                                          metadata,
                                          group.by,
                                          ident.1,
                                          ident.2 = NULL,
                                          covariates = NULL,
                                          design.formula = NULL,
                                          contrast.type = "simple",
                                          logfc.threshold = 0.25) {
  
  # Prepare metadata
  metadata$group <- factor(metadata$group)
  
  # Set conditions
  if (is.null(ident.2)) {
    metadata$condition <- ifelse(metadata$group == ident.1, "target", "reference")
  } else {
    # Only keep ident.1 and ident.2 samples
    keep_samples <- metadata$group %in% c(ident.1, ident.2)
    pb_matrix <- pb_matrix[, keep_samples]
    metadata <- metadata[keep_samples, ]
    metadata$condition <- ifelse(metadata$group == ident.1, "target", "reference")
  }
  
  metadata$condition <- factor(metadata$condition, levels = c("reference", "target"))
  
  # Create DGEList
  y <- DGEList(counts = round(pb_matrix),
               samples = metadata,
               group   = metadata$condition)
  
  # ② filterByExpr 호출 시 design 전달
  keep <- filterByExpr(y, design = model.matrix(~0 + condition, data = metadata))
  y <- y[keep, , keep.lib.sizes = FALSE]
  
  # Normalize
  y <- calcNormFactors(y)
  
  # Create design matrix
  if (!is.null(design.formula)) {
    # Use custom formula
    design <- model.matrix(design.formula, data = metadata)
  } else if (!is.null(covariates)) {
    # Build formula with covariates
    formula_str <- paste("~ 0 + condition", paste(covariates, collapse = " + "), sep = " + ")
    design <- model.matrix(as.formula(formula_str), data = metadata)
  } else {
    # Simple design
    design <- model.matrix(~ 0 + condition, data = metadata)
  }
  
  # Ensure we have condition columns in the design matrix
  if (!any(grepl("condition", colnames(design)))) {
    stop("Design matrix must include condition factor")
  }
  
  # Estimate dispersion
  y <- estimateDisp(y, design)
  
  # Fit model
  fit <- glmQLFit(y, design)
  
  # Make contrast based on contrast type
  if (contrast.type == "simple") {
    # Find the columns for conditions
    target_col <- grep("conditiontarget", colnames(design))
    ref_col <- grep("conditionreference", colnames(design))
    
    if (length(target_col) == 0 || length(ref_col) == 0) {
      # Try alternative naming
      contrast <- makeContrasts(
        conditiontarget - conditionreference, 
        levels = design
      )
    } else {
      contrast_vector <- numeric(ncol(design))
      contrast_vector[target_col] <- 1
      contrast_vector[ref_col] <- -1
      contrast <- contrast_vector
    }
  } else if (contrast.type == "interaction") {
    # For interaction terms, user should provide custom contrast
    stop("For interaction contrasts, please provide a custom design formula")
  }
  
  # Perform test
  qlf <- glmQLFTest(fit, contrast = contrast)
  
  # Get results
  de_results <- topTags(qlf, n = Inf)$table
  
  # Format results
  de_results$avg_log2FC <- de_results$logFC
  de_results$p_val <- de_results$PValue
  de_results$p_val_adj <- de_results$FDR
  de_results$avg_expr <- de_results$logCPM
  
  # Filter by log fold change
  de_results <- de_results[abs(de_results$avg_log2FC) > logfc.threshold, ]
  
  # Select columns
  de_results <- de_results[, c("avg_log2FC", "p_val", "p_val_adj", "avg_expr")]
  
  # Sort by adjusted p-value
  de_results <- de_results[order(de_results$p_val_adj), ]
  
  return(de_results)
}

#' Run advanced DESeq2 analysis on pseudobulk data
run_DESeq2_pseudobulk_advanced <- function(pb_matrix,
                                           metadata,
                                           group.by,
                                           ident.1,
                                           ident.2 = NULL,
                                           covariates = NULL,
                                           design.formula = NULL,
                                           logfc.threshold = 0.25) {
  
  # Prepare metadata
  metadata$group <- factor(metadata$group)
  
  # Set conditions
  if (is.null(ident.2)) {
    metadata$condition <- ifelse(metadata$group == ident.1, "target", "reference")
  } else {
    keep_samples <- metadata$group %in% c(ident.1, ident.2)
    pb_matrix <- pb_matrix[, keep_samples]
    metadata <- metadata[keep_samples, ]
    metadata$condition <- ifelse(metadata$group == ident.1, "target", "reference")
  }
  
  metadata$condition <- factor(metadata$condition, levels = c("reference", "target"))
  
  # Create design
  if (!is.null(design.formula)) {
    design <- design.formula
  } else if (!is.null(covariates)) {
    formula_str <- paste("~ condition", paste(covariates, collapse = " + "), sep = " + ")
    design <- as.formula(formula_str)
  } else {
    design <- ~ condition
  }
  
  # Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(
    countData = round(pb_matrix),
    colData = metadata,
    design = design
  )
  
  # Filter low counts
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep, ]
  
  # Run DESeq2
  dds <- DESeq(dds)
  
  # Get results
  res <- results(dds, contrast = c("condition", "target", "reference"))
  
  # Format results
  de_results <- as.data.frame(res)
  de_results$avg_log2FC <- de_results$log2FoldChange
  de_results$p_val <- de_results$pvalue
  de_results$p_val_adj <- de_results$padj
  
  # Filter
  de_results <- de_results[!is.na(de_results$p_val_adj), ]
  de_results <- de_results[abs(de_results$avg_log2FC) > logfc.threshold, ]
  
  # Select columns
  de_results <- de_results[, c("avg_log2FC", "p_val", "p_val_adj", "baseMean")]
  names(de_results)[4] <- "avg_expr"
  
  # Sort
  de_results <- de_results[order(de_results$p_val_adj), ]
  
  return(de_results)
}

#' Pseudobulk FindAllMarkers function
#'
#' @param seurat_obj Seurat object
#' @param group.by Column name in metadata to use for identity classes
#' @param sample.by Column name in metadata defining biological replicates/samples
#' @param covariates Character vector of additional metadata columns to include in the model
#' @param aggregate.by Additional column to create finer pseudobulk aggregations
#' @param min.cells Minimum number of cells per sample to include
#' @param method DE method to use: "DESeq2" or "edgeR"
#' @param design.formula Custom design formula
#' @param logfc.threshold Log fold change threshold
#' @param min.pct Minimum percentage of samples expressing the gene
#' @param only.pos Only return positive markers
#' @param assay Assay to use
#' @param slot Slot to use (counts recommended for pseudobulk)
#' @return Data frame with DE results for all clusters
FindAllMarkers_pseudobulk <- function(seurat_obj,
                                      group.by = "seurat_clusters",
                                      sample.by,
                                      covariates = NULL,
                                      aggregate.by = NULL,
                                      min.cells = 10,
                                      method = "DESeq2",
                                      design.formula = NULL,
                                      logfc.threshold = 0.25,
                                      min.pct = 0.1,
                                      only.pos = FALSE,
                                      assay = DefaultAssay(seurat_obj),
                                      slot = "counts") {
  
  # Get all unique identities
  Idents(seurat_obj) <- group.by
  all_idents <- levels(Idents(seurat_obj))
  
  # Run FindMarkers for each identity
  all_markers <- list()
  
  for (ident in all_idents) {
    message(paste0("Finding markers for ", group.by, " ", ident))
    
    tryCatch({
      markers <- FindMarkers_pseudobulk(
        seurat_obj = seurat_obj,
        ident.1 = ident,
        ident.2 = NULL,
        group.by = group.by,
        sample.by = sample.by,
        covariates = covariates,
        aggregate.by = aggregate.by,
        min.cells = min.cells,
        method = method,
        design.formula = design.formula,
        logfc.threshold = logfc.threshold,
        min.pct = min.pct,
        assay = assay,
        slot = slot
      )
      
      if (nrow(markers) > 0) {
        markers$gene <- rownames(markers)
        markers$cluster <- ident
        all_markers[[as.character(ident)]] <- markers
      }
    }, error = function(e) {
      warning(paste0("Failed to find markers for ", group.by, " ", ident, ": ", e$message))
    })
  }
  
  # Combine all results
  all_markers_df <- do.call(rbind, all_markers)
  rownames(all_markers_df) <- NULL
  
  # Filter for only positive markers if requested
  if (only.pos) {
    all_markers_df <- all_markers_df[all_markers_df$avg_log2FC > 0, ]
  }
  
  # Sort by cluster and p-value
  all_markers_df <- all_markers_df[order(all_markers_df$cluster, all_markers_df$p_val_adj), ]
  
  return(all_markers_df)
}

# Example usage:
# 1. Simple comparison with sample aggregation
# markers_cluster1 <- FindMarkers_pseudobulk(
#   seurat_obj = pbmc,
#   ident.1 = "1",
#   sample.by = "patient_id",
#   group.by = "seurat_clusters",
#   method = "edgeR"
# )
#
# 2. Comparison with finer aggregation by subclusters
# markers_with_subcluster <- FindMarkers_pseudobulk(
#   seurat_obj = pbmc,
#   ident.1 = "NK_cells",
#   ident.2 = "T_cells",
#   sample.by = "patient_id",
#   group.by = "cell_type",
#   aggregate.by = "subcluster",  # Creates separate pseudobulk for each patient-celltype-subcluster
#   method = "edgeR"
# )
#
# 3. Including batch effect as covariate
# markers_batch_corrected <- FindMarkers_pseudobulk(
#   seurat_obj = pbmc,
#   ident.1 = "1",
#   sample.by = "patient_id",
#   group.by = "seurat_clusters",
#   covariates = c("batch", "sex"),  # Include batch and sex in the model
#   method = "edgeR"
# )
#
# 4. Custom design formula for complex comparisons
# markers_interaction <- FindMarkers_pseudobulk(
#   seurat_obj = pbmc,
#   ident.1 = "treated",
#   ident.2 = "control",
#   sample.by = "patient_id",
#   group.by = "treatment",
#   covariates = c("timepoint"),
#   design.formula = ~ condition + timepoint + condition:timepoint,
#   method = "edgeR"
# )
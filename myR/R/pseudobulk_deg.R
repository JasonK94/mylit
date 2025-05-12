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
  # Need to handle potential differences in sample naming (e.g., if group_col source uses 'g' prefix)
  # Create a distinct mapping table first
  sample_group_map <- seurat_obj@meta.data %>%
    select(!!sym(sample_col), !!sym(group_col)) %>%
    distinct() %>%
    mutate(across(!!sym(sample_col), as.character)) # Ensure character type for matching
  
  # Check if mapping is one-to-one (one group per sample)
  if (any(duplicated(sample_group_map[[sample_col]]))) {
    warning("하나의 샘플(", sample_col, ")에 여러 그룹(", group_col, ")이 매핑됩니다. 첫 번째 그룹을 사용합니다.")
    sample_group_map <- sample_group_map %>% distinct(!!sym(sample_col), .keep_all = TRUE)
  }
  
  
  meta_pb <- meta_pb %>%
    mutate(patient = as.character(patient)) %>% # Ensure character for matching
    left_join(sample_group_map, by = setNames(sample_col, "patient")) %>% # Use setNames for dynamic join column name
    select(pb_sample_id, patient, ctype, !!sym(group_col))
  
  # Check if all group values were mapped
  if (any(is.na(meta_pb[[group_col]]))) {
    warning("일부 pseudo-bulk 샘플에 그룹 정보를 매핑하지 못했습니다. ",
            sample_col, " 컬럼의 값이 메타데이터와 pseudo-bulk 컬럼 이름 간에 일치하는지 확인하세요.")
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
                               target_cluster = NULL,
                               min_samples_per_group_cluster = 2,
                               ...,
                               verbose = TRUE) {
  
  analysis_level <- match.arg(analysis_level)
  dot_args <- list(...)
  
  # Check if prepared data is provided or needs to be generated
  if (all(c("pb", "meta", "dge", "design", "contrast_levels") %in% names(dot_args))) {
    if (verbose) message("Using pre-prepared pseudo-bulk data.")
    prep_data <- dot_args
    # Extract necessary components
    pb   <- prep_data$pb
    meta <- prep_data$meta
    dge  <- prep_data$dge # This dge is for the 'overall' case initially
    design <- prep_data$design
    contrast_levels <- prep_data$contrast_levels
    # Need to get group_col name from the prep_data structure if possible, or assume from contrast
    # This part is tricky without knowing how prep_data was made. Let's try to find group_col.
    group_col <- setdiff(colnames(meta), c("pb_sample_id", "patient", "ctype"))[1] # Best guess
    if(length(group_col) != 1 || !group_col %in% colnames(meta)) {
      warning("Could not automatically determine group_col from provided pre-prepared data. ",
              "Ensure 'meta' dataframe contains the grouping column.")
      # We might need to rely on contrast length or levels in design matrix columns
      group_col <- names(attr(design, "contrasts"))[1] # Another guess
      if(is.null(group_col)) stop("Cannot determine grouping column.")
    }
    
    
  } else if ("seurat_obj" %in% names(dot_args)) {
    if (verbose) message("Preparing pseudo-bulk data using prepare_pseudobulk_edgeR...")
    # Check if necessary args for prepare_pseudobulk_edgeR are present
    req_args_prep <- c("seurat_obj", "sample_col", "cluster_col", "group_col")
    if (!all(req_args_prep %in% names(dot_args))) {
      stop("Seurat 객체를 사용할 경우, 다음 인수가 필요합니다: ", paste(req_args_prep, collapse=", "))
    }
    # Prepare data using the helper function
    prep_data <- do.call(prepare_pseudobulk_edgeR, c(dot_args, list(verbose = verbose)))
    pb   <- prep_data$pb
    meta <- prep_data$meta
    dge  <- prep_data$dge
    design <- prep_data$design
    contrast_levels <- prep_data$contrast_levels
    group_col <- dot_args$group_col # Get group_col from input args
  } else {
    stop("Seurat 객체(seurat_obj) 또는 준비된 pseudo-bulk 데이터 리스트(pb, meta, dge, design 포함)를 제공해야 합니다.")
  }
  
  # Validate contrast
  if (missing(contrast)) stop("'contrast' 인수는 필수입니다.")
  expected_contrast_len <- ncol(design) - (!grepl("~ 0 +|~0+", deparse(attributes(design)$terms))) # Adjust based on intercept
  if (length(contrast) != expected_contrast_len && analysis_level == "overall") {
    warning("제공된 contrast 길이(", length(contrast), ")가 디자인 매트릭스 컬럼 수(", ncol(design), ", 인터셉트 제외 시 ", expected_contrast_len ,")와 다를 수 있습니다. ",
            "Contrast levels: ", paste(contrast_levels, collapse=", "), "; Design columns: ", paste(colnames(design), collapse=", "))
  }
  
  
  # --- Perform DEG based on analysis level ---
  
  results <- NULL
  
  if (analysis_level == "overall") {
    if (verbose) message("Performing 'overall' DEG analysis...")
    # Use the already prepared dge and design
    dge <- estimateDisp(dge, design)
    fit <- glmQLFit(dge, design)
    qlf <- glmQLFTest(fit, contrast = contrast)
    results <- topTags(qlf, n = Inf)$table %>%
      rownames_to_column("gene") %>%
      tibble()
    if (verbose) message("  Overall analysis complete. Found ", nrow(results), " genes.")
    
  } else if (analysis_level == "per_cluster" || analysis_level == "specific_cluster") {
    
    target_clusters <- if (analysis_level == "specific_cluster") {
      if (is.null(target_cluster)) stop("'specific_cluster' 분석 시 'target_cluster'를 지정해야 합니다.")
      if (!target_cluster %in% unique(meta$ctype)) stop("지정된 'target_cluster' (", target_cluster, ")가 메타데이터에 없습니다.")
      target_cluster
    } else {
      unique(meta$ctype)
    }
    
    if (verbose) message("Performing analysis for cluster(s): ", paste(target_clusters, collapse=", "))
    
    res_list <- list()
    
    for (cl in target_clusters) {
      if (verbose) message("  Processing cluster: ", cl)
      ix <- meta$ctype == cl
      if (sum(ix) == 0) {
        warning("클러스터 '", cl, "'에 해당하는 pseudo-bulk 샘플이 없습니다. 건너<0xEB>니다.")
        next
      }
      pb_sub <- pb[, ix, drop = FALSE]
      md_sub <- meta[ix, , drop = FALSE]
      
      # Check for sufficient samples per group within the cluster
      sample_counts <- md_sub %>% count(!!sym(group_col))
      if(any(sample_counts$n < min_samples_per_group_cluster)) {
        warning("클러스터 '", cl, "'는 그룹별 최소 샘플 수(", min_samples_per_group_cluster, ")를 만족하지 못합니다: ",
                paste(sample_counts[[group_col]], "=", sample_counts$n, collapse=", "), ". 건너<0xEB>니다.")
        next
      }
      # Ensure group is factor with all original levels for consistent design matrix structure
      md_sub[[group_col]] <- factor(md_sub[[group_col]], levels = contrast_levels)
      
      
      # 1) DGEList (subsetted)
      dge_sub <- DGEList(counts = pb_sub, group = md_sub[[group_col]], samples = md_sub)
      
      # 2) Filtering & Normalization (on subset)
      # Use design matrix appropriate for this subset comparison
      # Create design formula dynamically, preferring no intercept for within-cluster
      formula_sub_str <- paste("~", group_col) # Or use ~ 0 + group_col if preferred
      design_sub <- tryCatch({
        model.matrix(as.formula(formula_sub_str), data = md_sub)
      }, error = function(e) {
        warning("클러스터 '", cl, "'에 대한 디자인 매트릭스 생성 실패: ", e$message, ". 건너<0xEB>니다.")
        return(NULL) # Return NULL to skip this cluster
      })
      
      if(is.null(design_sub)) next # Skip if design matrix failed
      
      colnames(design_sub) <- make.names(colnames(design_sub)) # Clean names
      
      keep_sub <- filterByExpr(dge_sub, design = design_sub, group = md_sub[[group_col]], min.count = min_count) # Use design here
      if (verbose) message("    - Filtering for '", cl, "': Kept ", sum(keep_sub), " out of ", nrow(dge_sub), " genes.")
      if(sum(keep_sub) == 0) {
        warning("클러스터 '", cl, "'에서 filterByExpr 후 남은 유전자가 없습니다. 건너<0xEB>니다.")
        next
      }
      dge_sub <- dge_sub[keep_sub, , keep.lib.sizes = FALSE]
      dge_sub <- calcNormFactors(dge_sub)
      
      # Validate contrast for subset design matrix
      expected_contrast_len_sub <- ncol(design_sub) - (!grepl("~ 0 +|~0+", deparse(attributes(design_sub)$terms)))
      if (length(contrast) != expected_contrast_len_sub) {
        warning("클러스터 '", cl, "' 분석에서 contrast 길이(", length(contrast), ")가 디자인 매트릭스 컬럼 수(", ncol(design_sub), ", 인터셉트 고려 시 ", expected_contrast_len_sub,")와 다를 수 있습니다. ",
                "Design columns: ", paste(colnames(design_sub), collapse=", "))
        # Attempt to create contrast dynamically if possible (e.g., for treatment vs control)
        # This requires knowing which columns correspond to which groups, which can be complex.
        # For now, rely on user providing correct contrast based on levels and formula.
      }
      
      
      # 3) Dispersion Estimation & Fit (on subset)
      dge_sub <- tryCatch({
        estimateDisp(dge_sub, design_sub)
      }, warning = function(w){
        warning("Dispersion estimation 경고 (클러스터 ", cl,"): ", w$message)
        # Try robust estimation if default fails or gives warning (sometimes helps)
        message("    - Robust dispersion estimation 시도 중...")
        estimateDisp(dge_sub, design_sub, robust=TRUE)
      }, error = function(e){
        warning("Dispersion estimation 에러 (클러스터 ", cl,"): ", e$message, ". 건너<0xEB>니다.")
        return(NULL)
      })
      
      if(is.null(dge_sub)) next # Skip if dispersion failed
      
      fit_sub <- tryCatch({
        glmQLFit(dge_sub, design_sub)
      }, error = function(e){
        warning("glmQLFit 에러 (클러스터 ", cl, "): ", e$message, ". 건너<0xEB>니다.")
        return(NULL)
      })
      
      if(is.null(fit_sub)) next # Skip if fit failed
      
      # 4) QL F-Test (on subset)
      qlf_sub <- tryCatch({
        glmQLFTest(fit_sub, contrast = contrast)
      }, error = function(e){
        warning("glmQLFTest 에러 (클러스터 ", cl, "): ", e$message, ". 건너<0xEB>니다.")
        return(NULL)
      })
      
      if(is.null(qlf_sub)) next # Skip if test failed
      
      # 5) Results
      res_sub <- topTags(qlf_sub, n = Inf)$table %>%
        rownames_to_column("gene") %>%
        mutate(cluster = cl, .before = 1) %>%
        tibble()
      
      res_list[[cl]] <- res_sub
      if (verbose) message("    - Cluster '", cl, "' analysis complete. Found ", nrow(res_sub), " genes.")
      
    } # end for loop over clusters
    
    if (length(res_list) > 0) {
      results <- bind_rows(res_list)
    } else {
      warning("분석 수준 '", analysis_level, "'에 대한 결과를 생성하지 못했습니다.")
      results <- tibble() # Return empty tibble
    }
    
  } else {
    stop("알 수 없는 analysis_level 입니다: ", analysis_level)
  }
  
  return(results)
}

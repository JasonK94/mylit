# =============================================================================
# Pseudobulk Differential Expression Analysis - Refactored Functions
# =============================================================================
#
# 네이밍 컨벤션 (Naming Convention):
# =============================================================================
# - pb_prepare_*    : 분석 준비 함수
#                     예: pb_prepare_edgeR - edgeR 분석을 위한 데이터 준비
#
# - pb_matrix_*    : Pseudobulk 매트릭스 생성 및 조작
#                     예: pb_matrix_create - 유연한 집계 옵션을 가진 매트릭스 생성
#
# - pb_deg_*        : Differential Expression Analysis (메서드별 분리)
#                     예: pb_deg_edgeR, pb_deg_DESeq2, pb_deg_wilcox, pb_deg_ttest, 
#                         pb_deg_roc, pb_deg_linear
#
# - pb_deg_markers*: Seurat 스타일의 마커 찾기 함수
#                     예: pb_deg_markers - FindMarkers 스타일
#                         pb_deg_all_markers - FindAllMarkers 스타일
#
# - pb_deg_cluster* : 클러스터별 DEG 분석
#                     예: pb_deg_cluster - 클러스터 마커 찾기
#
# - pb_post_hoc_*   : 사후 분석 함수
#                     예: pb_post_hoc_slope - 그룹 간 기울기 비교
#
# - pb_utils_*      : 공통 유틸리티 함수 (내부 사용)
#                     예: pb_utils_map_sample_group, pb_utils_create_design, 
#                         pb_utils_check_sample_distribution
#
# 주요 변경사항:
# =============================================================================
# 1. prepare_pseudobulk_edgeR → pb_prepare_edgeR
# 2. run_pseudobulk_deg → pb_deg_edgeR (analysis_level 파라미터로 통합)
# 3. pseudobulk_linear_fit → pb_deg_linear
# 4. post_hoc_slope_comparison → pb_post_hoc_slope
# 5. cluster_pseudobulk_deg → pb_deg_cluster (향후 구현 예정)
# 6. pseudobulk_deg → pb_deg_* (메서드별로 분리: pb_deg_edgeR, pb_deg_DESeq2 등)
# 7. FindMarkers_pseudobulk → pb_deg_markers
# 8. FindAllMarkers_pseudobulk → pb_deg_all_markers
# 9. create_pseudobulk_matrix_advanced → pb_matrix_create (중복 제거)
# 10. 공통 로직을 pb_utils_* 함수로 모듈화
#
# 사용 예시:
# =============================================================================
# # edgeR 준비 및 분석
# prep_data <- pb_prepare_edgeR(seurat_obj, sample_col="sample", 
#                                cluster_col="cluster", group_col="group")
# results <- pb_deg_edgeR(prep_data, contrast=c(-1, 1), analysis_level="overall")
#
# # DESeq2 분석
# pb_data <- pb_matrix_create(seurat_obj, sample_col="sample", group.by="cluster")
# results <- pb_deg_DESeq2(pb_data$matrix, pb_data$metadata, group_col="group",
#                          ident.1="A", ident.2="B")
#
# # 선형 회귀 분석
# results <- pb_deg_linear(sobj, genes=c("GENE1", "GENE2"), 
#                          sample_col="sample", numeric_predictor="score")
# slope_comparison <- pb_post_hoc_slope(results)
# =============================================================================

# =============================================================================
# 공통 유틸리티 함수
# =============================================================================

#' Create sample-group mapping from Seurat metadata
#' @keywords internal
#' @importFrom dplyr select mutate group_by summarise filter distinct
#' @importFrom rlang sym
pb_utils_map_sample_group <- function(seurat_obj, sample_col, group_col, verbose = TRUE) {
  original_meta_subset <- seurat_obj@meta.data %>%
    dplyr::select(dplyr::all_of(c(sample_col, group_col))) %>%
    dplyr::mutate(!!rlang::sym(sample_col) := trimws(as.character(.data[[sample_col]])))
  
  # Apply g-prefixing logic for numeric-like sample IDs
  needs_g_prefix_flags <- grepl("^[0-9]+$", original_meta_subset[[sample_col]]) &
    !startsWith(original_meta_subset[[sample_col]], "g")
  
  if (any(needs_g_prefix_flags) && verbose) {
    message("Debug: Applying 'g' prefix to numeric-like values in '", sample_col, "'")
  }
  original_meta_subset[[sample_col]][needs_g_prefix_flags] <- paste0("g", original_meta_subset[[sample_col]][needs_g_prefix_flags])
  
  # Check for samples mapped to multiple groups
  sample_to_multiple_groups_check <- original_meta_subset %>%
    dplyr::group_by(!!rlang::sym(sample_col)) %>%
    dplyr::summarise(n_distinct_groups = dplyr::n_distinct(!!rlang::sym(group_col)), .groups = "drop") %>%
    dplyr::filter(.data$n_distinct_groups > 1)
  
  if (nrow(sample_to_multiple_groups_check) > 0 && verbose) {
    warning("Multiple distinct groups found for the same sample ID. ",
            "Affected sample(s): ", paste(utils::head(sample_to_multiple_groups_check[[sample_col]], 5), collapse=", "))
  }
  
  # Create unique mapping
  sample_group_map <- original_meta_subset %>%
    dplyr::select(!!rlang::sym(sample_col), !!rlang::sym(group_col)) %>%
    dplyr::distinct(!!rlang::sym(sample_col), .keep_all = TRUE)
  
  return(sample_group_map)
}

#' Create design matrix from formula and metadata
#' @keywords internal
#' @importFrom stats model.matrix
pb_utils_create_design <- function(formula_str, data, group_col = NULL) {
  if (is.null(formula_str)) {
    if (is.null(group_col)) stop("Either formula_str or group_col must be provided")
    formula_str <- paste("~", group_col)
  }
  
  design <- tryCatch({
    stats::model.matrix(as.formula(formula_str), data = data)
  }, error = function(e) {
    stop("Design matrix creation failed: ", e$message,
         "\n   Formula: ", formula_str,
         "\n   Metadata head:\n", paste(utils::capture.output(head(data)), collapse="\n"))
  })
  
  # Clean up column names
  colnames(design) <- make.names(colnames(design))
  return(design)
}

#' Check sample distribution per group
#' @keywords internal
#' @importFrom dplyr group_by summarise filter n_distinct
#' @importFrom rlang sym
pb_utils_check_sample_distribution <- function(meta, grouping_col, patient_col, min_reps = 2, verbose = TRUE) {
  sample_counts <- meta %>%
    dplyr::group_by(!!rlang::sym(grouping_col)) %>%
    dplyr::summarise(
      n_pseudo_bulks = dplyr::n(),
      n_unique_samples = dplyr::n_distinct(!!rlang::sym(patient_col)),
      .groups = "drop"
    )
  
  if (verbose) {
    message("    Sample distribution:")
    print(sample_counts)
  }
  
  if (nrow(sample_counts) < 2) {
    return(list(valid = FALSE, reason = "Less than two groups found"))
  }
  
  if (any(sample_counts$n_unique_samples < min_reps, na.rm = TRUE)) {
    problematic <- sample_counts %>%
      dplyr::filter(is.na(.data$n_unique_samples) | .data$n_unique_samples < min_reps)
    return(list(valid = FALSE, reason = "Insufficient samples in one or more groups", details = problematic))
  }
  
  return(list(valid = TRUE, counts = sample_counts))
}

# =============================================================================
# Pseudobulk Matrix Creation
# =============================================================================

#' Create pseudobulk expression matrix with flexible aggregation
#'
#' @param seurat_obj Seurat object
#' @param sample_col Column name defining samples
#' @param group.by Column name defining groups for comparison
#' @param aggregate.by Additional column for finer aggregation (optional)
#' @param covariates Additional metadata columns to include (optional)
#' @param assay Assay to use
#' @param slot Slot to use (typically "counts")
#' @param min.cells Minimum cells per pseudobulk sample
#' @return List with pseudobulk matrix and metadata
#'
#' @importFrom Seurat GetAssayData DefaultAssay
#' @importFrom dplyr distinct select all_of
#' @importFrom Matrix rowSums
#' @export
pb_matrix_create <- function(seurat_obj,
                             sample_col,
                             group.by,
                             aggregate.by = NULL,
                             covariates = NULL,
                             assay = NULL,
                             slot = "counts",
                             min.cells = 10) {
  
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(seurat_obj)
  }
  
  # Get expression matrix
  expr_matrix <- Seurat::GetAssayData(seurat_obj, assay = assay, slot = slot)
  
  # Prepare metadata columns
  meta_cols <- c(sample_col, group.by)
  if (!is.null(aggregate.by)) meta_cols <- c(meta_cols, aggregate.by)
  if (!is.null(covariates)) meta_cols <- c(meta_cols, covariates)
  
  # Get metadata
  metadata <- seurat_obj@meta.data[, meta_cols, drop = FALSE]
  metadata$cell_id <- colnames(seurat_obj)
  
  # Create aggregation ID
  if (!is.null(aggregate.by)) {
    metadata$pb_sample_id <- paste(metadata[[sample_col]], 
                                   metadata[[group.by]], 
                                   metadata[[aggregate.by]], 
                                   sep = "_")
  } else {
    metadata$pb_sample_id <- paste(metadata[[sample_col]], 
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
        sample = sample_info[[sample_col]],
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

# =============================================================================
# Preparation Functions
# =============================================================================

#' Prepare pseudobulk data for edgeR analysis
#'
#' @param seurat_obj Seurat object
#' @param assay Assay name (default: "SCT")
#' @param slot Slot name (default: "counts")
#' @param sample_col Sample identifier column
#' @param cluster_col Cluster identifier column
#' @param group_col Group comparison column
#' @param min_count Minimum count for filterByExpr (default: 10)
#' @param norm_method Normalization method (default: "TMM")
#' @param design_formula Design formula (optional)
#' @param verbose Verbose output (default: TRUE)
#' @return List with pb (matrix), meta (metadata), dge (DGEList), design, contrast_levels
#'
#' @importFrom Seurat AggregateExpression
#' @importFrom edgeR DGEList filterByExpr calcNormFactors
#' @importFrom dplyr %>% mutate select filter distinct left_join
#' @importFrom tidyr separate
#' @importFrom tibble tibble column_to_rownames
#' @export
pb_prepare_edgeR <- function(seurat_obj,
                             assay = "SCT",
                             slot = "counts",
                             sample_col,
                             cluster_col,
                             group_col,
                             min_count = 10,
                             norm_method = "TMM",
                             design_formula = NULL,
                             verbose = TRUE) {
  
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat package required")
  if (!requireNamespace("edgeR", quietly = TRUE)) stop("edgeR package required")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr package required")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("tidyr package required")
  if (!requireNamespace("tibble", quietly = TRUE)) stop("tibble package required")
  
  # Check input columns
  meta_cols <- colnames(seurat_obj@meta.data)
  required_cols <- c(sample_col, cluster_col, group_col)
  if (!all(required_cols %in% meta_cols)) {
    stop("Missing columns in metadata: ",
         paste(required_cols[!required_cols %in% meta_cols], collapse = ", "))
  }
  
  # Ensure group_col is a factor
  seurat_obj@meta.data[[group_col]] <- factor(seurat_obj@meta.data[[group_col]])
  
  # Aggregate Expression
  if (verbose) message("1. Aggregating expression to pseudo-bulk...")
  pb <- Seurat::AggregateExpression(
    seurat_obj,
    assays = assay,
    slot = slot,
    group.by = c(sample_col, cluster_col),
    return.seurat = FALSE
  )[[assay]]
  
  if (nrow(pb) == 0 || ncol(pb) == 0) {
    stop("AggregateExpression resulted in empty matrix")
  }
  
  # Create metadata for pseudo-bulk samples
  if (verbose) message("2. Creating metadata for pseudo-bulk samples...")
  col_example <- colnames(pb)[1]
  sep <- ifelse(grepl("_", col_example), "_", "-")
  
  meta_pb <- tryCatch({
    tibble::tibble(pb_sample_id = colnames(pb)) %>%
      tidyr::separate(pb_sample_id, into = c("patient", "ctype"), sep = sep, extra = "merge", remove = FALSE)
  }, error = function(e) {
    stop("Cannot separate pseudo-bulk column names: ", e$message)
  })
  
  # Map group information
  sample_group_map <- pb_utils_map_sample_group(seurat_obj, sample_col, group_col, verbose)
  
  join_key_map <- stats::setNames(sample_col, "patient")
  meta_pb <- meta_pb %>%
    dplyr::left_join(sample_group_map, by = join_key_map)
  
  if (sum(is.na(meta_pb[[group_col]])) > 0 && verbose) {
    message("Warning: Group mapping failed for some samples")
    if ("patient" %in% colnames(meta_pb)) {
      failed_patients <- meta_pb %>% 
        dplyr::filter(is.na(!!rlang::sym(group_col))) %>% 
        dplyr::pull(.data$patient) %>% 
        unique()
      message("Failed patient IDs: ", paste(head(failed_patients, 10), collapse=", "))
    }
  }
  
  # Ensure factor
  if (!is.factor(meta_pb[[group_col]])) {
    meta_pb[[group_col]] <- factor(meta_pb[[group_col]])
  }
  
  # Match row order with column order
  if ("pb_sample_id" %in% colnames(meta_pb)) {
    meta_pb <- meta_pb[match(colnames(pb), meta_pb$pb_sample_id), ]
    meta_pb <- tibble::column_to_rownames(meta_pb, var = "pb_sample_id")
  } else {
    # Fallback: match by first column if pb_sample_id not found
    meta_pb <- meta_pb[match(colnames(pb), meta_pb[[1]]), ]
  }
  
  # Prepare for edgeR
  if (verbose) message("3. Preparing DGEList and design matrix...")
  dge <- edgeR::DGEList(counts = pb, group = meta_pb[[group_col]], samples = meta_pb)
  
  # Filtering
  keep <- edgeR::filterByExpr(dge, group = meta_pb[[group_col]], min.count = min_count)
  if (verbose) message("   - Filtering: Kept ", sum(keep), " out of ", nrow(dge), " genes.")
  if (sum(keep) == 0) {
    warning("No genes left after filterByExpr")
  }
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  
  # Normalization
  dge <- edgeR::calcNormFactors(dge, method = norm_method)
  
  # Design Matrix
  contrast_levels <- levels(meta_pb[[group_col]])
  formula_str <- if (is.null(design_formula)) {
    paste("~", group_col)
  } else {
    gsub("group_col", group_col, deparse(design_formula))
  }
  
  design <- pb_utils_create_design(formula_str, meta_pb)
  
  if (verbose) message("4. Preparation complete.")
  
  return(list(
    pb = pb,
    meta = meta_pb,
    dge = dge,
    design = design,
    contrast_levels = contrast_levels
  ))
}

# =============================================================================
# Core DEG Analysis Functions
# =============================================================================

#' Perform edgeR-based DEG analysis on prepared pseudobulk data
#'
#' @param prepared_data Output from pb_prepare_edgeR
#' @param contrast Contrast vector for comparison
#' @param analysis_level Analysis level: "overall", "per_cluster", or "specific_cluster"
#' @param target_cluster Cluster name for specific_cluster analysis
#' @param min_samples_per_group Minimum samples per group per cluster
#' @param verbose Verbose output
#' @return Data frame with DEG results
#'
#' @importFrom edgeR estimateDisp glmQLFit glmQLFTest topTags filterByExpr calcNormFactors
#' @importFrom dplyr %>% filter bind_rows mutate rename
#' @importFrom tibble tibble rownames_to_column
#' @export
pb_deg_edgeR <- function(prepared_data,
                        contrast,
                        analysis_level = c("overall", "per_cluster", "specific_cluster"),
                        target_cluster = NULL,
                        min_samples_per_group = 2,
                        min_count = 10,
                        verbose = TRUE) {
  
  analysis_level <- match.arg(analysis_level)
  
  if (!requireNamespace("edgeR", quietly = TRUE)) stop("edgeR package required")
  
  pb <- prepared_data$pb
  meta <- prepared_data$meta
  dge <- prepared_data$dge
  design <- prepared_data$design
  contrast_levels <- prepared_data$contrast_levels
  group_col <- setdiff(colnames(meta), c("pb_sample_id", "patient", "ctype"))[1]
  
  if (missing(contrast)) stop("contrast argument is required")
  
  results <- NULL
  
  if (analysis_level == "overall") {
    if (verbose) message("Performing 'overall' DEG analysis...")
    dge_overall <- edgeR::estimateDisp(dge, design)
    fit_overall <- edgeR::glmQLFit(dge_overall, design)
    qlf_overall <- edgeR::glmQLFTest(fit_overall, contrast = contrast)
    results <- edgeR::topTags(qlf_overall, n = Inf)$table %>%
      rownames_to_column("gene") %>%
      tibble::tibble()
    
  } else if (analysis_level %in% c("per_cluster", "specific_cluster")) {
    all_meta_ctypes <- trimws(as.character(meta$ctype))
    unique_ctypes <- unique(all_meta_ctypes)
    
    if (analysis_level == "specific_cluster") {
      if (is.null(target_cluster)) {
        stop("target_cluster required for specific_cluster analysis")
      }
      target_normalized <- gsub("_", "-", trimws(as.character(target_cluster)))
      if (!target_normalized %in% unique_ctypes) {
        stop("target_cluster '", target_normalized, "' not found in data")
      }
      clusters_to_analyze <- target_normalized
    } else {
      clusters_to_analyze <- unique_ctypes
    }
    
    res_list <- list()
    for (cl_name in clusters_to_analyze) {
      if (verbose) message("  Processing cluster: ", cl_name)
      
      ix <- meta$ctype == cl_name
      if (sum(ix) == 0) {
        warning("No samples for cluster '", cl_name, "'. Skipping.")
        next
      }
      
      pb_sub <- pb[, ix, drop = FALSE]
      md_sub <- meta[ix, , drop = FALSE]
      
      sample_counts <- md_sub %>% dplyr::count(!!rlang::sym(group_col))
      if (any(sample_counts$n < min_samples_per_group)) {
        warning("Cluster '", cl_name, "' has insufficient samples. Skipping.")
        next
      }
      
      md_sub[[group_col]] <- factor(md_sub[[group_col]], levels = contrast_levels)
      
      dge_sub <- edgeR::DGEList(counts = pb_sub, group = md_sub[[group_col]], samples = md_sub)
      
      formula_sub_str <- paste("~", group_col)
      design_sub <- pb_utils_create_design(formula_sub_str, md_sub)
      
      keep_sub <- edgeR::filterByExpr(dge_sub, design = design_sub, group = md_sub[[group_col]], min.count = min_count)
      if (verbose) message("    - Filtering for '", cl_name, "': Kept ", sum(keep_sub), " genes.")
      if (sum(keep_sub) == 0) {
        warning("No genes left for cluster '", cl_name, "'. Skipping.")
        next
      }
      
      dge_sub <- dge_sub[keep_sub, , keep.lib.sizes = FALSE]
      dge_sub <- edgeR::calcNormFactors(dge_sub)
      
      dge_sub_disp <- tryCatch({
        edgeR::estimateDisp(dge_sub, design_sub)
      }, error = function(e) {
        warning("Dispersion estimation failed for cluster ", cl_name, ": ", e$message)
        return(NULL)
      })
      
      if (is.null(dge_sub_disp)) next
      
      fit_sub <- tryCatch({
        edgeR::glmQLFit(dge_sub_disp, design_sub)
      }, error = function(e) {
        warning("Model fitting failed for cluster ", cl_name, ": ", e$message)
        return(NULL)
      })
      
      if (is.null(fit_sub)) next
      
      qlf_sub <- tryCatch({
        edgeR::glmQLFTest(fit_sub, contrast = contrast)
      }, error = function(e) {
        warning("Test failed for cluster ", cl_name, ": ", e$message)
        return(NULL)
      })
      
      if (is.null(qlf_sub)) next
      
      res_sub <- edgeR::topTags(qlf_sub, n = Inf)$table %>%
        rownames_to_column("gene") %>%
        dplyr::mutate(cluster = cl_name, .before = 1) %>%
        tibble::tibble()
      
      res_list[[cl_name]] <- res_sub
    }
    
    if (length(res_list) > 0) {
      results <- dplyr::bind_rows(res_list)
    } else {
      results <- tibble::tibble()
    }
  }
  
  return(results)
}

#' Perform DESeq2-based DEG analysis
#'
#' @param pb_matrix Pseudobulk count matrix
#' @param metadata Sample metadata
#' @param group_col Group comparison column
#' @param ident.1 First group identifier
#' @param ident.2 Second group identifier (optional)
#' @param covariates Additional covariates (optional)
#' @param design.formula Custom design formula (optional)
#' @param logfc.threshold Log fold change threshold
#' @return Data frame with DEG results
#'
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results counts
#' @export
pb_deg_DESeq2 <- function(pb_matrix,
                          metadata,
                          group_col,
                          ident.1,
                          ident.2 = NULL,
                          covariates = NULL,
                          design.formula = NULL,
                          logfc.threshold = 0.25) {
  
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("DESeq2 package required")
  }
  
  # Prepare metadata
  metadata$group <- factor(metadata[[group_col]])
  
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
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(pb_matrix),
    colData = metadata,
    design = design
  )
  
  # Filter low counts
  keep <- rowSums(DESeq2::counts(dds)) >= 10
  dds <- dds[keep, ]
  
  # Run DESeq2
  dds <- DESeq2::DESeq(dds)
  
  # Get results
  res <- DESeq2::results(dds, contrast = c("condition", "target", "reference"))
  
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

#' Perform Wilcoxon rank-sum test on pseudobulk data
#'
#' @param pb_matrix Pseudobulk matrix (normalized)
#' @param metadata Sample metadata
#' @param group_col Group comparison column
#' @param ident.1 First group identifier
#' @param ident.2 Second group identifier
#' @param p_adjust_method P-value adjustment method
#' @return Data frame with DEG results
#'
#' @importFrom stats wilcox.test p.adjust
#' @export
pb_deg_wilcox <- function(pb_matrix,
                          metadata,
                          group_col,
                          ident.1,
                          ident.2 = NULL,
                          p_adjust_method = "BH") {
  
  # Normalize if needed
  if (requireNamespace("edgeR", quietly = TRUE)) {
    pb_norm <- edgeR::cpm(pb_matrix, log = TRUE, prior.count = 1)
  } else {
    pb_norm <- log2(sweep(pb_matrix, 2, colSums(pb_matrix)/1e6, FUN="*") + 1)
  }
  
  # Get groups
  if (is.null(ident.2)) {
    samples1 <- rownames(metadata)[metadata[[group_col]] == ident.1]
    samples2 <- rownames(metadata)[metadata[[group_col]] != ident.1]
  } else {
    samples1 <- rownames(metadata)[metadata[[group_col]] == ident.1]
    samples2 <- rownames(metadata)[metadata[[group_col]] == ident.2]
  }
  
  # Perform tests
  results_list <- lapply(rownames(pb_norm), function(g) {
    vec1 <- pb_norm[g, samples1]
    vec2 <- pb_norm[g, samples2]
    
    p_val <- tryCatch(stats::wilcox.test(vec1, vec2)$p.value, error = function(e) NA)
    avg_lfc <- mean(vec1, na.rm=TRUE) - mean(vec2, na.rm=TRUE)
    
    data.frame(
      gene = g,
      avg_log2FC = avg_lfc,
      p_val = p_val,
      stringsAsFactors = FALSE
    )
  })
  
  res_df <- do.call(rbind, results_list)
  res_df$p_val_adj <- stats::p.adjust(res_df$p_val, method = p_adjust_method)
  
  return(res_df)
}

# =============================================================================
# Linear Regression Analysis
# =============================================================================

#' Perform linear regression on pseudobulk data
#'
#' @param sobj Seurat object
#' @param genes Genes to analyze
#' @param sample_col Sample identifier column
#' @param numeric_predictor Numeric predictor variable
#' @param group_col Group column for stratified analysis (optional)
#' @param p_adjust_method P-value adjustment method
#' @param min_samples_per_group Minimum samples per group
#' @param min_distinct_predictor_values Minimum distinct predictor values
#' @return Data frame with regression results
#'
#' @importFrom Seurat GetAssayData DefaultAssay
#' @importFrom dplyr group_by summarise across all_of select distinct left_join filter
#' @importFrom stats lm coef summary.lm p.adjust as.formula na.omit
#' @export
pb_deg_linear <- function(sobj,
                          genes,
                          sample_col = "sample",
                          numeric_predictor = "severity_score",
                          group_col = NULL,
                          p_adjust_method = "BH",
                          min_samples_per_group = 3,
                          min_distinct_predictor_values = 2) {
  
  if (!inherits(sobj, "Seurat")) stop("sobj must be a Seurat object")
  if (!is.character(genes) || length(genes) == 0) stop("genes must be a character vector")
  if (!all(genes %in% rownames(sobj))) {
    missing_genes <- genes[!genes %in% rownames(sobj)]
    stop("Missing genes: ", paste(missing_genes, collapse=", "))
  }
  
  meta_cols <- colnames(sobj@meta.data)
  if (!sample_col %in% meta_cols) stop("sample_col not found in metadata")
  if (!numeric_predictor %in% meta_cols) stop("numeric_predictor not found in metadata")
  if (!is.null(group_col) && !group_col %in% meta_cols) stop("group_col not found in metadata")
  
  # Convert predictor to numeric if needed
  predictor_data <- sobj@meta.data[[numeric_predictor]]
  if (!is.numeric(predictor_data)) {
    predictor_data <- as.numeric(as.character(predictor_data))
    if (all(is.na(predictor_data))) stop("Cannot convert numeric_predictor to numeric")
    sobj@meta.data[[numeric_predictor]] <- predictor_data
  }
  
  # Get expression data
  assay <- Seurat::DefaultAssay(sobj)
  expr_data <- Seurat::GetAssayData(sobj, assay = assay, slot = "data")[genes, , drop = FALSE]
  
  # Create pseudobulk by averaging per sample
  expr_df_transposed <- as.data.frame(t(as.matrix(expr_data)))
  expr_df_transposed[[sample_col]] <- sobj@meta.data[rownames(expr_df_transposed), sample_col]
  
  avg_expr_all_samples <- expr_df_transposed %>%
    dplyr::group_by(!!rlang::sym(sample_col)) %>%
    dplyr::summarise(dplyr::across(dplyr::all_of(genes), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  
  # Merge with metadata
  cols_to_select <- c(sample_col, numeric_predictor)
  if (!is.null(group_col)) cols_to_select <- c(cols_to_select, group_col)
  
  meta_data_full <- sobj@meta.data %>%
    dplyr::select(dplyr::all_of(unique(cols_to_select))) %>%
    dplyr::distinct()
  
  df_merged_full <- dplyr::left_join(avg_expr_all_samples, meta_data_full, by = sample_col)
  
  # Run regression per gene (and per group if specified)
  valid_groups <- if (is.null(group_col)) {
    "all_samples"
  } else {
    unique(na.omit(df_merged_full[[group_col]]))
  }
  
  final_results_list <- list()
  
  for (grp_val in valid_groups) {
    current_data <- df_merged_full
    if (!is.null(group_col) && grp_val != "all_samples") {
      current_data <- df_merged_full %>% dplyr::filter(!!rlang::sym(group_col) == grp_val)
    }
    
    group_results_list <- lapply(genes, function(gene) {
      formula_str <- paste0("`", gene, "` ~ `", numeric_predictor, "`")
      model_data <- current_data[, c(gene, numeric_predictor, sample_col), drop = FALSE]
      model_data_complete <- stats::na.omit(model_data)
      
      n_samples <- dplyr::n_distinct(model_data_complete[[sample_col]])
      n_distinct_preds <- dplyr::n_distinct(model_data_complete[[numeric_predictor]])
      
      if (n_samples < min_samples_per_group || n_distinct_preds < min_distinct_predictor_values) {
        return(data.frame(
          gene = gene,
          intercept = NA_real_,
          slope = NA_real_,
          slope_se = NA_real_,
          n_samples_in_model = n_samples,
          p_value = NA_real_,
          r_squared = NA_real_,
          stringsAsFactors = FALSE
        ))
      }
      
      model <- stats::lm(as.formula(formula_str), data = model_data_complete)
      summary_model <- summary(model)
      coefs <- coef(summary_model)
      
      intercept_val <- if ("(Intercept)" %in% rownames(coefs)) coefs["(Intercept)", "Estimate"] else NA_real_
      slope_val <- if (nrow(coefs) >= 2) coefs[2, "Estimate"] else NA_real_
      slope_se_val <- if (nrow(coefs) >= 2) coefs[2, "Std. Error"] else NA_real_
      p_val <- if (nrow(coefs) >= 2) coefs[2, "Pr(>|t|)"] else NA_real_
      
      data.frame(
        gene = gene,
        intercept = intercept_val,
        slope = slope_val,
        slope_se = slope_se_val,
        n_samples_in_model = n_samples,
        p_value = p_val,
        r_squared = summary_model$r.squared,
        stringsAsFactors = FALSE
      )
    })
    
    group_results_df <- do.call(rbind, group_results_list)
    group_results_df$group <- grp_val
    final_results_list[[grp_val]] <- group_results_df
  }
  
  final_results_df <- do.call(rbind, final_results_list)
  
  # Adjust p-values
  valid_p_indices <- !is.na(final_results_df$p_value)
  if (any(valid_p_indices)) {
    final_results_df$adj_p_value <- NA_real_
    final_results_df$adj_p_value[valid_p_indices] <- stats::p.adjust(
      final_results_df$p_value[valid_p_indices], method = p_adjust_method
    )
  } else {
    final_results_df$adj_p_value <- NA_real_
  }
  
  cols_ordered <- c("group", "gene", "intercept", "slope", "slope_se", 
                    "n_samples_in_model", "p_value", "adj_p_value", "r_squared")
  final_results_df <- final_results_df[, intersect(cols_ordered, names(final_results_df)), drop = FALSE]
  
  return(final_results_df)
}

# =============================================================================
# Post-hoc Analysis
# =============================================================================

#' Post-hoc slope comparison between groups
#'
#' @param results_df Output from pb_deg_linear
#' @param gene_col Gene column name
#' @param group_col Group column name
#' @param slope_col Slope column name
#' @param se_col Standard error column name
#' @param n_samples_col Sample count column name
#' @param p_adjust_method P-value adjustment method
#' @param adjustment_scope Adjustment scope: "global" or "per_gene"
#' @return Data frame with comparison results
#'
#' @importFrom dplyr group_by do filter arrange mutate ungroup slice
#' @importFrom utils combn
#' @importFrom stats pt p.adjust
#' @export
pb_post_hoc_slope <- function(results_df,
                              gene_col = "gene",
                              group_col = "group",
                              slope_col = "slope",
                              se_col = "slope_se",
                              n_samples_col = "n_samples_in_model",
                              p_adjust_method = "BH",
                              adjustment_scope = "global") {
  
  required_cols <- c(gene_col, group_col, slope_col, se_col, n_samples_col)
  if (!all(required_cols %in% names(results_df))) {
    missing_cols <- required_cols[!required_cols %in% names(results_df)]
    stop("Missing columns: ", paste(missing_cols, collapse=", "))
  }
  
  if (!adjustment_scope %in% c("global", "per_gene")) {
    stop("adjustment_scope must be 'global' or 'per_gene'")
  }
  
  # Remove duplicates
  results_df_deduplicated <- results_df %>%
    dplyr::group_by(!!rlang::sym(gene_col), !!rlang::sym(group_col)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
  
  # Filter valid data
  results_df_complete <- results_df_deduplicated %>%
    dplyr::filter(
      !is.na(.data[[slope_col]]),
      !is.na(.data[[se_col]]),
      .data[[se_col]] > 0,
      !is.na(.data[[n_samples_col]]),
      .data[[n_samples_col]] >= 3
    )
  
  if (nrow(results_df_complete) == 0) {
    warning("No valid data for comparison")
    return(data.frame())
  }
  
  # Perform comparisons
  comparison_results <- results_df_complete %>%
    dplyr::group_by(!!rlang::sym(gene_col)) %>%
    dplyr::do({
      gene_data <- .
      if (nrow(gene_data) < 2) return(data.frame())
      
      group_combinations <- utils::combn(gene_data[[group_col]], 2)
      
      pair_results_list <- lapply(seq_len(ncol(group_combinations)), function(i) {
        g1 <- group_combinations[1, i]
        g2 <- group_combinations[2, i]
        
        d1 <- gene_data[gene_data[[group_col]] == g1, ]
        d2 <- gene_data[gene_data[[group_col]] == g2, ]
        
        b1 <- d1[[slope_col]]; se1 <- d1[[se_col]]; n1 <- d1[[n_samples_col]]
        b2 <- d2[[slope_col]]; se2 <- d2[[se_col]]; n2 <- d2[[n_samples_col]]
        
        df_reg1 <- n1 - 2
        df_reg2 <- n2 - 2
        
        t_stat <- (b1 - b2) / sqrt(se1^2 + se2^2)
        numerator_ws <- (se1^2 + se2^2)^2
        denominator_ws <- (se1^4 / df_reg1) + (se2^4 / df_reg2)
        df_ws <- numerator_ws / denominator_ws
        p_val <- 2 * stats::pt(abs(t_stat), df = df_ws, lower.tail = FALSE)
        
        data.frame(
          group1 = g1, group2 = g2,
          slope1 = b1, slope2 = b2,
          se1 = se1, se2 = se2,
          n1 = n1, n2 = n2,
          t_statistic = t_stat, df = df_ws, p_value_raw = p_val,
          stringsAsFactors = FALSE
        )
      })
      
      if (length(pair_results_list) > 0) {
        do.call(rbind, pair_results_list)
      } else {
        data.frame()
      }
    }) %>%
    dplyr::ungroup()
  
  # Adjust p-values
  if (adjustment_scope == "per_gene") {
    comparison_results <- comparison_results %>%
      dplyr::group_by(!!rlang::sym(gene_col)) %>%
      dplyr::mutate(adj_p_value = stats::p.adjust(.data$p_value_raw, method = p_adjust_method)) %>%
      dplyr::ungroup()
  } else {
    comparison_results <- comparison_results %>%
      dplyr::mutate(adj_p_value = stats::p.adjust(.data$p_value_raw, method = p_adjust_method))
  }
  
  comparison_results <- comparison_results %>%
    dplyr::arrange(!!rlang::sym(gene_col), 
                   desc(is.na(.data$adj_p_value)), 
                   .data$adj_p_value)
  
  return(comparison_results)
}

# =============================================================================
# High-level Interface Functions
# =============================================================================

#' Find markers using pseudobulk DEG (similar to Seurat FindMarkers)
#'
#' @param sobj Seurat object
#' @param ident.1 First identity
#' @param ident.2 Second identity (optional)
#' @param group.by Group column name
#' @param sample.by Sample column name
#' @param method DEG method: "edgeR" or "DESeq2"
#' @param covariates Additional covariates
#' @param aggregate.by Additional aggregation column
#' @param min.cells Minimum cells per pseudobulk
#' @param logfc.threshold Log fold change threshold
#' @param min.pct Minimum percentage threshold
#' @param assay Assay to use
#' @param slot Slot to use
#' @return Data frame with marker genes
#'
#' @importFrom Seurat Idents WhichCells subset DefaultAssay
#' @export
pb_deg_markers <- function(sobj,
                           ident.1,
                           ident.2 = NULL,
                           group.by = "seurat_clusters",
                           sample.by,
                           method = "DESeq2",
                           covariates = NULL,
                           aggregate.by = NULL,
                           min.cells = 10,
                           logfc.threshold = 0.25,
                           min.pct = 0.1,
                           assay = NULL,
                           slot = "counts") {
  
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(sobj)
  }
  
  # Set identities
  Seurat::Idents(sobj) <- group.by
  
  # Get cells for each group
  cells.1 <- Seurat::WhichCells(sobj, idents = ident.1)
  if (is.null(ident.2)) {
    cells.2 <- Seurat::WhichCells(sobj, idents = setdiff(levels(Seurat::Idents(sobj)), ident.1))
  } else {
    cells.2 <- Seurat::WhichCells(sobj, idents = ident.2)
  }
  
  # Subset object
  cells.use <- c(cells.1, cells.2)
  seurat_subset <- Seurat::subset(sobj, cells = cells.use)
  
  # Create pseudobulk matrix
  pb_data <- pb_matrix_create(
    seurat_obj = seurat_subset,
    sample_col = sample.by,
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
    de_results <- pb_deg_DESeq2(
      pb_matrix = pb_matrix,
      metadata = pb_metadata,
      group_col = "group",
      ident.1 = ident.1,
      ident.2 = ident.2,
      covariates = covariates,
      design.formula = NULL,
      logfc.threshold = logfc.threshold
    )
  } else if (method == "edgeR") {
    # For edgeR, use the advanced function from original file if needed
    # or implement similar logic here
    stop("edgeR method for pb_deg_markers not yet implemented. Use pb_prepare_edgeR and pb_deg_edgeR instead.")
  } else {
    stop("Method must be either 'DESeq2' or 'edgeR'")
  }
  
  return(de_results)
}

#' Find all markers using pseudobulk DEG (similar to Seurat FindAllMarkers)
#'
#' @param sobj Seurat object
#' @param group.by Group column name
#' @param sample.by Sample column name
#' @param method DEG method: "DESeq2" or "edgeR"
#' @param covariates Additional covariates
#' @param aggregate.by Additional aggregation column
#' @param min.cells Minimum cells per pseudobulk
#' @param logfc.threshold Log fold change threshold
#' @param min.pct Minimum percentage threshold
#' @param only.pos Only return positive markers
#' @param assay Assay to use
#' @param slot Slot to use
#' @return Data frame with marker genes for all clusters
#'
#' @export
pb_deg_all_markers <- function(sobj,
                                group.by = "seurat_clusters",
                                sample.by,
                                method = "DESeq2",
                                covariates = NULL,
                                aggregate.by = NULL,
                                min.cells = 10,
                                logfc.threshold = 0.25,
                                min.pct = 0.1,
                                only.pos = FALSE,
                                assay = NULL,
                                slot = "counts") {
  
  # Get all unique identities
  Seurat::Idents(sobj) <- group.by
  all_idents <- levels(Seurat::Idents(sobj))
  
  # Run FindMarkers for each identity
  all_markers <- list()
  
  for (ident in all_idents) {
    message(paste0("Finding markers for ", group.by, " ", ident))
    
    tryCatch({
      markers <- pb_deg_markers(
        sobj = sobj,
        ident.1 = ident,
        ident.2 = NULL,
        group.by = group.by,
        sample.by = sample.by,
        method = method,
        covariates = covariates,
        aggregate.by = aggregate.by,
        min.cells = min.cells,
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
  if (only.pos && "avg_log2FC" %in% colnames(all_markers_df)) {
    all_markers_df <- all_markers_df[all_markers_df$avg_log2FC > 0, ]
  }
  
  # Sort by cluster and p-value
  if ("cluster" %in% colnames(all_markers_df) && "p_val_adj" %in% colnames(all_markers_df)) {
    all_markers_df <- all_markers_df[order(all_markers_df$cluster, all_markers_df$p_val_adj), ]
  }
  
  return(all_markers_df)
}


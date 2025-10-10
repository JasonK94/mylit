# ===
# GEOMX SEURAT ANALYSIS SUITE - REFACTORED
# ===
# Purpose: Analyze paired GeoMx spatial transcriptomics data using linear 
#          mixed models to account for patient-level random effects
# 
# Data Structure: AOIs (Area of Illumination) treated as "cells" in Seurat
# Key Variables: 
#   - timepoint: "pre" vs "post" treatment
#   - response: "R" (Responder) vs "NR" (Non-Responder)  
#   - drug: "Infliximab", "Ustekinumab", "Vedolizumab"
# 

# Package Dependencies ---
library(dplyr)
library(tidyr)
library(glue)
library(Seurat)
library(lme4)
library(lmerTest)
library(emmeans)
library(parallel)
library(broom.mixed)
library(ggpubr)
library(pheatmap)
library(viridis)
library(purrr)
library(ggplot2)
library(tibble)
library(openxlsx)


# ===
# SECTION 1: CONFIGURATION & SETUP
# ===

#' Create Analysis Configuration Object
#'
#' Centralizes metadata column names for consistent reference throughout analysis.
#' This configuration object should be created once and passed to all functions.
#'
#' @param patient String. Column name for patient/subject identifier.
#' @param drug String. Column name for drug/treatment type.
#' @param timepoint String. Column name for timepoint (e.g., "pre"/"post").
#' @param ck String. Column name for stratification variable (e.g., CK status).
#' @param response String. Column name for treatment response classification.
#' @param aoi String. Column name for AOI (Area of Interest) unique identifier.
#'
#' @return A list containing standardized column name mappings.
#' @export
#'
#' @examples
#' config <- create_analysis_config(
#'   patient = "patient_id",
#'   drug = "treatment_drug",
#'   timepoint = "visit",
#'   response = "clinical_response"
#' )
create_analysis_config <- function(
    patient = "patient_id",
    drug = "drug",
    timepoint = "timepoint",
    ck = "ck_status",
    response = "response",
    aoi = "aoi_id"
) {
  config <- list(
    patient = patient,
    drug = drug,
    timepoint = timepoint,
    ck = ck,
    response = response,
    aoi = aoi
  )
  
  class(config) <- c("geomx_config", "list")
  return(config)
}

#' Validate Analysis Configuration
#'
#' Checks that all required columns exist in metadata and have valid values.
#'
#' @param metadata Data frame. Seurat metadata or standalone metadata.
#' @param config List. Configuration object from create_analysis_config().
#' @param required_cols Character vector. Which config columns are required.
#'
#' @return Invisibly returns TRUE if valid, otherwise throws informative error.
#' @keywords internal
validate_config <- function(metadata, config, 
                            required_cols = c("patient", "timepoint")) {
  
  # Check required columns exist in config
  missing_config <- setdiff(required_cols, names(config))
  if (length(missing_config) > 0) {
    stop(sprintf(
      "Config missing required fields: %s",
      paste(missing_config, collapse = ", ")
    ))
  }
  
  # Check columns exist in metadata
  config_values <- unlist(config[required_cols])
  missing_cols <- setdiff(config_values, colnames(metadata))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Metadata missing required columns: %s\nAvailable columns: %s",
      paste(missing_cols, collapse = ", "),
      paste(head(colnames(metadata), 10), collapse = ", ")
    ))
  }
  
  invisible(TRUE)
}


# ===
# SECTION 2: DATA LOADING & PREPARATION
# ===

#' Load and Prepare GeoMx Data for Seurat Analysis
#'
#' Reads count matrix and metadata from files (Excel or CSV), creates a Seurat
#' object, and adds derived metadata columns for convenience in downstream analysis.
#'
#' @param count_file String. Path to count data file (.xlsx or .csv). 
#'   If NULL, must provide count_matrix.
#' @param metadata_file String. Path to metadata file (.xlsx or .csv).
#'   If NULL, must provide metadata.
#' @param count_matrix Matrix. Raw count matrix (genes x samples). 
#'   Alternative to count_file.
#' @param metadata Data.frame. Sample metadata. Alternative to metadata_file.
#' @param config List. Configuration from create_analysis_config().
#' @param normalize_method String. One of "none", "log", or "quantile". 
#'   - "none": Assumes data is already normalized (e.g., Q3 normalized)
#'   - "log": Standard Seurat log-normalization
#'   - "quantile": For pre-normalized data, just converts to log scale
#' @param min_cells Integer. Minimum cells expressing a gene to keep it.
#' @param min_features Integer. Minimum features per cell to keep it.
#'
#' @return A Seurat object with the following added metadata columns:
#'   - condition_full: timepoint_drug_response (e.g., "pre_Infliximab_R")
#'   - condition_simple: timepoint_drug (e.g., "pre_Infliximab")
#'   - drug_response: drug_response (e.g., "Infliximab_R")
#'
#' @export
#'
#' @examples
#' sobj <- prepare_geomx_data(
#'   count_file = "counts.xlsx",
#'   metadata_file = "metadata.xlsx",
#'   config = config,
#'   normalize_method = "none"  # Already Q3 normalized
#' )
prepare_geomx_data <- function(
    count_file = NULL,
    metadata_file = NULL,
    count_matrix = NULL,
    metadata = NULL,
    config = create_analysis_config(),
    normalize_method = c("none", "log", "quantile"),
    min_cells = 0,
    min_features = 0
) {
  
  normalize_method <- match.arg(normalize_method)
  
  # ---
  # Load Data
  # ---
  
  # Load count matrix
  if (!is.null(count_file)) {
    message("Loading count matrix from file...")
    if (grepl("\\.xlsx$", count_file, ignore.case = TRUE)) {
      count_matrix <- openxlsx::read.xlsx(
        count_file, 
        sheet = 1, 
        rowNames = TRUE
      )
    } else if (grepl("\\.(csv|txt)$", count_file, ignore.case = TRUE)) {
      count_matrix <- read.csv(count_file, row.names = 1, check.names = FALSE)
    } else {
      stop("count_file must be .xlsx, .csv, or .txt")
    }
  }
  
  if (is.null(count_matrix)) {
    stop("Must provide either count_file or count_matrix")
  }
  
  # Load metadata
  if (!is.null(metadata_file)) {
    message("Loading metadata from file...")
    if (grepl("\\.xlsx$", metadata_file, ignore.case = TRUE)) {
      metadata <- openxlsx::read.xlsx(metadata_file, sheet = 1)
    } else if (grepl("\\.(csv|txt)$", metadata_file, ignore.case = TRUE)) {
      metadata <- read.csv(metadata_file, check.names = FALSE)
    } else {
      stop("metadata_file must be .xlsx, .csv, or .txt")
    }
  }
  
  if (is.null(metadata)) {
    stop("Must provide either metadata_file or metadata")
  }
  
  # ---
  # Validate & Align Data
  # ---
  
  # Set AOI as rownames
  if (!config$aoi %in% colnames(metadata)) {
    stop(sprintf("AOI column '%s' not found in metadata", config$aoi))
  }
  rownames(metadata) <- metadata[[config$aoi]]
  
  # Ensure samples match between count matrix and metadata
  common_samples <- intersect(colnames(count_matrix), rownames(metadata))
  
  if (length(common_samples) == 0) {
    stop("No overlapping samples between count matrix and metadata!")
  }
  
  if (length(common_samples) < ncol(count_matrix)) {
    warning(sprintf(
      "Dropping %d samples not found in metadata",
      ncol(count_matrix) - length(common_samples)
    ))
  }
  
  if (length(common_samples) < nrow(metadata)) {
    warning(sprintf(
      "Dropping %d metadata rows not found in count matrix",
      nrow(metadata) - length(common_samples)
    ))
  }
  
  # Align order
  count_matrix <- count_matrix[, common_samples]
  metadata <- metadata[common_samples, ]
  
  # Validate configuration
  validate_config(metadata, config, required_cols = c("patient", "timepoint"))
  
  # ---
  # Create Seurat Object
  # ---
  
  message(sprintf(
    "Creating Seurat object with %d genes x %d samples",
    nrow(count_matrix), ncol(count_matrix)
  ))
  
  seurat_obj <- CreateSeuratObject(
    counts = as.matrix(count_matrix),
    meta.data = metadata,
    min.cells = min_cells,
    min.features = min_features
  )
  
  # Apply normalization
  if (normalize_method == "log") {
    message("Applying log-normalization...")
    seurat_obj <- NormalizeData(
      seurat_obj, 
      normalization.method = "LogNormalize",
      scale.factor = 10000
    )
  } else if (normalize_method == "quantile") {
    message("Data assumed to be Q3 normalized, converting to log scale...")
    # For Q3-normalized data, the 'data' slot should be log2(Q3norm + 1)
    seurat_obj@assays$RNA@data <- log2(
      seurat_obj@assays$RNA@counts + 1
    )
  } else {
    message("No normalization applied (assuming pre-normalized data)")
    seurat_obj@assays$RNA@data <- seurat_obj@assays$RNA@counts
  }
  
  # ---
  # Add Derived Metadata
  # ---
  
  # Helper to safely paste metadata columns
  safe_paste <- function(..., sep = "_") {
    vals <- list(...)
    # Replace NA with "NA" string
    vals <- lapply(vals, function(x) {
      x <- as.character(x)
      x[is.na(x)] <- "NA"
      return(x)
    })
    do.call(paste, c(vals, list(sep = sep)))
  }
  
  # Create combination variables for easy subsetting
  if (all(c(config$timepoint, config$drug, config$response) %in% 
          colnames(metadata))) {
    seurat_obj$condition_full <- safe_paste(
      seurat_obj@meta.data[[config$timepoint]],
      seurat_obj@meta.data[[config$drug]],
      seurat_obj@meta.data[[config$response]]
    )
  }
  
  if (all(c(config$timepoint, config$drug) %in% colnames(metadata))) {
    seurat_obj$condition_simple <- safe_paste(
      seurat_obj@meta.data[[config$timepoint]],
      seurat_obj@meta.data[[config$drug]]
    )
  }
  
  if (all(c(config$drug, config$response) %in% colnames(metadata))) {
    seurat_obj$drug_response <- safe_paste(
      seurat_obj@meta.data[[config$drug]],
      seurat_obj@meta.data[[config$response]]
    )
  }
  
  message("✅ Seurat object created successfully")
  return(seurat_obj)
}


# ===
# SECTION 3: DATA VALIDATION - PARITY CHECKING
# ===

#' Diagnose Paired Sample Completeness
#'
#' For paired/longitudinal analyses, checks whether each patient (or group) has
#' complete and balanced representation of required conditions. For example,
#' ensures each patient has both "pre" and "post" timepoints with equal numbers
#' of samples, and no other timepoints.
#'
#' @param seurat_obj Seurat object.
#' @param config List. Configuration from create_analysis_config().
#' @param grouping_vars Character vector. Metadata columns defining groups
#'   (e.g., c("patient", "drug")). Each unique combination will be checked.
#' @param check_var String. The metadata column to check for parity
#'   (typically "timepoint").
#' @param required_values Character vector. Values that MUST be present in
#'   check_var for each group (e.g., c("pre", "post")).
#' @param allow_extra Logical. If FALSE, groups with additional values beyond
#'   required_values will fail validation.
#'
#' @return A list with components:
#'   \describe{
#'     \item{summary}{Data frame with one row per group, showing counts and pass/fail status}
#'     \item{passed_samples}{Character vector of sample IDs that passed validation}
#'     \item{failed_groups}{Character vector of group IDs that failed}
#'     \item{message}{Summary message}
#'   }
#'
#' @export
#'
#' @examples
#' # Check if each patient has both pre and post timepoints
#' parity_check <- diagnose_sample_parity(
#'   seurat_obj = sobj,
#'   config = config,
#'   grouping_vars = c("patient", "drug"),
#'   check_var = "timepoint",
#'   required_values = c("pre", "post"),
#'   allow_extra = FALSE
#' )
#' 
#' # Subset to valid samples
#' sobj_paired <- subset(sobj, cells = parity_check$passed_samples)
diagnose_sample_parity <- function(
    seurat_obj,
    config,
    grouping_vars = "patient",
    check_var = "timepoint",
    required_values = c("pre", "post"),
    allow_extra = FALSE
) {
  
  # Validate inputs
  metadata <- seurat_obj@meta.data
  
  # Map config names to actual column names
  grouping_cols <- sapply(grouping_vars, function(v) {
    if (v %in% names(config)) config[[v]] else v
  })
  
  check_col <- if (check_var %in% names(config)) {
    config[[check_var]]
  } else {
    check_var
  }
  
  missing_cols <- setdiff(
    c(grouping_cols, check_col), 
    colnames(metadata)
  )
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Columns not found in metadata: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }
  
  # ---
  # Create Diagnostic Table
  # ---
  
  # Create unique group ID
  meta_df <- metadata %>%
    tibble::rownames_to_column("sample_id") %>%
    mutate(
      group_id = do.call(paste, c(across(all_of(grouping_cols)), sep = "|")),
      check_value = as.character(.data[[check_col]])
    )
  
  # Count samples per group per check_value
  value_counts <- meta_df %>%
    count(group_id, check_value, name = "n_samples")
  
  # Create summary with one row per group
  summary_df <- meta_df %>%
    distinct(group_id, .keep_all = TRUE) %>%
    select(group_id, all_of(grouping_cols)) %>%
    left_join(
      # Total samples per group
      meta_df %>% count(group_id, name = "total_samples"),
      by = "group_id"
    )
  
  # Add columns for each required value
  for (val in required_values) {
    col_name <- paste0("n_", make.names(val))
    summary_df <- summary_df %>%
      left_join(
        value_counts %>%
          filter(check_value == val) %>%
          select(group_id, n_samples) %>%
          rename(!!col_name := n_samples),
        by = "group_id"
      ) %>%
      mutate(!!col_name := replace_na(.data[[col_name]], 0))
  }
  
  # Calculate derived columns
  required_cols <- paste0("n_", make.names(required_values))
  
  summary_df <- summary_df %>%
    rowwise() %>%
    mutate(
      # Total samples with required values
      n_required = sum(c_across(all_of(required_cols))),
      # Samples with other values
      n_other = total_samples - n_required,
      # Has all required values present
      has_all_values = all(c_across(all_of(required_cols)) > 0),
      # All required values have equal counts (balanced)
      is_balanced = length(unique(c_across(all_of(required_cols)))) == 1,
      # Pass criteria
      passes = has_all_values & is_balanced & (allow_extra | n_other == 0)
    ) %>%
    ungroup()
  
  # ---
  # Identify Passing Samples
  # ---
  
  passed_groups <- summary_df %>%
    filter(passes) %>%
    pull(group_id)
  
  passed_samples <- meta_df %>%
    filter(group_id %in% passed_groups) %>%
    pull(sample_id)
  
  failed_groups <- summary_df %>%
    filter(!passes) %>%
    pull(group_id)
  
  # ---
  # Generate Summary Message
  # ---
  
  n_total <- nrow(summary_df)
  n_passed <- length(passed_groups)
  n_failed <- length(failed_groups)
  
  if (n_failed == 0) {
    msg <- sprintf(
      "✅ All %d groups passed validation",
      n_total
    )
  } else {
    msg <- sprintf(
      "⚠️  %d/%d groups failed validation:\n%s",
      n_failed, n_total,
      paste(head(failed_groups, 5), collapse = "\n")
    )
    if (n_failed > 5) {
      msg <- paste0(msg, sprintf("\n... and %d more", n_failed - 5))
    }
  }
  
  message(msg)
  message(sprintf(
    "Passed samples: %d/%d (%.1f%%)",
    length(passed_samples),
    nrow(metadata),
    100 * length(passed_samples) / nrow(metadata)
  ))
  
  # ---
  # Return Results
  # ---
  
  return(list(
    summary = summary_df,
    passed_samples = passed_samples,
    failed_groups = failed_groups,
    message = msg
  ))
}


# ===
# SECTION 4: GENE SCREENING
# ===

#' Quick Gene Screening Using Wilcoxon Tests
#'
#' Performs rapid statistical screening to identify candidate genes for 
#' computationally intensive LMM analysis. Uses Seurat's FindAllMarkers
#' for efficiency.
#'
#' @param seurat_obj Seurat object.
#' @param config List. Configuration from create_analysis_config().
#' @param grouping_var String. Which variable to compare ("timepoint", "drug",
#'   "response", or name of a derived column like "condition_full").
#' @param subset_var String. Optional. Metadata column to subset by.
#' @param subset_value Value to filter subset_var by (e.g., "Tumor" for tissue type).
#' @param min_pct Numeric. Minimum fraction of samples expressing gene in either group.
#' @param logfc_threshold Numeric. Minimum log fold-change threshold.
#' @param top_n Integer. Number of top genes to return.
#' @param test_method String. Statistical test method (default "wilcox").
#'
#' @return A list containing:
#'   \describe{
#'     \item{all_markers}{Full results data frame from FindAllMarkers}
#'     \item{top_genes}{Character vector of top N gene names}
#'     \item{subset_obj}{The Seurat object used (potentially subsetted)}
#'   }
#'
#' @export
#'
#' @examples
#' # Screen for genes changing between timepoints
#' screening <- screen_genes(
#'   seurat_obj = sobj,
#'   config = config,
#'   grouping_var = "timepoint",
#'   top_n = 500
#' )
#' 
#' # Use top genes for LMM analysis
#' lmm_results <- run_lmm_analysis(
#'   seurat_obj = sobj,
#'   genes = screening$top_genes,
#'   config = config
#' )
screen_genes <- function(
    seurat_obj,
    config,
    grouping_var = "timepoint",
    subset_var = NULL,
    subset_value = NULL,
    min_pct = 0.1,
    logfc_threshold = 0.1,
    top_n = 1000,
    test_method = "wilcox"
) {
  
  # ---
  # Subset Data if Requested
  # ---
  
  obj <- seurat_obj
  
  if (!is.null(subset_var) && !is.null(subset_value)) {
    # Map config name to column name
    subset_col <- if (subset_var %in% names(config)) {
      config[[subset_var]]
    } else {
      subset_var
    }
    
    if (!subset_col %in% colnames(obj@meta.data)) {
      stop(sprintf("Subset column '%s' not found in metadata", subset_col))
    }
    
    message(sprintf("Subsetting to %s = %s", subset_col, subset_value))
    obj <- subset(obj, subset = .data[[subset_col]] == subset_value)
    
    if (ncol(obj) == 0) {
      stop("No samples remain after subsetting!")
    }
  }
  
  # ---
  # Set Identity for Comparison
  # ---
  
  # Map config name to column name
  group_col <- if (grouping_var %in% names(config)) {
    config[[grouping_var]]
  } else {
    grouping_var
  }
  
  if (!group_col %in% colnames(obj@meta.data)) {
    stop(sprintf("Grouping column '%s' not found in metadata", group_col))
  }
  
  # Set identity (avoids modifying original object)
  Idents(obj) <- obj@meta.data[[group_col]]
  
  n_groups <- length(unique(Idents(obj)))
  message(sprintf(
    "Screening genes across %d groups in '%s'",
    n_groups, group_col
  ))
  
  # ---
  # Find Markers
  # ---
  
  markers <- FindAllMarkers(
    obj,
    test.use = test_method,
    min.pct = min_pct,
    logfc.threshold = logfc_threshold,
    only.pos = FALSE,
    return.thresh = 1,  # Return all genes
    verbose = FALSE
  )
  
  if (nrow(markers) == 0) {
    warning("No significant markers found with current thresholds!")
    return(list(
      all_markers = markers,
      top_genes = character(0),
      subset_obj = obj
    ))
  }
  
  # ---
  # Select Top Genes
  # ---
  
  # Rank by adjusted p-value, then by absolute log fold-change
  top_genes <- markers %>%
    arrange(p_val_adj, desc(abs(avg_log2FC))) %>%
    head(top_n) %>%
    pull(gene) %>%
    unique()
  
  message(sprintf("Selected %d candidate genes for LMM analysis", length(top_genes)))
  
  return(list(
    all_markers = markers,
    top_genes = top_genes,
    subset_obj = obj
  ))
}


# ===
# SECTION 5: LINEAR MIXED MODEL ANALYSIS
# ===

#' Fit Linear Mixed Model for Single Gene
#'
#' Internal helper function that fits an lmer model for one gene. Handles
#' formula construction, factor releveling, and error catching.
#'
#' @param gene String. Gene name.
#' @param expr_vector Numeric vector. Expression values for the gene.
#' @param metadata Data frame. Sample metadata.
#' @param formula_str String. Complete model formula (e.g., 
#'   "Expression ~ drug * timepoint + (1|patient)").
#' @param factor_cols Character vector. Columns to convert to factors.
#' @param reference_levels Named list. Factor reference levels 
#'   (e.g., list(timepoint = "pre", drug = "Infliximab")).
#'
#' @return List with model results or error information.
#' @keywords internal
fit_single_gene_lmm <- function(
    gene,
    expr_vector,
    metadata,
    formula_str,
    factor_cols = NULL,
    reference_levels = NULL
) {
  
  # Prepare data frame
  df <- data.frame(
    Expression = expr_vector,
    metadata,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  
  # Convert to factors and set reference levels
  if (!is.null(factor_cols)) {
    for (col in factor_cols) {
      if (col %in% colnames(df)) {
        df[[col]] <- factor(df[[col]])
        
        # Set reference level if specified
        if (!is.null(reference_levels) && col %in% names(reference_levels)) {
          ref_level <- reference_levels[[col]]
          if (ref_level %in% levels(df[[col]])) {
            df[[col]] <- relevel(df[[col]], ref = ref_level)
          }
        }
      }
    }
  }
  
  # Check for sufficient variation
  if (length(unique(df$Expression)) < 2) {
    return(list(
      gene = gene,
      converged = FALSE,
      error = "No variation in expression",
      formula = formula_str
    ))
  }
  
  # Fit model
  tryCatch({
    model <- lmer(
      as.formula(formula_str),
      data = df,
      REML = FALSE,
      control = lmerControl(optimizer = "bobyqa")
    )
    
    # Check convergence
    if (length(model@optinfo$conv$lme4$messages) > 0) {
      warning(sprintf(
        "Gene %s: %s",
        gene,
        model@optinfo$conv$lme4$messages
      ))
    }
    
    # Extract results
    coef_summary <- summary(model)$coefficients
    effects_df <- as.data.frame(coef_summary) %>%
      tibble::rownames_to_column("term") %>%
      rename(
        estimate = Estimate,
        std_error = `Std. Error`,
        df = df,
        t_value = `t value`,
        p_value = `Pr(>|t|)`
      )
    
    return(list(
      gene = gene,
      model = model,
      effects = effects_df,
      anova = anova(model),
      converged = TRUE,
      singular = isSingular(model),
      formula = formula_str
    ))
    
  }, error = function(e) {
    return(list(
      gene = gene,
      converged = FALSE,
      error = as.character(e$message),
      formula = formula_str
    ))
  })
}


#' Run Linear Mixed Model Analysis for Multiple Genes
#'
#' Main function for LMM analysis. Fits mixed-effects models to account for
#' patient-level random effects while testing fixed effects of treatment,
#' timepoint, and their interactions.
#'
#' @param seurat_obj Seurat object.
#' @param genes Character vector. Gene names to analyze. If NULL, uses all genes.
#' @param config List. Configuration from create_analysis_config().
#' @param fixed_effects Character vector. Fixed effect terms 
#'   (e.g., c("drug", "timepoint", "response")).
#' @param interactions Character vector. Interaction terms 
#'   (e.g., c("drug:timepoint", "drug:response:timepoint")).
#' @param random_effects String. Random effects formula (default: "(1|patient)").
#' @param reference_levels Named list. Reference levels for factors 
#'   (e.g., list(timepoint = "pre")).
#' @param n_cores Integer. Number of CPU cores for parallel processing.
#' @param verbose Logical. Print progress messages.
#'
#' @return A list with components:
#'   \describe{
#'     \item{results}{List of model results for each gene (includes lmerMod objects)}
#'     \item{summary}{Tidy data frame with all model coefficients}
#'     \item{converged_count}{Number of successfully fitted models}
#'     \item{total_count}{Total number of genes attempted}
#'     \item{formula}{The model formula used}
#'   }
#'
#' @export
#'
#' @examples
#' # Full interaction model
#' lmm_res <- run_lmm_analysis(
#'   seurat_obj = sobj,
#'   genes = top_genes,
#'   config = config,
#'   fixed_effects = c("drug", "timepoint", "response"),
#'   interactions = c("drug:timepoint", "drug:response", 
#'                    "timepoint:response", "drug:timepoint:response"),
#'   reference_levels = list(timepoint = "pre", response = "NR")
#' )
#' 
#' # Simple pre-post comparison by drug
#' lmm_res <- run_lmm_analysis(
#'   seurat_obj = sobj,
#'   genes = top_genes,
#'   config = config,
#'   fixed_effects = c("drug", "timepoint"),
#'   interactions = c("drug:timepoint")
#' )
run_lmm_analysis <- function(
    seurat_obj,
    genes = NULL,
    config = create_analysis_config(),
    fixed_effects = c("drug", "timepoint", "response"),
    interactions = c("drug:timepoint:response"),
    random_effects = "(1|patient)",
    reference_levels = list(timepoint = "pre"),
    n_cores = parallel::detectCores() - 1,
    verbose = TRUE
) {
  
  # ---
  # Input Validation
  # ---
  
  if (is.null(genes)) {
    genes <- rownames(seurat_obj)
    if (length(genes) > 100) {
      warning(sprintf(
        "No genes specified. Using first 100 genes as test. (Total: %d)",
        length(genes)
      ))
      genes <- genes[1:100]
    }
  }
  
  # Map config names to column names
  map_term <- function(term) {
    for (config_name in names(config)) {
      pattern <- paste0("\\b", config_name, "\\b")
      term <- gsub(pattern, config[[config_name]], term)
    }
    return(term)
  }
  
  fixed_mapped <- sapply(fixed_effects, map_term)
  interactions_mapped <- sapply(interactions, map_term)
  random_mapped <- map_term(random_effects)
  
  # Build formula
  formula_str <- paste(
    "Expression ~",
    paste(c(fixed_mapped, interactions_mapped), collapse = " + "),
    "+",
    random_mapped
  )
  
  if (verbose) {
    message("Model formula: ", formula_str)
  }
  
  # Identify factor columns
  all_terms <- c(fixed_mapped, 
                 unlist(strsplit(interactions_mapped, ":")),
                 gsub(".*\\|(.*)\\)", "\\1", random_mapped))
  factor_cols <- unique(all_terms[all_terms %in% colnames(seurat_obj@meta.data)])
  
  # Map reference levels
  ref_levels_mapped <- list()
  if (!is.null(reference_levels)) {
    for (name in names(reference_levels)) {
      col <- if (name %in% names(config)) config[[name]] else name
      ref_levels_mapped[[col]] <- reference_levels[[name]]
    }
  }
  
  # ---
  # Extract Data
  # ---
  
  expr_matrix <- GetAssayData(seurat_obj, slot = "data")
  metadata <- seurat_obj@meta.data
  
  # Filter to genes present in matrix
  genes <- intersect(genes, rownames(expr_matrix))
  if (length(genes) == 0) {
    stop("No valid genes found in expression matrix!")
  }
  
  if (verbose) {
    message(sprintf(
      "Fitting LMMs for %d genes using %d cores...",
      length(genes), n_cores
    ))
    message(sprintf("Samples: %d", ncol(expr_matrix)))
  }
  
  # ---
  # Fit Models
  # ---
  
  fit_wrapper <- function(gene) {
    fit_single_gene_lmm(
      gene = gene,
      expr_vector = as.numeric(expr_matrix[gene, ]),
      metadata = metadata,
      formula_str = formula_str,
      factor_cols = factor_cols,
      reference_levels = ref_levels_mapped
    )
  }
  
  if (n_cores > 1) {
    # Parallel execution
    cl <- makeCluster(n_cores)
    on.exit(stopCluster(cl), add = TRUE)
    
    clusterEvalQ(cl, {
      library(lme4)
      library(lmerTest)
      library(dplyr)
      library(tibble)
    })
    
    clusterExport(
      cl,
      c("fit_single_gene_lmm", "metadata", "expr_matrix",
        "formula_str", "factor_cols", "ref_levels_mapped"),
      envir = environment()
    )
    
    results <- parLapply(cl, genes, fit_wrapper)
    
  } else {
    # Sequential execution
    results <- lapply(genes, fit_wrapper)
  }
  
  names(results) <- genes
  
  # ---
  # Summarize Results
  # ---
  
  converged <- results[sapply(results, function(x) x$converged)]
  failed <- results[!sapply(results, function(x) x$converged)]
  singular <- converged[sapply(converged, function(x) 
    isTRUE(x$singular))]
  
  if (verbose) {
    message(sprintf(
      "\n✅ Converged: %d/%d (%.1f%%)",
      length(converged), length(results),
      100 * length(converged) / length(results)
    ))
    
    if (length(singular) > 0) {
      message(sprintf(
        "⚠️  Singular fit warning: %d models",
        length(singular)
      ))
    }
    
    if (length(failed) > 0) {
      message(sprintf(
        "❌ Failed: %d models",
        length(failed)
      ))
      
      # Show error summary
      error_types <- table(sapply(failed, function(x) x$error))
      for (err in names(error_types)) {
        message(sprintf("   - %s: %d", err, error_types[err]))
      }
    }
  }
  
  # Create summary table
  if (length(converged) > 0) {
    summary_df <- purrr::map_dfr(converged, "effects", .id = "gene") %>%
      group_by(term) %>%
      mutate(
        p_adj = p.adjust(p_value, method = "BH"),
        significant = p_adj < 0.05
      ) %>%
      ungroup() %>%
      arrange(term, p_adj)
  } else {
    warning("No models converged successfully!")
    summary_df <- NULL
  }
  
  # ---
  # Return Results
  # ---
  
  return(list(
    results = results,
    summary = summary_df,
    converged_count = length(converged),
    total_count = length(results),
    formula = formula_str,
    singular_count = length(singular),
    failed_count = length(failed)
  ))
}


# ===
# SECTION 6: POST-HOC ANALYSIS & GENE SELECTION
# ===

#' Extract Genes with Significant Effects
#'
#' Filters LMM summary table to genes with significant effects for specific terms.
#'
#' @param lmm_summary Data frame. The $summary component from run_lmm_analysis().
#' @param term_pattern String. Regex pattern to match model terms 
#'   (e.g., "drug.*:.*timepoint" for drug-by-timepoint interactions).
#' @param p_threshold Numeric. Adjusted p-value threshold (default 0.05).
#' @param effect_threshold Numeric. Minimum absolute effect size (default 0).
#' @param top_n Integer. Maximum number of genes to return.
#' @param rank_by String. How to rank genes: "effect" or "pvalue".
#'
#' @return Data frame of significant genes with their statistics.
#' @export
#'
#' @examples
#' # Find genes with drug-by-timepoint interaction
#' interaction_genes <- extract_significant_genes(
#'   lmm_summary = lmm_res$summary,
#'   term_pattern = "drug.*:.*timepoint",
#'   top_n = 50
#' )
#' 
#' # Find genes with main effect of response
#' response_genes <- extract_significant_genes(
#'   lmm_summary = lmm_res$summary,
#'   term_pattern = "^response",  # ^ for exact start
#'   top_n = 100
#' )
extract_significant_genes <- function(
    lmm_summary,
    term_pattern,
    p_threshold = 0.05,
    effect_threshold = 0,
    top_n = 50,
    rank_by = c("effect", "pvalue")
) {
  
  rank_by <- match.arg(rank_by)
  
  # Filter by term pattern
  filtered <- lmm_summary %>%
    filter(grepl(term_pattern, term, perl = TRUE))
  
  if (nrow(filtered) == 0) {
    warning(sprintf("No terms matched pattern: %s", term_pattern))
    return(data.frame())
  }
  
  # Filter by significance thresholds
  sig_genes <- filtered %>%
    filter(
      p_adj < p_threshold,
      abs(estimate) > effect_threshold
    )
  
  if (nrow(sig_genes) == 0) {
    message("No genes met significance thresholds")
    return(data.frame())
  }
  
  # Summarize per gene (take strongest effect if multiple terms)
  gene_summary <- sig_genes %>%
    group_by(gene) %>%
    summarise(
      n_significant_terms = n(),
      max_abs_effect = max(abs(estimate)),
      min_p_adj = min(p_adj),
      terms = paste(term, collapse = "; "),
      .groups = "drop"
    )
  
  # Rank genes
  if (rank_by == "effect") {
    gene_summary <- gene_summary %>%
      arrange(desc(max_abs_effect), min_p_adj)
  } else {
    gene_summary <- gene_summary %>%
      arrange(min_p_adj, desc(max_abs_effect))
  }
  
  # Select top N
  result <- gene_summary %>%
    head(top_n)
  
  message(sprintf(
    "Found %d significant genes (showing top %d)",
    nrow(gene_summary), nrow(result)
  ))
  
  return(result)
}


#' Find Genes with Drug-Specific Response Patterns
#'
#' Identifies genes where treatment response differs by drug type.
#'
#' @param lmm_summary Data frame from run_lmm_analysis()$summary.
#' @param config List. Configuration object.
#' @param focus_drug String. Optional. Specific drug to focus on.
#' @param ... Additional arguments passed to extract_significant_genes().
#'
#' @return Data frame of genes with drug-specific response patterns.
#' @export
find_drug_response_genes <- function(
    lmm_summary,
    config,
    focus_drug = NULL,
    ...
) {
  
  # Build pattern for drug:response interaction
  drug_col <- config$drug
  response_col <- config$response
  
  pattern <- sprintf("(%s.*:.*%s|%s.*:.*%s)",
                     drug_col, response_col,
                     response_col, drug_col)
  
  if (!is.null(focus_drug)) {
    pattern <- sprintf("(%s).*%s", pattern, focus_drug)
  }
  
  extract_significant_genes(
    lmm_summary = lmm_summary,
    term_pattern = pattern,
    ...
  )
}


#' Find Genes with Temporal Response Patterns
#'
#' Identifies genes where expression changes over time depend on treatment response.
#'
#' @param lmm_summary Data frame from run_lmm_analysis()$summary.
#' @param config List. Configuration object.
#' @param ... Additional arguments passed to extract_significant_genes().
#'
#' @return Data frame of genes with time-by-response interactions.
#' @export
find_temporal_response_genes <- function(
    lmm_summary,
    config,
    ...
) {
  
  timepoint_col <- config$timepoint
  response_col <- config$response
  
  pattern <- sprintf("(%s.*:.*%s|%s.*:.*%s)",
                     timepoint_col, response_col,
                     response_col, timepoint_col)
  
  extract_significant_genes(
    lmm_summary = lmm_summary,
    term_pattern = pattern,
    ...
  )
}


# ===
# SECTION 7: VISUALIZATION
# ===

#' Create Volcano Plot from LMM Results
#'
#' Visualizes effect sizes and p-values for a specific model term.
#'
#' @param lmm_summary Data frame from run_lmm_analysis()$summary.
#' @param term_pattern String. Regex to filter model terms.
#' @param effect_threshold Numeric. Threshold for "large effect".
#' @param p_threshold Numeric. P-value significance threshold.
#' @param label_top Integer. Number of top genes to label.
#' @param title String. Plot title.
#'
#' @return ggplot object.
#' @export
#'
#' @examples
#' plot_volcano(
#'   lmm_summary = lmm_res$summary,
#'   term_pattern = "drug.*:.*timepoint",
#'   title = "Drug-by-Time Interaction Effects"
#' )
plot_volcano <- function(
    lmm_summary,
    term_pattern = NULL,
    effect_threshold = 0.5,
    p_threshold = 0.05,
    label_top = 10,
    title = "Volcano Plot"
) {
  
  # Filter data
  plot_data <- lmm_summary
  if (!is.null(term_pattern)) {
    plot_data <- plot_data %>%
      filter(grepl(term_pattern, term, perl = TRUE))
  }
  
  if (nrow(plot_data) == 0) {
    stop("No data after filtering by term_pattern")
  }
  
  # Prepare plotting variables
  plot_data <- plot_data %>%
    mutate(
      log_p = -log10(p_value),
      category = case_when(
        abs(estimate) > effect_threshold & p_adj < p_threshold ~ "Significant",
        abs(estimate) > effect_threshold ~ "Large effect only",
        p_adj < p_threshold ~ "Significant only",
        TRUE ~ "Not significant"
      ),
      category = factor(
        category,
        levels = c("Significant", "Large effect only", 
                   "Significant only", "Not significant")
      )
    )
  
  # Identify top genes to label
  top_genes <- plot_data %>%
    filter(category == "Significant") %>%
    arrange(p_adj) %>%
    head(label_top)
  
  # Create plot
  p <- ggplot(plot_data, aes(x = estimate, y = log_p, color = category)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_hline(
      yintercept = -log10(p_threshold),
      linetype = "dashed",
      color = "gray40"
    ) +
    geom_vline(
      xintercept = c(-effect_threshold, effect_threshold),
      linetype = "dashed",
      color = "gray40"
    ) +
    scale_color_manual(
      values = c(
        "Significant" = "#d62728",
        "Large effect only" = "#ff7f0e",
        "Significant only" = "#1f77b4",
        "Not significant" = "gray70"
      )
    ) +
    labs(
      title = title,
      x = "Effect Size (Estimate)",
      y = "-log10(p-value)",
      color = "Category"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold")
    )
  
  # Add gene labels if any
  if (nrow(top_genes) > 0) {
    p <- p +
      ggrepel::geom_text_repel(
        data = top_genes,
        aes(label = gene),
        size = 3,
        max.overlaps = 20,
        color = "black"
      )
  }
  
  return(p)
}


#' Plot Expression Trajectories with Mixed Model Fits
#'
#' Visualizes individual patient trajectories and estimated marginal means
#' from a fitted LMM. This is the recommended way to plot results after
#' running run_lmm_analysis().
#'
#' @param lmm_result List. Result for a single gene from run_lmm_analysis()$results.
#' @param gene String. Gene name (for titles).
#' @param seurat_obj Seurat object (for raw data).
#' @param config List. Configuration object.
#' @param facet_by String. Optional metadata column to facet by.
#'
#' @return List with three ggplot objects: emmeans, trajectories, and contrasts.
#' @export
#'
#' @examples
#' # Plot a specific gene's results
#' gene_plots <- plot_lmm_results(
#'   lmm_result = lmm_res$results[["CD8A"]],
#'   gene = "CD8A",
#'   seurat_obj = sobj,
#'   config = config
#' )
#' 
#' # Display plots
#' gene_plots$emmeans
#' gene_plots$trajectories
#' gene_plots$contrasts
plot_lmm_results <- function(
    lmm_result,
    gene,
    seurat_obj,
    config,
    facet_by = NULL
) {
  
  if (!lmm_result$converged) {
    stop(sprintf("Model for %s did not converge: %s", 
                 gene, lmm_result$error))
  }
  
  model <- lmm_result$model
  
  # ---
  # 1. Estimated Marginal Means
  # ---
  
  # Determine model terms
  terms <- names(fixef(model))
  has_drug <- any(grepl(config$drug, terms))
  has_time <- any(grepl(config$timepoint, terms))
  has_response <- any(grepl(config$response, terms))
  
  # Build emmeans formula
  if (has_drug && has_time) {
    emm_formula <- as.formula(sprintf("~ %s * %s", 
                                      config$timepoint, config$drug))
    emm <- emmeans(model, emm_formula)
  } else if (has_time) {
    emm_formula <- as.formula(sprintf("~ %s", config$timepoint))
    emm <- emmeans(model, emm_formula)
  } else {
    stop("Model must include timepoint for trajectory plotting")
  }
  
  emm_df <- as.data.frame(confint(emm))
  
  # Plot EMMs
  p_emm <- ggplot(emm_df, aes_string(
    x = config$timepoint,
    y = "emmean",
    group = if(has_drug) config$drug else "1",
    color = if(has_drug) config$drug else NULL
  )) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    geom_errorbar(
      aes(ymin = lower.CL, ymax = upper.CL),
      width = 0.1,
      linewidth = 0.8
    ) +
    labs(
      title = sprintf("%s: Estimated Marginal Means", gene),
      subtitle = "Model-adjusted means ± 95% CI",
      x = "Timepoint",
      y = "Expression (EMM)",
      color = if(has_drug) "Drug" else NULL
    ) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  
  # ---
  # 2. Individual Trajectories
  # ---
  
  # Extract raw data
  expr_data <- GetAssayData(seurat_obj, slot = "data")[gene, ]
  plot_df <- data.frame(
    expression = as.numeric(expr_data),
    seurat_obj@meta.data
  )
  
  # Plot trajectories
  p_traj <- ggplot(plot_df, aes_string(
    x = config$timepoint,
    y = "expression",
    group = config$patient,
    color = if(has_drug) config$drug else NULL
  )) +
    geom_line(alpha = 0.4) +
    geom_point(alpha = 0.6, size = 2) +
    labs(
      title = sprintf("%s: Individual Trajectories", gene),
      subtitle = "Raw data, one line per patient",
      x = "Timepoint",
      y = "Expression",
      color = if(has_drug) "Drug" else NULL
    ) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  
  if (has_drug) {
    p_traj <- p_traj + facet_wrap(as.formula(paste("~", config$drug)))
  }
  
  # ---
  # 3. Pre-Post Contrasts
  # ---
  
  if (has_drug) {
    em_by_drug <- emmeans(model, as.formula(sprintf(
      "~ %s | %s", config$timepoint, config$drug
    )))
    contrasts <- contrast(em_by_drug, method = "revpairwise", by = config$drug)
  } else {
    contrasts <- contrast(emm, method = "revpairwise")
  }
  
  contrast_df <- as.data.frame(confint(contrasts))
  
  p_contrast <- ggplot(contrast_df, aes(
    x = if(has_drug) .data[[config$drug]] else "Overall",
    y = estimate,
    ymin = lower.CL,
    ymax = upper.CL,
    color = if(has_drug) .data[[config$drug]] else NULL
  )) +
    geom_pointrange(size = 1, linewidth = 1.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    labs(
      title = sprintf("%s: Change from Pre to Post", gene),
      subtitle = "Estimated contrast ± 95% CI",
      x = if(has_drug) "Drug" else "",
      y = "Δ Expression (Post - Pre)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "none"
    )
  
  # ---
  # Return Results
  # ---
  
  return(list(
    emmeans = p_emm,
    trajectories = p_traj,
    contrasts = p_contrast,
    emm_data = emm_df,
    contrast_data = contrast_df
  ))
}


#' Create Boxplot for Pre-Post Comparison
#'
#' Simple visualization for a single gene's expression across timepoints.
#'
#' @param seurat_obj Seurat object.
#' @param gene String. Gene name.
#' @param config List. Configuration object.
#' @param split_by String. Optional metadata column to facet by.
#' @param add_stats Logical. Add statistical comparison.
#'
#' @return ggplot object.
#' @export
plot_gene_boxplot <- function(
    seurat_obj,
    gene,
    config,
    split_by = NULL,
    add_stats = TRUE
) {
  
  # Extract data
  expr_data <- GetAssayData(seurat_obj, slot = "data")[gene, ]
  plot_df <- data.frame(
    expression = as.numeric(expr_data),
    seurat_obj@meta.data
  )
  
  # Create base plot
  p <- ggplot(plot_df, aes_string(
    x = config$timepoint,
    y = "expression",
    fill = config$timepoint
  )) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.4, size = 1.5) +
    labs(
      title = sprintf("Expression: %s", gene),
      x = "Timepoint",
      y = "Normalized Expression",
      fill = "Timepoint"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none")
  
  # Add statistical comparison
  if (add_stats) {
    p <- p + stat_compare_means(
      paired = TRUE,
      method = "wilcox.test",
      label = "p.format"
    )
  }
  
  # Add faceting if requested
  if (!is.null(split_by)) {
    if (split_by %in% names(config)) {
      split_by <- config[[split_by]]
    }
    p <- p + facet_wrap(as.formula(paste("~", split_by)))
  }
  
  return(p)
}


# ===
# SECTION 8: ALTERNATIVE ANALYSIS - DELTA MATRIX APPROACH
# ===

#' Create Delta (Change Score) Matrix
#'
#' Alternative analysis approach that calculates per-patient change scores
#' (post - pre) for each gene. Useful for PCA, clustering, or when paired
#' analysis assumptions are violated.
#'
#' @param seurat_obj Seurat object.
#' @param config List. Configuration object.
#' @param assay String. Assay name (default "RNA").
#' @param slot String. Data slot (default "data").
#' @param subset_var String. Optional variable to subset by before calculation.
#' @param subset_value Value to filter subset_var by.
#' @param aggregate_fun Function to aggregate multiple AOIs per patient-timepoint.
#'
#' @return List with delta_matrix (genes x patients) and patient_metadata.
#' @export
#'
#' @examples
#' # Calculate change scores
#' delta_data <- create_delta_matrix(
#'   seurat_obj = sobj,
#'   config = config,
#'   subset_var = "ck",
#'   subset_value = "Positive"
#' )
#' 
#' # Use for downstream analysis
#' delta_sobj <- CreateSeuratObject(
#'   counts = delta_data$delta_matrix,
#'   meta.data = delta_data$patient_metadata
#' )
create_delta_matrix <- function(
    seurat_obj,
    config,
    assay = "RNA",
    slot = "data",
    subset_var = NULL,
    subset_value = NULL,
    aggregate_fun = function(x) mean(x, na.rm = TRUE)
) {
  
  # Extract data
  expr_matrix <- GetAssayData(seurat_obj, assay = assay, slot = slot)
  metadata <- seurat_obj@meta.data
  
  # Subset if requested
  if (!is.null(subset_var) && !is.null(subset_value)) {
    subset_col <- if (subset_var %in% names(config)) {
      config[[subset_var]]
    } else {
      subset_var
    }
    
    keep_samples <- metadata[[subset_col]] == subset_value
    expr_matrix <- expr_matrix[, keep_samples, drop = FALSE]
    metadata <- metadata[keep_samples, , drop = FALSE]
    
    message(sprintf(
      "Subsetted to %d samples where %s = %s",
      sum(keep_samples), subset_col, subset_value
    ))
  }
  
  # Validate required columns
  required_cols <- c(config$patient, config$timepoint, config$drug)
  missing <- setdiff(required_cols, colnames(metadata))
  if (length(missing) > 0) {
    stop(sprintf("Missing columns: %s", paste(missing, collapse = ", ")))
  }
  
  # Convert to long format
  long_data <- as.data.frame(t(expr_matrix)) %>%
    tibble::rownames_to_column("sample_id") %>%
    left_join(
      metadata %>%
        select(all_of(required_cols)) %>%
        tibble::rownames_to_column("sample_id"),
      by = "sample_id"
    ) %>%
    pivot_longer(
      cols = -c(sample_id, all_of(required_cols)),
      names_to = "gene",
      values_to = "expression"
    )
  
  # Get timepoint levels
  timepoint_levels <- unique(long_data[[config$timepoint]])
  if (length(timepoint_levels) != 2) {
    stop(sprintf(
      "Expected 2 timepoints, found %d: %s",
      length(timepoint_levels),
      paste(timepoint_levels, collapse = ", ")
    ))
  }
  
  # Assume first is "pre" and second is "post" (or alphabetical order)
  timepoint_levels <- sort(timepoint_levels)
  message(sprintf(
    "Computing delta: %s - %s",
    timepoint_levels[2], timepoint_levels[1]
  ))
  
  # Calculate delta per patient-drug-gene
  delta_long <- long_data %>%
    group_by(
      .data[[config$patient]],
      .data[[config$drug]],
      gene
    ) %>%
    summarise(
      pre_expr = aggregate_fun(
        expression[.data[[config$timepoint]] == timepoint_levels[1]]
      ),
      post_expr = aggregate_fun(
        expression[.data[[config$timepoint]] == timepoint_levels[2]]
      ),
      delta = post_expr - pre_expr,
      .groups = "drop"
    ) %>%
    filter(!is.na(delta))  # Remove patients missing either timepoint
  
  # Create patient ID for matrix columns
  delta_long <- delta_long %>%
    mutate(
      patient_id = paste(
        .data[[config$patient]],
        .data[[config$drug]],
        sep = "|"
      )
    )
  
  # Convert to wide matrix format
  delta_matrix <- delta_long %>%
    select(gene, patient_id, delta) %>%
    pivot_wider(
      names_from = patient_id,
      values_from = delta,
      values_fill = NA
    ) %>%
    column_to_rownames("gene") %>%
    as.matrix()
  
  # Create patient metadata
  patient_metadata <- delta_long %>%
    distinct(patient_id, .data[[config$patient]], .data[[config$drug]]) %>%
    column_to_rownames("patient_id")
  
  # Add response info if available
  if (config$response %in% colnames(metadata)) {
    response_map <- metadata %>%
      distinct(.data[[config$patient]], .data[[config$response]])
    
    patient_metadata <- patient_metadata %>%
      tibble::rownames_to_column("patient_id") %>%
      left_join(response_map, by = config$patient) %>%
      column_to_rownames("patient_id")
  }
  
  message(sprintf(
    "Created delta matrix: %d genes x %d patients",
    nrow(delta_matrix), ncol(delta_matrix)
  ))
  
  return(list(
    delta_matrix = delta_matrix,
    patient_metadata = patient_metadata,
    delta_long = delta_long
  ))
}


# ===
# SECTION 9: RESULTS EXPORT & REPORTING
# ===

#' Export LMM Results to Excel
#'
#' Saves LMM analysis results to a multi-sheet Excel file for easy sharing.
#'
#' @param lmm_results List from run_lmm_analysis().
#' @param output_file String. Output Excel file path.
#' @param include_model_objects Logical. Save model objects to RDS file.
#' @param top_genes_n Integer. Number of top genes to highlight in summary sheet.
#'
#' @return Invisibly returns the output file path.
#' @export
#'
#' @examples
#' export_lmm_results(
#'   lmm_results = lmm_res,
#'   output_file = "lmm_analysis_results.xlsx",
#'   top_genes_n = 100
#' )
export_lmm_results <- function(
    lmm_results,
    output_file,
    include_model_objects = TRUE,
    top_genes_n = 100
) {
  
  # Create workbook
  wb <- createWorkbook()
  
  # ---
  # Sheet 1: Analysis Summary
  # ---
  
  summary_sheet <- data.frame(
    Metric = c(
      "Total Genes Analyzed",
      "Successfully Converged",
      "Failed to Converge",
      "Singular Fit Warnings",
      "Model Formula"
    ),
    Value = c(
      lmm_results$total_count,
      lmm_results$converged_count,
      lmm_results$failed_count,
      lmm_results$singular_count,
      lmm_results$formula
    )
  )
  
  addWorksheet(wb, "Summary")
  writeData(wb, "Summary", summary_sheet)
  
  # ---
  # Sheet 2: All Coefficients
  # ---
  
  if (!is.null(lmm_results$summary)) {
    addWorksheet(wb, "All_Coefficients")
    writeData(wb, "All_Coefficients", lmm_results$summary)
    
    # Add conditional formatting for significant results
    sig_rows <- which(lmm_results$summary$p_adj < 0.05) + 1  # +1 for header
    if (length(sig_rows) > 0) {
      addStyle(
        wb, "All_Coefficients",
        style = createStyle(fgFill = "#FFF4E6"),
        rows = sig_rows,
        cols = 1:ncol(lmm_results$summary),
        gridExpand = TRUE
      )
    }
  }
  
  # ---
  # Sheet 3: Top Significant Genes
  # ---
  
  if (!is.null(lmm_results$summary)) {
    top_sig <- lmm_results$summary %>%
      filter(significant) %>%
      group_by(gene) %>%
      summarise(
        n_sig_terms = n(),
        min_p_adj = min(p_adj),
        max_abs_effect = max(abs(estimate)),
        .groups = "drop"
      ) %>%
      arrange(min_p_adj) %>%
      head(top_genes_n)
    
    addWorksheet(wb, "Top_Genes")
    writeData(wb, "Top_Genes", top_sig)
  }
  
  # ---
  # Sheet 4: Failed Models
  # ---
  
  failed_results <- lmm_results$results[
    !sapply(lmm_results$results, function(x) x$converged)
  ]
  
  if (length(failed_results) > 0) {
    failed_df <- data.frame(
      gene = names(failed_results),
      error = sapply(failed_results, function(x) x$error),
      stringsAsFactors = FALSE
    )
    
    addWorksheet(wb, "Failed_Models")
    writeData(wb, "Failed_Models", failed_df)
  }
  
  # Save workbook
  saveWorkbook(wb, output_file, overwrite = TRUE)
  message(sprintf("✅ Results exported to: %s", output_file))
  
  # ---
  # Save Model Objects
  # ---
  
  if (include_model_objects) {
    rds_file <- sub("\\.[^.]+$", "_models.rds", output_file)
    saveRDS(lmm_results$results, rds_file)
    message(sprintf("✅ Model objects saved to: %s", rds_file))
  }
  
  invisible(output_file)
}


#' Generate HTML Report for LMM Analysis
#'
#' Creates an HTML summary report with key visualizations and tables.
#'
#' @param lmm_results List from run_lmm_analysis().
#' @param seurat_obj Seurat object (needed for plotting).
#' @param config List. Configuration object.
#' @param output_file String. Output HTML file path.
#' @param top_genes Integer. Number of top genes to visualize.
#'
#' @return Invisibly returns the output file path.
#' @export
generate_lmm_report <- function(
    lmm_results,
    seurat_obj,
    config,
    output_file = "lmm_analysis_report.html",
    top_genes = 10
) {
  
  # This is a simplified version - in practice you'd use rmarkdown
  # to generate a more comprehensive report
  
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("Package 'rmarkdown' needed for report generation")
  }
  
  message("Generating HTML report...")
  
  # Create temporary Rmd file
  rmd_content <- sprintf('
---
title: "Linear Mixed Model Analysis Report"
date: "`r Sys.Date()`"
output: 
  html_document:
    toc: true
    toc_float: true
    theme: cosmo
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)
library(ggplot2)
library(dplyr)
library(knitr)
```

## Analysis Summary

- **Total genes analyzed**: %d
- **Successfully converged**: %d (%.1f%%)
- **Model formula**: `%s`

## Convergence Statistics

```{r}
data.frame(
  Status = c("Converged", "Failed", "Singular Fit"),
  Count = c(%d, %d, %d)
) %%>%% kable()
```

## Top Significant Genes

```{r}
lmm_results$summary %%>%%
  filter(significant) %%>%%
  group_by(gene) %%>%%
  summarise(
    n_sig_terms = n(),
    min_p_adj = min(p_adj),
    max_effect = max(abs(estimate))
  ) %%>%%
  arrange(min_p_adj) %%>%%
  head(%d) %%>%%
  kable(digits = 4)
```
',
                         lmm_results$total_count,
                         lmm_results$converged_count,
                         100 * lmm_results$converged_count / lmm_results$total_count,
                         lmm_results$formula,
                         lmm_results$converged_count,
                         lmm_results$failed_count,
                         lmm_results$singular_count,
                         top_genes
  )
  
  # Write and render
  rmd_file <- tempfile(fileext = ".Rmd")
  writeLines(rmd_content, rmd_file)
  
  rmarkdown::render(
    rmd_file,
    output_file = output_file,
    quiet = TRUE
  )
  
  message(sprintf("✅ Report generated: %s", output_file))
  invisible(output_file)
}


# ===
# SECTION 10: UTILITY FUNCTIONS
# ===

#' Print Analysis Configuration
#'
#' @param config List from create_analysis_config().
#' @export
print.geomx_config <- function(config) {
  cat("GeoMx Analysis Configuration\n")
  cat("============================\n\n")
  
  for (name in names(config)) {
    cat(sprintf("%-12s: %s\n", name, config[[name]]))
  }
  
  invisible(config)
}


#' Quick Summary of Seurat Object for GeoMx Data
#'
#' @param seurat_obj Seurat object.
#' @param config List. Configuration object.
#' @export
summarize_geomx_object <- function(seurat_obj, config = NULL) {
  
  cat("\nGeoMx Seurat Object Summary\n")
  cat("===========================\n\n")
  
  # Basic dimensions
  cat(sprintf("Genes: %d\n", nrow(seurat_obj)))
  cat(sprintf("Samples (AOIs): %d\n", ncol(seurat_obj)))
  
  # Metadata summary
  meta <- seurat_obj@meta.data
  
  if (!is.null(config)) {
    cat("\nMetadata Summary:\n")
    
    for (var in c("patient", "drug", "timepoint", "response")) {
      if (var %in% names(config)) {
        col <- config[[var]]
        if (col %in% colnames(meta)) {
          vals <- table(meta[[col]])
          cat(sprintf("\n%s (%s):\n", var, col))
          print(vals)
        }
      }
    }
  }
  
  # Assay info
  cat("\n\nAssays:", paste(names(seurat_obj@assays), collapse = ", "), "\n")
  
  invisible(seurat_obj)
}


#' Extract Model Coefficients for Specific Genes
#'
#' Convenience function to get detailed coefficient info for genes of interest.
#'
#' @param lmm_results List from run_lmm_analysis().
#' @param genes Character vector of gene names.
#'
#' @return Data frame of coefficients for specified genes.
#' @export
extract_gene_coefficients <- function(lmm_results, genes) {
  
  lmm_results$summary %>%
    filter(gene %in% genes) %>%
    arrange(gene, term)
}


#' Compare Screening vs LMM Results
#'
#' Diagnostic function to see how well initial screening predicts LMM significance.
#'
#' @param screening_results List from screen_genes().
#' @param lmm_results List from run_lmm_analysis().
#' @param term_pattern String. Which LMM terms to focus on.
#'
#' @return Data frame comparing both methods.
#' @export
compare_screening_vs_lmm <- function(
    screening_results,
    lmm_results,
    term_pattern = NULL
) {
  
  # Get screening p-values
  screen_df <- screening_results$all_markers %>%
    select(gene, cluster, p_val_adj_screen = p_val_adj, 
           avg_log2FC_screen = avg_log2FC)
  
  # Get LMM results
  lmm_df <- lmm_results$summary
  
  if (!is.null(term_pattern)) {
    lmm_df <- lmm_df %>%
      filter(grepl(term_pattern, term))
  }
  
  lmm_df <- lmm_df %>%
    group_by(gene) %>%
    summarise(
      min_p_adj_lmm = min(p_adj),
      max_effect_lmm = max(abs(estimate)),
      .groups = "drop"
    )
  
  # Merge
  comparison <- screen_df %>%
    left_join(lmm_df, by = "gene") %>%
    mutate(
      sig_screen = p_val_adj_screen < 0.05,
      sig_lmm = min_p_adj_lmm < 0.05,
      concordance = case_when(
        sig_screen & sig_lmm ~ "Both significant",
        sig_screen & !sig_lmm ~ "Screen only",
        !sig_screen & sig_lmm ~ "LMM only",
        TRUE ~ "Neither significant"
      )
    )
  
  # Print summary
  cat("\nConcordance Summary:\n")
  print(table(comparison$concordance))
  
  return(comparison)
}


# ===
# SECTION 11: EXAMPLE WORKFLOW FUNCTION
# ===

#' Complete LMM Analysis Workflow
#'
#' End-to-end analysis pipeline combining screening and LMM fitting.
#' This is a convenience wrapper for the most common analysis workflow.
#'
#' @param seurat_obj Seurat object.
#' @param config List from create_analysis_config().
#' @param screening_params List of parameters for screen_genes().
#' @param lmm_params List of parameters for run_lmm_analysis().
#' @param check_parity Logical. Run parity diagnostics first.
#' @param export_results Logical. Export results to Excel.
#' @param output_prefix String. Prefix for output files.
#'
#' @return List with screening results, LMM results, and file paths.
#' @export
#'
#' @examples
#' # Complete workflow
#' results <- run_complete_analysis(
#'   seurat_obj = sobj,
#'   config = config,
#'   screening_params = list(
#'     grouping_var = "timepoint",
#'     top_n = 500
#'   ),
#'   lmm_params = list(
#'     fixed_effects = c("drug", "timepoint", "response"),
#'     interactions = c("drug:timepoint:response")
#'   ),
#'   output_prefix = "IBD_analysis"
#' )
run_complete_analysis <- function(
    seurat_obj,
    config,
    screening_params = list(),
    lmm_params = list(),
    check_parity = TRUE,
    export_results = TRUE,
    output_prefix = "analysis"
) {
  
  cat("\n")
  cat("╔════════════════════════════════════════════════════════╗\n")
  cat("║   GeoMx Linear Mixed Model Analysis Pipeline           ║\n")
  cat("╚════════════════════════════════════════════════════════╝\n\n")
  
  results <- list()
  
  # ---
  # Step 1: Parity Check
  # ---
  
  if (check_parity) {
    cat("STEP 1: Checking sample parity...\n")
    cat("─────────────────────────────────────\n")
    
    parity <- diagnose_sample_parity(
      seurat_obj = seurat_obj,
      config = config,
      grouping_vars = c("patient", "drug"),
      check_var = "timepoint",
      required_values = c("pre", "post")
    )
    
    results$parity <- parity
    
    if (length(parity$failed_groups) > 0) {
      warning("Some groups failed parity check. Consider subsetting.")
    }
    
    cat("\n")
  }
  
  # ---
  # Step 2: Gene Screening
  # ---
  
  cat("STEP 2: Screening candidate genes...\n")
  cat("─────────────────────────────────────\n")
  
  screening_args <- c(
    list(seurat_obj = seurat_obj, config = config),
    screening_params
  )
  
  screening <- do.call(screen_genes, screening_args)
  results$screening <- screening
  
  cat(sprintf("\n✓ Selected %d candidate genes\n\n", 
              length(screening$top_genes)))
  
  # ---
  # Step 3: LMM Analysis
  # ---
  
  cat("STEP 3: Fitting linear mixed models...\n")
  cat("─────────────────────────────────────\n")
  
  lmm_args <- c(
    list(
      seurat_obj = seurat_obj,
      genes = screening$top_genes,
      config = config
    ),
    lmm_params
  )
  
  lmm_results <- do.call(run_lmm_analysis, lmm_args)
  results$lmm <- lmm_results
  
  cat("\n")
  
  # ---
  # Step 4: Export Results
  # ---
  
  if (export_results && !is.null(lmm_results$summary)) {
    cat("STEP 4: Exporting results...\n")
    cat("─────────────────────────────────────\n")
    
    excel_file <- paste0(output_prefix, "_lmm_results.xlsx")
    export_lmm_results(
      lmm_results = lmm_results,
      output_file = excel_file
    )
    
    results$output_files <- list(excel = excel_file)
    
    cat("\n")
  }
  
  # ---
  # Final Summary
  # ---
  
  cat("╔════════════════════════════════════════════════════════╗\n")
  cat("║   Analysis Complete!                                   ║\n")
  cat("╚════════════════════════════════════════════════════════╝\n\n")
  
  cat("Summary:\n")
  cat(sprintf("  • Genes screened: %d\n", 
              nrow(screening$all_markers)))
  cat(sprintf("  • Genes modeled: %d\n", 
              lmm_results$total_count))
  cat(sprintf("  • Models converged: %d (%.1f%%)\n",
              lmm_results$converged_count,
              100 * lmm_results$converged_count / lmm_results$total_count))
  
  if (!is.null(lmm_results$summary)) {
    n_sig <- sum(lmm_results$summary$significant)
    cat(sprintf("  • Significant effects: %d\n", n_sig))
  }
  
  cat("\n")
  
  return(results)
}


# ===
# END OF REFACTORED CODE
# ===

# Usage Example:
# 
# # 1. Setup
# config <- create_analysis_config(
#   patient = "patient_id",
#   drug = "treatment",
#   timepoint = "visit",
#   response = "clinical_response"
# )
#
# # 2. Load data
# sobj <- prepare_geomx_data(
#   count_file = "counts.xlsx",
#   metadata_file = "metadata.xlsx",
#   config = config,
#   normalize_method = "none"
# )
#
# # 3. Run complete analysis
# results <- run_complete_analysis(
#   seurat_obj = sobj,
#   config = config,
#   output_prefix = "IBD_biologics"
# )
#
# # 4. Extract and visualize top genes
# top_genes <- extract_significant_genes(
#   lmm_summary = results$lmm$summary,
#   term_pattern = "drug.*:.*timepoint"
# )
#
# # 5. Plot specific gene
# plots <- plot_lmm_results(
#   lmm_result = results$lmm$results[["CD8A"]],
#   gene = "CD8A",
#   seurat_obj = sobj,
#   config = config
# )
# plots$emmeans
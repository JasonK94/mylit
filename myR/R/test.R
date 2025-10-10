# ---
# title: "GeoMxSeurat: A package for analyzing GeoMx data using Seurat and LMM"
# author: "Refactored by Gemini"
# date: "2025-10-10"
# ---

# Pacakge Dependencies
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


#
# SECTION 1: DATA PREPARATION & VALIDATION
#

#' @title Diagnose Metadata Parity for Paired Analysis
#' @description Checks if groups in metadata have balanced and exclusive representation
#'   of specified conditions. For example, ensures each patient has both 'pre' and
#'   'post' samples and no other timepoints.
#'
#' @param sobj A Seurat object.
#' @param patient A bare column name for the patient identifier.
#' @param treatment A bare column name for the treatment/drug identifier.
#' @param timepoint A bare column name for the timepoint identifier.
#' @param tissue A bare column name for the tissue identifier.
#' @param check_key The metadata column name (as a string) to check for parity
#'   (e.g., "timepoint").
#' @param check_values A character vector of values that must all be present
#'   and represented equally within each group (e.g., `c("pre", "post")`).
#'
#' @return A list containing:
#'   \itemize{
#'     \item `diag_table`: A summary table with diagnostics for each group.
#'     \item `filtered_meta`: The metadata filtered to include only groups that pass the checks.
#'     \item `message`: A summary message indicating success or listing failing groups.
#'   }
#' @export
diagnosis_parity_legacy <- function(sobj, patient, treatment, timepoint, tissue,
                             check_key, check_values) {
  # Create a unique ID for each group and prepare the check variable
  meta <- sobj@meta.data %>%
    mutate(
      unique_id = paste0({{patient}}, {{treatment}}, {{timepoint}}, {{tissue}}),
      .val_chr  = as.character(.data[[check_key]])
    )
  
  # Filter for groups that contain all the required `check_values`
  filt <- meta %>%
    group_by(unique_id) %>%
    filter(all(check_values %in% .val_chr)) %>%
    ungroup()
  
  # Create a diagnostic table
  totals <- filt %>%
    count(unique_id, name = "total_n")
  
  counts_sel <- filt %>%
    filter(.val_chr %in% check_values) %>%
    count(unique_id, .val_chr) %>%
    tidyr::pivot_wider(names_from = .val_chr, values_from = n, values_fill = 0)
  
  diag <- totals %>%
    left_join(counts_sel, by = "unique_id") %>%
    mutate(across(all_of(intersect(check_values, names(.))), ~replace_na(.x, 0))) %>%
    rowwise() %>%
    mutate(
      n_selected = sum(c_across(all_of(check_values)), na.rm = TRUE),
      n_other    = total_n - n_selected,
      has_all    = all(c_across(all_of(check_values)) > 0, na.rm = TRUE),
      # Check if all selected value counts are identical
      parity_eq  = {
        vals <- c_across(all_of(check_values))
        vals <- vals[!is.na(vals)]
        length(vals) > 0 && length(unique(vals)) == 1
      },
      ok = has_all && parity_eq && n_other == 0
    ) %>%
    ungroup()
  
  # Generate a summary message for groups that did not pass
  no_pass_ids <- diag %>% filter(!ok) %>% pull(unique_id)
  msg <- if (length(no_pass_ids) == 0) {
    "✅ All groups passed parity & exclusivity checks."
  } else {
    glue::glue("Checkout {paste(no_pass_ids, collapse = ', ')}")
  }
  
  list(
    diag_table = diag,
    filtered_meta = filt,
    message = msg
  )
}


#' @title Detailed Metadata Parity Diagnosis
#' @description A more robust version of `diagnosis_parity`. It dynamically
#'   handles column names and provides detailed output tables for inspection.
#'
#' @param sobj A Seurat object.
#' @param patient,treatment,timepoint,tissue Optional. Column names (as strings)
#'   to define unique groups.
#' @param check_key The metadata column name (as a string) to check.
#' @param check_values A character vector of values to check for parity.
#'
#' @return A list containing:
#'   \itemize{
#'     \item `detailed`: The original metadata with added diagnostic columns.
#'     \item `diag_table`: A one-row-per-group summary table.
#'     \item `message`: A summary message.
#'     \item `paired_ids`: A vector of `SegmentDisplayName`s for samples that passed the checks.
#'   }
#' @export
diagnosis_parity <- function(sobj,
                                      patient = NULL, treatment = NULL, timepoint = NULL, tissue = NULL,
                                      check_key, check_values) {
  # Helper to create unique IDs from specified columns
  .make_unique_id <- function(meta, cols) {
    cols <- cols[!cols %in% check_key]
    cols <- cols[!sapply(cols, is.null)]
    if (length(cols) == 0) {
      meta$unique_id <- seq_len(nrow(meta))
    } else {
      meta <- meta %>%
        mutate(across(all_of(cols), ~ tidyr::replace_na(as.character(.x), ""))) %>%
        mutate(unique_id = do.call(paste0, .[cols]))
    }
    return(meta)
  }
  
  meta <- sobj@meta.data %>%
    .make_unique_id(cols = c(patient, treatment, timepoint, tissue)) %>%
    mutate(.val_chr = as.character(.data[[check_key]]))
  
  # Dynamically create names for count columns, e.g., "pre" -> "n_pre"
  safe_counts_names <- paste0("n_", make.names(check_values, unique = TRUE))
  
  # Calculate counts and diagnostic columns without using rowwise for efficiency
  diag_table <- meta %>%
    count(unique_id, .val_chr) %>%
    complete(unique_id, .val_chr = check_values, fill = list(n = 0)) %>%
    pivot_wider(names_from = .val_chr, values_from = n, values_fill = 0) %>%
    rename_with(~ paste0("n_", make.names(.x)), .cols = all_of(check_values)) %>%
    left_join(meta %>% count(unique_id, name = "total_n"), by = "unique_id") %>%
    mutate(across(everything(), ~ replace_na(.x, 0))) %>%
    mutate(
      n_selected = purrr::reduce(across(all_of(safe_counts_names)), `+`),
      n_other = total_n - n_selected,
      has_all = if_all(all_of(safe_counts_names), ~ .x > 0),
      parity_eq = purrr::pmap_lgl(across(all_of(safe_counts_names)), ~ length(unique(c(...))) == 1),
      ok = has_all & parity_eq & n_other == 0
    ) %>%
    dplyr::select(unique_id, total_n, all_of(safe_counts_names), n_selected, n_other, has_all, parity_eq, ok)
  
  # Join diagnostics back to the full metadata
  detailed <- meta %>% left_join(diag_table, by = "unique_id")
  
  # Generate message and identify passing sample IDs
  no_pass_ids <- diag_table %>% filter(!ok) %>% pull(unique_id)
  msg <- if (length(no_pass_ids) == 0) {
    "✅ All groups passed parity & exclusivity checks."
  } else {
    sprintf("Checkout %s", paste(no_pass_ids, collapse = ", "))
  }
  
  message(msg)
  message("Access results with `$detailed` or `$diag_table`.")
  
  ok_id <- detailed %>%
    filter(ok) %>%
    pull(SegmentDisplayName)
  
  list(
    detailed = detailed,
    diag_table = diag_table,
    message = msg,
    paired_ids = ok_id
  )
}


# SECTION 2: ANALYSIS CONFIGURATION & DATA LOADING

#' @title Create Analysis Configuration
#' @description Creates a list to manage metadata column names consistently
#'   throughout the analysis.
#'
#' @param patient String, column name for patient ID.
#' @param drug String, column name for drug/treatment.
#' @param timepoint String, column name for timepoint (e.g., pre/post).
#' @param ck String, column name for a stratification variable (e.g., CK status).
#' @param response String, column name for treatment response.
#' @param aoi String, column name for AOI (Area of Interest) ID.
#'
#' @return A named list containing the column names.
#' @export
create_analysis_config <- function(
    patient = "patient_id",
    drug = "drug",
    timepoint = "timepoint",
    ck = "ck_status",
    response = "response",
    aoi = "aoi_id"
) {
  list(
    patient = patient,
    drug = drug,
    timepoint = timepoint,
    ck = ck,
    response = response,
    aoi = aoi
  )
}

#' @title Prepare GeoMx Data for Seurat Analysis
#' @description Reads count and metadata files (Excel or CSV) and creates a
#'   Seurat object, ready for analysis.
#'
#' @param count_file Path to the count data file (.xlsx or .csv).
#' @param metadata_file Path to the metadata file (.xlsx or .csv).
#' @param count_matrix A raw count matrix (can be provided instead of `count_file`).
#' @param metadata A metadata data frame (can be provided instead of `metadata_file`).
#' @param config An analysis configuration object from `create_analysis_config()`.
#' @param q3_normalize Logical. If FALSE, performs standard log-normalization.
#'   Assumes data is already Q3 normalized if TRUE.
#'
#' @return A Seurat object with combined variables added to metadata.
#' @export
prepare_geomx_data <- function(count_file = NULL,
                               metadata_file = NULL,
                               count_matrix = NULL,
                               metadata = NULL,
                               config = create_analysis_config(),
                               q3_normalize = TRUE) {
  
  if (!is.null(count_file)) {
    if (grepl("\\.xlsx$", count_file)) {
      count_matrix <- openxlsx::read.xlsx(count_file, sheet = 1, rowNames = TRUE)
    } else {
      count_matrix <- read.csv(count_file, row.names = 1)
    }
  }
  
  if (!is.null(metadata_file)) {
    if (grepl("\\.xlsx$", metadata_file)) {
      metadata <- openxlsx::read.xlsx(metadata_file, sheet = 2)
    } else {
      metadata <- read.csv(metadata_file)
    }
  }
  
  # Match metadata rows to count matrix columns
  rownames(metadata) <- metadata[[config$aoi]]
  
  seurat_obj <- CreateSeuratObject(
    counts = as.matrix(count_matrix[, rownames(metadata)]), # Ensure order and presence
    meta.data = metadata,
    min.cells = 0,
    min.features = 0
  )
  
  if (!q3_normalize) {
    seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize")
  }
  
  # Create convenience combination variables
  seurat_obj$condition_full <- paste(
    seurat_obj@meta.data[[config$timepoint]],
    seurat_obj@meta.data[[config$drug]],
    seurat_obj@meta.data[[config$response]],
    sep = "_"
  )
  seurat_obj$condition_simple <- paste(
    seurat_obj@meta.data[[config$timepoint]],
    seurat_obj@meta.data[[config$drug]],
    sep = "_"
  )
  
  return(seurat_obj)
}


#
# SECTION 3: GENE SCREENING & LMM ANALYSIS
#

#' @title Quick Gene Screening using Wilcoxon Test
#' @description Performs a rapid screening of genes using `FindAllMarkers` to
#'   identify candidate genes for the more computationally intensive LMM analysis.
#'
#' @param seurat_obj A Seurat object.
#' @param config An analysis configuration object.
#' @param comparison_type One of "pre_post", "drug", "response", or "combined".
#' @param ck_subset A value to subset the data by from the `ck` column in `config`.
#' @param min_pct Minimum percentage of cells expressing the gene in either group.
#' @param logfc_threshold Log-fold change threshold.
#' @param top_n The number of top genes to return.
#' @param use_adj Logical. If TRUE, ranks genes by adjusted p-value; otherwise, by raw p-value.
#'
#' @return A list containing `markers` (all results), `top_genes` (a character
#'   vector of top gene names), and the `seurat_obj` (potentially subsetted).
#' @export
quick_screen_genes <- function(seurat_obj,
                               config = create_analysis_config(),
                               comparison_type = "pre_post",
                               ck_subset = NULL,
                               min_pct = 0.1,
                               logfc_threshold = 0.25,
                               top_n = 1000,
                               use_adj = FALSE) {
  
  if (!is.null(ck_subset)) {
    seurat_obj <- subset(seurat_obj,
                         subset = !!sym(config$ck) == ck_subset)
  }
  
  # Set identity for comparison
  ident_group <- switch(comparison_type,
                        "pre_post" = config$timepoint,
                        "drug" = config$drug,
                        "response" = config$response,
                        "combined" = "condition_full",
                        stop("Invalid comparison_type")
  )
  Idents(seurat_obj) <- seurat_obj@meta.data[[ident_group]]
  
  # Find markers
  markers <- FindAllMarkers(
    seurat_obj,
    min.pct = min_pct,
    logfc.threshold = logfc_threshold,
    only.pos = FALSE,
    return.thresh = 1 # Return all genes
  )
  
  # Select top genes based on p-value
  p_val_col <- if (use_adj) "p_val_adj" else "p_val"
  top_genes <- markers %>%
    arrange(.data[[p_val_col]]) %>%
    head(top_n) %>%
    pull(gene) %>%
    unique()
  
  return(list(
    markers = markers,
    top_genes = top_genes,
    seurat_obj = seurat_obj
  ))
}


#' @title Fit a Linear Mixed Model for a Single Gene
#' @description An internal helper function to fit an `lmer` model for one gene.
#'
#' @param gene_expr A numeric vector of expression values for a single gene.
#' @param metadata The metadata data frame for the samples.
#' @param config An analysis configuration object.
#' @param formula_components A list specifying `fixed`, `interactions`, and `random`
#'   parts of the formula. If NULL, a default full model is used.
#' @param use_config_names Logical. If TRUE, maps short names in `formula_components`
#'   (e.g., "drug") to full names in `config` (e.g., "drug_name_in_meta").
#'
#' @return A list containing the model, effects summary, ANOVA table, and convergence status.
#' @keywords internal
fit_lmm_single_gene <- function(gene_expr,
                                metadata,
                                config = create_analysis_config(),
                                formula_str=NULL,
                                formula_components = NULL,
                                use_config_names = TRUE) {
  
  df <- cbind(data.frame(Expression = gene_expr), metadata)
  
  # Generate model formula
  if(is.null(formula_str)){
    if (is.null(formula_components)) {
      # Default formula using all interactions
      formula_str <- glue::glue(
        "Expression ~ {config$drug}*{config$timepoint}*{config$response} + (1|{config$patient})"
      )
    } else {
      # Custom formula from components
      fixed_effects <- formula_components$fixed
      interactions <- formula_components$interactions
      random_effects <- formula_components$random
      
      if (use_config_names) {
        # Map generic terms like "drug" to the actual column names in config
        for (name in names(config)) {
          fixed_effects <- gsub(paste0("\\b", name, "\\b"), config[[name]], fixed_effects)
          interactions <- gsub(paste0("\\b", name, "\\b"), config[[name]], interactions)
          random_effects <- gsub(paste0("\\b", name, "\\b"), config[[name]], random_effects)
        }
      }
      formula_str <- paste(
        "Expression ~",
        paste(c(fixed_effects, interactions), collapse = " + "),
        "+", random_effects
      )
    }
  }else{
    formula_str=as.formula(formula_str)
    cat(paste0("formula is... ", formula_str))
  }
  
  tryCatch({
    # Set reference levels for factors to ensure consistent contrasts
    for (col in c(config$drug, config$response, config$timepoint)) {
      if (col %in% names(df) && is.factor(df[[col]])) {
        df[[col]] <- relevel(df[[col]], ref = levels(df[[col]])[1])
      } else if (col %in% names(df)) {
        df[[col]] <- factor(df[[col]])
        df[[col]] <- relevel(df[[col]], ref = levels(df[[col]])[1])
      }
    }
    
    model <- lmer(as.formula(formula_str), data = df, REML = FALSE)
    coef_summary <- summary(model)$coefficients
    effects <- as.data.frame(coef_summary) %>%
      tibble::rownames_to_column("term") %>%
      rename(
        estimate = Estimate, std_error = `Std. Error`,
        t_value = `t value`, p_value = `Pr(>|t|)`
      )
    
    return(list(
      model = model,
      effects = effects,
      anova = anova(model),
      converged = TRUE,
      formula = formula_str
    ))
  }, error = function(e) {
    return(list(
      converged = FALSE,
      error = e$message,
      formula = formula_str
    ))
  })
}

#' @title Summarize LMM Results
#' @description An internal helper to combine results from multiple LMMs and
#'   calculate adjusted p-values.
#'
#' @param lmm_results A list of results from `fit_lmm_single_gene`.
#' @param config An analysis configuration object.
#'
#' @return A tidy data frame summarizing all model effects.
#' @keywords internal
summarize_lmm_results <- function(lmm_results, config) {
  converged <- lmm_results[sapply(lmm_results, `[[`, "converged")]
  if (length(converged) == 0) {
    warning("No models converged!")
    return(NULL)
  }
  
  all_effects <- purrr::map_dfr(converged, "effects", .id = "gene")
  
  # Adjust p-values for each term across all genes
  all_effects <- all_effects %>%
    group_by(term) %>%
    mutate(
      p_adj = p.adjust(p_value, method = "BH"),
      significant = p_adj < 0.05
    ) %>%
    ungroup()
  
  return(all_effects)
}


#' @title Run Linear Mixed Models for Multiple Genes
#' @description Applies `fit_lmm_single_gene` to a list of genes in parallel.
#'   This is the main workhorse function for the LMM analysis.
#'
#' @param seurat_obj A Seurat object.
#' @param genes A character vector of gene names to analyze.
#' @param config An analysis configuration object.
#' @param formula_components A list specifying the model formula. See `fit_lmm_single_gene`.
#' @param use_config_names Logical. See `fit_lmm_single_gene`.
#' @param n_cores The number of CPU cores to use for parallel processing.
#' @param verbose Logical. If TRUE, prints progress messages.
#'
#' @return A list containing:
#'   \itemize{
#'     \item `raw_results`: A list where each element is the complete result
#'           from `fit_lmm_single_gene` for one gene (including the `lmerMod` object).
#'     \item `summary`: A tidy data frame summarizing all model effects, created
#'           by `summarize_lmm_results`.
#'     \item `converged_genes`: The number of successfully fitted models.
#'     \item `total_genes`: The total number of genes attempted.
#'   }
#' @export
run_lmm_multiple_genes <- function(seurat_obj,
                                   genes = NULL,
                                   config = create_analysis_config(),
                                   formula_str = NULL,
                                   formula_components = NULL,
                                   use_config_names = TRUE,
                                   n_cores = parallel::detectCores() - 1,
                                   verbose = TRUE) {
  
  if (is.null(genes)) {
    genes <- rownames(seurat_obj)[1:100]
    warning("No genes specified. Using first 100 genes as a test.")
  }
  
  expr_matrix <- GetAssayData(seurat_obj, slot = "data")
  metadata <- seurat_obj@meta.data
  
  if (verbose) {
    message(sprintf("Running LMM for %d genes using %d cores...", length(genes), n_cores))
  }
  
  # run_fun은 부모 함수의 formula_str을 그대로 사용하면 됩니다.
  run_fun <- function(gene) {
    if (gene %in% rownames(expr_matrix)) {
      return(fit_lmm_single_gene(
        gene_expr = as.numeric(expr_matrix[gene, ]),
        metadata = metadata,
        config = config,
        formula_str = formula_str, # as.formula 제거하고 그냥 전달
        formula_components = formula_components,
        use_config_names = use_config_names
      ))
    } else {
      return(list(gene = gene, converged = FALSE, error = "Gene not found"))
    }
  }
  
  # if/else 문을 단순화했습니다. run_fun 내부 로직은 formula_str이 NULL이든 아니든 동일합니다.
  if (n_cores > 1) {
    cl <- makeCluster(n_cores)
    clusterEvalQ(cl, { library(lme4); library(lmerTest) })
    
    # <<-- ⭐️ FIX 1: 여기에 "formula_str"를 추가합니다. -->>
    clusterExport(cl, c("fit_lmm_single_gene", "config", "metadata", "expr_matrix",
                        "formula_components", "use_config_names", "formula_str"), 
                  envir = environment())
    
    results <- parLapply(cl, genes, run_fun)
    stopCluster(cl)
  } else {
    results <- lapply(genes, run_fun)
  }
  
  names(results) <- genes
  summary_df <- summarize_lmm_results(results, config)
  
  return(list(
    raw_results = results,
    summary = summary_df,
    converged_genes = sum(sapply(results, `[[`, "converged")),
    total_genes = length(genes)
  ))
}


#
# SECTION 4: POST-HOC ANALYSIS & INTERPRETATION
#

#' @title Find Genes with Differential Response to Treatment
#' @description From a LMM summary, this function identifies genes where the
#'   treatment response (e.g., Responder vs. Non-Responder) differs by drug.
#'
#' @param lmm_summary The `$summary` data frame from `run_lmm_multiple_genes`.
#' @param config An analysis configuration object.
#' @param drug_name An optional drug name to focus the analysis on.
#' @param top_n The number of top genes to return.
#'
#' @return A data frame of top genes ranked by effect size.
#' @export
find_response_differential_genes <- function(lmm_summary,
                                             config,
                                             drug_name = NULL,
                                             top_n = 50) {
  
  interaction_pattern <- paste0(config$drug, ".*:", ".*", config$response)
  if (!is.null(drug_name)) {
    interaction_pattern <- paste0("(", interaction_pattern, ").*", drug_name)
  }
  
  drug_terms <- lmm_summary %>%
    filter(grepl(interaction_pattern, term, perl = TRUE))
  
  top_genes <- drug_terms %>%
    group_by(gene) %>%
    summarize(
      max_effect = max(abs(estimate)),
      min_p_adj = min(p_adj),
      .groups = "drop"
    ) %>%
    arrange(desc(max_effect)) %>%
    head(top_n)
  
  return(top_genes)
}

#' @title Find Drug-Specific Genes
#' @description Identifies genes whose expression is most strongly associated with
#'   specific drugs, based on main effects or interactions from the LMM summary.
#'
#' @param lmm_summary The `$summary` data frame from `run_lmm_multiple_genes`.
#' @param config An analysis configuration object.
#' @param top_n The number of top genes to return.
#'
#' @return A data frame of top genes ranked by effect size.
#' @export
find_drug_specific_genes <- function(lmm_summary,
                                     config,
                                     top_n = 50) {
  
  drug_terms <- lmm_summary %>%
    filter(grepl(config$drug, term) & !grepl(":", term)) # Main effect of drug
  
  top_genes <- drug_terms %>%
    group_by(gene) %>%
    summarize(
      max_effect = max(abs(estimate)),
      min_p_adj = min(p_adj),
      .groups = "drop"
    ) %>%
    arrange(desc(max_effect)) %>%
    head(top_n)
  
  return(top_genes)
}

#' @title Find Genes with Strong Interaction Effects
#' @description Identifies genes with the largest interaction effects from the LMM summary.
#'
#' @param lmm_summary The `$summary` data frame from `run_lmm_multiple_genes`.
#' @param config An analysis configuration object.
#' @param interaction_type One of "drug_response", "drug_time", or "triple".
#' @param top_n The number of top genes to return.
#'
#' @return A data frame of top genes ranked by interaction strength.
#' @export
find_interaction_genes <- function(lmm_summary,
                                   config,
                                   interaction_type = "drug_response",
                                   top_n = 50) {
  
  pattern <- switch(interaction_type,
                    "drug_response" = paste0(config$drug, ".*:", ".*", config$response),
                    "drug_time" = paste0(config$drug, ".*:", ".*", config$timepoint),
                    "triple" = paste0(config$drug, ".*:", ".*", config$timepoint, ".*:", ".*", config$response),
                    stop("Invalid interaction_type")
  )
  
  interaction_genes <- lmm_summary %>%
    filter(grepl(pattern, term)) %>%
    group_by(gene) %>%
    summarize(
      interaction_strength = max(abs(estimate)),
      min_p_adj = min(p_adj),
      .groups = "drop"
    ) %>%
    arrange(desc(interaction_strength)) %>%
    head(top_n)
  
  return(interaction_genes)
}


#
# SECTION 5: VISUALIZATION
#

#' @title Plot Volcano from LMM Results
#' @description A generalized volcano plot function for LMM summary tables.
#'
#' @param lmm_summary The `$summary` data frame from `run_lmm_multiple_genes`.
#' @param x_col String, column name for the x-axis (effect size).
#' @param y_col String, column name for the y-axis (p-value).
#' @param filter_col String, column name to filter by.
#' @param filter_pattern A regex pattern to select specific terms for plotting.
#' @param title Plot title.
#' @param effect_threshold Threshold for significant effect size.
#' @param p_threshold P-value threshold.
#' @param resize_x Logical, automatically resize x-axis to the 95th percentile of effects.
#'
#' @return A ggplot object.
#' @export
plot_volcano <- function(lmm_summary,
                         x_col = "estimate",
                         y_col = "p_value",
                         filter_col = "term",
                         filter_pattern = NULL,
                         title = "Volcano Plot",
                         effect_threshold = 0.5,
                         p_threshold = 0.05,
                         resize_x = TRUE) {
  
  plot_data <- lmm_summary
  if (!is.null(filter_pattern)) {
    plot_data <- plot_data %>%
      filter(grepl(filter_pattern, .data[[filter_col]]))
  }
  
  plot_data <- plot_data %>%
    mutate(
      effect_size = .data[[x_col]],
      pval = .data[[y_col]],
      log_p = -log10(pval),
      category = case_when(
        abs(effect_size) > effect_threshold & pval < p_threshold ~ "Significant",
        abs(effect_size) > effect_threshold ~ "Large effect",
        pval < p_threshold ~ "Small effect",
        TRUE ~ "Not significant"
      )
    )
  
  p <- ggplot(plot_data, aes(x = effect_size, y = log_p, color = category)) +
    geom_point(alpha = 0.6) +
    geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "gray") +
    geom_vline(xintercept = c(-effect_threshold, effect_threshold), linetype = "dashed", color = "gray") +
    scale_color_manual(values = c("Significant" = "red", "Large effect" = "orange",
                                  "Small effect" = "blue", "Not significant" = "gray")) +
    labs(title = title, x = "Effect Size (Estimate)", y = "-log10(p-value)") +
    theme_bw()
  
  if (resize_x) {
    xlim_val <- quantile(abs(plot_data$effect_size), 0.95, na.rm = TRUE)
    p <- p + coord_cartesian(xlim = c(-xlim_val, xlim_val))
  }
  return(p)
}


#' @title Plot Pre-Post Changes with Boxplots
#' @description Visualize gene expression changes between timepoints.
#'
#' @param seurat_obj A Seurat object.
#' @param gene The name of the gene to plot.
#' @param config An analysis configuration object.
#' @param split_by A character vector of metadata columns to facet the plot by.
#'
#' @return A ggplot object.
#' @export
plot_pre_post_boxplot <- function(seurat_obj,
                                  gene,
                                  config,
                                  split_by = NULL) {
  
  df <- data.frame(
    Expression = GetAssayData(seurat_obj, slot = "data")[gene, ],
    seurat_obj@meta.data
  )
  
  p <- ggplot(df, aes_string(x = config$timepoint, y = "Expression")) +
    geom_boxplot(aes_string(fill = config$timepoint), alpha = 0.7, outlier.shape = NA) +
    geom_point(position = position_jitterdodge(), alpha = 0.5) +
    stat_compare_means(paired = TRUE, method = "wilcox.test") +
    labs(title = paste("Expression Changes:", gene), y = "Normalized Expression") +
    theme_bw()
  
  if (!is.null(split_by)) {
    if (length(split_by) == 1) {
      p <- p + facet_wrap(as.formula(paste("~", split_by)))
    } else if (length(split_by) > 1) {
      p <- p + facet_grid(as.formula(paste(split_by[1], "~", split_by[2])))
    }
  }
  return(p)
}

#' @title Plot LMM Interaction Effects from a Fitted Model
#' @description **Recommended Method.** Takes a pre-computed `lmer` model object
#'   and generates plots for estimated marginal means, individual patient trajectories,
#'   and the delta (post-pre) contrast. This is highly efficient for plotting
#'   many genes after running `run_lmm_multiple_genes`.
#'
#' @param lmm_fit An `lmerMod` object, typically from `$raw_results[[gene]]$model`.
#' @param gene The name of the gene (for plot titles).
#' @param sobj A Seurat object, needed to get raw data for the spaghetti plot.
#' @param config An analysis configuration object.
#'
#' @return A list of ggplot objects: `plot_emmeans`, `plot_spaghetti`, `plot_delta`.
#' @export
plot_lmm_interaction <- function(lmm_fit, gene, sobj, config) {
  
  # 1. Estimated Marginal Means Plot
  emm <- emmeans(lmm_fit, ~ time * drug)
  emm_df <- as.data.frame(confint(emm))
  
  p1 <- ggplot(emm_df, aes(x = time, y = emmean, group = drug, color = drug)) +
    geom_point(size = 3) +
    geom_line(linewidth = 1) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1) +
    labs(
      title = sprintf("%s: time × drug (Estimated Marginal Means ± 95%% CI)", gene),
      x = "Time", y = "Estimated Mean Expression"
    ) +
    theme_bw()
  
  # 2. Individual Trajectory (Spaghetti) Plot
  df_raw <- data.frame(
    y = GetAssayData(sobj, slot = "data")[gene, ],
    sobj@meta.data
  ) %>%
    rename(time = !!sym(config$timepoint),
           drug = !!sym(config$drug),
           patient = !!sym(config$patient))
  
  p2 <- ggplot(df_raw, aes(x = time, y = y, group = patient, color = drug)) +
    geom_line(alpha = 0.35) +
    geom_point(size = 2, alpha = 0.7) +
    facet_wrap(~drug) +
    labs(
      title = sprintf("%s: Individual pre→post trajectories (raw data)", gene),
      x = "Time", y = "Expression"
    ) +
    theme_bw()
  
  # 3. Delta (post - pre) Contrast Plot
  em_time_by_drug <- emmeans(lmm_fit, ~ time | drug)
  deltas <- contrast(em_time_by_drug, method = "revpairwise", by = "drug")
  deltas_df <- as.data.frame(confint(deltas))
  
  p3 <- ggplot(deltas_df, aes(x = drug, y = estimate, ymin = lower.CL, ymax = upper.CL, color = drug)) +
    geom_pointrange() +
    geom_hline(yintercept = 0, linetype = 2) +
    labs(
      title = sprintf("%s: Δ(post-pre) by drug (EMMeans Contrast)", gene),
      x = "Drug", y = expression(Delta ~ "(post - pre)")
    ) +
    theme_bw() +
    theme(legend.position = "none")
  
  list(
    plot_emmeans = p1,
    plot_spaghetti = p2,
    plot_delta = p3,
    emm_data = emm_df,
    delta_data = deltas_df
  )
}

#' @title Plot LMM Interaction Effects Directly from Seurat Object
#' @description **Legacy/Convenience Method.** Fits a new `lmer` model for a
#'   single gene and generates interaction plots.
#' @note This is inefficient for analyzing many genes as it re-runs the model
#'   every time. For batch analysis, please use the `run_lmm_multiple_genes` ->
#'   `plot_lmm_interaction` workflow.
#'
#' @param sobj A Seurat object.
#' @param gene The gene to analyze.
#' @param config An analysis configuration object.
#' @param nested Logical, if TRUE, fits `(1 | drug/patient)`, else `(1 | patient)`.
#' @return A list of plots and data frames, similar to `plot_lmm_interaction`.
#' @export
plot_interaction_for_gene_direct <- function(sobj,
                                             gene,
                                             config,
                                             nested = FALSE) {
  
  df <- data.frame(
    y = GetAssayData(sobj, slot = "data")[gene, ],
    sobj@meta.data
  ) %>%
    rename(
      time = !!sym(config$timepoint),
      drug = !!sym(config$drug),
      patient = !!sym(config$patient)
    ) %>%
    mutate(
      time = factor(time),
      drug = factor(drug),
      patient = factor(patient)
    )
  
  # Fit model on the fly
  formula_str <- if (nested) {
    "y ~ time * drug + (1 | drug/patient)"
  } else {
    "y ~ time * drug + (1 | patient)"
  }
  fit <- lmer(as.formula(formula_str), data = df, REML = FALSE)
  
  # Re-use the efficient plotting function
  return(plot_lmm_interaction(fit, gene, sobj, config))
}


#
# SECTION 6: ALTERNATIVE DELTA-MATRIX ANALYSIS
#

#' @title Create a Delta Matrix by Patient
#' @description An alternative analysis approach. Calculates the delta (post - pre)
#'   expression value for each gene for each patient, creating a new "delta matrix".
#'   This matrix can then be used for downstream analyses like PCA or clustering
#'   to find patterns in patient-level responses.
#'
#' @param sobj A Seurat object.
#' @param time_col,drug_col,patient_col,segment_col Strings for metadata columns.
#' @param time_levels A two-element character vector, e.g., `c("pre", "post")`.
#' @param subset_comp An optional value to subset by from the `segment_col`.
#' @param agg_fun A function to aggregate expression values if a patient has
#'   multiple AOIs for a timepoint (e.g., `mean` or `median`).
#'
#' @return A list containing `delta_mat` (genes x samples matrix) and
#'   `sample_meta` (metadata for the new samples).
#' @export
make_delta_matrix_by_patient <- function(
    sobj,
    assay = "RNA", layer = "data",
    time_col = "treatment", drug_col = "drug", patient_col = "emrid", segment_col = "ck",
    time_levels = c("pre", "post"), subset_comp = NULL,
    agg_fun = function(v) mean(v, na.rm = TRUE)
){
  X <- GetAssayData(sobj, assay = assay, slot = layer)
  meta <- sobj@meta.data
  
  if (!is.null(subset_comp)) {
    keep <- meta[[segment_col]] == subset_comp
    X <- X[, keep, drop = FALSE]
    meta <- meta[keep, , drop = FALSE]
  }
  
  # Prepare data in a long format
  long_data <- as.data.frame(t(X)) %>%
    tibble::rownames_to_column("sample_id") %>%
    left_join(
      meta %>%
        select(all_of(c(patient_col, drug_col, segment_col, time_col))) %>%
        tibble::rownames_to_column("sample_id"),
      by = "sample_id"
    ) %>%
    pivot_longer(
      cols = -c(sample_id, !!sym(patient_col), !!sym(drug_col), !!sym(segment_col), !!sym(time_col)),
      names_to = "gene",
      values_to = "expression"
    )
  
  # Calculate delta values
  delta_long <- long_data %>%
    group_by(!!sym(patient_col), !!sym(drug_col), !!sym(segment_col), gene) %>%
    # Aggregate multiple AOIs per timepoint if they exist
    summarise(
      expr = agg_fun(expression[!!sym(time_col) == time_levels[2]]) -
        agg_fun(expression[!!sym(time_col) == time_levels[1]]),
      .groups = "drop"
    ) %>%
    rename(delta = expr)
  
  # Create the final delta matrix and sample metadata
  delta_long <- delta_long %>%
    mutate(sample_id_new = paste(.data[[patient_col]], .data[[drug_col]], .data[[segment_col]], sep = "|"))
  
  delta_mat <- delta_long %>%
    dplyr::select(gene, sample_id_new, delta) %>%
    tidyr::pivot_wider(names_from = sample_id_new, values_from = delta) %>%
    tibble::column_to_rownames("gene") %>%
    as.matrix()
  
  sample_meta <- delta_long %>%
    distinct(sample_id_new, .data[[patient_col]], .data[[drug_col]], .data[[segment_col]]) %>%
    tibble::column_to_rownames("sample_id_new")
  
  list(delta_mat = delta_mat, sample_meta = sample_meta)
}







####################################################
####################################################
####################################################

# parity checker
```{r}
library(dplyr)
library(tidyr)
library(glue)

diagnosis_parity <- function(sobj, patient, treatment, timepoint, tissue,
                             check_key, check_values) {
  # 1) 메타데이터 + unique_id 만들기
  meta <- sobj@meta.data %>%
    mutate(
      unique_id = paste0({{patient}}, {{treatment}}, {{timepoint}}, {{tissue}}),
      .val_chr  = as.character(.data[[check_key]])  # factor 안전 처리
    )
  
  # 2) 필터: 각 unique_id 그룹에 check_values가 "모두" 실제 관측되었는지
  filt <- meta %>%
    group_by(unique_id) %>%
    filter(all(check_values %in% .val_chr)) %>%
    ungroup()
  
  # 3) 진단 테이블: 각 unique_id별로 check_values의 개수와 기타(n_other) 계산
  totals <- filt %>%
    count(unique_id, name = "total_n")
  
  counts_sel <- filt %>%
    filter(.val_chr %in% check_values) %>%
    count(unique_id, .val_chr) %>%
    tidyr::pivot_wider(names_from = .val_chr, values_from = n, values_fill = 0)
  
  diag <- totals %>%
    left_join(counts_sel, by = "unique_id") %>%
    # 혹시라도 존재하지 않는 컬럼이 생기면 0으로 채움
    mutate(across(all_of(intersect(check_values, names(.))), ~replace_na(.x, 0))) %>%
    # n_other = 전체 - 선택값들의 합
    rowwise() %>%
    mutate(
      n_selected = sum(c_across(all_of(check_values))[!is.na(c_across(all_of(check_values)))]),
      n_other    = total_n - n_selected,
      # has_all: 모든 check_values가 최소 1개 이상
      has_all    = all(c_across(all_of(check_values)) > 0, na.rm = TRUE),
      # parity: 선택한 값들의 개수가 서로 모두 동일(= pre/post가 같음, 3개면 모두 같음)
      parity_eq  = { 
        vals <- c_across(all_of(check_values))
        vals <- vals[!is.na(vals)]
        length(vals) > 0 && length(unique(vals)) == 1
      },
      ok = has_all && parity_eq && n_other == 0
    ) %>%
    ungroup()
  
  # 4) 패스 못 한 id 메시지
  no_pass_ids <- diag %>% filter(!ok) %>% pull(unique_id)
  msg <- if (length(no_pass_ids) == 0) {
    "✅ All groups passed parity & exclusivity checks."
  } else {
    glue("checkout {paste(no_pass_ids, collapse = ', ')}")
  }
  
  # 5) 결과 반환(진단표 + 메시지 + 필터된 데이터)
  list(
    diag_table = diag,
    filtered_meta = filt,
    message = msg
  )
}

diagnosis_parity_detailed <- function(sobj,
                                      patient=NULL, treatment=NULL, timepoint=NULL, tissue=NULL,   # 문자열 컬럼명 허용
                                      check_key, check_values) {
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
      n_selected = reduce(across(all_of(safe_counts_names)), `+`),
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
  print(msg)
  print("View($detailed)")
  print("View($diag_table)")
  ok_id=detailed%>%
    filter(ok)%>%
    pull(SegmentDisplayName)
  
  list(
    detailed    = detailed,    # 원본 행 + n_pre/n_post/... + ok 등
    diag_table  = diag_table,  # id별 한 줄 요약
    message     = msg,
    paired_ids = ok_id
  )
}
```


#6 ver 11 (main) QSG, RLMG -> LMM, etc. (Claude)


```{r}
library(tidyverse)
library(Seurat)
library(lme4)
library(lmerTest)
library(emmeans)
library(parallel)
library(broom.mixed)
library(ggpubr)
library(pheatmap)
library(viridis)

# ========================================
# PART 1: 데이터 준비 및 설정
# ========================================

#' 분석 설정 클래스 - 모든 변수명을 여기서 관리
create_analysis_config <- function(
    patient = "patient_id",
    drug = "drug", 
    timepoint = "timepoint",
    ck = "ck_status",
    response = "response",
    aoi = "aoi_id"
) {
  list(
    patient = patient,
    drug = drug,
    timepoint = timepoint,
    ck = ck,
    response = response,
    aoi = aoi
  )
}

#' GeoMx 데이터 준비 (Excel/CSV에서 읽기)
#' @param count_file 발현값 파일 경로
#' @param metadata_file 메타데이터 파일 경로  
#' @param config 변수명 설정
prepare_geomx_data <- function(count_file = NULL, 
                               metadata_file = NULL,
                               count_matrix = NULL,
                               metadata = NULL,
                               config = create_analysis_config(),
                               q3_normalize = TRUE) {
  
  # 파일에서 읽기 또는 직접 입력
  if (!is.null(count_file)) {
    if (grepl("\\.xlsx$", count_file)) {
      count_matrix <- openxlsx::read.xlsx(count_file, sheet = 1, rowNames = TRUE)
    } else {
      count_matrix <- read.csv(count_file, row.names = 1)
    }
  }
  
  if (!is.null(metadata_file)) {
    if (grepl("\\.xlsx$", metadata_file)) {
      metadata <- openxlsx::read.xlsx(metadata_file, sheet = 2)
    } else {
      metadata <- read.csv(metadata_file)
    }
  }
  
  # AOI 이름 매칭
  rownames(metadata) <- metadata[[config$aoi]]
  
  # Seurat 객체 생성
  seurat_obj <- CreateSeuratObject(
    counts = as.matrix(count_matrix),
    meta.data = metadata,
    min.cells = 0,
    min.features = 0
  )
  
  # Q3 정규화는 이미 되어있다고 가정
  if (!q3_normalize) {
    seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize")
  }
  
  # 조합 변수 생성 (빠른 스크리닝용)
  seurat_obj$condition_full <- paste0(
    seurat_obj@meta.data[[config$timepoint]], "_",
    seurat_obj@meta.data[[config$drug]], "_", 
    seurat_obj@meta.data[[config$response]]
  )
  
  seurat_obj$condition_simple <- paste0(
    seurat_obj@meta.data[[config$timepoint]], "_",
    seurat_obj@meta.data[[config$drug]]
  )
  
  return(seurat_obj)
}

# ========================================
# PART 2: 빠른 스크리닝 (Wilcoxon)
# ========================================

#' Step 1: 빠른 유전자 스크리닝
#' @param seurat_obj Seurat 객체
#' @param config 변수명 설정
#' @param comparison_type "pre_post", "drug", "response", "combined" 중 선택
#' @param use_adj p_val_adj 사용 여부 (FALSE면 p_val 사용)
quick_screen_genes <- function(seurat_obj,
                               config = create_analysis_config(),
                               comparison_type = "pre_post",
                               ck_subset = NULL,
                               min_pct = 0.1,
                               logfc_threshold = 0.25,
                               top_n = 1000,
                               use_adj = FALSE) {
  
  # CK 서브셋 필터링
  if (!is.null(ck_subset)) {
    seurat_obj <- subset(seurat_obj, 
                         subset = !!sym(config$ck) == ck_subset)
  }
  
  # 비교 그룹 설정
  if (comparison_type == "pre_post") {
    Idents(seurat_obj) <- seurat_obj@meta.data[[config$timepoint]]
  } else if (comparison_type == "drug") {
    Idents(seurat_obj) <- seurat_obj@meta.data[[config$drug]]
  } else if (comparison_type == "response") {
    Idents(seurat_obj) <- seurat_obj@meta.data[[config$response]]
  } else if (comparison_type == "combined") {
    Idents(seurat_obj) <- seurat_obj$condition_full
  }
  
  # FindAllMarkers 실행
  markers <- FindAllMarkers(
    seurat_obj,
    min.pct = min_pct,
    logfc.threshold = logfc_threshold,
    only.pos = FALSE,
    return.thresh = 1  # 모든 유전자 반환
  )
  
  # p-value 선택 및 상위 유전자 선택
  if (use_adj) {
    top_genes <- markers %>%
      arrange(p_val_adj) %>%
      head(top_n) %>%
      pull(gene) %>%
      unique()
  } else {
    top_genes <- markers %>%
      arrange(p_val) %>%
      head(top_n) %>%
      pull(gene) %>%
      unique()
  }
  
  return(list(
    markers = markers,
    top_genes = top_genes,
    seurat_obj = seurat_obj
  ))
}

# ========================================
# PART 3: Linear Mixed Model 분석
# ========================================

#' Step 2: 단일 유전자 LMM
fit_lmm_single_gene <- function(gene_expr,
                                metadata,
                                config = create_analysis_config(),
                                formula_components = NULL,
                                use_config_names = TRUE) {
  
  # 데이터프레임 생성
  df <- cbind(
    data.frame(Expression = gene_expr),
    metadata
  )
  
  # Formula 생성
  if (is.null(formula_components)) {
    # 기본값: config 변수명 사용
    fixed_effects <- c(config$drug, config$timepoint, config$response)
    interactions <- c(
      paste(config$drug, config$timepoint, sep = ":"),
      paste(config$drug, config$response, sep = ":"),
      paste(config$timepoint, config$response, sep = ":"),
      paste(config$drug, config$timepoint, config$response, sep = ":")
    )
    random <- paste0("(1|", config$patient, ")")
    
    formula_str <- paste0(
      "Expression ~ ",
      paste(c(fixed_effects, interactions), collapse = " + "),
      " + ", random
    )
  } else {
    # formula_components 제공 시
    if (use_config_names) {
      # config 매핑 사용
      mapping <- list(
        "patient" = config$patient,
        "drug" = config$drug,
        "timepoint" = config$timepoint,
        "response" = config$response,
        "ck" = config$ck
      )
      
      # fixed effects 매핑
      fixed_mapped <- sapply(formula_components$fixed, function(x) {
        if (x %in% names(mapping)) mapping[[x]] else x
      })
      
      # interactions 매핑
      if (!is.null(formula_components$interactions)) {
        interactions_mapped <- sapply(formula_components$interactions, function(x) {
          parts <- strsplit(x, ":")[[1]]
          mapped_parts <- sapply(parts, function(p) {
            if (p %in% names(mapping)) mapping[[p]] else p
          })
          paste(mapped_parts, collapse = ":")
        })
      } else {
        interactions_mapped <- NULL
      }
      
      # random effects 매핑
      random_mapped <- formula_components$random
      for (key in names(mapping)) {
        random_mapped <- gsub(key, mapping[[key]], random_mapped)
      }
      
      formula_str <- paste0(
        "Expression ~ ",
        paste(c(fixed_mapped, interactions_mapped), collapse = " + "),
        " + ", random_mapped
      )
    } else {
      # 직접 사용 (매핑 없이)
      fixed_effects <- formula_components$fixed
      interactions <- formula_components$interactions
      formula_str <- paste0(
        "Expression ~ ",
        paste(c(fixed_effects, interactions), collapse = " + "),
        " + ", formula_components$random
      )
    }
  }
  
  # 모델 적합
  tryCatch({
    # Factor 레벨 재정렬 (reference level 설정)
    if (config$drug %in% names(df)) {
      df[[config$drug]] <- relevel(as.factor(df[[config$drug]]), ref = levels(as.factor(df[[config$drug]]))[1])
    }
    if (config$response %in% names(df)) {
      df[[config$response]] <- relevel(as.factor(df[[config$response]]), ref = levels(as.factor(df[[config$response]]))[1])
    }
    
    model <- lmer(as.formula(formula_str), data = df, REML = FALSE)
    
    # 결과 추출
    coef_summary <- summary(model)$coefficients
    anova_result <- anova(model)
    
    # 모든 약물 대비 추출 (reference level 포함)
    all_drugs <- unique(df[[config$drug]])
    ref_drug <- levels(as.factor(df[[config$drug]]))[1]
    
    # 주요 효과 계산
    effects <- data.frame(
      term = rownames(coef_summary),
      estimate = coef_summary[, "Estimate"],
      std_error = coef_summary[, "Std. Error"],
      t_value = coef_summary[, "t value"],
      p_value = coef_summary[, "Pr(>|t|)"],
      ref_drug = ref_drug
    )
    
    return(list(
      model = model,
      effects = effects,
      anova = anova_result,
      converged = TRUE,
      formula = formula_str,
      all_drugs = all_drugs,
      ref_drug = ref_drug
    ))
    
  }, error = function(e) {
    return(list(
      converged = FALSE,
      error = e$message,
      formula = formula_str
    ))
  })
}

#' Step 3: 다중 유전자 LMM (병렬처리)
run_lmm_multiple_genes <- function(seurat_obj,
                                   genes = NULL,
                                   config = create_analysis_config(),
                                   formula_components = NULL,
                                   use_config_names = TRUE,
                                   n_cores = 16,
                                   verbose = TRUE) {
  
  if (is.null(genes)) {
    genes <- rownames(seurat_obj)[1:100]  # 기본값: 상위 100개
    warning("No genes specified. Using first 100 genes.")
  }
  
  # 발현 매트릭스와 메타데이터 추출
  expr_matrix <- GetAssayData(seurat_obj, slot = "data")
  metadata <- seurat_obj@meta.data
  
  # 진행상황 표시
  if (verbose) {
    message(sprintf("Running LMM for %d genes using %d cores...", 
                    length(genes), n_cores))
  }
  
  # 병렬 처리
  if (n_cores > 1) {
    cl <- makeCluster(n_cores)
    clusterEvalQ(cl, {
      library(lme4)
      library(lmerTest)
    })
    clusterExport(cl, c("fit_lmm_single_gene", "config", "metadata", 
                        "expr_matrix", "formula_components", "use_config_names"),
                  envir = environment())
    
    results <- parLapply(cl, genes, function(gene) {
      if (gene %in% rownames(expr_matrix)) {
        result <- fit_lmm_single_gene(
          gene_expr = as.numeric(expr_matrix[gene, ]),
          metadata = metadata,
          config = config,
          formula_components = formula_components,
          use_config_names = use_config_names
        )
        result$gene <- gene
        return(result)
      } else {
        return(list(gene = gene, converged = FALSE, 
                    error = "Gene not found"))
      }
    })
    
    stopCluster(cl)
  } else {
    # 순차 처리
    results <- lapply(genes, function(gene) {
      if (verbose && which(genes == gene) %% 100 == 0) {
        message(sprintf("  Processing gene %d/%d", 
                        which(genes == gene), length(genes)))
      }
      
      if (gene %in% rownames(expr_matrix)) {
        result <- fit_lmm_single_gene(
          gene_expr = as.numeric(expr_matrix[gene, ]),
          metadata = metadata,
          config = config,
          formula_components = formula_components
        )
        result$gene <- gene
        return(result)
      } else {
        return(list(gene = gene, converged = FALSE, 
                    error = "Gene not found"))
      }
    })
  }
  
  names(results) <- genes
  
  # 결과 요약
  summary_df <- summarize_lmm_results(results, config)
  
  return(list(
    raw_results = results,
    summary = summary_df,
    converged_genes = sum(sapply(results, function(x) x$converged)),
    total_genes = length(genes)
  ))
}

#' LMM 결과 요약; RLMG 내부에서 잘 작동하므로 건드리지 않아도 됨. 이 결과가 $summary로 나간다.
summarize_lmm_results <- function(lmm_results, config) {
  # 수렴한 모델만 처리
  converged <- lmm_results[sapply(lmm_results, function(x) x$converged)]
  
  if (length(converged) == 0) {
    warning("No models converged!")
    return(NULL)
  }
  
  # 모든 효과 수집
  all_effects <- do.call(rbind, lapply(names(converged), function(gene) {
    effects <- converged[[gene]]$effects
    effects$gene <- gene
    return(effects)
  }))
  
  # p-value 보정
  all_effects <- all_effects %>%
    group_by(term) %>%
    mutate(
      p_adj = p.adjust(p_value, method = "BH"),
      significant = p_adj < 0.05
    ) %>%
    ungroup()
  
  # Reference drug 정보 추가
  if ("ref_drug" %in% names(converged[[1]]$effects)) {
    ref_drugs <- unique(sapply(converged, function(x) x$ref_drug))
    attr(all_effects, "ref_drug") <- ref_drugs[1]
  }
  
  return(all_effects)
}

# ========================================
# PART 4: 특정 질문 분석
# ========================================

#' 질문 iv) 약물별 R/NR 효과가 큰 유전자 (reference drug 처리)
find_response_differential_genes <- function(lmm_summary,
                                             config,
                                             drug_name = NULL,
                                             top_n = 50) {
  
  # Reference drug 확인
  ref_drug <- attr(lmm_summary, "ref_drug")
  if (is.null(ref_drug) && "ref_drug" %in% names(lmm_summary)) {
    ref_drug <- unique(lmm_summary$ref_drug)[1]
  }
  
  if (!is.null(drug_name)) {
    if (drug_name == ref_drug) {
      # Reference drug인 경우: 다른 약물과의 대비 효과 확인
      message(sprintf("%s is the reference level. Showing contrasts with other drugs.", drug_name))
      
      # 다른 약물들의 main effect와 interaction 수집
      drug_terms <- lmm_summary %>%
        filter(grepl(config$drug, term) & grepl(config$response, term))
      
    } else {
      # Non-reference drug: 직접 효과 확인
      interaction_pattern <- paste0(drug_name, ".*:", ".*", config$response, "|",
                                    config$response, ".*:", ".*", drug_name)
      drug_terms <- lmm_summary %>%
        filter(grepl(interaction_pattern, term))
    }
  } else {
    # 모든 약물의 response interaction
    drug_terms <- lmm_summary %>%
      filter(grepl(paste0(config$drug, ".*:", ".*", config$response), term))
  }
  
  top_genes <- drug_terms %>%
    group_by(gene) %>%
    dplyr::summarize(
      max_effect = max(abs(estimate)),
      min_p = min(p_value),
      min_p_adj = min(p_adj),
      n_sig = sum(significant)
    ) %>%
    arrange(desc(max_effect)) %>%
    head(top_n)
  
  return(top_genes)
}

#' 질문 v) 약물 특이적 유전자
find_drug_specific_genes <- function(lmm_summary,
                                     config,
                                     response_subset = NULL,
                                     top_n = 50) {
  
  # Drug 주효과 및 상호작용
  drug_terms <- lmm_summary %>%
    filter(grepl(config$drug, term))
  
  if (!is.null(response_subset)) {
    drug_terms <- drug_terms %>%
      filter(grepl(response_subset, term))
  }
  
  top_genes <- drug_terms %>%
    group_by(gene) %>%
    dplyr::summarize(
      mean_effect = mean(abs(estimate)),
      max_effect = max(abs(estimate)),
      n_significant = sum(significant)
    ) %>%
    arrange(desc(max_effect)) %>%
    head(top_n)
  
  return(top_genes)
}

#' 질문 vi) Drug x Response 상호작용
find_interaction_genes <- function(lmm_summary,
                                   config,
                                   interaction_type = "drug_response",
                                   top_n = 50) {
  
  if (interaction_type == "drug_response") {
    pattern <- paste0(config$drug, ".*:", ".*", config$response)
  } else if (interaction_type == "drug_time") {
    pattern <- paste0(config$drug, ".*:", ".*", config$timepoint)
  } else if (interaction_type == "triple") {
    pattern <- paste0(config$drug, ".*:", ".*", 
                      config$timepoint, ".*:", ".*", config$response)
  }
  
  interaction_genes <- lmm_summary %>%
    filter(grepl(pattern, term)) %>%
    group_by(gene) %>%
    dplyr::summarize(
      interaction_strength = max(abs(estimate)),
      min_p_adj = min(p_adj),
      n_sig_terms = sum(significant)
    ) %>%
    arrange(desc(interaction_strength)) %>%
    head(top_n)
  
  return(interaction_genes)
}

# ========================================
# PART 5: 시각화
# ========================================

#' Pre-Post 변화 박스플롯
plot_pre_post_boxplot <- function(seurat_obj,
                                  gene,
                                  config,
                                  split_by = NULL) {
  
  # 데이터 준비
  df <- data.frame(
    Expression = GetAssayData(seurat_obj)[gene, ],
    seurat_obj@meta.data
  )
  
  # 기본 플롯
  p <- ggplot(df, aes_string(x = config$timepoint, y = "Expression")) +
    geom_boxplot(aes_string(fill = config$timepoint), alpha = 0.7) +
    geom_point(aes_string(fill = config$timepoint),
               position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.2), alpha = 0.5) +
    stat_compare_means(paired = TRUE, method = "wilcox.test") +
    labs(title = paste("Expression Changes:", gene),
         y = "Normalized Expression",
         x = "Timepoint") +
    theme_bw()
  
  # 분할 옵션
  if (!is.null(split_by)) {
    if (length(split_by) == 1) {
      p <- p + facet_wrap(as.formula(paste("~", split_by)))
    } else if (length(split_by) == 2) {
      p <- p + facet_grid(as.formula(paste(split_by[1], "~", split_by[2])))
    }
  }
  
  return(p)
}

#' 효과 크기 화산도
plot_volcano_legacy <- function(lmm_summary,
                                term_pattern = NULL,
                                title = "Volcano Plot",
                                effect_threshold = 0.5,
                                p_threshold = 0.05,
                                use_adj=FALSE) {
  
  plot_data <- lmm_summary
  
  if (!is.null(term_pattern)) {
    plot_data <- plot_data %>%
      filter(grepl(term_pattern, term))
  }
  
  if(use_adj){y_label="-log10(adj_p-value)"}else{y_label="-log10(p-value)"}
  
  
  plot_data <- plot_data %>%
    mutate(
      log_p = -log10(p_value),
      category = case_when(
        abs(estimate) > effect_threshold & p_adj < p_threshold & use_adj ~ "Significant",
        abs(estimate) > effect_threshold & p_value < p_threshold & !use_adj ~ "Significant",
        abs(estimate) > effect_threshold ~ "Large effect",
        p_adj  < p_threshold & use_adj ~ "Small effect",
        p_value  < p_threshold& !use_adj ~ "Small effect",
        TRUE ~ "Not significant"
      )
    )
  xlim_val=quantile(abs(plot_data$estimate),0.95,na.rm=TRUE)
  ggplot(plot_data, aes(x = estimate, y = log_p, color = category)) +
    geom_point(alpha = 0.6) +
    geom_hline(yintercept = -log10(p_threshold),
               linetype = "dashed", color = "gray") +
    geom_vline(xintercept = c(-effect_threshold, effect_threshold),
               linetype = "dashed", color = "gray") +
    scale_color_manual(values = c("Significant" = "red",
                                  "Large effect" = "orange",
                                  "Small effect" = "blue",
                                  "Not significant" = "gray")) +
    labs(title = title,
         x = "Effect Size (Estimate)",
         y = y_label) +
    theme_bw()+
    coord_cartesian(xlim = c(-xlim_val, xlim_val))
}

# ========================================
# PART 6: 완전한 워크플로우
# ========================================

#' 전체 분석 워크플로우
run_complete_workflow <- function(
    count_data,      # 매트릭스 또는 파일 경로
    metadata,        # 데이터프레임 또는 파일 경로
    config = create_analysis_config(),
    output_dir = "./results",
    n_cores = 16
) {
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ===== STEP 1: 데이터 준비 =====
  message("Step 1: Preparing data...")
  
  if (is.character(count_data)) {
    seurat_obj <- prepare_geomx_data(
      count_file = count_data,
      metadata_file = metadata,
      config = config
    )
  } else {
    seurat_obj <- prepare_geomx_data(
      count_matrix = count_data,
      metadata = metadata,
      config = config
    )
  }
  
  message(sprintf("  Loaded %d genes x %d AOIs", 
                  nrow(seurat_obj), ncol(seurat_obj)))
  
  # ===== STEP 2: CK별 분리 =====
  message("Step 2: Separating by CK status...")
  
  ck_values <- unique(seurat_obj@meta.data[[config$ck]])
  message(sprintf("  Found CK groups: %s", paste(ck_values, collapse = ", ")))
  
  results_by_ck <- list()
  
  for (ck_val in ck_values) {
    message(sprintf("\n=== Analyzing %s ===", ck_val))
    
    # CK 서브셋
    seurat_ck <- subset(seurat_obj, 
                        subset = !!sym(config$ck) == ck_val)
    
    # ===== STEP 3: 빠른 스크리닝 =====
    message("Step 3: Quick screening with Wilcoxon test...")
    
    screen_results <- quick_screen_genes(
      seurat_ck,
      config = config,
      comparison_type = "combined",
      top_n = 1000
    )
    
    message(sprintf("  Selected top %d genes", 
                    length(screen_results$top_genes)))
    
    # ===== STEP 4: LMM 분석 =====
    message("Step 4: Running Linear Mixed Models...")
    
    # 간단한 모델부터 시작 (config 키워드 사용)
    simple_formula <- list(
      fixed = c("drug", "timepoint", "response"),
      interactions = c("drug:timepoint"),
      random = "(1|patient)"
    )
    
    lmm_results <- run_lmm_multiple_genes(
      seurat_ck,
      genes = screen_results$top_genes,
      config = config,
      formula_components = simple_formula,
      use_config_names = TRUE,  # config 매핑 사용
      n_cores = n_cores
    )
    
    message(sprintf("  %d/%d models converged", 
                    lmm_results$converged_genes,
                    lmm_results$total_genes))
    
    # ===== STEP 5: 특정 분석 =====
    message("Step 5: Analyzing specific questions...")
    
    # 약물별 분석
    drugs <- unique(seurat_ck@meta.data[[config$drug]])
    drug_response_effects <- list()
    
    for (drug in drugs) {
      drug_response_effects[[drug]] <- find_response_differential_genes(
        lmm_results$summary,
        config,
        drug_name = drug,
        top_n = 30
      )
    }
    
    # Drug-specific genes
    drug_specific <- find_drug_specific_genes(
      lmm_results$summary,
      config,
      top_n = 50
    )
    
    # Interaction genes
    interaction_genes <- find_interaction_genes(
      lmm_results$summary,
      config,
      interaction_type = "drug_response",
      top_n = 50
    )
    
    # ===== STEP 6: 결과 저장 =====
    results_by_ck[[ck_val]] <- list(
      seurat_obj = seurat_ck,
      screening = screen_results,
      lmm = lmm_results,
      drug_response_effects = drug_response_effects,
      drug_specific = drug_specific,
      interaction_genes = interaction_genes
    )
    
    # CSV 저장
    write.csv(drug_specific,
              file.path(output_dir, paste0("drug_specific_", ck_val, ".csv")),
              row.names = FALSE)
    
    write.csv(interaction_genes,
              file.path(output_dir, paste0("interaction_genes_", ck_val, ".csv")),
              row.names = FALSE)
    
    # ===== STEP 7: 시각화 =====
    message("Step 6: Creating visualizations...")
    
    # Top 5 유전자 플롯
    top_genes_to_plot <- drug_specific$gene[1:min(5, nrow(drug_specific))]
    
    pdf(file.path(output_dir, paste0("plots_", ck_val, ".pdf")), 
        width = 12, height = 8)
    
    for (gene in top_genes_to_plot) {
      p <- plot_pre_post_boxplot(
        seurat_ck,
        gene = gene,
        config = config,
        split_by = c(config$drug, config$response)
      )
      print(p)
    }
    
    # 화산도
    p_volcano <- plot_volcano(
      lmm_results$summary,
      term_pattern = config$drug,
      title = paste("Drug Effects -", ck_val)
    )
    print(p_volcano)
    
    dev.off()
  }
  
  # 전체 결과 저장
  saveRDS(results_by_ck, 
          file.path(output_dir, "complete_results.rds"))
  
  message("\n=== Analysis Complete ===")
  message(sprintf("Results saved to: %s", output_dir))
  
  return(results_by_ck)
}



```
## patching
```{r}
plot_volcano <- function(lmm_summary,
                         x_col="estimate",
                         y_col="p_value",
                         filter_col="term",
                         filter_pattern = NULL,
                         title = "Volcano Plot",
                         effect_threshold = 0.5,
                         p_threshold = 0.05,
                         resize_x=TRUE,
                         resize_y=TRUE) {
  #ver2: generalized
  plot_data <- lmm_summary
  
  if (!is.null(filter_pattern)) {
    plot_data <- plot_data %>%
      filter(grepl(filter_pattern, term))
  }
  if (!is.null(x_col)){plot_data=plot_data%>%mutate(effect_size= !!sym(x_col))}
  if (!is.null(y_col)){plot_data=plot_data%>%mutate(p_value= !!sym(y_col))}
  
  # plot label definition
  x_label="Effect Size"
  y_label="Significance"
  
  # plot resizing
  if(resize_x){xlim_val=quantile(abs(plot_data$effect_size),0.95,na.rm=TRUE)}
  if(resize_y){ylim_val=quantile(abs(plot_data$effect_size),0.95,na.rm=TRUE)}
  
  # gene significance label arrange
  plot_data <- plot_data %>%
    mutate(
      log_p = -log10(p_value),
      category = case_when(
        abs(effect_size) > effect_threshold & p_value < p_threshold ~ "Significant",
        abs(effect_size) > effect_threshold ~ "Large effect",
        p_value  < p_threshold ~ "Small effect",
        TRUE ~ "Not significant"
      )
    )
  ggplot(plot_data, aes(x = effect_size, y = log_p, color = category)) +
    geom_point(alpha = 0.6) +
    geom_hline(yintercept = -log10(p_threshold), 
               linetype = "dashed", color = "gray") +
    geom_vline(xintercept = c(-effect_threshold, effect_threshold), 
               linetype = "dashed", color = "gray") +
    scale_color_manual(values = c("Significant" = "red",
                                  "Large effect" = "orange", 
                                  "Small effect" = "blue",
                                  "Not significant" = "gray")) +
    labs(title = title,
         x = x_label,
         y = y_label) +
    theme_bw()+
    coord_cartesian(xlim = c(-xlim_val, xlim_val))
}
```


# GPT Epithelial .... plotting
##11 function
```{r}
plot_interaction_for_gene_legacy <- function(sobj, gene, assay="RNA", layer="data",
                                             drug_col="drug", nested = FALSE, # ref_drug=NULL,
                                             time_col="treatment", time_levels=c("pre","post"),time_labels=c("pre","post"),
                                             segment_col="ck" ,subset_comp = NULL, 
                                             patient_col="emrid" #,optimizer="bobyqa"
) {
  # segment_col은 GeoMx를 상정.
  # subset_col은, segment_col에서 어떤 value를 가진 데이터만 subset할 건지,
  # nested는 drug/patient로 nesting하여 분석할 건지 (lmer ~ (1|drug/patient))
  expr=GetAssayData(sobj, assay=assay, layer=layer)
  meta=sobj@meta.data %>%
    mutate(
      time=factor(.data[[time_col]], levels=time_levels,labels=time_labels),
      drug=factor(.data[[drug_col]]),
      patient=factor(.data[[patient_col]]),
      comp=factor(.data[[segment_col]])
    )
  # 선택적 compartment 필터
  if (!is.null(subset_comp)) {
    keep <- meta$comp == subset_comp
  } else {
    keep <- rep(TRUE, nrow(meta))
  }
  y <- as.numeric(expr[gene, keep])
  df <- cbind.data.frame(y = y, meta[keep, , drop = FALSE])
  
  # 변동성/결측 체크
  if (all(is.na(y)) || sd(y, na.rm = TRUE) == 0) {
    stop(sprintf("Gene %s has no variance or all NA in the selected subset.", gene))
  }
  
  # 모델: nested = TRUE면 patient가 drug에 nested
  # (1|patient)도 괜찮지만 엄밀히는 patient가 drug에 속하므로 (1|drug/patient) 권장
  if (nested) {
    fit <- lmer(y ~ time * drug + (1 | drug/patient), data = df, REML = FALSE)
  } else {
    fit <- lmer(y ~ time * drug + (1 | patient), data = df, REML = FALSE)
  }
  
  # 2-1) emmeans로 시간*약제 조합의 추정 평균과 CI
  emm <- emmeans(fit, ~ time * drug)
  emm_df <- as.data.frame(emm) %>%
    mutate(time = factor(time, levels = c("pre","post")))
  
  p1 <- ggplot(emm_df, aes(x = time, y = emmean, group = drug, color = drug)) +
    geom_point(size = 3) +
    geom_line(size = 1) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1) +
    labs(
      title = sprintf("%s: time × drug (estimated marginal means ± 95%% CI)", gene),
      x = "Time", y = "Estimated mean expression"
    ) +
    theme_bw()
  
  # 2-2) 환자별 pre→post 스파게티 (실측치)
  # (개별 변화의 분산/이질성을 직관적으로 확인)
  df_plot <- df %>%
    dplyr::select(y, time, drug, patient) %>%
    mutate(time = factor(time, levels = c("pre","post")))
  
  p2 <- ggplot(df_plot, aes(x = time, y = y, group = patient, color = drug)) +
    geom_line(alpha = 0.35) +
    geom_point(size = 2, alpha = 0.7) +
    labs(
      title = sprintf("%s: individual pre→post (raw, by patient)", gene),
      x = "Time", y = "Expression (Q3/log scale)"
    ) +
    theme_bw()
  
  # 2-3) 약제별 Δ(post-pre) 비교 (emmeans 대비/대조)
  # 각 약제 내 post-pre contrast
  delta_by_drug <- contrast(emm, method = "revpairwise", by = "drug", adjust = "none") %>%
    # revpairwise(pre,post) → (post - pre)로 보려면 pairs(~ time | drug) 사용 권장
    # 아래는 명시적 방법:
    pairs(by = "drug", adjust = "none")
  # 더 직관적인 방식:
  deltas <- contrast(emmeans(fit, ~ time | drug), method = list(delta = c(-1, 1)))
  deltas_df <- as.data.frame(deltas)
  
  p3 <- ggplot(deltas_df, aes(x = drug, y = estimate, ymin = lower.CL, ymax = upper.CL, color = drug)) +
    geom_pointrange(position = position_dodge(width = 0.2)) +
    geom_hline(yintercept = 0, linetype = 2) +
    labs(
      title = sprintf("%s: Δ(post-pre) by drug (EMMean contrast ± 95%% CI)", gene),
      x = "Drug", y = "Δ(post - pre)"
    ) +
    theme_bw() +
    theme(legend.position = "none")
  
  # 결과와 객체 반환
  list(
    fit = fit,
    emm = emm_df,
    deltas = deltas_df,
    plot_emmeans = p1,
    plot_spaghetti = p2,
    plot_delta = p3
  )
}
```


## 13 function v2

📌 plot_interaction_for_gene v2.0
변경/추가 사항
ref_drug 인자 추가 → 기준 drug 레벨 지정 가능 (NULL이면 factor 첫 레벨).
optimizer 인자 추가 → lmer optimizer 선택 가능 (기본 "bobyqa").
factor 수준 체크 강화 (time, patient 최소 2수준).
confint() 사용하여 항상 lower.CL, upper.CL 포함.
plot_delta_heatmap 추가: 여러 유전자 결과(deltas_df)를 모아 heatmap 요약 가능.
deltas 데이터프레임에 gene 이름을 붙여서 list에 반환.
plot_delta_heatmap는 ComplexHeatmap 또는 ggplot2::geom_tile() 기반. (여기선 ggplot2 버전으로 제공)
```{r}
plot_interaction_for_gene <- function(
    sobj, gene, assay = "RNA", layer = "data",
    drug_col = "drug", nested = FALSE, ref_drug = NULL, 
    time_col = "treatment", time_levels = c("pre","post"), time_labels = c("pre","post"),
    segment_col = "ck",subset_comp = NULL,
    patient_col = "emrid",
    optimizer = "bobyqa"
) {
  # ==== Data prep ====
  expr <- GetAssayData(sobj, assay = assay, slot = layer)
  meta <- sobj@meta.data %>%
    mutate(
      time    = factor(.data[[time_col]], levels = time_levels, labels = time_labels),
      drug    = factor(.data[[drug_col]]),
      patient = factor(.data[[patient_col]])
    )
  if (!is.null(segment_col)){
    meta=meta%>%mutate(
      comp=.data[[segment_col]]
    )
  }else{subset_comp=NULL} #현재 상태로는 comp가 있을지 여부가 갈리므로 후속 분석에 안 좋을 가능성이 크다. 현재는 문제없음.
  
  if (!is.null(subset_comp)) {
    keep <- meta$comp == subset_comp
  } else {
    keep <- rep(TRUE, nrow(meta))
  }
  
  if (!is.null(ref_drug)) meta$drug <- relevel(meta$drug, ref = ref_drug)
  
  y <- as.numeric(expr[gene, keep])
  df <- cbind.data.frame(y = y, meta[keep, , drop = FALSE])
  
  # ==== Checks ====
  if (all(is.na(y)) || sd(y, na.rm = TRUE) == 0) {
    stop(sprintf("Gene %s has no variance or all NA in this subset.", gene))
  }
  if (nlevels(df$time) < 2) stop("time has <2 levels in this subset.")
  if (nlevels(df$patient) < 2) stop("patient has <2 levels in this subset.")
  
  # ==== Model ====
  ctrl <- lmerControl(optimizer = optimizer, calc.derivs = TRUE,
                      check.conv.singular = "ignore")
  if (nested) {
    fit <- lmer(y ~ time * drug + (1 | drug/patient),
                data = df, REML = FALSE, control = ctrl)
  } else {
    fit <- lmer(y ~ time * drug + (1 | patient),
                data = df, REML = FALSE, control = ctrl)
  }
  
  # ==== emmeans ====
  emm <- emmeans(fit, ~ time * drug)
  emm_df <- as.data.frame(confint(emm)) %>%
    mutate(time = factor(time, levels = time_levels, labels = time_labels))
  
  p1 <- ggplot(emm_df, aes(x = time, y = emmean, group = drug, color = drug)) +
    geom_point(size = 3) +
    geom_line(size = 1) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1) +
    labs(
      title = sprintf("%s: time × drug (EMMeans ±95%% CI)", gene),
      x = "Time", y = "Estimated mean expression"
    ) +
    theme_bw()
  
  # ==== spaghetti ====
  df_plot <- df %>%
    dplyr::select(y, time, drug, patient) %>%
    mutate(time = factor(time, levels = time_levels, labels = time_labels))
  
  p2 <- ggplot(df_plot, aes(x = time, y = y, group = interaction(patient, drug), color = drug)) +
    geom_line(alpha = 0.35) +
    geom_point(size = 2, alpha = 0.7) +
    labs(
      title = sprintf("%s: individual pre→post (raw)", gene),
      x = "Time", y = "Expression (Q3/log scale)"
    ) +
    theme_bw()
  
  # ==== delta (post-pre) ====
  em_time_by_drug <- emmeans(fit, ~ time | drug)
  deltas <- contrast(em_time_by_drug, method = list(delta = c(-1, 1)))
  deltas_df <- as.data.frame(confint(deltas)) %>%
    mutate(gene = gene)
  
  p3 <- ggplot(deltas_df, aes(x = drug, y = estimate,
                              ymin = lower.CL, ymax = upper.CL, color = drug)) +
    geom_pointrange(position = position_dodge(width = 0.2)) +
    geom_hline(yintercept = 0, linetype = 2) +
    labs(
      title = sprintf("%s: Δ(post-pre) by drug (EMMeans ±95%% CI)", gene),
      x = "Drug", y = expression(Delta~"(post - pre)")
    ) +
    theme_bw() +
    theme(legend.position = "none")
  
  # Return
  list(
    df = df,
    fit = fit,
    emm = emm_df,
    deltas = deltas_df,
    plot_emmeans = p1,
    plot_spaghetti = p2,
    plot_delta = p3
  )
}

plot_delta_heatmap <- function(results_list) {
  # results_list: list of outputs from plot_interaction_for_gene (여러 gene)
  all_deltas <- purrr::map_dfr(results_list, "deltas")
  
  # drug x gene 매트릭스
  mat <- all_deltas %>%
    dplyr::select(gene, drug, estimate) %>%
    tidyr::pivot_wider(names_from = drug, values_from = estimate)
  
  df_long <- all_deltas
  
  p <- ggplot(df_long, aes(x = drug, y = gene, fill = estimate)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
    labs(title = "Δ(post-pre) heatmap by drug",
         x = "Drug", y = "Gene", fill = expression(Delta)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  list(matrix = mat, plot = p, raw = df_long)
}
```

## 14 MDMP: delta calculator - patient by patient (cell-wise)
```{r}
make_delta_matrix_by_patient <- function(
    sobj,
    assay="RNA", layer="data",
    time_col="treatment", drug_col="drug", patient_col="emrid", segment_col="ck",
    time_levels=c("pre","post"), subset_comp=NULL,
    agg_fun = function(v) mean(v, na.rm=TRUE)
){
  X <- GetAssayData(sobj, assay=assay, layer=layer) # genes x AOIs
  meta <- sobj@meta.data %>%
    mutate(
      time    = factor(.data[[time_col]], levels=time_levels),
      drug    = factor(.data[[drug_col]]),
      patient = factor(.data[[patient_col]]),
      comp    = factor(.data[[segment_col]])
    )
  stopifnot(identical(colnames(X), rownames(meta)))
  
  if (!is.null(subset_comp)) {
    keep <- meta$comp == subset_comp
    X <- X[, keep, drop=FALSE]
    meta <- meta[keep, , drop=FALSE]
  }
  
  # meta column filtering
  gene_names <- rownames(X)
  
  run=0#
  while(length(which(gene_names%in%names(meta)))>0){
    num_meta_in_gene=which(names(meta)%in%gene_names)
    names(meta)[num_meta_in_gene]=paste0(names(meta)[num_meta_in_gene],"_re")
    run=run+1#
    print(run)#
    if(run>5)break#
  }
  # AOI별 표현을 long으로
  long <- as.data.frame(t(X)) %>%
    tibble::rownames_to_column("AOI") %>%
    cbind(meta, .)
  
  # pre/post를 같은 환자/약제/comp 안에서 집계 후 Δ = post − pre
  # 메모리 절약을 위해 gene-wise로 루프
  make_delta_one_gene <- function(g){
    dd <- long %>%
      dplyr::select(patient, drug, comp, time, !!sym(g)) %>%
      group_by(patient, drug, comp, time) %>%
      summarise(expr = agg_fun(.data[[g]]), .groups="drop") %>%
      tidyr::pivot_wider(names_from = time, values_from = expr) %>%
      mutate(delta = .data[[time_levels[2]]] - .data[[time_levels[1]]]) %>%
      dplyr::select(patient, drug, comp, delta)
    dd$gene <- g
    dd
  }
  
  delta_long <- purrr::map_dfr(gene_names, make_delta_one_gene)
  
  # 샘플 정의: 환자*약제*(comp)
  delta_long <- delta_long %>%
    mutate(sample_id = if (!is.null(subset_comp)) 
      paste(patient, drug, sep="|") else
        paste(patient, drug, comp, sep="|"))
  
  # gene x sample 행렬
  delta_mat <- delta_long %>%
    dplyr::select(gene, sample_id, delta) %>%
    tidyr::pivot_wider(names_from = sample_id, values_from = delta) %>%
    tibble::column_to_rownames("gene") %>%
    as.matrix()
  
  # 메타(샘플) 테이블
  sample_meta <- delta_long %>%
    distinct(sample_id, patient, drug, comp) %>%
    tibble::column_to_rownames("sample_id")
  
  list(delta_mat = delta_mat, sample_meta = sample_meta)
}

# # 예시: CK만 사용
# out <- make_delta_matrix_by_patient(
#   sobj = data_seurat, subset_comp = TRUE,
#   time_col="treatment", drug_col="drug", patient_col="emrid", segment_col="ck",
#   time_levels=c("pre","post")
# )
# 
# # Seurat로 올려서 HVG 선택/스케일/PCA/클러스터링
# delta_mat <- out$delta_mat                       # genes x samples
# sample_meta <- out$sample_meta                   # samples x {patient,drug,comp}
# 
# # Δ 값은 음수/양수 포함 → counts 슬롯 대신 data 슬롯에 넣는 게 안전
# obj_delta <- CreateSeuratObject(counts = matrix(0, nrow=nrow(delta_mat), ncol=ncol(delta_mat)))
# obj_delta[["DELTA"]] <- CreateAssayObject(counts = matrix(0, nrow=nrow(delta_mat), ncol=ncol(delta_mat)))
# DefaultAssay(obj_delta) <- "DELTA"
# obj_delta[["DELTA"]]@data <- delta_mat
# obj_delta <- AddMetaData(obj_delta, sample_meta[ colnames(obj_delta), , drop=FALSE ])
# 
# # HVG는 Δ-변이로 선별 (원하시면 전체 사용도 가능)
# obj_delta <- FindVariableFeatures(obj_delta, selection.method="vst", nfeatures=3000, assay="DELTA")
# 
# # Scale & PCA & UMAP
# obj_delta <- ScaleData(obj_delta, features=VariableFeatures(obj_delta), assay="DELTA", clip.max=3)
# obj_delta <- RunPCA(obj_delta, features=VariableFeatures(obj_delta), assay="DELTA")
# obj_delta <- FindNeighbors(obj_delta, dims=1:20)
# obj_delta <- FindClusters(obj_delta, resolution=0.6)
# obj_delta <- RunUMAP(obj_delta, dims=1:20)
```
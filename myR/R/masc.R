#' MASC pipeline orchestrator
#'
#' @description
#' Wraps the MASC (Mixed-effects Association testing for Single Cells) workflow
#' into reusable steps with configurable saving, force-run, and plotting behaviour.
#' MASC tests whether a specified covariate influences the membership of single cells
#' in any of multiple cellular subsets while accounting for technical confounds and
#' biological variation using logistic mixed-effects models.
#'
#' The implementation is based on the original MASC package from
#' https://github.com/immunogenomics/masc and adapted to work seamlessly with
#' Seurat objects and integrate into the myR package ecosystem.
#'
#' @param seurat_obj A `Seurat` object already loaded in memory. If `NULL`,
#'   `seurat_qs_path` must be provided.
#' @param seurat_qs_path Optional path to a `.qs` file containing a `Seurat`
#'   object. Used when `seurat_obj` is `NULL`.
#' @param cluster_var Character string; metadata column name containing cluster
#'   assignments for each cell.
#' @param contrast_var Character string; metadata column name containing the
#'   variable to be tested for association with cluster abundance (e.g., "status",
#'   "treatment"). Must be a factor.
#' @param random_effects Optional character vector; column names in metadata to be
#'   modeled as random effects (e.g., c("donor_id", "batch")).
#' @param fixed_effects Optional character vector; column names in metadata to be
#'   modeled as fixed effects (e.g., c("age", "sex")).
#' @param save Logical; if `TRUE` (default) intermediate objects are written to
#'   disk using `qs::qsave()`.
#' @param output_dir,prefix,suffix Configure save location. When `suffix` is
#'   `NULL`, a numeric suffix (`"", "_1", ...`) is auto-generated to avoid
#'   overwriting existing files.
#' @param force_run Logical scalar or named logical vector controlling whether
#'   each step (`preparation`, `analysis`) is recomputed even if cached files
#'   are present.
#' @param plotting Logical; if `TRUE` (default), generates visualization plots.
#' @param save_models Logical; if `TRUE`, saves the mixed-effects model objects
#'   for each cluster. Default is `FALSE`.
#' @param adjust_pvalue Logical; if `TRUE` (default), applies FDR correction to
#'   p-values using `p.adjust(method = "fdr")`.
#' @param verbose Logical; emits progress messages when `TRUE`.
#'
#' @return A list with:
#'   - `masc_results`: data frame with association p-values, OR, and CI for each cluster
#'   - `models`: (optional) list of model objects for each cluster
#'   - `plots`: (optional) list of ggplot objects; `NULL` when `plotting = FALSE`
#'
#' @references
#' Fonseka et al. (2018) MASC: Mixed-effects Association testing for Single Cells.
#' https://github.com/immunogenomics/masc
#'
#' @export
run_masc_pipeline <- function(
    seurat_obj = NULL,
    seurat_qs_path = NULL,
    cluster_var,
    contrast_var,
    random_effects = NULL,
    fixed_effects = NULL,
    save = TRUE,
    output_dir = file.path(tempdir(), "masc"),
    prefix = "masc",
    suffix = NULL,
    force_run = FALSE,
    plotting = TRUE,
    save_models = FALSE,
    adjust_pvalue = TRUE,
    verbose = TRUE
) {
    .masc_require_packages()

    # Resolve Seurat object
    input_seurat_obj <- seurat_obj
    seurat_obj <- .masc_prepare_seurat(
        seurat_obj = input_seurat_obj,
        seurat_qs_path = seurat_qs_path,
        verbose = verbose
    )

    force_flags <- .masc_normalize_force_flags(force_run)
    paths <- .masc_resolve_paths(
        output_dir = output_dir,
        prefix = prefix,
        suffix = suffix,
        save = save
    )

    if (save) {
        dir.create(paths$output_dir, recursive = TRUE, showWarnings = FALSE)
    }

    # Validate metadata columns
    .masc_validate_metadata(
        seurat_obj = seurat_obj,
        cluster_var = cluster_var,
        contrast_var = contrast_var,
        random_effects = random_effects,
        fixed_effects = fixed_effects
    )

    # ---- Step 1: Prepare data ----
    has_cached_data <- save && file.exists(paths$files$data)
    if (!force_flags["preparation"] && has_cached_data) {
        if (verbose) {
            if (verbose) cat(sprintf("Loading cached data from %s\n", paths$files$data))
        }
        masc_data <- qs::qread(paths$files$data)
    } else {
        if (verbose) cat("Preparing data for MASC analysis.\n")
        masc_data <- .masc_prepare_data(
            seurat_obj = seurat_obj,
            cluster_var = cluster_var,
            contrast_var = contrast_var,
            random_effects = random_effects,
            fixed_effects = fixed_effects,
            verbose = verbose
        )
        if (save) {
            qs::qsave(masc_data, paths$files$data)
            if (verbose) cat(sprintf("Saved data to %s\n", paths$files$data))
        }
    }

    # ---- Step 2: Run MASC analysis ----
    has_cached_results <- save && file.exists(paths$files$results)
    if (!force_flags["analysis"] && has_cached_results) {
        if (verbose) {
            if (verbose) cat(sprintf("Loading cached MASC results from %s\n", paths$files$results))
        }
        masc_results <- qs::qread(paths$files$results)
        cluster_models <- NULL
        if (save_models && file.exists(paths$files$models)) {
            cluster_models <- qs::qread(paths$files$models)
        }
    } else {
        if (verbose) cat("Running MASC analysis.\n")
        masc_bundle <- .masc_run_analysis(
            dataset = masc_data$dataset,
            cluster = masc_data$cluster,
            contrast = contrast_var,
            random_effects = random_effects,
            fixed_effects = fixed_effects,
            save_models = save_models,
            save_model_dir = if (save_models) paths$output_dir else NULL,
            verbose = verbose
        )
        masc_results <- masc_bundle$results
        cluster_models <- masc_bundle$models

        # Apply FDR correction if requested
        if (adjust_pvalue && "model.pvalue" %in% names(masc_results)) {
            masc_results$model.pvalue.fdr <- p.adjust(masc_results$model.pvalue, method = "fdr")
            if (verbose) {
                n_sig <- sum(masc_results$model.pvalue.fdr < 0.05, na.rm = TRUE)
                if (verbose) cat(sprintf("FDR correction applied: %d clusters with FDR < 0.05\n", n_sig))
            }
        }

        if (save) {
            qs::qsave(masc_results, paths$files$results)
            if (verbose) cat(sprintf("Saved MASC results to %s\n", paths$files$results))
            if (save_models && !is.null(cluster_models)) {
                qs::qsave(cluster_models, paths$files$models)
                if (verbose) cat(sprintf("Saved model objects to %s\n", paths$files$models))
            }
        }
    }

    # ---- Step 3: Generate plots ----
    plots <- NULL
    if (plotting) {
        if (verbose) cat("Generating MASC plots.\n")
        plots <- .masc_plot_bundle(
            masc_results = masc_results,
            cluster_var = cluster_var,
            contrast_var = contrast_var,
            save = save,
            save_path = paths$files$plots,
            verbose = verbose
        )
    }

    # Return results
    result_list <- list(
        masc_results = masc_results
    )
    if (save_models && !is.null(cluster_models)) {
        result_list$models <- cluster_models
    }
    if (!is.null(plots)) {
        result_list$plots <- plots
    }

    result_list
}

# ==============================================================================
# Helper Functions
# ==============================================================================

.masc_require_packages <- function() {
    pkgs <- c("Seurat", "lme4", "Matrix", "qs", "cli", "ggplot2", "dplyr")
    missing <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
    if (length(missing) > 0) {
        stop("Missing required packages: ", paste(missing, collapse = ", "),
             ". Install with: install.packages(c('", paste(missing, collapse = "', '"), "'))")
    }
}

.masc_prepare_seurat <- function(seurat_obj, seurat_qs_path, verbose) {
    if (!is.null(seurat_obj)) {
        if (!inherits(seurat_obj, "Seurat")) {
            stop("`seurat_obj` must be a Seurat object.")
        }
        return(seurat_obj)
    }
    if (!is.null(seurat_qs_path)) {
        if (!file.exists(seurat_qs_path)) {
            stop(sprintf("File not found: %s", seurat_qs_path))
        }
        if (verbose) cat(sprintf("Loading Seurat object from %s\n", seurat_qs_path))
        return(qs::qread(seurat_qs_path))
    }
    stop("Either `seurat_obj` or `seurat_qs_path` must be provided.")
}

.masc_normalize_force_flags <- function(force_run) {
    # Initialize flags
    flags <- c(preparation = FALSE, analysis = FALSE)
    
    if (is.logical(force_run) && length(force_run) == 1) {
        flags["preparation"] <- force_run
        flags["analysis"] <- force_run
    } else if (is.logical(force_run) && length(force_run) > 1 && !is.null(names(force_run))) {
        if ("preparation" %in% names(force_run)) flags["preparation"] <- force_run["preparation"]
        if ("analysis" %in% names(force_run)) flags["analysis"] <- force_run["analysis"]
        # Handle default/all if present
        if ("all" %in% names(force_run)) {
            flags["preparation"] <- force_run["all"]
            flags["analysis"] <- force_run["all"]
        }
    }
    return(flags)
}

.masc_resolve_paths <- function(output_dir, prefix, suffix, save) {
    if (!save) {
        return(list(
            output_dir = output_dir,
            files = list()
        ))
    }

    base_paths <- list(
        data = file.path(output_dir, paste0(prefix, "_data.qs")),
        results = file.path(output_dir, paste0(prefix, "_results.qs")),
        models = file.path(output_dir, paste0(prefix, "_models.qs")),
        plots = file.path(output_dir, paste0(prefix, "_plots.qs"))
    )

    # Auto-generate suffix if needed
    if (is.null(suffix)) {
        suffix <- ""
        counter <- 1
        while (any(file.exists(unlist(base_paths)))) {
            suffix <- paste0("_", counter)
            base_paths <- list(
                data = file.path(output_dir, paste0(prefix, suffix, "_data.qs")),
                results = file.path(output_dir, paste0(prefix, suffix, "_results.qs")),
                models = file.path(output_dir, paste0(prefix, suffix, "_models.qs")),
                plots = file.path(output_dir, paste0(prefix, suffix, "_plots.qs"))
            )
            counter <- counter + 1
        }
    } else {
        suffix_str <- if (suffix == "") "" else paste0("_", suffix)
        base_paths <- list(
            data = file.path(output_dir, paste0(prefix, suffix_str, "_data.qs")),
            results = file.path(output_dir, paste0(prefix, suffix_str, "_results.qs")),
            models = file.path(output_dir, paste0(prefix, suffix_str, "_models.qs")),
            plots = file.path(output_dir, paste0(prefix, suffix_str, "_plots.qs"))
        )
    }

    list(
        output_dir = output_dir,
        suffix = suffix,
        files = base_paths
    )
}

.masc_validate_metadata <- function(seurat_obj, cluster_var, contrast_var, random_effects, fixed_effects) {
    meta_cols <- colnames(seurat_obj@meta.data)
    required_vars <- c(cluster_var, contrast_var)
    if (!is.null(random_effects)) required_vars <- c(required_vars, random_effects)
    if (!is.null(fixed_effects)) required_vars <- c(required_vars, fixed_effects)

    missing <- setdiff(required_vars, meta_cols)
    if (length(missing) > 0) {
        stop("Missing metadata columns: ", paste(missing, collapse = ", "))
    }

    # Check that contrast_var is a factor
    if (!is.factor(seurat_obj@meta.data[[contrast_var]])) {
        warning(sprintf("Converting '%s' to factor. Please ensure levels are ordered correctly.",
                       contrast_var), immediate. = TRUE)
        seurat_obj@meta.data[[contrast_var]] <- as.factor(seurat_obj@meta.data[[contrast_var]])
    }

    # Check factor levels
    contrast_levels <- levels(seurat_obj@meta.data[[contrast_var]])
    if (length(contrast_levels) < 2) {
        stop(sprintf("Contrast variable '%s' must have at least 2 levels, found: %d",
                    contrast_var, length(contrast_levels)))
    }

    # Check cluster assignment
    cluster_values <- unique(seurat_obj@meta.data[[cluster_var]])
    if (length(cluster_values) < 2) {
        stop(sprintf("Cluster variable '%s' must have at least 2 unique values, found: %d",
                    cluster_var, length(cluster_values)))
    }
}

.masc_prepare_data <- function(seurat_obj, cluster_var, contrast_var, random_effects, fixed_effects, verbose) {
    # Extract metadata
    dataset <- as.data.frame(seurat_obj@meta.data)
    cluster <- dataset[[cluster_var]]

    # Ensure cluster is character
    cluster <- as.character(cluster)

    # Ensure contrast is factor
    if (!is.factor(dataset[[contrast_var]])) {
        dataset[[contrast_var]] <- as.factor(dataset[[contrast_var]])
    }

    # Verify required columns
    all_vars <- c(contrast_var)
    if (!is.null(random_effects)) all_vars <- c(all_vars, random_effects)
    if (!is.null(fixed_effects)) all_vars <- c(all_vars, fixed_effects)

    # Remove any rows with NA in required variables
    complete_cases <- complete.cases(dataset[, all_vars, drop = FALSE])
    if (sum(!complete_cases) > 0 && verbose) {
        warning(sprintf("Removing %d rows with missing values in required variables", sum(!complete_cases)))
    }
    dataset <- dataset[complete_cases, all_vars, drop = FALSE]
    cluster <- cluster[complete_cases]

    # Convert random effects to factors if they are character
    if (!is.null(random_effects)) {
        for (var in random_effects) {
            if (var %in% colnames(dataset)) {
                if (is.character(dataset[[var]]) || is.numeric(dataset[[var]])) {
                    dataset[[var]] <- as.factor(dataset[[var]])
                }
            }
        }
    }

    # Convert fixed effects - numeric should stay numeric, character should become factor
    if (!is.null(fixed_effects)) {
        for (var in fixed_effects) {
            if (var %in% colnames(dataset)) {
                # Skip if already properly formatted
                if (is.numeric(dataset[[var]])) {
                    next  # Keep numeric variables as numeric
                }
                
                if (is.character(dataset[[var]]) || is.logical(dataset[[var]])) {
                    # Convert to factor, handling NA
                    var_values <- dataset[[var]]
                    var_factor <- as.factor(var_values)
                    
                    # Check if too many levels (might be problematic)
                    n_levels <- length(levels(var_factor))
                    if (n_levels > 20 && verbose) {
                        warning(sprintf("Variable '%s' has %d levels, which may cause model fitting issues", var, n_levels))
                    }
                    
                    dataset[[var]] <- var_factor
                    if (verbose && n_levels <= 20) {
                        if (verbose) cat(sprintf("Converted '%s' to factor with %d levels\n", var, n_levels))
                    }
                } else if (!is.factor(dataset[[var]])) {
                    # Try to convert to numeric
                    dataset[[var]] <- tryCatch({
                        as.numeric(dataset[[var]])
                    }, warning = function(w) {
                        as.factor(dataset[[var]])
                    }, error = function(e) {
                        as.factor(dataset[[var]])
                    })
                }
            }
        }
    }

    if (verbose) {
        n_cells <- nrow(dataset)
        n_clusters <- length(unique(cluster))
        n_contrast_levels <- length(levels(dataset[[contrast_var]]))
        if (verbose) cat(sprintf("Prepared data: %d cells, %d clusters, %d contrast levels\n", n_cells, n_clusters, n_contrast_levels))
    }

    list(dataset = dataset, cluster = cluster)
}

.masc_run_analysis <- function(dataset, cluster, contrast, random_effects, fixed_effects,
                                save_models, save_model_dir, verbose) {
    # This is the core MASC function, adapted from the original implementation
    # Generate design matrix from cluster assignments
    cluster <- as.character(cluster)
    
    # Ensure cluster and dataset have same length
    if (length(cluster) != nrow(dataset)) {
        stop(sprintf("Length mismatch: cluster has %d elements, dataset has %d rows", 
                     length(cluster), nrow(dataset)))
    }
    
    # Create cluster data frame for model.matrix
    cluster_df <- data.frame(cluster = cluster, stringsAsFactors = FALSE)
    
    if (verbose) {
        cat(sprintf("Creating design matrix for %d clusters\n", length(unique(cluster))))
    }
    
    designmat <- tryCatch({
        mm <- model.matrix(~ cluster + 0, data = cluster_df)
        if (verbose) {
            cat(sprintf("Design matrix created: %d x %d\n", nrow(mm), ncol(mm)))
        }
        mm
    }, error = function(e) {
        stop(sprintf("Failed to create design matrix: %s", e$message))
    })
    
    # Clean column names - remove "cluster" prefix if present
    original_colnames <- colnames(designmat)
    colnames(designmat) <- gsub("^cluster", "", colnames(designmat))
    
    if (verbose && any(original_colnames != colnames(designmat))) {
        cat(sprintf("Cleaned column names: %s -> %s\n", 
                       paste(head(original_colnames, 2), collapse = ", "),
                       paste(head(colnames(designmat), 2), collapse = ", ")))
    }
    
    dataset <- cbind(designmat, dataset)

    # Create output list to hold results - use colnames instead of attributes
    cluster_names <- colnames(designmat)
    res <- vector(mode = "list", length = length(cluster_names))
    names(res) <- cluster_names

    # Create model formulas
    if (!is.null(fixed_effects) && !is.null(random_effects)) {
        model_rhs <- paste0(c(paste0(fixed_effects, collapse = " + "),
                              paste0("(1|", random_effects, ")", collapse = " + ")),
                            collapse = " + ")
        if (verbose) {
            cat(sprintf("Using model: cluster ~ %s + %s\n", contrast, model_rhs))
        }
    } else if (!is.null(fixed_effects) && is.null(random_effects)) {
        model_rhs <- paste0(fixed_effects, collapse = " + ")
        if (verbose) {
            cat(sprintf("Using model: cluster ~ %s + %s\n", contrast, model_rhs))
        }
        stop("MASC requires at least one random effect term")
    } else if (is.null(fixed_effects) && !is.null(random_effects)) {
        model_rhs <- paste0("(1|", random_effects, ")", collapse = " + ")
        if (verbose) {
            cat(sprintf("Using model: cluster ~ %s + %s\n", contrast, model_rhs))
        }
    } else {
        model_rhs <- "1"
        if (verbose) {
            cat(sprintf("Using model: cluster ~ %s + %s\n", contrast, model_rhs))
        }
        stop("MASC requires at least one random effect term")
    }

    # Initialize list to store model objects
    cluster_models <- vector(mode = "list", length = length(cluster_names))
    names(cluster_models) <- cluster_names

    # Run nested mixed-effects models for each cluster
    for (i in seq_along(cluster_names)) {
        test_cluster <- cluster_names[i]
        if (verbose) {
            cat(sprintf("Creating logistic mixed models for %s\n", test_cluster))
        }

        # Ensure test_cluster column exists in dataset
        if (!test_cluster %in% colnames(dataset)) {
            if (verbose) {
                warning(sprintf("Cluster column '%s' not found in dataset. Available: %s", 
                               test_cluster, paste(head(colnames(dataset), 10), collapse = ", ")))
            }
            cluster_models[[i]] <- NULL
            next
        }
        
        null_fm_str <- paste0(c(paste0(test_cluster, " ~ 1 + "),
                               model_rhs), collapse = "")
        full_fm_str <- paste0(c(paste0(test_cluster, " ~ ", contrast, " + "),
                               model_rhs), collapse = "")
        
        if (verbose) {
            cat(sprintf("Formula for %s - Null: %s\n", test_cluster, null_fm_str))
        }
        
        null_fm <- tryCatch({
            as.formula(null_fm_str)
        }, error = function(e) {
            stop(sprintf("Failed to create null formula for %s: %s", test_cluster, e$message))
        })
        
        full_fm <- tryCatch({
            as.formula(full_fm_str)
        }, error = function(e) {
            stop(sprintf("Failed to create full formula for %s: %s", test_cluster, e$message))
        })

        # Run null and full mixed-effects models
        null_model <- tryCatch({
            result <- lme4::glmer(formula = null_fm, data = dataset,
                                 family = binomial, nAGQ = 1, verbose = 0,
                                 control = lme4::glmerControl(optimizer = "bobyqa"))
            result
        }, error = function(e) {
            if (verbose) {
                warning(sprintf("Null model failed for cluster %s: %s", test_cluster, e$message))
            }
            NULL
        }, warning = function(w) {
            # Suppress convergence warnings, but continue
            if (verbose && grepl("convergence|singular", w$message, ignore.case = TRUE)) {
                # Just continue - model may still be usable
            }
            NULL
        })

        full_model <- tryCatch({
            result <- lme4::glmer(formula = full_fm, data = dataset,
                                 family = binomial, nAGQ = 1, verbose = 0,
                                 control = lme4::glmerControl(optimizer = "bobyqa"))
            result
        }, error = function(e) {
            if (verbose) {
                warning(sprintf("Full model failed for cluster %s: %s", test_cluster, e$message))
            }
            NULL
        }, warning = function(w) {
            # Suppress convergence warnings
            if (verbose && grepl("convergence|singular", w$message, ignore.case = TRUE)) {
                # Just continue
            }
            NULL
        })

        if (!is.null(null_model) && !is.null(full_model)) {
            # Perform LRT
            model_lrt <- NULL
            tryCatch({
                model_lrt <- anova(null_model, full_model)
            }, error = function(e) {
                if (verbose) {
                    warning(sprintf("LRT failed for cluster %s: %s", test_cluster, e$message))
                }
                model_lrt <<- NULL
            })
            
            if (is.null(model_lrt)) {
                cluster_models[[i]] <- NULL
                next
            }
            
            # Calculate confidence intervals for contrast term beta
            contrast_levels <- levels(dataset[[contrast]])
            if (length(contrast_levels) < 2) {
                if (verbose) {
                    warning(sprintf("Contrast variable has < 2 levels for cluster %s", test_cluster))
                }
                cluster_models[[i]] <- NULL
                next
            }
            
            contrast_lvl2 <- paste0(contrast, contrast_levels[2])
            
            # Check if contrast term exists in model
            coef_names <- names(lme4::fixef(full_model))
            if (!contrast_lvl2 %in% coef_names) {
                if (verbose) {
                    warning(sprintf("Contrast term %s not found in model for cluster %s", contrast_lvl2, test_cluster))
                }
                contrast_ci <- matrix(NA, nrow = 1, ncol = 2, 
                                     dimnames = list(contrast_lvl2, c("2.5 %", "97.5 %")))
            } else {
                contrast_ci <- tryCatch({
                    lme4::confint.merMod(full_model, method = "Wald", parm = contrast_lvl2)
                }, error = function(e) {
                    if (verbose) {
                        warning(sprintf("CI calculation failed for %s: %s", test_cluster, e$message))
                    }
                    matrix(NA, nrow = 1, ncol = 2, 
                          dimnames = list(contrast_lvl2, c("2.5 %", "97.5 %")))
                })
            }

            cluster_models[[i]]$null_model <- null_model
            cluster_models[[i]]$full_model <- full_model
            cluster_models[[i]]$model_lrt <- model_lrt
            cluster_models[[i]]$confint <- contrast_ci
        } else {
            cluster_models[[i]] <- NULL
        }
    }

    # Organize results into output dataframe
    output <- data.frame(
        cluster = cluster_names,
        size = colSums(designmat),
        stringsAsFactors = FALSE
    )

    contrast_levels <- levels(dataset[[contrast]])
    contrast_lvl2 <- paste0(contrast, contrast_levels[2])

    output$model.pvalue <- sapply(cluster_models, function(x) {
        if (is.null(x) || is.null(x$model_lrt)) return(NA_real_)
        x$model_lrt[["Pr(>Chisq)"]][2]
    })

    output[[paste(contrast_lvl2, "OR", sep = ".")]] <- sapply(cluster_models, function(x) {
        if (is.null(x) || is.null(x$full_model)) return(NA_real_)
        coef_val <- lme4::fixef(x$full_model)[[contrast_lvl2]]
        if (is.na(coef_val)) return(NA_real_)
        exp(coef_val)
    })

    output[[paste(contrast_lvl2, "OR", "95pct.ci.lower", sep = ".")]] <- sapply(cluster_models, function(x) {
        if (is.null(x) || is.null(x$confint)) return(NA_real_)
        ci_val <- x$confint[contrast_lvl2, "2.5 %"]
        if (is.na(ci_val)) return(NA_real_)
        exp(ci_val)
    })

    output[[paste(contrast_lvl2, "OR", "95pct.ci.upper", sep = ".")]] <- sapply(cluster_models, function(x) {
        if (is.null(x) || is.null(x$confint)) return(NA_real_)
        ci_val <- x$confint[contrast_lvl2, "97.5 %"]
        if (is.na(ci_val)) return(NA_real_)
        exp(ci_val)
    })

    list(results = output, models = if (save_models) cluster_models else NULL)
}

.masc_plot_bundle <- function(masc_results, cluster_var, contrast_var, save, save_path, verbose) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        if (verbose) warning("ggplot2 not available; skipping plots.")
        return(NULL)
    }

    plots <- list()

    # Plot 1: OR with confidence intervals (forest plot style)
    if (any(grepl("\\.OR$", names(masc_results)))) {
        or_col <- grep("\\.OR$", names(masc_results), value = TRUE)[1]
        ci_lower_col <- grep("ci\\.lower", names(masc_results), value = TRUE)[1]
        ci_upper_col <- grep("ci\\.upper", names(masc_results), value = TRUE)[1]

        if (!is.null(or_col) && !is.null(ci_lower_col) && !is.null(ci_upper_col)) {
            plot_df <- masc_results[, c("cluster", or_col, ci_lower_col, ci_upper_col)]
            names(plot_df) <- c("cluster", "OR", "CI_lower", "CI_upper")
            plot_df <- plot_df[!is.na(plot_df$OR), ]

            plots$or_forest <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$cluster, y = .data$OR)) +
                ggplot2::geom_point() +
                ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$CI_lower, ymax = .data$CI_upper),
                                      width = 0.2) +
                ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
                ggplot2::coord_flip() +
                ggplot2::labs(
                    title = "MASC Results: Odds Ratios",
                    x = "Cluster",
                    y = "Odds Ratio (95% CI)"
                ) +
                ggplot2::theme_minimal()
        }
    }

    # Plot 2: P-value bar plot
    if ("model.pvalue" %in% names(masc_results)) {
        plot_df <- masc_results[, c("cluster", "model.pvalue")]
        plot_df <- plot_df[order(plot_df$model.pvalue), ]
        plot_df$cluster <- factor(plot_df$cluster, levels = plot_df$cluster)

        plots$pvalue_bar <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$cluster, y = -log10(.data$model.pvalue))) +
            ggplot2::geom_bar(stat = "identity") +
            ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
            ggplot2::coord_flip() +
            ggplot2::labs(
                title = "MASC Results: P-values",
                x = "Cluster",
                y = "-log10(p-value)"
            ) +
            ggplot2::theme_minimal()
    }

    if (save && length(plots) > 0) {
        qs::qsave(plots, save_path)
        if (verbose) cat(sprintf("Saved plots to %s\n", save_path))
    }

    if (length(plots) == 0) return(NULL)
    plots
}

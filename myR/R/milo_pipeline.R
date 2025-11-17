#' Milo pipeline orchestrator
#'
#' @description
#' Wraps the Milo differential abundance workflow into reusable steps with
#' configurable saving, force-run, and plotting behaviour. The implementation is
#' distilled from the `st/KDW_251110.Rmd` notebook (lines 3476-3868) and
#' incorporates lessons learned during manual debugging (see `context_Korean.md`).
#'
#' @param seurat_obj A `Seurat` object already loaded in memory. If `NULL`,
#'   `seurat_qs_path` must be provided.
#' @param seurat_qs_path Optional path to a `.qs` file containing a `Seurat`
#'   object. Used when `seurat_obj` is `NULL`.
#' @param patient_var,cluster_var,target_var,batch_var Metadata column names
#'   required for Milo DA testing.
#' @param graph_reduction Reduction slot used to build Milo neighbourhoods.
#' @param layout_reduction Reduction slot used to populate the Milo UMAP layout.
#' @param k,d,prop Parameters forwarded to `buildGraph()` / `makeNhoods()`.
#' @param alpha Default significance threshold used across plots.
#' @param cell_metric Metric used to colour cells in UMAP (supports
#'   `"SpatialFDR"`, `"PValue"`, `"logFC_percentile"`).
#' @param beeswarm_metric Metric added to `da_results` for colouring the
#'   beeswarm plot. Same options as `cell_metric`.
#' @param beeswarm_alpha Threshold passed to `plotDAbeeswarm()`.
#' @param fdr_breaks,fdr_labels Bin specification for colouring cells by metric.
#' @param target_levels Optional character vector that fixes the factor order of
#'   `target_var` (e.g. `c("low", "high")`). The first level becomes the
#'   reference (logFC < 0 ⇒ enrichment there) and the second level becomes the
#'   test (logFC > 0 ⇒ enrichment there).
#' @param target_include Optional character vector. When supplied, only cells
#'   whose `target_var` matches one of the values are kept before running Milo.
#' @param save Logical; if `TRUE` (default) intermediate objects are written to
#'   disk using `qs::qsave()`.
#' @param output_dir,prefix,suffix Configure save location. When `suffix` is
#'   `NULL`, a numeric suffix (`"", "_1", ...`) is auto-generated to avoid
#'   overwriting existing files.
#' @param cache_files Optional named list that overrides auto-generated cache
#'   paths. Supported names: `nhoods`, `distances`, `da_milo`, `da_results`,
#'   `plots`. This is useful when you want to pin a specific `.qs` file instead
#'   of relying on the default suffix-based naming.
#' @param force_run Logical scalar or named logical vector controlling whether
#'   each step (`nhoods`, `distances`, `testing`) is recomputed even if cached
#'   files are present.
#' @param plotting Logical; if `FALSE`, plotting helpers are skipped.
#' @param max_cells Optional integer. When set and lower than the total number
#'   of cells, a random subset is used to accelerate development runs.
#' @param seed Integer passed to `set.seed()` before sampling or Milo routines.
#' @param verbose Logical; emits progress messages when `TRUE`.
#'
#' @return A list with:
#'   - `milo`: Milo object with `metadata(milo)$commands[['.milo_run_da']]`
#'     capturing the parameters, cache files, and comparison labels used.
#'   - `da_results`: data frame returned by `miloR::testNhoods()` plus
#'     `comparison_reference`, `comparison_test`, `comparison_label`,
#'     `enriched_in` columns that spell out the logFC direction.
#'   - `plots`: list of ggplot objects (`NULL` when `plotting = FALSE`),
#'     cached as `.qs` when `save = TRUE`.
#'
#' @export
run_milo_pipeline <- function(
    seurat_obj = NULL,
    seurat_qs_path = NULL,
    patient_var,
    cluster_var,
    target_var,
    batch_var,
    graph_reduction = "integrated.scvi",
    layout_reduction = "umap.scvi",
    k = 30,
    d = 30,
    prop = 0.1,
    alpha = 0.1,
    cell_metric = c("SpatialFDR", "PValue", "logFC_percentile"),
    beeswarm_metric = c("SpatialFDR", "PValue", "logFC_percentile"),
    beeswarm_alpha = alpha,
    fdr_breaks = c(0, 0.1, 0.2, 0.3, 0.5, 1),
    fdr_labels = c("< 0.1", "0.1 - 0.2", "0.2 - 0.3", "0.3 - 0.5", ">= 0.5"),
    target_levels = NULL,
    target_include = NULL,
    save = TRUE,
    output_dir = file.path(tempdir(), "milo"),
    prefix = "milo",
    suffix = NULL,
    cache_files = NULL,
    force_run = FALSE,
    plotting = TRUE,
    max_cells = NULL,
    seed = 1L,
    verbose = TRUE
) {
    .milo_require_packages()

    input_seurat_obj <- seurat_obj
    seurat_supplier <- local({
        prepared <- NULL
        function() {
            if (!is.null(prepared)) return(prepared)
            prepared <<- .milo_prepare_seurat(
                seurat_obj = input_seurat_obj,
                seurat_qs_path = seurat_qs_path,
                max_cells = max_cells,
                seed = seed,
                verbose = verbose
            )
            prepared
        }
    })
    seurat_obj <- if (!is.null(input_seurat_obj)) seurat_supplier() else NULL
    set.seed(seed)

    force_flags <- .milo_normalize_force_flags(force_run)
    paths <- .milo_resolve_paths(
        output_dir = output_dir,
        prefix = prefix,
        suffix = suffix,
        save = save
    )
    paths <- .milo_apply_cache_overrides(paths, cache_files, verbose)

    if (save) {
        dir.create(paths$output_dir, recursive = TRUE, showWarnings = FALSE)
    }

    # Check cache status for all steps
    cache_status <- list(
        nhoods = save && !force_flags["nhoods"] && file.exists(paths$files$nhoods),
        distances = save && !force_flags["distances"] && file.exists(paths$files$distances),
        testing = save && !force_flags["testing"] && file.exists(paths$files$da_milo) && file.exists(paths$files$da_results)
    )
    
    if (verbose && save) {
        cache_msg <- c("Cache status:")
        for (step in names(cache_status)) {
            status <- if (cache_status[[step]]) "✓" else "✗"
            cache_msg <- c(cache_msg, paste0("  ", step, ": ", status))
        }
        cli::cli_inform(cache_msg)
    }

    # If all steps are cached, load final results directly
    if (cache_status$testing) {
        if (verbose) {
            cli::cli_inform(c("step" = "All steps cached. Loading final results from {.path {paths$files$da_results}}"))
        }
        milo <- qs::qread(paths$files$da_milo)
        da_results <- qs::qread(paths$files$da_results)
    } else {
        if (!is.null(seurat_obj)) {
            if (verbose) cli::cli_inform(c("Preparing metadata for Milo conversion."))
            seurat_obj <- .milo_prepare_metadata(seurat_obj, patient_var, cluster_var, target_var, batch_var)
        }

        # ---- Step 1: Build Milo object with neighbourhoods ----
        d_eff <- NULL  # Will be set after building nhoods
        if (cache_status$nhoods) {
            if (verbose) {
                cli::cli_inform(c("step" = "Loading cached Milo object with neighbourhoods from {.path {paths$files$nhoods}}"))
            }
            milo <- qs::qread(paths$files$nhoods)
            # Try to recover d_eff from the cached object
            if ("GRAPH" %in% SingleCellExperiment::reducedDimNames(milo)) {
                d_eff <- ncol(SingleCellExperiment::reducedDim(milo, "GRAPH"))
            }
        } else {
            if (verbose) cli::cli_inform(c("step" = "Building Milo neighbourhoods."))
            seurat_obj <- seurat_supplier()
            if (is.null(seurat_obj)) {
                stop("Seurat object required to build Milo neighbourhoods. Provide `seurat_obj` or `seurat_qs_path`.")
            }
            seurat_obj <- .milo_prepare_metadata(seurat_obj, patient_var, cluster_var, target_var, batch_var)
            environment(seurat_supplier)$prepared <- seurat_obj
            build_result <- .milo_build_nhoods(
                seurat_obj = seurat_obj,
                graph_reduction = graph_reduction,
                k = k,
                d = d,
                prop = prop,
                verbose = verbose
            )
            milo <- build_result$milo
            d_eff <- build_result$d_eff
            if (save) {
                qs::qsave(milo, paths$files$nhoods)
                if (verbose) {
                    cli::cli_inform(c("cache" = "Saved Milo object to {.path {paths$files$nhoods}}"))
                }
            }
        }
        
        # Ensure d_eff is available
        if (is.null(d_eff)) {
            if ("GRAPH" %in% SingleCellExperiment::reducedDimNames(milo)) {
                d_eff <- ncol(SingleCellExperiment::reducedDim(milo, "GRAPH"))
            } else {
                # Fallback: use the provided d parameter
                d_eff <- d
                if (verbose) {
                    cli::cli_warn(paste0("Could not determine effective dimension from Milo object. Using provided d=", d, "."))
                }
            }
        }

        # ---- Step 2: Distances ----
        if (cache_status$distances) {
            if (verbose) {
                cli::cli_inform(c("step" = "Loading cached Milo object with neighbourhood distances from {.path {paths$files$distances}}"))
            }
            milo <- qs::qread(paths$files$distances)
        } else {
            if (verbose) {
                cli::cli_inform(c("step" = "Calculating Milo neighbourhood distances (this may take a while)."))
            }
            milo <- miloR::calcNhoodDistance(milo, d = d_eff, reduced.dim = "GRAPH")
            if (save) {
                qs::qsave(milo, paths$files$distances)
                if (verbose) cli::cli_inform(c("cache" = "Saved Milo object to {.path {paths$files$distances}}"))
            }
        }

        # ---- Step 3: DA testing ----
        da_results <- NULL
        if (cache_status$testing) {
            if (verbose) {
                cli::cli_inform(c("step" = "Loading cached Milo DA results from {.path {paths$files$da_results}}"))
            }
            milo <- qs::qread(paths$files$da_milo)
            da_results <- qs::qread(paths$files$da_results)
        } else {
            if (verbose) cli::cli_inform(c("step" = "Running Milo differential abundance test."))
            da_bundle <- .milo_run_da(
                milo = milo,
                patient_var = patient_var,
                target_var = target_var,
                batch_var = batch_var,
                cluster_var = cluster_var,
                target_levels = target_levels,
                target_include = target_include,
                cache_files_used = paths$files,
                verbose = verbose
            )
            milo <- da_bundle$milo
            da_results <- da_bundle$da_results
            if (save) {
                qs::qsave(milo, paths$files$da_milo)
                qs::qsave(da_results, paths$files$da_results)
                if (verbose) {
                    cli::cli_inform(c("cache" = "Saved Milo object and DA results to {.path {paths$files$da_milo}} / {.path {paths$files$da_results}}"))
                }
            }
        }
    }

    plots <- NULL
    if (plotting) {
        if (verbose) cli::cli_inform(c("step" = "Preparing Milo plots."))
        if (!"UMAP" %in% SingleCellExperiment::reducedDimNames(milo)) {
            seurat_obj <- seurat_supplier()
        }
        beeswarm_alpha <- .milo_normalize_alpha(beeswarm_alpha, name = "beeswarm_alpha")

        plots <- .milo_plot_bundle(
            milo = milo,
            da_results = da_results,
            seurat_obj = seurat_obj,
            layout_reduction = layout_reduction,
            cluster_var = cluster_var,
            alpha = alpha,
            cell_metric = match.arg(cell_metric),
            beeswarm_metric = match.arg(beeswarm_metric),
            beeswarm_alpha = beeswarm_alpha,
            fdr_breaks = fdr_breaks,
            fdr_labels = fdr_labels,
            save = save,
            save_path = paths$files$plots,
            verbose = verbose
        )
    }

    list(
        milo = milo,
        da_results = da_results,
        plots = plots
    )
}

.milo_require_packages <- function() {
    pkgs <- c("Seurat", "SingleCellExperiment", "miloR", "Matrix", "qs", "cli", "ggplot2", "patchwork", "scater")
    missing <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
    if (length(missing) > 0) {
        stop("Missing required packages: ", paste(missing, collapse = ", "))
    }
}

.milo_set_nhood_counts <- function(milo, counts) {
    milo@nhoodCounts <- counts
    milo
}

.milo_log_command <- function(milo, command_name, args_list) {
    if (is.null(args_list$timestamp)) {
        args_list$timestamp <- format(Sys.time(), tz = "UTC", usetz = TRUE)
    }
    meta <- S4Vectors::metadata(milo)
    if (is.null(meta) || !is.list(meta)) {
        meta <- list()
    }
    commands <- meta$commands
    if (is.null(commands) || !is.list(commands)) {
        commands <- list()
    }
    commands[[command_name]] <- args_list
    meta$commands <- commands
    S4Vectors::metadata(milo) <- meta
    milo
}

.milo_prepare_metadata <- function(sobj, patient_var, cluster_var, target_var, batch_var) {
    meta_vars <- c(patient_var, cluster_var, target_var, batch_var)
    missing <- meta_vars[!meta_vars %in% colnames(sobj@meta.data)]
    if (length(missing) > 0) {
        stop("Missing metadata columns in Seurat object: ", paste(missing, collapse = ", "))
    }
    for (var in meta_vars) {
        sobj@meta.data[[var]] <- as.character(sobj@meta.data[[var]])
    }
    sobj
}

.milo_build_nhoods <- function(seurat_obj, graph_reduction, k, d, prop, verbose = TRUE) {
    if (!graph_reduction %in% Seurat::Reductions(seurat_obj)) {
        stop(sprintf("Reduction '%s' not found in Seurat object.", graph_reduction))
    }

    G_emb <- Seurat::Embeddings(seurat_obj, reduction = graph_reduction)
    n_cells <- nrow(G_emb)
    n_dims <- ncol(G_emb)
    
    if (n_cells == 0 || n_dims == 0) {
        stop(sprintf("Reduction '%s' has invalid dimensions: %d cells x %d dimensions", 
                     graph_reduction, n_cells, n_dims))
    }
    
    # Validate k parameter
    if (k < 2) {
        stop(sprintf("k (%d) must be at least 2 for building k-nearest neighbor graph.", k))
    }
    if (k >= n_cells) {
        warning(sprintf("k (%d) is >= number of cells (%d). Setting k to %d.", 
                       k, n_cells, max(2, n_cells - 1)))
        k <- max(2, n_cells - 1)
    }
    
    d_eff <- min(n_dims, d)
    if (d_eff < 2) {
        stop(sprintf("Effective dimension (%d) is too small. Reduction '%s' has %d dimensions, but at least 2 are required.", 
                     d_eff, graph_reduction, n_dims))
    }
    
    if (verbose) {
        if (d_eff < d) {
            cli::cli_inform(c("info" = paste0("Using effective dimension ", d_eff, " (reduction has ", n_dims, " dimensions, requested d=", d, ").")))
        } else {
            cli::cli_inform(c("info" = paste0("Using effective dimension ", d_eff, " (requested d=", d, ").")))
        }
    }
    
    G_use <- G_emb[, seq_len(d_eff), drop = FALSE]

    sce <- Seurat::as.SingleCellExperiment(seurat_obj, assay = Seurat::DefaultAssay(seurat_obj))
    SingleCellExperiment::reducedDim(sce, "GRAPH") <- as.matrix(G_use)

    milo <- miloR::Milo(sce)
    # Use d_eff instead of d to match the actual embedding dimensions
    milo <- miloR::buildGraph(milo, k = k, d = d_eff, reduced.dim = "GRAPH")
    miloR::makeNhoods(milo, prop = prop, k = k, refined = TRUE, reduced_dims = "GRAPH")
    
    list(milo = milo, d_eff = d_eff)
}

.milo_run_da <- function(
    milo,
    patient_var,
    target_var,
    batch_var,
    cluster_var,
    target_levels = NULL,
    target_include = NULL,
    cache_files_used = NULL,
    verbose = TRUE
) {
    if (verbose) {
        cli::cli_inform(c("info" = "Running countCells..."))
    }
    
    meta_df <- as.data.frame(SingleCellExperiment::colData(milo))

    target_include_effective <- NULL
    target_mask <- rep(TRUE, nrow(meta_df))
    if (!is.null(target_include)) {
        include_vals <- intersect(target_include, unique(meta_df[[target_var]]))
        if (!length(include_vals)) {
            stop(
                sprintf(
                    "target_include (%s) does not match any values in '%s'.",
                    paste(target_include, collapse = ", "),
                    target_var
                )
            )
        }
        dropped <- setdiff(target_include, include_vals)
        if (length(dropped) > 0 && verbose) {
            cli::cli_warn(c(
                "warning" = paste0(
                    "Dropped target levels not present in data: ",
                    paste(dropped, collapse = ", ")
                )
            ))
        }
        target_mask <- meta_df[[target_var]] %in% include_vals
        if (!any(target_mask)) {
            stop("No cells left after applying target_include filter.")
        }
        target_include_effective <- include_vals
    }
    milo <- miloR::countCells(
        milo,
        samples = patient_var,
        meta.data = meta_df
    )

    if (verbose) {
        cli::cli_inform(c("info" = "Creating sample design matrix..."))
    }
    
    keep_cols <- c(patient_var, target_var, batch_var)
    missing_cols <- setdiff(keep_cols, colnames(meta_df))
    if (length(missing_cols) > 0) {
        stop("Missing required metadata columns: ", paste(missing_cols, collapse = ", "))
    }

    sample_ids <- colnames(miloR::nhoodCounts(milo))
    if (length(sample_ids) == 0) {
        stop("countCells returned zero samples. Check patient_var input.")
    }

    duplicated_samples <- any(duplicated(meta_df[[patient_var]]))
    if (duplicated_samples && verbose) {
        cli::cli_warn(c("warning" = paste0("patient_var '", patient_var, "' has duplicates. Using first occurrence for each sample.")))
    }

    sample_design_list <- lapply(sample_ids, function(sid) {
        rows <- which(meta_df[[patient_var]] == sid)
        if (length(rows) == 0) {
            stop(sprintf("Sample ID '%s' not found in metadata column '%s'.", sid, patient_var))
        }
        row_data <- meta_df[rows[1], keep_cols, drop = FALSE]
        row_data$sample_id <- sid
        row_data
    })

    sample_design <- do.call(rbind, sample_design_list)
    rownames(sample_design) <- sample_design$sample_id
    sample_design$sample_id <- NULL

    if (!is.null(target_include_effective)) {
        sample_keep <- sample_design[[target_var]] %in% target_include_effective
        sample_design <- sample_design[sample_keep, , drop = FALSE]
        if (!nrow(sample_design)) {
            stop("No samples remain after applying target_include filter.")
        }
        counts <- miloR::nhoodCounts(milo)
        keep_cols <- colnames(counts) %in% rownames(sample_design)
        if (!any(keep_cols)) {
            stop("No matching samples found between nhoodCounts and filtered target levels.")
        }
        counts <- counts[, keep_cols, drop = FALSE]
        milo <- .milo_set_nhood_counts(milo, counts)
    }

    sample_ids <- colnames(miloR::nhoodCounts(milo))

    sample_design[[batch_var]] <- factor(sample_design[[batch_var]])

    sample_design[[target_var]] <- factor(sample_design[[target_var]])
    target_levels_used <- levels(sample_design[[target_var]])

    if (!is.null(target_levels)) {
        present_levels <- intersect(target_levels, target_levels_used)
        if (!length(present_levels)) {
            stop(
                sprintf(
                    "target_levels (%s) not found in '%s'.",
                    paste(target_levels, collapse = ", "),
                    target_var
                )
            )
        }
        missing_levels <- setdiff(target_levels, present_levels)
        if (length(missing_levels) > 0 && verbose) {
            cli::cli_warn(c(
                "warning" = paste0(
                    "Levels not present in design and ignored: ",
                    paste(missing_levels, collapse = ", ")
                )
            ))
        }
        other_levels <- setdiff(target_levels_used, present_levels)
        sample_design[[target_var]] <- factor(
            sample_design[[target_var]],
            levels = c(present_levels, other_levels)
        )
    } else {
        sample_design[[target_var]] <- droplevels(sample_design[[target_var]])
    }

    target_levels_used <- levels(sample_design[[target_var]])
    if (length(target_levels_used) < 2) {
        stop("target_var must have at least two levels after filtering.")
    }
    if (length(target_levels_used) > 2 && verbose) {
        cli::cli_warn(c(
            "warning" = paste0(
                "target_var has ",
                length(target_levels_used),
                " levels; logFC reflects contrasts relative to the reference level '",
                target_levels_used[1],
                "'. Consider setting target_include/target_levels for binary comparisons."
            )
        ))
    }
    comparison_reference <- target_levels_used[1]
    comparison_test <- if (length(target_levels_used) >= 2) target_levels_used[2] else NA_character_

    if (verbose) {
        cli::cli_inform(c("info" = "Running testNhoods..."))
    }
    
    formula_str <- sprintf("~ %s + %s", batch_var, target_var)
    fdr_weighting_used <- "graph-overlap"
    da_results <- tryCatch({
        miloR::testNhoods(
            milo,
            design = as.formula(formula_str),
            design.df = sample_design,
            reduced.dim = "GRAPH",
            fdr.weighting = "graph-overlap"
        )
    }, error = function(e) {
        if (grepl("median\\(f75\\)|missing value|NA", conditionMessage(e), ignore.case = TRUE)) {
            if (verbose) {
                cli::cli_warn("TMM normalization failed. Trying with fdr.weighting='k-distance'...")
            }
            fdr_weighting_used <<- "k-distance"
            result_k <- tryCatch({
                miloR::testNhoods(
                    milo,
                    design = as.formula(formula_str),
                    design.df = sample_design,
                    reduced.dim = "GRAPH",
                    fdr.weighting = "k-distance"
                )
            }, error = function(e2) {
                if (verbose) {
                    cli::cli_warn("Still failing. Trying without fdr.weighting...")
                }
                fdr_weighting_used <<- "none"
                miloR::testNhoods(
                    milo,
                    design = as.formula(formula_str),
                    design.df = sample_design,
                    reduced.dim = "GRAPH"
                )
            })
            return(result_k)
        } else {
            stop(e)
        }
    })

    if (verbose) {
        cli::cli_inform(c("info" = "Adding cluster annotations..."))
    }
    
    nhood_matrix <- miloR::nhoods(milo)
    cell_clusters <- SingleCellExperiment::colData(milo)[[cluster_var]]
    majority_cluster <- vapply(
        seq_len(ncol(nhood_matrix)),
        function(i) {
            idx <- as.logical(nhood_matrix[, i])
            labs <- cell_clusters[idx]
            if (!length(labs)) return(NA_character_)
            tbl <- table(labs)
            names(tbl)[which.max(tbl)]
        },
        FUN.VALUE = character(1)
    )
    da_results$major_cluster <- majority_cluster[da_results$Nhood]

    if (!is.null(comparison_reference) && !is.null(comparison_test)) {
        da_results$comparison_reference <- comparison_reference
        da_results$comparison_test <- comparison_test
        da_results$comparison_label <- paste(comparison_test, "vs", comparison_reference)
        da_results$enriched_in <- ifelse(
            is.na(da_results$logFC),
            NA_character_,
            ifelse(da_results$logFC > 0, comparison_test, comparison_reference)
        )
    } else {
        da_results$comparison_reference <- NA_character_
        da_results$comparison_test <- NA_character_
        da_results$comparison_label <- NA_character_
        da_results$enriched_in <- NA_character_
    }

    milo <- .milo_log_command(
        milo,
        ".milo_run_da",
        list(
            patient_var = patient_var,
            target_var = target_var,
            batch_var = batch_var,
            cluster_var = cluster_var,
            reduced_dim = "GRAPH",
            fdr_weighting = fdr_weighting_used,
            sample_count = length(sample_ids),
            duplicates_in_patient_var = duplicated_samples,
            target_include = target_include_effective,
            target_levels_requested = target_levels,
            target_levels_used = target_levels_used,
            comparison_reference = comparison_reference,
            comparison_test = comparison_test,
            cache_files = cache_files_used
        )
    )

    list(milo = milo, da_results = da_results)
}

.milo_plot_bundle <- function(
    milo,
    da_results,
    seurat_obj,
    layout_reduction,
    cluster_var,
    alpha,
    cell_metric,
    beeswarm_metric,
    beeswarm_alpha,
    fdr_breaks,
    fdr_labels,
    save,
    save_path,
    verbose
) {
    milo <- .milo_ensure_nhood_graph(milo, verbose = verbose)
    milo <- .milo_attach_umap(milo, seurat_obj, layout_reduction, verbose = verbose)

    da_results <- .milo_prepare_metric_column(da_results, metric = beeswarm_metric)
    milo <- .milo_attach_cell_metric(
        milo = milo,
        da_results = da_results,
        metric = cell_metric,
        breaks = fdr_breaks,
        labels = fdr_labels,
        verbose = verbose
    )

    beeswarm_payload <- .milo_prepare_beeswarm_inputs(
        da_results = da_results,
        preferred_metric = beeswarm_metric,
        fallback_metric = "PValue",
        alpha = beeswarm_alpha,
        verbose = verbose
    )

    graph_plot <- miloR::plotNhoodGraphDA(milo, da_results, alpha = alpha) +
        ggplot2::ggtitle("DA Results (logFC)")

    beeswarm_plot <- miloR::plotDAbeeswarm(
        beeswarm_payload$results,
        group.by = "major_cluster",
        alpha = beeswarm_payload$alpha
    ) + ggplot2::ggtitle(beeswarm_payload$title)

    umap_plot <- scater::plotReducedDim(
        milo,
        dimred = "UMAP",
        colour_by = "milo_metric_bin"
    ) +
        ggplot2::scale_colour_manual(values = .milo_metric_palette(fdr_labels), name = "Cell metric") +
        ggplot2::ggtitle(sprintf("Cells by %s", cell_metric))

    combined <- patchwork::wrap_plots(umap_plot, graph_plot)

    plots <- list(
        graph = graph_plot,
        beeswarm = beeswarm_plot,
        umap = umap_plot,
        combined = combined
    )

    if (save) {
        qs::qsave(plots, file = save_path)
        if (verbose) cli::cli_inform(c("cache" = "Saved plot bundle to {.path {save_path}}"))
    }

    plots
}

.milo_prepare_metric_column <- function(da_results, metric) {
    metric <- match.arg(metric, c("SpatialFDR", "PValue", "logFC_percentile"))
    if (metric == "logFC_percentile") {
        if (!"logFC" %in% names(da_results)) {
            stop("`da_results` lacks 'logFC' needed for percentile metric.")
        }
        ranks <- rank(abs(da_results$logFC), na.last = "keep", ties.method = "average")
        da_results$logFC_percentile <- ranks / max(ranks, na.rm = TRUE)
    }
    da_results
}

.milo_attach_cell_metric <- function(milo, da_results, metric, breaks, labels, verbose) {
    metric <- match.arg(metric, c("SpatialFDR", "PValue", "logFC_percentile"))
    if (!metric %in% names(da_results)) {
        stop(sprintf("Metric '%s' not available in `da_results`.", metric))
    }

    metric_values <- da_results[[metric]]
    metric_values[is.na(metric_values)] <- max(breaks, na.rm = TRUE)

    N_matrix <- miloR::nhoods(milo)
    if (verbose) {
        dims <- if (!is.null(dim(N_matrix))) paste(dim(N_matrix), collapse = " x ") else length(N_matrix)
        cli::cli_inform(c("debug" = paste0("Computing cell-level metric from nhood matrix (class: ",
                                           paste(class(N_matrix), collapse = ", "),
                                           "; size: ", dims, ").")))
    }

    nz_idx <- Matrix::summary(N_matrix)
    nz_idx <- data.frame(i = nz_idx$i, j = nz_idx$j, stringsAsFactors = FALSE)
    if (nrow(nz_idx) == 0) {
        warning("Neighbourhood matrix has no memberships; defaulting all cell metrics to max(breaks).", call. = FALSE)
        cell_metric <- rep(max(breaks, na.rm = TRUE), nrow(N_matrix))
    } else {
        nz_idx$metric <- metric_values[nz_idx$j]
        min_metric_per_cell <- tapply(nz_idx$metric, nz_idx$i, min, na.rm = TRUE)
        cell_metric <- rep(max(breaks, na.rm = TRUE), nrow(N_matrix))
        cell_metric[as.integer(names(min_metric_per_cell))] <- min_metric_per_cell
    }

    metric_bins <- cut(
        cell_metric,
        breaks = breaks,
        labels = labels,
        right = TRUE,
        include.lowest = TRUE
    )

    # Update colData using SummarizedExperiment setter
    cd <- SummarizedExperiment::colData(milo)
    cd$milo_metric_value <- cell_metric
    cd$milo_metric_bin <- metric_bins
    SummarizedExperiment::colData(milo) <- cd

    if (verbose) {
        bin_table <- table(metric_bins, useNA = "ifany")
        summary_text <- paste(sprintf("%s=%d", names(bin_table), bin_table), collapse = ", ")
        cli::cli_inform(paste("Cell metric bins:", summary_text))
    }
    milo
}

.milo_prepare_beeswarm_inputs <- function(da_results, preferred_metric, fallback_metric, alpha, verbose) {
    preferred_metric <- match.arg(preferred_metric, c("SpatialFDR", "PValue", "logFC_percentile"))
    fallback_metric <- match.arg(fallback_metric, c("SpatialFDR", "PValue", "logFC_percentile"))

    colour_by <- preferred_metric
    metric_vals <- da_results[[preferred_metric]]
    metric_numeric <- is.numeric(metric_vals)

    needs_fallback <- !metric_numeric ||
        all(is.na(metric_vals)) ||
        (preferred_metric %in% c("SpatialFDR", "PValue") && all(metric_vals >= alpha, na.rm = TRUE))

    if (needs_fallback) {
        if (verbose) {
            cli::cli_warn(c("warning" = paste0("Beeswarm metric '", preferred_metric,
                                               "' unsuitable (non-numeric, all NA, or all ≥ alpha). Falling back to '",
                                               fallback_metric, "'.")))
        }
        if (fallback_metric == "logFC_percentile" && !"logFC_percentile" %in% names(da_results)) {
            da_results <- .milo_prepare_metric_column(da_results, metric = "logFC_percentile")
        }
        colour_by <- fallback_metric
        metric_vals <- da_results[[fallback_metric]]
    }

    list(
        results = da_results,
        colour_by = colour_by,
        alpha = if (colour_by %in% c("SpatialFDR", "PValue")) alpha else 1,
        title = if (colour_by == preferred_metric) "DA Results by Cell Type"
        else sprintf("DA Results by Cell Type (%s fallback)", colour_by)
    )
}

.milo_prepare_seurat <- function(seurat_obj, seurat_qs_path, max_cells, seed, verbose) {
    if (!is.null(seurat_obj)) {
        return(.milo_subset_seurat(seurat_obj, max_cells, seed, verbose))
    }
    if (!is.null(seurat_qs_path)) {
        if (verbose) cli::cli_inform(c("Loading Seurat object from {.path {seurat_qs_path}}"))
        sobj <- qs::qread(seurat_qs_path)
        return(.milo_subset_seurat(sobj, max_cells, seed, verbose))
    }
    if (verbose) cli::cli_warn("No Seurat object or path supplied; downstream steps may rely on cached reductions.")
    NULL
}

.milo_subset_seurat <- function(seurat_obj, max_cells, seed, verbose) {
    if (!inherits(seurat_obj, "Seurat")) {
        stop("`seurat_obj` must be a Seurat object.")
    }
    if (!is.null(max_cells) && max_cells > 0 && max_cells < ncol(seurat_obj)) {
        set.seed(seed)
        if (verbose) {
            cli::cli_inform(c("Subsampling cells" = "{ncol(seurat_obj)} → {max_cells} for development run."))
        }
        keep_cells <- sample(colnames(seurat_obj), size = max_cells)
        seurat_obj <- subset(seurat_obj, cells = keep_cells)
    }
    seurat_obj
}

.milo_normalize_alpha <- function(alpha, name) {
    alpha_num <- suppressWarnings(as.numeric(alpha))
    if (is.na(alpha_num)) {
        stop(sprintf("`%s` must be numeric; received value that cannot be coerced.", name))
    }
    alpha_num
}

.milo_metric_palette <- function(labels) {
    base_palette <- c(
        "< 0.1" = "red",
        "0.1 - 0.2" = "orange",
        "0.2 - 0.3" = "yellow",
        "0.3 - 0.5" = "grey60",
        ">= 0.5" = "grey90"
    )
    mapped <- base_palette[match(labels, names(base_palette))]
    mapped[is.na(mapped)] <- "grey80"
    names(mapped) <- labels
    mapped
}

.milo_ensure_nhood_graph <- function(milo, verbose) {
    graph_obj <- miloR::nhoodGraph(milo)
    needs_build <- is.null(graph_obj) || (is.list(graph_obj) && length(graph_obj) == 0)
    if (needs_build) {
        if (verbose) cli::cli_inform(c("step" = "Running buildNhoodGraph() before plotting."))
        milo <- miloR::buildNhoodGraph(milo)
    }
    milo
}

.milo_attach_umap <- function(milo, seurat_obj, layout_reduction, verbose) {
    if (!"UMAP" %in% SingleCellExperiment::reducedDimNames(milo)) {
        if (is.null(seurat_obj)) {
            if (verbose) {
                cli::cli_warn("Reduction 'UMAP' missing on Milo object and no Seurat object supplied; skipping UMAP attachment.")
            }
            return(milo)
        }
        if (layout_reduction %in% Seurat::Reductions(seurat_obj)) {
            if (verbose) cli::cli_inform(c("step" = paste0("Attaching UMAP layout '", layout_reduction, "' to Milo object.")))
            umap_emb <- Seurat::Embeddings(seurat_obj, reduction = layout_reduction)
            umap_emb <- as.matrix(umap_emb)
            SingleCellExperiment::reducedDim(milo, "UMAP") <- umap_emb
        } else if (verbose) {
            cli::cli_warn(paste0("Reduction '", layout_reduction, "' not found in Seurat object; UMAP plots may fail."))
        }
    }
    milo
}

.milo_normalize_force_flags <- function(force_run) {
    steps <- c("nhoods", "distances", "testing")
    if (length(force_run) == 1L) {
        flags <- rep(as.logical(force_run), length(steps))
        names(flags) <- steps
        return(flags)
    }
    if (is.null(names(force_run))) {
        stop("Named logical vector expected for `force_run` when length > 1. Use names: ", paste(steps, collapse = ", "))
    }
    flags <- rep(FALSE, length(steps))
    names(flags) <- steps
    overlap <- intersect(names(force_run), steps)
    flags[overlap] <- as.logical(force_run[overlap])
    flags
}

.milo_resolve_paths <- function(output_dir, prefix, suffix, save) {
    files <- list()
    if (!save) {
        return(list(output_dir = output_dir, files = files, suffix = suffix))
    }

    base_names <- c(
        nhoods = "01_nhoods_built.qs",
        distances = "02_distances_calculated.qs",
        da_milo = "03_tested.qs",
        da_results = "03_da_results.qs",
        plots = "04_plots.qs"
    )

    suffix_info <- .milo_pick_suffix(output_dir, prefix, suffix, base_names)
    suffix_str <- suffix_info$suffix_str

    files <- lapply(base_names, function(base) {
        file.path(output_dir, paste0(prefix, "_", .milo_strip_ext(base), suffix_str, ".", tools::file_ext(base)))
    })

    names(files) <- names(base_names)
    list(output_dir = output_dir, files = files, suffix = suffix_info$suffix)
}

.milo_apply_cache_overrides <- function(paths, cache_files, verbose) {
    if (is.null(cache_files) || !length(cache_files)) {
        return(paths)
    }
    if (is.null(paths$files)) {
        paths$files <- list()
    }
    valid <- intersect(names(paths$files), names(cache_files))
    for (nm in valid) {
        override <- cache_files[[nm]]
        if (is.null(override) || identical(override, "")) next
        paths$files[[nm]] <- override
        if (isTRUE(verbose)) {
            cli::cli_inform(c("cache" = paste0("Using manual cache for ", nm, ": ", override)))
        }
    }
    paths
}

.milo_strip_ext <- function(filename) sub(paste0(".", tools::file_ext(filename), "$"), "", filename)

.milo_pick_suffix <- function(output_dir, prefix, suffix, base_names) {
    make_paths <- function(sfx) {
        suffix_str <- if (is.null(sfx) || identical(sfx, "")) "" else paste0("_", sfx)
        vapply(
            base_names,
            function(base) {
                file.path(output_dir, paste0(prefix, "_", .milo_strip_ext(base), suffix_str, ".", tools::file_ext(base)))
            },
            character(1)
        )
    }

    if (!is.null(suffix)) {
        suffix_str <- if (identical(suffix, "") || is.null(suffix)) "" else paste0("_", suffix)
        return(list(suffix = as.character(suffix), suffix_str = suffix_str, paths = make_paths(suffix)))
    }

    base_paths <- make_paths(NULL)
    if (any(file.exists(base_paths))) {
        return(list(suffix = NULL, suffix_str = "", paths = base_paths))
    }
    return(list(suffix = NULL, suffix_str = "", paths = base_paths))
}

#' Test cluster-wise logFC bias in Milo DA results
#'
#' @description
#' Tests whether logFC values are systematically biased (enriched/depleted) per cluster,
#' accounting for non-independence between neighborhoods using block permutation,
#' correlation-adjusted t-test, or mixed-effects models.
#'
#' @param da_results Data frame from `miloR::testNhoods()` with columns `Nhood`, `logFC`, `PValue`, etc.
#' @param milo Milo object containing nhood membership matrix and graph.
#' @param cluster_col Character string; column name in `da_results` containing cluster identities (default: `"major_cluster"`).
#' @param block_method Character; method for creating blocks: `"sample"` (extract from cell names), `"community"` (graph-based), or `"none"` (no blocking).
#' @param test_methods Character vector; which tests to run: `"permutation"`, `"neff"`, `"lmm"`, `"ashr"`.
#' @param n_perm Integer; number of permutations for block permutation test (default: 2000).
#' @param max_nhoods Optional integer; if provided, randomly sample this many neighborhoods for faster testing.
#' @param block_var Optional character; column name in `colData(milo)` to use for blocking (e.g., `patient_var`, `batch_var`). If `NULL`, attempts to extract from cell names.
#' @param seed Integer; random seed for reproducibility (default: 1).
#' @param verbose Logical; emit progress messages.
#'
#' @return Data frame with columns:
#'   - `cluster`: cluster identity
#'   - `mean_logFC`: mean logFC per cluster
#'   - `p_perm`: p-value from block permutation (if `"permutation"` in `test_methods`)
#'   - `p_neff`: p-value from correlation-adjusted t-test (if `"neff"` in `test_methods`)
#'   - `p_lmm`: p-value from mixed-effects model (if `"lmm"` in `test_methods`)
#'   - `eb_mean`: Empirical Bayes posterior mean (if `"ashr"` in `test_methods`)
#'
#' @importFrom stringr str_extract
#' @importFrom Matrix sparse.model.matrix
#' @importFrom stats qnorm pt
#' @export
test_cluster_logfc_bias <- function(
    da_results,
    milo,
    cluster_col = "major_cluster",
    block_method = c("sample", "community", "none"),
    test_methods = c("permutation", "neff"),
    n_perm = 2000L,
    max_nhoods = NULL,
    block_var = NULL,
    seed = 1L,
    verbose = TRUE
) {
    if (!cluster_col %in% names(da_results)) {
        stop(sprintf("Column '%s' not found in `da_results`.", cluster_col))
    }

    block_method <- match.arg(block_method)
    test_methods <- match.arg(test_methods, several.ok = TRUE)

    # Downsample if requested
    if (!is.null(max_nhoods) && max_nhoods < nrow(da_results)) {
        if (verbose) cli::cli_inform("Downsampling to {max_nhoods} neighborhoods for faster testing...")
        set.seed(seed)
        keep_idx <- sample.int(nrow(da_results), size = max_nhoods)
        da_results <- da_results[keep_idx, , drop = FALSE]
    }

    # Restore or create block variable
    if (block_method == "sample") {
        block_id <- NULL
        
        # Method 1: Use block_var from colData if provided (preferred)
        if (!is.null(block_var) && block_var %in% names(SingleCellExperiment::colData(milo))) {
            cell_meta <- as.data.frame(SingleCellExperiment::colData(milo))
            cell_block_id <- as.character(cell_meta[[block_var]])
            names(cell_block_id) <- rownames(cell_meta)
            
            W <- miloR::nhoods(milo)
            cells <- rownames(W)
            if (is.null(cells)) {
                cells <- rownames(SingleCellExperiment::colData(milo))
            }
            block_id <- cell_block_id[cells]
            
            if (any(is.na(block_id)) || length(unique(block_id)) < 2) {
                if (verbose) cli::cli_warn("block_var has insufficient unique values; falling back to cell name extraction.")
                block_id <- NULL
            }
        }
        
        # Method 2: Extract from cell names (fallback)
        if (is.null(block_id) || all(is.na(block_id))) {
            W <- miloR::nhoods(milo)
            cells <- rownames(W)
            if (is.null(cells)) {
                cells <- rownames(SingleCellExperiment::colData(milo))
            }
            
            # Try pattern "S\\d+$" (e.g., "S1", "S2")
            block_id <- stringr::str_extract(cells, "S\\d+$")
            if (all(is.na(block_id))) {
                # Try first 8 characters (common in some naming schemes)
                block_id <- substr(cells, 1, 8)
            }
            if (all(is.na(block_id))) {
                # Last resort: use first part before any separator
                block_id <- sub("[-_].*$", "", cells)
            }
        }
        
        # Validate extracted block IDs
        if (all(is.na(block_id)) || length(unique(block_id)) < 2) {
            if (verbose) cli::cli_warn("Could not extract block_id from cell names or metadata; falling back to community blocks.")
            block_method <- "community"
        } else {
            S <- Matrix::sparse.model.matrix(~0 + factor(block_id))
            W <- miloR::nhoods(milo)
            # S is cells × blocks, W is cells × nhoods
            # t(S) is blocks × cells, so t(S) %*% W is blocks × nhoods
            C_block <- Matrix::t(S) %*% W
            top_block <- apply(C_block, 2, function(v) {
                idx <- which.max(v)
                if (length(idx) && idx > 0) colnames(C_block)[idx] else NA_character_
            })
            da_results$block_id <- top_block[da_results$Nhood]
            if (verbose && any(is.na(da_results$block_id))) {
                cli::cli_warn("Some neighborhoods lack block_id; falling back to community blocks.")
                block_method <- "community"
            }
        }
    }

    if (block_method == "community" || (block_method == "sample" && any(is.na(da_results$block_id)))) {
        G <- miloR::nhoodGraph(milo)
        if (is.null(G) || length(G) == 0) {
            if (verbose) cli::cli_warn("nhoodGraph not available; building it.")
            milo <- miloR::buildNhoodGraph(milo)
            G <- miloR::nhoodGraph(milo)
        }
        if (requireNamespace("igraph", quietly = TRUE)) {
            comm <- igraph::cluster_louvain(G)$membership
            da_results$block_id <- as.character(comm[da_results$Nhood])
        } else {
            if (verbose) cli::cli_warn("Package 'igraph' not available; skipping community blocks.")
            da_results$block_id <- "all"
        }
    }

    if (block_method == "none") {
        da_results$block_id <- "all"
    }

    # Remove rows with missing cluster or block
    da_results <- da_results[!is.na(da_results[[cluster_col]]) & !is.na(da_results$block_id), ]

    # Group by cluster and compute statistics
    cluster_stats <- da_results %>%
        dplyr::group_by(.data[[cluster_col]]) %>%
        dplyr::summarise(
            mean_logFC = mean(.data$logFC, na.rm = TRUE),
            n_nhoods = dplyr::n(),
            .groups = "drop"
        )

    # Block permutation test
    if ("permutation" %in% test_methods) {
        if (verbose) cli::cli_inform("Running block permutation tests...")
        perm_test_block <- function(df, block_var = "block_id", n = n_perm) {
            block_summary <- df %>%
                dplyr::group_by(.data[[block_var]]) %>%
                dplyr::summarise(
                    mean_logFC = mean(.data$logFC, na.rm = TRUE),
                    weight = dplyr::n(),
                    .groups = "drop"
                )
            if (nrow(block_summary) < 2) {
                return(NA_real_)
            }
            obs <- sum(block_summary$mean_logFC * block_summary$weight) / sum(block_summary$weight)
            perm_means <- replicate(n, {
                perm_idx <- sample.int(nrow(block_summary))
                sum(block_summary$mean_logFC[perm_idx] * block_summary$weight) / sum(block_summary$weight)
            })
            (sum(abs(perm_means) >= abs(obs), na.rm = TRUE) + 1) / (n + 1)
        }

        cluster_perm <- da_results %>%
            dplyr::group_by(.data[[cluster_col]]) %>%
            dplyr::summarise(
                p_perm = perm_test_block(dplyr::cur_data(), block_var = "block_id", n = n_perm),
                .groups = "drop"
            )
        cluster_stats <- dplyr::left_join(cluster_stats, cluster_perm, by = cluster_col)
    }

    # Correlation-adjusted t-test (neff correction)
    if ("neff" %in% test_methods) {
        if (verbose) cli::cli_inform("Computing correlation-adjusted t-tests...")
        G <- tryCatch(miloR::nhoodGraph(milo), warning = function(w) NULL, error = function(e) NULL)
        if (is.null(G) || length(G) == 0) {
            if (verbose) cli::cli_inform("Building nhoodGraph...")
            milo <- miloR::buildNhoodGraph(milo)
            G <- miloR::nhoodGraph(milo)
        }
        if (requireNamespace("igraph", quietly = TRUE)) {
            adj_mat <- igraph::as_adjacency_matrix(G, sparse = TRUE)
            nhood_idx <- match(da_results$Nhood, as.numeric(igraph::V(G)$name))
            valid_idx <- !is.na(nhood_idx)
            if (sum(valid_idx) > 0) {
                adj_sub <- adj_mat[nhood_idx[valid_idx], nhood_idx[valid_idx], drop = FALSE]
                eigen_vals <- eigen(as.matrix(adj_sub), only.values = TRUE)$values
                neff <- sum(eigen_vals > 0.1)
                neff <- max(neff, 1)
            } else {
                neff <- nrow(da_results)
            }
        } else {
            neff <- nrow(da_results)
        }

        cluster_neff <- da_results %>%
            dplyr::group_by(.data[[cluster_col]]) %>%
            dplyr::summarise(
                mean_lfc = mean(.data$logFC, na.rm = TRUE),
                sd_lfc = sd(.data$logFC, na.rm = TRUE),
                n = dplyr::n(),
                .groups = "drop"
            ) %>%
            dplyr::mutate(
                neff_cluster = pmax(2, round(neff * (.data$n / nrow(da_results)))),
                sd_lfc = pmax(.data$sd_lfc, 1e-10),
                t_stat = .data$mean_lfc / (.data$sd_lfc / sqrt(.data$neff_cluster)),
                p_neff = ifelse(is.finite(.data$t_stat) & .data$neff_cluster > 1,
                    2 * (1 - pt(abs(.data$t_stat), df = .data$neff_cluster - 1)),
                    NA_real_)
            ) %>%
            dplyr::select(.data[[cluster_col]], .data$p_neff)
        cluster_stats <- dplyr::left_join(cluster_stats, cluster_neff, by = cluster_col)
    }

    # Mixed-effects model
    if ("lmm" %in% test_methods) {
        if (requireNamespace("lme4", quietly = TRUE)) {
            if (verbose) cli::cli_inform("Fitting mixed-effects models...")
            cluster_lmm <- da_results %>%
                dplyr::group_by(.data[[cluster_col]]) %>%
                dplyr::summarise(
                    p_lmm = tryCatch({
                        mod <- lme4::lmer(logFC ~ 1 + (1 | block_id), data = dplyr::cur_data())
                        coefs <- summary(mod)$coefficients
                        if (nrow(coefs) > 0) {
                            coefs[1, "Pr(>|t|)"]
                        } else {
                            NA_real_
                        }
                    }, error = function(e) NA_real_),
                    .groups = "drop"
                )
            cluster_stats <- dplyr::left_join(cluster_stats, cluster_lmm, by = cluster_col)
        } else {
            if (verbose) cli::cli_warn("Package 'lme4' not available; skipping LMM tests.")
        }
    }

    # Empirical Bayes (ashr)
    if ("ashr" %in% test_methods) {
        if (requireNamespace("ashr", quietly = TRUE)) {
            if (verbose) cli::cli_inform("Computing Empirical Bayes estimates...")
            da_results <- da_results %>%
                dplyr::mutate(
                    z = qnorm(1 - .data$PValue / 2, lower.tail = FALSE) * sign(.data$logFC),
                    se = abs(.data$logFC) / pmax(abs(.data$z), 1e-6)
                )
            cluster_ashr <- da_results %>%
                dplyr::group_by(.data[[cluster_col]]) %>%
                dplyr::summarise(
                    eb_mean = tryCatch({
                        fit <- ashr::ash(.data$logFC, sebetahat = .data$se, method = "fdr")
                        mean(fit$result$PosteriorMean, na.rm = TRUE)
                    }, error = function(e) NA_real_),
                    .groups = "drop"
                )
            cluster_stats <- dplyr::left_join(cluster_stats, cluster_ashr, by = cluster_col)
        } else {
            if (verbose) cli::cli_warn("Package 'ashr' not available; skipping Empirical Bayes estimates.")
        }
    }

    cluster_stats
}


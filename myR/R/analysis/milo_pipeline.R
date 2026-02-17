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
#' @param save Logical; if `TRUE` (default) intermediate objects are written to
#'   disk using `qs::qsave()`.
#' @param output_dir,prefix,suffix Configure save location. When `suffix` is
#'   `NULL`, a numeric suffix (`"", "_1", ...`) is auto-generated to avoid
#'   overwriting existing files.
#' @param force_run Logical scalar or named logical vector controlling whether
#'   each step (`nhoods`, `distances`, `testing`) is recomputed even if cached
#'   files are present.
#' @param plotting Logical; if `FALSE`, plotting helpers are skipped.
#' @param max_cells Optional integer. When set and lower than the total number
#'   of cells, a random subset is used to accelerate development runs.
#' @param seed Integer passed to `set.seed()` before sampling or Milo routines.
#' @param verbose Logical; emits progress messages when `TRUE`.
#'
#' @return A list with `milo` (the Milo object), `da_results` (data frame), and
#'   `plots` (list of ggplot objects; `NULL` when `plotting = FALSE`).
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
    save = TRUE,
    output_dir = file.path(tempdir(), "milo"),
    prefix = "milo",
    suffix = NULL,
    force_run = FALSE,
    plotting = TRUE,
    max_cells = NULL,
    seed = 1L,
    verbose = TRUE
) {
    .milo_require_packages()

    if (is.null(seurat_obj)) {
        if (is.null(seurat_qs_path)) {
            stop("Provide either `seurat_obj` or `seurat_qs_path`.")
        }
        if (verbose) cli::cli_inform(c("Loading Seurat object from {.path {seurat_qs_path}}"))
        seurat_obj <- qs::qread(seurat_qs_path)
    }

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

    set.seed(seed)

    force_flags <- .milo_normalize_force_flags(force_run)
    paths <- .milo_resolve_paths(
        output_dir = output_dir,
        prefix = prefix,
        suffix = suffix,
        save = save
    )

    if (save) {
        dir.create(paths$output_dir, recursive = TRUE, showWarnings = FALSE)
    }

    if (verbose) cli::cli_inform(c("Preparing metadata for Milo conversion."))
    seurat_obj <- .milo_prepare_metadata(seurat_obj, patient_var, cluster_var, target_var, batch_var)

    # ---- Step 1: Build Milo object with neighbourhoods ----
    if (!force_flags["nhoods"] && save && file.exists(paths$files$nhoods)) {
        if (verbose) {
            cli::cli_inform(c("step" = "Loading cached Milo object with neighbourhoods from {.path {paths$files$nhoods}}"))
        }
        milo <- qs::qread(paths$files$nhoods)
    } else {
        if (verbose) cli::cli_inform(c("step" = "Building Milo neighbourhoods."))
        milo <- .milo_build_nhoods(
            seurat_obj = seurat_obj,
            graph_reduction = graph_reduction,
            k = k,
            d = d,
            prop = prop
        )
        if (save) {
            qs::qsave(milo, paths$files$nhoods)
            if (verbose) cli::cli_inform(c("cache" = "Saved Milo object to {.path {paths$files$nhoods}}"))
        }
    }

    # ---- Step 2: Distances ----
    if (!force_flags["distances"] && save && file.exists(paths$files$distances)) {
        if (verbose) {
            cli::cli_inform(c("step" = "Loading cached Milo object with neighbourhood distances from {.path {paths$files$distances}}"))
        }
        milo <- qs::qread(paths$files$distances)
    } else {
        if (verbose) cli::cli_inform(c("step" = "Calculating Milo neighbourhood distances (this may take a while)."))
        milo <- miloR::calcNhoodDistance(milo, d = d, reduced.dim = "GRAPH")
        if (save) {
            qs::qsave(milo, paths$files$distances)
            if (verbose) cli::cli_inform(c("cache" = "Saved Milo object to {.path {paths$files$distances}}"))
        }
    }

    # ---- Step 3: DA testing ----
    da_results <- NULL
    has_cached_da <- save && file.exists(paths$files$da_milo) && file.exists(paths$files$da_results)
    if (!force_flags["testing"] && has_cached_da) {
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
            cluster_var = cluster_var
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

    plots <- NULL
    if (plotting) {
        if (verbose) cli::cli_inform(c("step" = "Preparing Milo plots."))
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
    pkgs <- c("Seurat", "SingleCellExperiment", "miloR", "Matrix", "qs", "cli", "ggplot2", "patchwork")
    missing <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
    if (length(missing) > 0) {
        stop("Missing required packages: ", paste(missing, collapse = ", "))
    }
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

.milo_build_nhoods <- function(seurat_obj, graph_reduction, k, d, prop) {
    if (!graph_reduction %in% Seurat::Reductions(seurat_obj)) {
        stop(sprintf("Reduction '%s' not found in Seurat object.", graph_reduction))
    }

    G_emb <- Seurat::Embeddings(seurat_obj, reduction = graph_reduction)
    d_eff <- min(ncol(G_emb), d)
    G_use <- G_emb[, seq_len(d_eff), drop = FALSE]

    sce <- Seurat::as.SingleCellExperiment(seurat_obj, assay = Seurat::DefaultAssay(seurat_obj))
    SingleCellExperiment::reducedDim(sce, "GRAPH") <- as.matrix(G_use)

    milo <- miloR::Milo(sce)
    milo <- miloR::buildGraph(milo, k = k, d = d, reduced.dim = "GRAPH")
    miloR::makeNhoods(milo, prop = prop, k = k, refined = TRUE, reduced_dims = "GRAPH")
}

.milo_run_da <- function(milo, patient_var, target_var, batch_var, cluster_var, reduced_dim = "GRAPH") {
    milo <- miloR::countCells(
        milo,
        samples = patient_var,
        meta.data = as.data.frame(SingleCellExperiment::colData(milo))
    )

    meta_df <- as.data.frame(SingleCellExperiment::colData(milo))
    keep_cols <- c(patient_var, target_var, batch_var)
    patient_info <- unique(meta_df[, keep_cols, drop = FALSE])
    rownames(patient_info) <- patient_info[[patient_var]]

    sample_design <- patient_info[colnames(miloR::nhoodCounts(milo)), , drop = FALSE]
    sample_design[[target_var]] <- factor(sample_design[[target_var]])
    sample_design[[batch_var]] <- factor(sample_design[[batch_var]])

    formula_str <- sprintf("~ %s + %s", batch_var, target_var)
    da_results <- miloR::testNhoods(
        milo,
        design = as.formula(formula_str),
        design.df = sample_design,
        reduced.dim = reduced_dim
    )

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

    graph_plot <- miloR::plotNhoodGraphDA(milo, da_results, alpha = alpha) +
        ggplot2::ggtitle("DA Results (logFC)")

    beeswarm_plot <- miloR::plotDAbeeswarm(
        da_results,
        group.by = "major_cluster",
        alpha = beeswarm_alpha
    ) + ggplot2::ggtitle("DA Results by Cell Type")

    umap_plot <- miloR::plotUMAP(milo, colour_by = "milo_metric_bin") +
        ggplot2::scale_color_manual(values = .milo_metric_palette(fdr_labels), name = "Cell metric") +
        ggplot2::ggtitle(sprintf("Cells by %s", cell_metric))

    combined <- patchwork::wrap_plots(umap_plot, graph_plot)

    plots <- list(
        graph = graph_plot,
        beeswarm = beeswarm_plot,
        umap = umap_plot,
        combined = combined
    )

    if (save) {
        saveRDS(plots, file = save_path)
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
    N_summary <- summary(N_matrix)
    N_summary$metric <- metric_values[N_summary$j]

    min_metric_per_cell <- tapply(N_summary$metric, N_summary$i, min, na.rm = TRUE)
    cell_metric <- rep(max(breaks, na.rm = TRUE), nrow(N_matrix))
    cell_metric[as.integer(names(min_metric_per_cell))] <- min_metric_per_cell

    metric_bins <- cut(
        cell_metric,
        breaks = breaks,
        labels = labels,
        right = TRUE,
        include.lowest = TRUE
    )

    SingleCellExperiment::colData(milo)$milo_metric_value <- cell_metric
    SingleCellExperiment::colData(milo)$milo_metric_bin <- metric_bins

    if (verbose) {
        bin_table <- table(metric_bins, useNA = "ifany")
        summary_text <- paste(sprintf("%s=%d", names(bin_table), bin_table), collapse = ", ")
        cli::cli_inform(paste("Cell metric bins:", summary_text))
    }
    milo
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
    if (is.null(miloR::nhoodGraph(milo))) {
        if (verbose) cli::cli_inform(c("step" = "Running buildNhoodGraph() before plotting."))
        milo <- miloR::buildNhoodGraph(milo)
    }
    milo
}

.milo_attach_umap <- function(milo, seurat_obj, layout_reduction, verbose) {
    if (!"UMAP" %in% SingleCellExperiment::reducedDimNames(milo)) {
        if (layout_reduction %in% Seurat::Reductions(seurat_obj)) {
            if (verbose) cli::cli_inform(c("step" = paste0("Attaching UMAP layout '", layout_reduction, "' to Milo object.")))
            umap_emb <- Seurat::Embeddings(seurat_obj, reduction = layout_reduction)
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
        return(list(output_dir = output_dir, files = files))
    }

    base_names <- c(
        nhoods = "01_nhoods_built.qs",
        distances = "02_distances_calculated.qs",
        da_milo = "03_tested.qs",
        da_results = "03_da_results.qs",
        plots = "04_plots.rds"
    )

    resolved_suffix <- .milo_pick_suffix(output_dir, prefix, suffix, base_names)
    suffix_str <- if (is.null(resolved_suffix)) "" else paste0("_", resolved_suffix)

    files <- lapply(base_names, function(base) {
        file.path(output_dir, paste0(prefix, "_", .milo_strip_ext(base), suffix_str, ".", tools::file_ext(base)))
    })

    names(files) <- names(base_names)
    list(output_dir = output_dir, files = files)
}

.milo_strip_ext <- function(filename) sub(paste0(".", tools::file_ext(filename), "$"), "", filename)

.milo_pick_suffix <- function(output_dir, prefix, suffix, base_names) {
    if (!is.null(suffix)) return(as.character(suffix))
    idx <- 0
    repeat {
        current_suffix <- if (idx == 0) NULL else idx
        suffix_str <- if (is.null(current_suffix)) "" else paste0("_", current_suffix)
        candidate_paths <- vapply(
            base_names,
            function(base) {
                file.path(output_dir, paste0(prefix, "_", .milo_strip_ext(base), suffix_str, ".", tools::file_ext(base)))
            },
            character(1)
        )
        if (!any(file.exists(candidate_paths))) {
            return(current_suffix)
        }
        idx <- idx + 1
    }
}


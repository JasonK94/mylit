#' Build Cluster-Level Marker Signatures for Module Scoring
#'
#' Convert a `FindAllMarkers`-style data.frame into a list of gene vectors that
#' can be consumed by `Seurat::AddModuleScore()`. Either supply a Seurat object
#' (the default, used to infer a background gene universe) or pass
#' `background_genes` explicitly.
#'
#' @param markers_df Data frame with at least `gene`, `cluster` and one of
#'   `avg_log2FC`/`avg_logFC` plus one of `p_val_adj`/`p_val`.
#' @param seurat_obj Optional Seurat object used to infer available features.
#' @param background_genes Optional character vector restricting the signature
#'   universe; takes precedence over `seurat_obj` if provided.
#' @param assay Assay to pull features from when `seurat_obj` is supplied. Falls
#'   back to `DefaultAssay(seurat_obj)`.
#' @param top_n Number of genes to keep per cluster (sorted by logFC).
#' @param p_adj_cut Adjusted p-value threshold.
#' @param min_genes Minimum genes required to retain a cluster signature.
#' @param verbose Logical; emit progress messages.
#'
#' @return A list with `features_list`, `map`, and `settings`.
#' @export
make_cluster_signatures <- function(markers_df,
                                    seurat_obj = NULL,
                                    background_genes = NULL,
                                    assay = NULL,
                                    top_n = 20,
                                    p_adj_cut = 0.05,
                                    min_genes = 5,
                                    verbose = FALSE) {
  .msg <- function(...) if (isTRUE(verbose)) cat(...)

  stopifnot(is.data.frame(markers_df))
  if (!all(c("gene", "cluster") %in% names(markers_df))) {
    stop("`markers_df` must contain columns: gene, cluster")
  }

  logfc_col <- if ("avg_log2FC" %in% names(markers_df)) {
    "avg_log2FC"
  } else if ("avg_logFC" %in% names(markers_df)) {
    "avg_logFC"
  } else {
    stop("`markers_df` needs avg_log2FC or avg_logFC")
  }
  p_col <- if ("p_val_adj" %in% names(markers_df)) {
    "p_val_adj"
  } else if ("p_val" %in% names(markers_df)) {
    "p_val"
  } else {
    stop("`markers_df` needs p_val_adj or p_val")
  }

  if (is.null(background_genes)) {
    if (is.null(seurat_obj)) {
      stop("Provide either `background_genes` or `seurat_obj`.")
    }
    if (is.null(assay)) {
      assay <- Seurat::DefaultAssay(seurat_obj)
    }
    assay_names <- tryCatch({
      ax <- Seurat::Assays(seurat_obj)
      if (is.character(ax)) ax else names(seurat_obj@assays)
    }, error = function(e) names(seurat_obj@assays))
    if (!(assay %in% assay_names)) {
      stop(sprintf("Assay '%s' not found. Available: %s",
                   assay, paste(assay_names, collapse = ", ")))
    }
    background_genes <- tryCatch(rownames(Seurat::GetAssay(seurat_obj, assay)),
                                 error = function(e) NULL)
    if (is.null(background_genes) || !length(background_genes)) {
      background_genes <- tryCatch(
        rownames(Seurat::GetAssayData(seurat_obj, assay = assay, layer = "counts")),
        error = function(e) NULL
      )
    }
    if (is.null(background_genes) || !length(background_genes)) {
      background_genes <- tryCatch(
        rownames(Seurat::GetAssayData(seurat_obj, assay = assay, layer = "data")),
        error = function(e) NULL
      )
    }
    if (is.null(background_genes) || !length(background_genes)) {
      stop("Could not infer feature names from `seurat_obj`; ",
           "provide `background_genes` explicitly.")
    }
  }
  background_genes <- as.character(background_genes)

  if (verbose) {
    .msg(sprintf("[make_cluster_signatures] Using %d background genes\n",
                 length(background_genes)))
  }

  markers_df$gene <- as.character(markers_df$gene)
  sig_df <- markers_df |>
    dplyr::filter(.data[[p_col]] < p_adj_cut) |>
    dplyr::group_by(.data$cluster) |>
    dplyr::slice_max(n = top_n, order_by = .data[[logfc_col]]) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      cluster = as.character(.data$cluster),
      gene = as.character(.data$gene)
    ) |>
    dplyr::filter(!is.na(.data$gene), .data$gene %in% background_genes) |>
    dplyr::distinct(.data$cluster, .data$gene, .keep_all = TRUE)

  if (!nrow(sig_df)) {
    stop("No marker genes passed filters / present in `background_genes`.")
  }

  features_by_cluster <- split(sig_df$gene, sig_df$cluster)
  keep <- lengths(features_by_cluster) >= min_genes
  features_by_cluster <- features_by_cluster[keep]
  if (!length(features_by_cluster)) {
    stop("No modules left after `min_genes` filter.")
  }

  features_list <- unname(features_by_cluster)
  map_tbl <- tibble::tibble(
    module_idx = seq_along(features_list),
    cluster = names(features_by_cluster),
    n_genes = lengths(features_by_cluster)
  )

  list(
    features_list = features_list,
    map = map_tbl,
    settings = list(
      assay = assay,
      top_n = top_n,
      p_adj_cut = p_adj_cut,
      min_genes = min_genes
    )
  )
}

#' Closed Log-Ratio Transform
#'
#' Apply a CLR transform (with pseudo-count) to compositional matrices.
#'
#' @param X Numeric matrix (rows = samples, columns = components).
#' @param eps Pseudo-count added before log transform.
#'
#' @return Matrix on CLR scale.
#' @export
clr_transform <- function(X, eps = 1e-6) {
  stopifnot(is.matrix(X))
  Xp <- X + eps
  logX <- log(Xp)
  sweep(logX, 1, rowMeans(logX), FUN = "-")
}

#' Compute Patient-Level Cluster Frequencies
#'
#' Collapse cell-level metadata into patient-by-cluster frequency matrices
#' (optionally CLR transformed).
#'
#' @param seurat_obj Seurat object containing at least `sample_col` and
#'   `cluster_col` in `@meta.data`.
#' @param sample_col Column with sample/patient identifiers.
#' @param cluster_col Column with cluster identities.
#' @param transform Either `"clr"` (default) or `"none"`.
#' @param eps Pseudo-count used when `transform = "clr"`.
#'
#' @return Numeric matrix with patients as rows and clusters as columns.
#' @export
compute_patient_cluster_frequency <- function(seurat_obj,
                                              sample_col = "hos_no",
                                              cluster_col = "anno3.scvi",
                                              transform = c("clr", "none"),
                                              eps = 1e-6) {
  transform <- match.arg(transform)
  meta <- seurat_obj@meta.data
  required <- c(sample_col, cluster_col)
  missing_cols <- setdiff(required, colnames(meta))
  if (length(missing_cols)) {
    stop("Missing columns in metadata: ", paste(missing_cols, collapse = ", "))
  }

  sample_sym <- rlang::sym(sample_col)
  cluster_sym <- rlang::sym(cluster_col)

  freq_tbl <- meta |>
    dplyr::count(!!sample_sym, !!cluster_sym, name = "n") |>
    dplyr::group_by(!!sample_sym) |>
    dplyr::mutate(freq = .data$n / sum(.data$n)) |>
    dplyr::ungroup() |>
    tidyr::pivot_wider(
      names_from = !!cluster_sym,
      values_from = .data$freq,
      values_fill = 0
    ) |>
    dplyr::arrange(!!sample_sym)

  freq_mat <- freq_tbl |>
    dplyr::select(-dplyr::all_of(sample_col)) |>
    as.matrix()
  rownames(freq_mat) <- freq_tbl[[sample_col]]

  if (identical(transform, "clr")) {
    freq_mat <- clr_transform(freq_mat, eps = eps)
  }

  freq_mat
}

#' Compute Patient-Level Signature Matrix from Module Scores
#'
#' Run `Seurat::AddModuleScore()` and aggregate the resulting signatures to
#' patient×cluster means.
#'
#' @param seurat_obj Seurat object with metadata columns `sample_col` and
#'   `cluster_col`.
#' @param markers_df Data frame supplied to `make_cluster_signatures()`.
#' @inheritParams make_cluster_signatures
#' @param sample_col Column with sample/patient identifiers.
#' @param cluster_col Column with cluster identities.
#' @param assay Assay passed to `AddModuleScore` (defaults to `DefaultAssay`).
#' @param score_prefix Prefix passed to `AddModuleScore::name`.
#'
#' @return List with `matrix`, `seurat_obj`, `map`, and `score_columns`.
#' @export
compute_patient_signature_matrix <- function(seurat_obj,
                                             markers_df,
                                             background_genes = NULL,
                                             assay = NULL,
                                             sample_col = "hos_no",
                                             cluster_col = "anno3.scvi",
                                             top_n = 20,
                                             p_adj_cut = 0.05,
                                             min_genes = 5,
                                             score_prefix = "sig_",
                                             verbose = FALSE) {
  stopifnot(inherits(seurat_obj, "Seurat"))
  meta <- seurat_obj@meta.data
  required <- c(sample_col, cluster_col)
  missing_cols <- setdiff(required, colnames(meta))
  if (length(missing_cols)) {
    stop("Missing columns in metadata: ", paste(missing_cols, collapse = ", "))
  }

  sig_res <- make_cluster_signatures(
    markers_df = markers_df,
    seurat_obj = if (is.null(background_genes)) seurat_obj else NULL,
    background_genes = background_genes,
    assay = assay,
    top_n = top_n,
    p_adj_cut = p_adj_cut,
    min_genes = min_genes,
    verbose = verbose
  )
  features_list <- sig_res$features_list

  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(seurat_obj)
  }

  if (isTRUE(verbose)) {
    cat("[compute_patient_signature_matrix] Running AddModuleScore...\n")
  }
  seurat_obj <- Seurat::AddModuleScore(
    object = seurat_obj,
    features = features_list,
    name = score_prefix,
    assay = assay
  )

  meta <- seurat_obj@meta.data
  score_cols <- grep(paste0("^", score_prefix), colnames(meta), value = TRUE)
  if (!length(score_cols)) {
    stop("No module score columns found with prefix: ", score_prefix)
  }

  sample_sym <- rlang::sym(sample_col)
  cluster_sym <- rlang::sym(cluster_col)

  agg_tbl <- meta |>
    dplyr::select(!!sample_sym, !!cluster_sym, dplyr::all_of(score_cols)) |>
    dplyr::group_by(!!sample_sym, !!cluster_sym) |>
    dplyr::summarise(
      dplyr::across(dplyr::all_of(score_cols), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(score_cols),
      names_to = "signature",
      values_to = "value"
    ) |>
    dplyr::mutate(
      feature = paste(.data[[cluster_col]], .data$signature, sep = "__")
    ) |>
    tidyr::pivot_wider(
      id_cols = !!sample_sym,
      names_from = .data$feature,
      values_from = .data$value,
      values_fill = 0
    ) |>
    dplyr::arrange(!!sample_sym)

  sig_mat <- agg_tbl |>
    dplyr::select(-dplyr::all_of(sample_col)) |>
    as.matrix()
  rownames(sig_mat) <- agg_tbl[[sample_col]]

  list(
    matrix = sig_mat,
    seurat_obj = seurat_obj,
    map = sig_res$map,
    score_columns = score_cols
  )
}

#' Aggregate Latent Embeddings by Patient and Cluster
#'
#' @param seurat_obj Seurat object with requested reduction stored in
#'   `Embeddings(seurat_obj, reduction)`.
#' @param reduction Reduction name (default `"integrated.scvi"`).
#' @param dims Optional integer vector selecting latent dimensions.
#' @param sample_col Column with sample/patient identifiers.
#' @param cluster_col Column with cluster identities.
#' @param aggregator Function applied within each patient×cluster cell (default
#'   `mean`).
#'
#' @return Numeric matrix with patients as rows and latent-dim__cluster columns.
#' @export
compute_patient_latent_matrix <- function(seurat_obj,
                                          reduction = "integrated.scvi",
                                          dims = NULL,
                                          sample_col = "hos_no",
                                          cluster_col = "anno3.scvi",
                                          aggregator = mean) {
  stopifnot(inherits(seurat_obj, "Seurat"))
  meta <- seurat_obj@meta.data
  required <- c(sample_col, cluster_col)
  missing_cols <- setdiff(required, colnames(meta))
  if (length(missing_cols)) {
    stop("Missing columns in metadata: ", paste(missing_cols, collapse = ", "))
  }

  latent <- Seurat::Embeddings(seurat_obj, reduction = reduction)
  if (!is.null(dims)) {
    latent <- latent[, dims, drop = FALSE]
  }

  latent_cols <- colnames(latent)
  latent_df <- cbind(
    meta[, c(sample_col, cluster_col), drop = FALSE],
    as.data.frame(latent, stringsAsFactors = FALSE)
  )
  latent_df[[sample_col]] <- as.character(latent_df[[sample_col]])
  latent_df[[cluster_col]] <- as.character(latent_df[[cluster_col]])

  sample_sym <- rlang::sym(sample_col)
  cluster_sym <- rlang::sym(cluster_col)

  latent_tbl <- latent_df |>
    dplyr::group_by(!!sample_sym, !!cluster_sym) |>
    dplyr::summarise(
      dplyr::across(dplyr::all_of(latent_cols), ~ aggregator(.x)),
      .groups = "drop"
    ) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(latent_cols),
      names_to = "latent",
      values_to = "value"
    ) |>
    dplyr::mutate(feature = paste(.data$latent, .data[[cluster_col]], sep = "__")) |>
    tidyr::pivot_wider(
      id_cols = !!sample_sym,
      names_from = .data$feature,
      values_from = .data$value,
      values_fill = 0
    ) |>
    dplyr::arrange(!!sample_sym)

  latent_mat <- latent_tbl |>
    dplyr::select(-dplyr::all_of(sample_col)) |>
    as.matrix()
  rownames(latent_mat) <- latent_tbl[[sample_col]]

  latent_mat
}

#' Scale a Feature Block with Optional Weighting
#'
#' @param M Numeric matrix.
#' @param weight Scalar applied after scaling.
#' @param center,scale Logical flags passed to `scale()`.
#' @param na_value Replacement value for `NA`s after scaling.
#'
#' @return Weighted, scaled matrix.
#' @export
scale_feature_block <- function(M,
                                weight = 1,
                                center = TRUE,
                                scale = TRUE,
                                na_value = 0) {
  stopifnot(is.matrix(M))
  if (center || scale) {
    M <- base::scale(M, center = center, scale = scale)
  }
  M[is.na(M)] <- na_value
  weight * M
}

#' Combine Multiple Patient-Level Feature Views
#'
#' @param views Named list of matrices (rows = patients).
#' @param weights Named numeric vector with same names as `views` (default 1).
#' @param intersect Logical; if TRUE (default) keep only patients present in all
#'   views.
#'
#' @return List with `matrix`, `ids`, and `view_columns`.
#' @export
combine_patient_feature_views <- function(views,
                                          weights = NULL,
                                          intersect = TRUE) {
  if (!length(views)) {
    stop("`views` must contain at least one matrix.")
  }
  is_null <- vapply(views, is.null, logical(1))
  views <- views[!is_null]
  if (!length(views)) {
    stop("All views were NULL; nothing to combine.")
  }

  if (is.null(weights)) {
    weights <- rep(1, length(views))
  }
  if (is.null(names(weights))) {
    names(weights) <- names(views)
  }
  weights <- weights[names(views)]
  weights[is.na(weights)] <- 1

  if (isTRUE(intersect)) {
    ids <- Reduce(intersect, lapply(views, rownames))
  } else {
    ids <- unique(unlist(lapply(views, rownames)))
  }
  if (!length(ids)) {
    stop("No overlapping patient IDs across views.")
  }

  aligned <- lapply(seq_along(views), function(i) {
    view <- views[[i]]
    w <- weights[[i]]
    if (!all(ids %in% rownames(view))) {
      missing_ids <- setdiff(ids, rownames(view))
      filler <- matrix(0, nrow = length(missing_ids), ncol = ncol(view),
                       dimnames = list(missing_ids, colnames(view)))
      view <- rbind(view, filler)
    }
    view <- view[ids, , drop = FALSE]
    scale_feature_block(view, weight = w, center = TRUE, scale = TRUE)
  })

  combined <- do.call(cbind, aligned)
  colnames(combined) <- make.unique(colnames(combined), sep = "_")

  list(
    matrix = combined,
    ids = ids,
    view_columns = lapply(aligned, colnames)
  )
}

#' Optional Batch Correction Wrapper
#'
#' @param X Matrix (rows = patients).
#' @param meta Metadata data frame with at least `sample_col`.
#' @param sample_col Column mapping rows in `X` to metadata.
#' @param batch_var Optional batch column; if NULL, no correction performed.
#'
#' @return Batch-corrected matrix.
#' @export
remove_batch_effect_if <- function(X,
                                   meta,
                                   sample_col = "hos_no",
                                   batch_var = NULL) {
  if (is.null(batch_var)) {
    return(X)
  }
  if (!batch_var %in% colnames(meta)) {
    warning("`batch_var` not found in metadata; skipping batch correction.")
    return(X)
  }
  sample_sym <- rlang::sym(sample_col)
  batch_sym <- rlang::sym(batch_var)
  batch_df <- meta |>
    dplyr::select(!!sample_sym, !!batch_sym) |>
    dplyr::distinct() |>
    dplyr::mutate(.order = match(.data[[sample_col]], rownames(X))) |>
    dplyr::filter(!is.na(.data$.order)) |>
    dplyr::arrange(.data$.order) |>
    dplyr::select(-.data$.order)

  if (!nrow(batch_df)) {
    warning("No batch annotations matched rownames; skipping correction.")
    return(X)
  }

  batch_vec <- batch_df[[batch_var]]
  corrected <- limma::removeBatchEffect(t(X), batch = batch_vec)
  t(corrected)
}

#' Choose a Safe Number of Principal Components
#'
#' @param X Numeric matrix.
#' @param max_k Maximum PCs to compute.
#' @param fraction_of_rank Fraction of rank used as cap.
#' @param extra_headroom Padding PCs above elbow.
#' @param elbow Optional vector of eigenvalues.
#'
#' @return Integer number of PCs.
#' @export
choose_pca_k <- function(X,
                         max_k = 50,
                         fraction_of_rank = 0.3,
                         extra_headroom = 5,
                         elbow = NULL) {
  n <- nrow(X)
  p <- ncol(X)
  r <- min(n, p)
  k_cap <- max(1L, min(max_k, r - 1L))

  if (is.null(elbow)) {
    k_probe <- min(20L, k_cap)
    probe <- tryCatch(
      irlba::prcomp_irlba(X, n = k_probe, center = TRUE, scale. = FALSE),
      error = function(e) NULL
    )
    if (!is.null(probe)) {
      s2 <- (probe$sdev)^2
      dif <- diff(s2)
      elbow_idx <- which.max(-dif)
      k_elbow <- max(3L, elbow_idx + extra_headroom)
    } else {
      k_elbow <- 10L
    }
  } else {
    dif <- diff(elbow)
    elbow_idx <- which.max(-dif)
    k_elbow <- max(3L, elbow_idx + extra_headroom)
  }

  k_frac <- floor(fraction_of_rank * r)
  max(3L, min(k_cap, max(k_elbow, k_frac)))
}

#' Run PCA with a Safe Backend Choice
#'
#' @param X Numeric matrix (rows = samples).
#' @param k Number of PCs to compute.
#' @param center,scale. Passed to PCA function.
#' @param verbose Logical.
#'
#' @return List containing PCA results and the engine used.
#' @export
run_pca_safely <- function(X, k, center = TRUE, scale. = FALSE, verbose = FALSE) {
  r <- min(nrow(X), ncol(X))
  frac <- k / r
  if (isTRUE(verbose)) {
    cat(sprintf("[run_pca_safely] r=%d, k=%d (%.1f%% of rank)\n", r, k, 100 * frac))
  }
  use_full <- frac >= 0.35
  if (use_full) {
    if (isTRUE(verbose)) {
      cat("[run_pca_safely] Using stats::prcomp\n")
    }
    pca <- stats::prcomp(X, center = center, scale. = scale.)
    keep <- seq_len(min(k, ncol(pca$x)))
    pca$x <- pca$x[, keep, drop = FALSE]
    pca$rotation <- pca$rotation[, keep, drop = FALSE]
    pca$sdev <- pca$sdev[keep]
    pca$engine <- "prcomp"
  } else {
    if (isTRUE(verbose)) {
      cat("[run_pca_safely] Using irlba::prcomp_irlba\n")
    }
    pca <- irlba::prcomp_irlba(X, n = k, center = center, scale. = scale.)
    pca$engine <- "prcomp_irlba"
  }
  pca
}

#' Patient-Level Dimensionality Reduction Pipeline
#'
#' Assemble multi-view patient features (cluster frequencies, marker signatures,
#' latent embeddings), perform optional batch correction, and compute PCA+UMAP.
#'
#' @param seurat_obj Seurat object containing columns `sample_col`,
#'   `cluster_col`, optional `batch_var`, and optional metadata for plotting
#'   (e.g. `g3`, `nih_change`).
#' @param markers_df Optional `FindAllMarkers` result; required when
#'   `include_signatures = TRUE`.
#' @param include_frequency Logical; include frequency view (default TRUE).
#' @param include_signatures Logical; include marker signatures (default TRUE
#'   when `markers_df` is supplied).
#' @param include_latent Logical; include latent reduction view (default TRUE).
#' @param sample_col,cluster_col Metadata columns for patient/cluster.
#' @param batch_var Batch column used in `removeBatchEffect`; set `NULL` to skip.
#' @param frequency_weight,signature_weight,latent_weight Numeric weights applied
#'   to each view during block scaling.
#' @param reduction Reduction name for latent embeddings.
#' @param latent_dims Optional integer vector selecting latent dimensions.
#' @param latent_aggregator Function applied when aggregating latent embeddings.
#' @param max_k,max_umap_dims,n_neighbors,min_dist PCA/UMAP parameters.
#' @param extra_plot_cols Additional metadata columns to join onto the output
#'   plotting data.frame.
#' @param verbose Logical; emit progress messages.
#' @param ... Forwarded to `compute_patient_signature_matrix()`.
#'
#' @return List with embedding, plot_df, pca, combined_features, view_columns.
#' @export
patient_dimensionality_reduction <- function(seurat_obj,
                                             markers_df = NULL,
                                             include_frequency = TRUE,
                                             include_signatures = !is.null(markers_df),
                                             include_latent = TRUE,
                                             sample_col = "hos_no",
                                             cluster_col = "anno3.scvi",
                                             batch_var = "GEM",
                                             frequency_weight = 1.5,
                                             signature_weight = 1.0,
                                             latent_weight = 1.0,
                                             reduction = "integrated.scvi",
                                             latent_dims = NULL,
                                             latent_aggregator = mean,
                                             max_k = 50,
                                             max_umap_dims = 30,
                                             n_neighbors = 15,
                                             min_dist = 0.2,
                                             extra_plot_cols = c("g3", "nih_change"),
                                             verbose = FALSE,
                                             ...) {
  stopifnot(inherits(seurat_obj, "Seurat"))

  meta <- seurat_obj@meta.data
  required <- c(sample_col, cluster_col)
  missing_cols <- setdiff(required, colnames(meta))
  if (length(missing_cols)) {
    stop("Missing columns in metadata: ", paste(missing_cols, collapse = ", "))
  }
  sample_sym <- rlang::sym(sample_col)

  views <- list()

  if (isTRUE(include_frequency)) {
    if (isTRUE(verbose)) {
      cat("[patient_dimensionality_reduction] Computing frequency view...\n")
    }
    views$frequency <- compute_patient_cluster_frequency(
      seurat_obj = seurat_obj,
      sample_col = sample_col,
      cluster_col = cluster_col,
      transform = "clr"
    )
  }

  sig_info <- NULL
  if (isTRUE(include_signatures)) {
    if (is.null(markers_df)) {
      stop("`markers_df` is required when `include_signatures = TRUE`.")
    }
    if (isTRUE(verbose)) {
      cat("[patient_dimensionality_reduction] Computing signature view...\n")
    }
    sig_info <- compute_patient_signature_matrix(
      seurat_obj = seurat_obj,
      markers_df = markers_df,
      sample_col = sample_col,
      cluster_col = cluster_col,
      verbose = verbose,
      ...
    )
    seurat_obj <- sig_info$seurat_obj
    views$signatures <- sig_info$matrix
  }

  if (isTRUE(include_latent)) {
    if (isTRUE(verbose)) {
      cat("[patient_dimensionality_reduction] Computing latent view...\n")
    }
    views$latent <- compute_patient_latent_matrix(
      seurat_obj = seurat_obj,
      reduction = reduction,
      dims = latent_dims,
      sample_col = sample_col,
      cluster_col = cluster_col,
      aggregator = latent_aggregator
    )
  }

  weights <- c(
    frequency = frequency_weight,
    signatures = signature_weight,
    latent = latent_weight
  )
  combined <- combine_patient_feature_views(
    views = views,
    weights = weights,
    intersect = TRUE
  )

  X <- combined$matrix
  ids <- combined$ids

  if (isTRUE(verbose)) {
    cat("[patient_dimensionality_reduction] Applying batch correction...\n")
  }
  X_bcorr <- remove_batch_effect_if(
    X = X,
    meta = seurat_obj@meta.data,
    sample_col = sample_col,
    batch_var = batch_var
  )

  k <- choose_pca_k(X_bcorr, max_k = max_k)
  if (isTRUE(verbose)) {
    cat(sprintf("[patient_dimensionality_reduction] Chosen k = %d\n", k))
  }
  pca <- run_pca_safely(
    X = X_bcorr,
    k = k,
    center = TRUE,
    scale. = FALSE,
    verbose = verbose
  )

  umap_dims <- min(max_umap_dims, ncol(pca$x))
  if (isTRUE(verbose)) {
    cat(sprintf("[patient_dimensionality_reduction] UMAP on %d PCs\n", umap_dims))
  }
  emb <- uwot::umap(
    X = pca$x[, seq_len(umap_dims), drop = FALSE],
    n_neighbors = n_neighbors,
    min_dist = min_dist,
    metric = "euclidean"
  )
  rownames(emb) <- ids

  plot_df <- data.frame(
    UMAP1 = emb[, 1],
    UMAP2 = emb[, 2],
    stringsAsFactors = FALSE
  )
  plot_df[[sample_col]] <- as.character(ids)
  plot_df <- dplyr::relocate(plot_df, !!sample_sym, .before = dplyr::everything())

  meta_cols <- unique(c(sample_col, extra_plot_cols, batch_var))
  meta_cols <- intersect(meta_cols, colnames(seurat_obj@meta.data))
  if (length(meta_cols)) {
    meta_join <- seurat_obj@meta.data[, meta_cols, drop = FALSE] |>
      dplyr::distinct(!!sample_sym, .keep_all = TRUE)
    meta_join[[sample_col]] <- as.character(meta_join[[sample_col]])
    plot_df <- dplyr::left_join(plot_df, meta_join, by = sample_col)
  }

  list(
    embedding = emb,
    plot_df = plot_df,
    pca = pca,
    combined_features = X_bcorr,
    view_columns = combined$view_columns,
    signature_info = sig_info,
    settings = list(
      sample_col = sample_col,
      cluster_col = cluster_col,
      batch_var = batch_var,
      weights = weights,
      reduction = reduction,
      latent_dims = latent_dims,
      n_neighbors = n_neighbors,
      min_dist = min_dist,
      max_k = max_k,
      max_umap_dims = max_umap_dims
    )
  )
}

#' Plot patient-level UMAP embedding
#'
#' Generate a lightweight ggplot scatter of the patient-level embedding
#' returned by `patient_dimensionality_reduction()`.
#'
#' @param plot_df Data frame containing at least `UMAP1`, `UMAP2`, and a
#'   patient identifier column (default `hos_no`).
#' @param color_by Column name in `plot_df` used for point colouring
#'   (default: `"g3"`).
#' @param label Logical; if TRUE, attempt to label points using `label_col`
#'   (requires `ggrepel`).
#' @param label_col Column name for point labels (default: `"hos_no"`).
#' @param title Optional plot title. If NULL, a default title is generated.
#'
#' @return A `ggplot` object.
#' @export
plot_patient_umap <- function(plot_df,
                              color_by = "g3",
                              label = FALSE,
                              label_col = "hos_no",
                              title = NULL) {
  required_cols <- c("UMAP1", "UMAP2")
  missing <- setdiff(required_cols, colnames(plot_df))
  if (length(missing)) {
    stop("`plot_df` is missing required columns: ", paste(missing, collapse = ", "))
  }
  if (!color_by %in% colnames(plot_df)) {
    stop("Column `", color_by, "` not found in `plot_df`.")
  }

  if (is.null(title)) {
    title <- paste("Patient UMAP coloured by", color_by)
  }

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$UMAP1, y = .data$UMAP2)) +
    ggplot2::geom_point(ggplot2::aes(color = .data[[color_by]]), size = 2.2, alpha = 0.85) +
    ggplot2::labs(
      x = "UMAP1",
      y = "UMAP2",
      color = color_by,
      title = title
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_line(color = "grey90"),
      legend.position = "right"
    )

  if (label && label_col %in% colnames(plot_df)) {
    if (!requireNamespace("ggrepel", quietly = TRUE)) {
      warning("Package 'ggrepel' not installed; skipping labels.", call. = FALSE)
    } else {
      p <- p + ggrepel::geom_text_repel(
        ggplot2::aes(label = .data[[label_col]]),
        size = 3,
        max.overlaps = Inf,
        min.segment.length = 0
      )
    }
  }

  p
}



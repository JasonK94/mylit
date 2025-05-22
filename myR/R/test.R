.get_feature_vector <- function(seu, feature, assay = "RNA", slot = "scale.data") {
  # Return numeric vector (cells) for either a gene (row in assay) or metadata col
  if (feature %in% colnames(seu@meta.data)) {
    v <- seu@meta.data[[feature]]
  } else if (feature %in% rownames(seu[[assay]])) {
    v <- GetAssayData(seu, assay = assay, slot = slot)[feature, ]
  } else {
    stop(feature, " not found in meta.data nor in assay rows.")
  }
  return(as.numeric(v))
}

.get_feature_vec <- function(obj, feature, assay = "RNA", slot = "scale.data") {
  if (feature %in% rownames(obj)) {
    return(as.numeric(GetAssayData(obj, assay = assay, slot = slot)[feature, ]))
  } else if (feature %in% colnames(obj@meta.data)) {
    return(obj@meta.data[[feature]])
  } else {
    stop(sprintf("Feature '%s' not found in genes or meta.data!", feature))
  }
}

# ===== 1. Scatter plot across cells ===========================================
#' Scatter plot of two cell‑level features (gene expression or metadata).
#'
#' @param seu Seurat object.
#' @param feature_x,feature_y Names of the x/y variables.
#' @param assay_x,assay_y Assays to pull gene expression from (ignored if metadata).
#' @param slot_x,slot_y Slots (layer) in assay to use (e.g. "scale.data", "data").
#' @param color_by Optional feature used for point colouring.
#' @param color_assay,color_slot Assay/slot for colour feature if gene.
#' @param regression One of "lm", "loess", or "lasso".
#' @param alpha Point transparency.
#' @param lambda Optional lambda for lasso. If NULL and regression=="lasso" performs cv.glmnet.
#' @param extra_predictors Optional character vector of additional cell‑level predictors for lasso.
#' @return ggplot2 object.
scatter_smooth_cells <- function(
    seu,
    feature_x,
    feature_y,
    assay_x = "RNA", slot_x = "scale.data",
    assay_y = "RNA", slot_y = "scale.data",
    color_by = NULL, color_assay = "RNA", color_slot = "scale.data",
    regression = c("lm", "loess", "lasso"),
    alpha = 0.5,
    lambda = NULL,
    extra_predictors = NULL) {
  
  regression <- match.arg(regression)
  df <- tibble(
    x = .get_feature_vector(seu, feature_x, assay_x, slot_x),
    y = .get_feature_vector(seu, feature_y, assay_y, slot_y))
  
  if (!is.null(color_by)) {
    df$col <- .get_feature_vector(seu, color_by, color_assay, color_slot)
  }
  
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_point(aes(colour = if (!is.null(color_by)) col else NULL), alpha = alpha) +
    labs(x = feature_x, y = feature_y, colour = color_by)
  
  # Smoothed line
  if (regression %in% c("lm", "loess")) {
    method <- ifelse(regression == "lm", "lm", "loess")
    p <- p + geom_smooth(method = method, se = FALSE, colour = "black")
  } else if (regression == "lasso") {
    # Build model matrix (at least 2 columns required by glmnet)
    X <- df[, c("x", extra_predictors), drop = FALSE] %>% as.matrix()
    if (ncol(X) < 2) {
      stop("glmnet needs >=2 predictors. Provide extra_predictors or use lm/loess.")
    }
    fit <- if (is.null(lambda)) {
      cv.glmnet(X, df$y, alpha = 1)
    } else {
      glmnet(X, df$y, alpha = 1, lambda = lambda)
    }
    beta <- as.numeric(coef(fit, s = if (is.null(lambda)) fit$lambda.min else lambda))
    # Predicted line (only varying x, others fixed at mean)
    x_seq <- seq(min(df$x), max(df$x), length.out = 100)
    Xnew <- cbind(x_seq, matrix(colMeans(X[, -1, drop = FALSE]), nrow = 100, ncol = ncol(X) - 1, byrow = TRUE))
    y_pred <- cbind(1, Xnew) %*% beta
    pred_df <- tibble(x = x_seq, y = as.numeric(y_pred))
    p <- p + geom_line(data = pred_df, colour = "black")
  }
  p + theme_classic()
}

# ===== 3. Scatter plot across genes ===========================================
#' Gene‑level scatter plot (each dot = gene).
#' Works with a tidy data.frame where each row is a gene, or generates such a
#' data.frame from a Seurat object by summarising across cells.
#'
#' @param data Seurat object or data.frame.
#' @param feature_x,feature_y Character names of columns / gene metrics.
#' @param color_by Optional column for colouring.
#' @param summarise_fun Function to turn cell matrix -> per‑gene metric if Seurat.
#'        Default = rowMeans of "scale.data".
#' @param regression "lm", "loess", or "lasso" (see scatter_smooth_cells).
scatter_smooth_genes <- function(
    data,
    feature_x,
    feature_y,
    color_by = NULL,
    regression = c("lm", "loess", "lasso"),
    summarise_fun = NULL,
    lambda = NULL,
    extra_predictors = NULL) {
  
  regression <- match.arg(regression)
  if (inherits(data, "Seurat")) {
    if (is.null(summarise_fun)) {
      summarise_fun <- function(mat) rowMeans(mat)
    }
    mat <- GetAssayData(data, slot = "scale.data")
    df <- tibble(
      gene = rownames(mat),
      x    = summarise_fun(mat[feature_x, , drop = FALSE]),
      y    = summarise_fun(mat[feature_y, , drop = FALSE])
    )
    if (!is.null(color_by)) {
      df$col <- summarise_fun(mat[color_by, , drop = FALSE])
    }
  } else {
    df <- as_tibble(data) %>% rename(x = !!feature_x, y = !!feature_y)
    if (!is.null(color_by)) df$col <- data[[color_by]]
  }
  
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_point(aes(colour = if (!is.null(color_by)) col else NULL), alpha = 0.7) +
    labs(x = feature_x, y = feature_y, colour = color_by)
  
  if (regression %in% c("lm", "loess")) {
    p <- p + geom_smooth(method = regression, se = FALSE, colour = "black")
  } else {
    X <- df[, c("x", extra_predictors), drop = FALSE] %>% as.matrix()
    if (ncol(X) < 2) {
      stop("glmnet needs >=2 predictors. Provide extra_predictors or use lm/loess.")
    }
    fit <- if (is.null(lambda)) cv.glmnet(X, df$y, alpha = 1) else glmnet(X, df$y, alpha = 1, lambda = lambda)
    beta <- as.numeric(coef(fit, s = if (is.null(lambda)) fit$lambda.min else lambda))
    x_seq <- seq(min(df$x), max(df$x), length.out = 100)
    Xnew <- cbind(x_seq, matrix(colMeans(X[, -1, drop = FALSE]), nrow = 100, ncol = ncol(X) - 1, byrow = TRUE))
    y_pred <- cbind(1, Xnew) %*% beta
    p <- p + geom_line(data = tibble(x = x_seq, y = y_pred), colour = "black")
  }
  p + theme_classic()
}

# ===== 4. Correlation between major feature & others ==========================
#' Compute correlation (Pearson/Spearman) of each feature with a major feature.
#' Handles optional grouping & splitting.
#'
#' @return A list (if split.by provided) or tibble.
corr_with_major <- function(
    data,
    major_feature,
    features,
    method = c("pearson", "spearman"),
    group.by = NULL,
    split.by = NULL,
    assay = "RNA", slot = "scale.data") {
  
  method <- match.arg(method)
  get_vec <- function(obj, feat) {
    if (inherits(obj, "Seurat")) .get_feature_vector(obj, feat, assay, slot) else obj[[feat]]
  }
  cor_fun <- function(v1, v2) {
    test <- cor.test(v1, v2, method = method)
    tibble(cor = test$estimate[[1]], p.value = test$p.value, n = length(v1))
  }
  
  if (is.null(split.by)) {
    df <- if (is.null(group.by)) {
      # No grouping, use raw cells
      map_dfr(features, ~ cor_fun(get_vec(data, major_feature), get_vec(data, .x)) %>% mutate(feature = .x), .id = NULL)
    } else {
      # Aggregate by group.by, then correlate across groups
      if (inherits(data, "Seurat")) {
        meta <- data@meta.data
        meta$..row <- rownames(meta)
      } else {
        meta <- data
      }
      meta <- meta %>% group_by(across(all_of(group.by))) %>% mutate(group_id = cur_group_id())
      agg <- map_dfc(c(major_feature, features), function(f) {
        tapply(get_vec(data, f), meta$group_id, mean, na.rm = TRUE)
      }) %>% as.data.frame()
      colnames(agg) <- c(major_feature, features)
      map_dfr(features, ~ cor_fun(agg[[major_feature]], agg[[.x]]) %>% mutate(feature = .x))
    }
    return(df)
  } else {
    splits <- if (inherits(data, "Seurat")) data@meta.data[[split.by]] else data[[split.by]]
    if (length(unique(splits)) > 10) warning("split.by has >10 categories; plots may be crowded.")
    res <- split(seq_along(splits), splits) %>% map(function(idx) {
      dsub <- if (inherits(data, "Seurat")) subset(data, cells = rownames(data)[idx]) else data[idx, ]
      corr_with_major(dsub, major_feature, features, method, group.by, split.by = NULL, assay, slot) %>% mutate(split = unique(splits[idx])[1])
    })
    return(res)
  }
}

# ===== 5. Full pairwise correlation ==========================================
#' Pairwise correlation among a set of features.
#' Optionally grouped or split (same semantics as corr_with_major).
#'
#' @return A list of correlation matrices (if split) or a single tibble of pairs.
corr_pairwise <- function(
    data,
    features,
    method = c("pearson", "spearman"),
    group.by = NULL,
    split.by = NULL,
    assay = "RNA", slot = "scale.data") {
  
  method <- match.arg(method)
  get_vec <- function(obj, feat) if (inherits(obj, "Seurat")) .get_feature_vector(obj, feat, assay, slot) else obj[[feat]]
  
  calc_mat <- function(indices = NULL) {
    vals <- lapply(features, function(f) get_vec(data, f)[indices])
    mat <- do.call(cbind, vals)
    colnames(mat) <- features
    cor(mat, method = method, use = "pairwise.complete.obs")
  }
  
  if (is.null(split.by)) {
    mat <- if (is.null(group.by)) {
      calc_mat()
    } else {
      # group‑aggregate means first
      splits <- if (inherits(data, "Seurat")) data@meta.data[[group.by]] else data[[group.by]]
      t(apply(sapply(features, function(f) tapply(get_vec(data, f), splits, mean, na.rm = TRUE)), 1, as.numeric)) %>%
        cor(method = method)
    }
    return(mat)
  } else {
    cat_vec <- if (inherits(data, "Seurat")) data@meta.data[[split.by]] else data[[split.by]]
    split(seq_along(cat_vec), cat_vec) %>% map(~ calc_mat(.x))
  }
}

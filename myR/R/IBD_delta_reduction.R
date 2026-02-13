#' Plot matched delta vectors in a 2D embedding space (e.g., PCA/UMAP)
#'
#' @param sobj Seurat object
#' @param match_id character vector of meta.data columns that uniquely match level1 vs level2 (e.g. c("emrid","ck_plot","drugno"))
#' @param compare_var meta.data column used to define the two levels to compare (default "treatment")
#' @param levels length-2 character vector, c(level1, level2) (e.g. c("pre","post"))
#' @param shape meta.data column mapped to point shape (default "drug")
#' @param color meta.data column mapped to point color (default "response")
#' @param reduction Seurat reduction name (default "pca")
#' @param axes either character vector of embedding column names (e.g. c("PC_1","PC_2")) or numeric vector of dims (e.g. c(1,2))
#' @param alpha_level1 alpha for level1 points (e.g. pre) (default 0.5)
#' @param alpha_arrows alpha for arrows (default 1; set <1 for transparency)
#' @param arrow_linetype linetype for arrows (default "dashed")
#' @param mean_group character vector of columns to group for mean arrows; default uses c(shape, color) if present plus other stable ids
#' @param draw_individual_arrows draw arrows for each matched pair (default TRUE)
#' @param draw_mean_arrows draw mean arrows for groups (default TRUE)
#' @param point_size point size for individual points
#' @param mean_point_size point size for mean post points
#' @return ggplot object
plot_matched_delta2d <- function(
  sobj,
  match_id = c("emrid", "ck_plot", "drugno"),
  compare_var = "treatment",
  levels = c("pre", "post"),
  shape = "drug",
  color = "response",
  reduction = "pca",
  axes = c(1, 2),
  alpha_level1 = 0.5,
  alpha_arrows = 1,
  arrow_linetype = "dashed",
  mean_group = NULL,
  draw_individual_arrows = TRUE,
  draw_mean_arrows = TRUE,
  point_size = 2,
  mean_point_size = 4
) {
    stopifnot(length(levels) == 2)
    if (!inherits(sobj, "Seurat")) stop("sobj must be a Seurat object.")
    if (!reduction %in% names(sobj@reductions)) stop(sprintf("Reduction '%s' not found in sobj@reductions.", reduction))

    md <- sobj@meta.data
    needed_md <- unique(c(match_id, compare_var, shape, color))
    missing_md <- setdiff(needed_md, colnames(md))
    if (length(missing_md) > 0) stop("Missing meta.data columns: ", paste(missing_md, collapse = ", "))

    emb_mat <- Seurat::Embeddings(sobj, reduction = reduction)
    if (is.numeric(axes)) {
        if (length(axes) != 2) stop("If axes is numeric, it must have length 2 (e.g. c(1,2)).")
        if (max(axes) > ncol(emb_mat)) stop("axes index exceeds available embedding dims.")
        xname <- colnames(emb_mat)[axes[1]]
        yname <- colnames(emb_mat)[axes[2]]
    } else {
        if (length(axes) != 2) stop("If axes is character, it must have length 2 (e.g. c('PC_1','PC_2')).")
        if (!all(axes %in% colnames(emb_mat))) stop("Specified axes not found in embedding columns: ", paste(axes, collapse = ", "))
        xname <- axes[1]
        yname <- axes[2]
    }

    emb <- as.data.frame(emb_mat[, c(xname, yname), drop = FALSE])
    emb$cell <- rownames(emb)

    md2 <- md %>%
        tibble::rownames_to_column("cell") %>%
        dplyr::select(dplyr::all_of(c("cell", match_id, compare_var, shape, color)))

    df <- dplyr::left_join(md2, emb, by = "cell") %>%
        dplyr::rename(x = !!xname, y = !!yname)

    # Split by levels
    level1 <- levels[1]
    level2 <- levels[2]

    pre <- df %>%
        dplyr::filter(.data[[compare_var]] == level1) %>%
        dplyr::select(dplyr::all_of(c(match_id, shape, color)), x_pre = x, y_pre = y)

    post <- df %>%
        dplyr::filter(.data[[compare_var]] == level2) %>%
        dplyr::select(dplyr::all_of(c(match_id, shape, color)), x_post = x, y_post = y)

    # Join matched pairs
    by_keys <- unique(c(match_id, shape, color))
    delta <- dplyr::inner_join(pre, post, by = by_keys)

    if (nrow(delta) == 0) {
        stop("No matched pairs found. Check match_id/compare_var/levels and that matches are 1:1.")
    }

    # Mean arrows grouping (default: drug_ck_r == shape + color; keep ck_plot if present to avoid mixing)
    if (is.null(mean_group)) {
        mean_group <- unique(c(shape, color, intersect(match_id, c("ck_plot"))))
        mean_group <- mean_group[!is.na(mean_group) & nzchar(mean_group)]
    }
    missing_mg <- setdiff(mean_group, colnames(delta))
    if (length(missing_mg) > 0) stop("mean_group columns not found: ", paste(missing_mg, collapse = ", "))

    delta_mean <- delta %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(mean_group))) %>%
        dplyr::summarise(
            x_pre = mean(x_pre, na.rm = TRUE),
            y_pre = mean(y_pre, na.rm = TRUE),
            x_post = mean(x_post, na.rm = TRUE),
            y_post = mean(y_post, na.rm = TRUE),
            n_pairs = dplyr::n(),
            .groups = "drop"
        )

    # Build plot
    p <- ggplot2::ggplot()

    if (draw_individual_arrows) {
        p <- p +
            ggplot2::geom_segment(
                data = delta,
                ggplot2::aes(x = x_pre, y = y_pre, xend = x_post, yend = y_post, color = .data[[color]]),
                alpha = alpha_arrows,
                linewidth = 0.5,
                linetype = arrow_linetype,
                arrow = ggplot2::arrow(length = grid::unit(0.14, "cm"))
            )
    }

    # level1 points (transparent)
    p <- p +
        ggplot2::geom_point(
            data = delta,
            ggplot2::aes(x = x_pre, y = y_pre, color = .data[[color]], shape = .data[[shape]]),
            alpha = alpha_level1,
            size = point_size
        ) +
        ggplot2::geom_point(
            data = delta,
            ggplot2::aes(x = x_post, y = y_post, color = .data[[color]], shape = .data[[shape]]),
            alpha = 0.95,
            size = point_size
        )

    if (draw_mean_arrows) {
        # Mean arrows (same color mapping; shape mapping uses mean_group's shape if present)
        has_shape_in_mean <- shape %in% mean_group
        has_color_in_mean <- color %in% mean_group

        p <- p +
            ggplot2::geom_segment(
                data = delta_mean,
                ggplot2::aes(
                    x = x_pre, y = y_pre, xend = x_post, yend = y_post,
                    color = if (has_color_in_mean) .data[[color]] else NULL
                ),
                linewidth = 1.2,
                linetype = "solid",
                arrow = ggplot2::arrow(length = grid::unit(0.22, "cm"))
            )

        p <- p +
            ggplot2::geom_point(
                data = delta_mean,
                ggplot2::aes(
                    x = x_post, y = y_post,
                    color = if (has_color_in_mean) .data[[color]] else NULL,
                    shape = if (has_shape_in_mean) .data[[shape]] else NULL
                ),
                size = mean_point_size,
                alpha = 1
            )
    }

    p +
        ggplot2::coord_equal() +
        ggplot2::theme_classic() +
        ggplot2::labs(
            x = xname, y = yname,
            color = color, shape = shape,
            title = sprintf("%s → %s delta in %s (%s vs %s)", level1, level2, reduction, compare_var, compare_var)
        )
}

#' Plot matched delta vectors from the origin (compare direction & magnitude only)
#'
#' - Individual arrows: (0,0) -> (x_post - x_pre, y_post - y_pre)
#' - Mean arrows (by mean_group): (0,0) -> (mean delta)
#' - Point(s): only level2 endpoints are plotted at delta coordinates (optional but default yes)
#'
#' @param sobj Seurat object
#' @param match_id meta.data columns used to 1:1 match (e.g. c("emrid","ck_plot","drugno"))
#' @param compare_var column defining two levels (default "treatment")
#' @param levels c(level1, level2) (e.g. c("pre","post"))
#' @param shape meta.data column mapped to arrowhead shape AND endpoint point shape (default "drug")
#' @param color meta.data column mapped to color (default "response")
#' @param reduction reduction name (default "pca")
#' @param axes numeric dims (e.g. c(1,2)) or character names (e.g. c("PC_1","PC_2"))
#' @param alpha_arrows alpha for individual arrows (default 0.3)
#' @param alpha_mean_arrow alpha for mean arrows (default 0.7)
#' @param arrow_linetype linetype for arrows (default "dashed")
#' @param mean_group grouping columns for mean arrows; default c(shape,color, intersect(match_id,"ck_plot"))
#' @param draw_individual_arrows whether to draw individual arrows (default TRUE)
#' @param draw_mean_arrows whether to draw mean arrows (default TRUE)
#' @param point_size endpoint (level2 delta endpoint) point size
#' @param mean_point_size mean endpoint point size
#' @param arrow_size line width for individual arrows
#' @param mean_arrow_size line width for mean arrows
#' @return ggplot object
plot_matched_delta2d_origin <- function(
  sobj,
  match_id = c("emrid", "ck_plot", "drugno"),
  compare_var = "treatment",
  levels = c("pre", "post"),
  shape = "drug",
  color = "response",
  reduction = "pca",
  axes = c(1, 2),
  alpha_arrows = 0.3,
  alpha_mean_arrow = 0.7,
  arrow_linetype = "dashed",
  mean_group = NULL,
  draw_individual_arrows = TRUE,
  draw_mean_arrows = TRUE,
  point_size = 2,
  mean_point_size = 4,
  arrow_size = 0.5,
  mean_arrow_size = 1.2
) {
    require(Seurat)
    require(dplyr)
    require(tibble)
    require(ggplot2)
    require(grid)

    stopifnot(length(levels) == 2)
    if (!inherits(sobj, "Seurat")) stop("sobj must be a Seurat object.")
    if (!reduction %in% names(sobj@reductions)) stop(sprintf("Reduction '%s' not found in sobj@reductions.", reduction))

    md <- sobj@meta.data
    needed_md <- unique(c(match_id, compare_var, shape, color))
    missing_md <- setdiff(needed_md, colnames(md))
    if (length(missing_md) > 0) stop("Missing meta.data columns: ", paste(missing_md, collapse = ", "))

    emb_mat <- Seurat::Embeddings(sobj, reduction = reduction)

    # Resolve axes -> xname/yname
    if (is.numeric(axes)) {
        if (length(axes) != 2) stop("If axes is numeric, it must have length 2 (e.g. c(1,2)).")
        if (max(axes) > ncol(emb_mat)) stop("axes index exceeds available embedding dims.")
        xname <- colnames(emb_mat)[axes[1]]
        yname <- colnames(emb_mat)[axes[2]]
    } else {
        if (length(axes) != 2) stop("If axes is character, it must have length 2 (e.g. c('PC_1','PC_2')).")
        if (!all(axes %in% colnames(emb_mat))) stop("Specified axes not found in embedding columns: ", paste(axes, collapse = ", "))
        xname <- axes[1]
        yname <- axes[2]
    }

    emb <- as.data.frame(emb_mat[, c(xname, yname), drop = FALSE])
    emb$cell <- rownames(emb)

    md2 <- md %>%
        tibble::rownames_to_column("cell") %>%
        dplyr::select(dplyr::all_of(c("cell", match_id, compare_var, shape, color)))

    df <- dplyr::left_join(md2, emb, by = "cell") %>%
        dplyr::rename(x = !!xname, y = !!yname)

    lvl1 <- levels[1]
    lvl2 <- levels[2]

    level1_df <- df %>%
        dplyr::filter(.data[[compare_var]] == lvl1) %>%
        dplyr::select(dplyr::all_of(c(match_id, shape, color)), x1 = x, y1 = y)

    level2_df <- df %>%
        dplyr::filter(.data[[compare_var]] == lvl2) %>%
        dplyr::select(dplyr::all_of(c(match_id, shape, color)), x2 = x, y2 = y)

    # Join matched pairs
    by_keys <- unique(c(match_id, shape, color))
    delta <- dplyr::inner_join(level1_df, level2_df, by = by_keys)

    if (nrow(delta) == 0) stop("No matched pairs found. Check match_id/compare_var/levels and that matches are 1:1.")

    # Compute deltas (endpoints in origin-centered plot)
    delta <- delta %>%
        dplyr::mutate(
            dx = x2 - x1,
            dy = y2 - y1
        )

    # Default mean grouping: "drug_ck_r" analog = shape + color (+ ck_plot if present)
    if (is.null(mean_group)) {
        mean_group <- unique(c(shape, color, intersect(match_id, "ck_plot")))
        mean_group <- mean_group[!is.na(mean_group) & nzchar(mean_group)]
    }
    missing_mg <- setdiff(mean_group, colnames(delta))
    if (length(missing_mg) > 0) stop("mean_group columns not found: ", paste(missing_mg, collapse = ", "))

    delta_mean <- delta %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(mean_group))) %>%
        dplyr::summarise(
            dx = mean(dx, na.rm = TRUE),
            dy = mean(dy, na.rm = TRUE),
            n_pairs = dplyr::n(),
            .groups = "drop"
        )

    # ---- Arrowheads that vary by 'shape' ----
    # ggplot2::geom_segment's arrow() cannot map "type"/"ends" per-row.
    # Workaround: use ggplot2::geom_spoke (which supports per-row "arrow" via after_scale? no),
    # or split data by shape levels and add one geom_segment per level with a fixed arrow spec.
    #
    # We'll implement the split-by-shape approach:
    shape_levels <- sort(unique(as.character(delta[[shape]])))
    shape_levels_mean <- if (shape %in% mean_group) sort(unique(as.character(delta_mean[[shape]]))) else character(0)

    # Define a simple mapping from shape level -> arrowhead style
    # You can customize this mapping as you like.
    arrow_style_for <- function(s) {
        # Example heuristic: alternate closed/open arrowheads
        # Return list(length_cm, type)
        idx <- match(s, shape_levels)
        if (is.na(idx)) idx <- 1
        list(len = 0.14, type = if (idx %% 2 == 1) "closed" else "open")
    }
    arrow_style_for_mean <- function(s) {
        idx <- match(s, shape_levels_mean)
        if (is.na(idx)) idx <- 1
        list(len = 0.22, type = if (idx %% 2 == 1) "closed" else "open")
    }

    p <- ggplot2::ggplot() +
        ggplot2::coord_equal() +
        ggplot2::theme_classic() +
        ggplot2::labs(
            x = paste0("Δ", xname), y = paste0("Δ", yname),
            color = color, shape = shape,
            title = sprintf("%s → %s delta from origin (%s)", lvl1, lvl2, reduction)
        )

    # Individual arrows: add one layer per shape level to vary arrowhead type
    if (draw_individual_arrows) {
        for (s in shape_levels) {
            ds <- delta %>% dplyr::filter(.data[[shape]] == s)
            st <- arrow_style_for(s)
            p <- p +
                ggplot2::geom_segment(
                    data = ds,
                    ggplot2::aes(x = 0, y = 0, xend = dx, yend = dy, color = .data[[color]]),
                    alpha = alpha_arrows,
                    linewidth = arrow_size,
                    linetype = arrow_linetype,
                    arrow = ggplot2::arrow(length = grid::unit(st$len, "cm"), type = st$type)
                )
        }
    }

    # Endpoints (level2 deltas): shape is mapped by the same 'shape' variable
    # These are points at (dx,dy)
    p <- p +
        ggplot2::geom_point(
            data = delta,
            ggplot2::aes(x = dx, y = dy, color = .data[[color]], shape = .data[[shape]]),
            alpha = 0.95,
            size = point_size
        )

    # Mean arrows + mean endpoints
    if (draw_mean_arrows && nrow(delta_mean) > 0) {
        if (shape %in% mean_group) {
            for (s in shape_levels_mean) {
                dm <- delta_mean %>% dplyr::filter(.data[[shape]] == s)
                st <- arrow_style_for_mean(s)
                p <- p +
                    ggplot2::geom_segment(
                        data = dm,
                        ggplot2::aes(x = 0, y = 0, xend = dx, yend = dy, color = .data[[color]]),
                        alpha = alpha_mean_arrow,
                        linewidth = mean_arrow_size,
                        linetype = "solid",
                        arrow = ggplot2::arrow(length = grid::unit(st$len, "cm"), type = st$type)
                    ) +
                    ggplot2::geom_point(
                        data = dm,
                        ggplot2::aes(x = dx, y = dy, color = .data[[color]], shape = .data[[shape]]),
                        alpha = 1,
                        size = mean_point_size
                    )
            }
        } else {
            # If shape isn't in mean_group, draw mean arrows without per-shape arrowheads
            p <- p +
                ggplot2::geom_segment(
                    data = delta_mean,
                    ggplot2::aes(x = 0, y = 0, xend = dx, yend = dy, color = .data[[color]]),
                    alpha = alpha_mean_arrow,
                    linewidth = mean_arrow_size,
                    linetype = "solid",
                    arrow = ggplot2::arrow(length = grid::unit(0.22, "cm"), type = "closed")
                ) +
                ggplot2::geom_point(
                    data = delta_mean,
                    ggplot2::aes(x = dx, y = dy, color = .data[[color]]),
                    alpha = 1,
                    size = mean_point_size
                )
        }
    }

    p
}



reduction_transfer_anchors <- function(
  reference, query,
  reference.reduction = "pca",
  query.assay = DefaultAssay(query),
  reference.assay = DefaultAssay(reference),
  features = NULL,
  dims = 1:30,
  projected.reduction.name = "pca_refproj"
) {
    require(Seurat)

    if (is.null(features)) {
        # 가능한 넓게: 둘의 공통 유전자
        features <- intersect(rownames(reference[[reference.assay]]), rownames(query[[query.assay]]))
    }

    anchors <- FindTransferAnchors(
        reference = reference,
        query = query,
        reference.assay = reference.assay,
        query.assay = query.assay,
        features = features,
        reference.reduction = reference.reduction,
        reduction = "pcaproject",
        dims = dims
    )

    # MapQuery는 projected PCA embedding을 query에 넣어줌 (DimReduc로)
    query2 <- MapQuery(
        anchorset = anchors,
        reference = reference,
        query = query,
        reference.reduction = reference.reduction,
        reduction.model = NULL
    )

    # MapQuery가 만든 reductions 중 projected pca가 무엇인지 버전에 따라 이름이 달라질 수 있음.
    # 필요하면 names(query2@reductions) 확인해서 rename 해도 됨.
    query2
}

plot_centroid_delta2d <- function(
  sobj,
  group_vars = c("cell_type", "drug", "response"),
  compare_var = "treatment",
  levels = c("pre", "post"),
  reduction = "pca",
  axes = c(1, 2), # 또는 c("PC_1","PC_2")
  color = "response",
  shape = "drug",
  alpha_level1 = 0.5,
  arrow_linetype = "dashed",
  arrow_alpha = 1,
  point_size = 3,
  mean_point_size = 4
) {
    require(dplyr)
    require(ggplot2)
    require(tibble)
    require(Seurat)
    require(grid)

    stopifnot(length(levels) == 2)
    if (!reduction %in% names(sobj@reductions)) stop("reduction 없음: ", reduction)

    md <- sobj@meta.data
    need <- unique(c(group_vars, compare_var, color, shape))
    miss <- setdiff(need, colnames(md))
    if (length(miss)) stop("meta.data에 없는 컬럼: ", paste(miss, collapse = ", "))

    emb <- Embeddings(sobj, reduction = reduction)
    if (is.numeric(axes)) {
        xname <- colnames(emb)[axes[1]]
        yname <- colnames(emb)[axes[2]]
    } else {
        xname <- axes[1]
        yname <- axes[2]
    }

    df <- md %>%
        tibble::rownames_to_column("cell") %>%
        dplyr::select(cell, dplyr::all_of(need)) %>%
        dplyr::left_join(
            as.data.frame(emb[, c(xname, yname), drop = FALSE]) %>%
                tibble::rownames_to_column("cell") %>%
                dplyr::rename(x = !!xname, y = !!yname),
            by = "cell"
        )

    # centroid 계산: (group_vars + compare_var)별 평균좌표
    cent <- df %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(c(group_vars, compare_var)))) %>%
        dplyr::summarise(
            x = mean(x, na.rm = TRUE),
            y = mean(y, na.rm = TRUE),
            n_cells = dplyr::n(),
            .groups = "drop"
        )

    # pre/post centroid 매칭
    lvl1 <- levels[1]
    lvl2 <- levels[2]
    c1 <- cent %>%
        dplyr::filter(.data[[compare_var]] == lvl1) %>%
        dplyr::select(dplyr::all_of(group_vars), x1 = x, y1 = y, n1 = n_cells)
    c2 <- cent %>%
        dplyr::filter(.data[[compare_var]] == lvl2) %>%
        dplyr::select(dplyr::all_of(group_vars), x2 = x, y2 = y, n2 = n_cells)

    seg <- dplyr::inner_join(c1, c2, by = group_vars)
    if (nrow(seg) == 0) stop("매칭되는 centroid 쌍이 없음. group_vars/levels를 확인.")

    # plot
    ggplot() +
        geom_segment(
            data = seg,
            aes(x = x1, y = y1, xend = x2, yend = y2, color = .data[[color]]),
            alpha = arrow_alpha,
            linetype = arrow_linetype,
            linewidth = 1.1,
            arrow = arrow(length = unit(0.22, "cm"))
        ) +
        geom_point(
            data = seg,
            aes(x = x1, y = y1, color = .data[[color]], shape = .data[[shape]]),
            alpha = alpha_level1,
            size = point_size
        ) +
        geom_point(
            data = seg,
            aes(x = x2, y = y2, color = .data[[color]], shape = .data[[shape]]),
            alpha = 0.95,
            size = mean_point_size
        ) +
        coord_equal() +
        theme_classic() +
        labs(
            x = xname, y = yname, color = color, shape = shape,
            title = paste0("Centroid delta: ", lvl1, " → ", lvl2, " (", reduction, ")")
        )
}




run_pc_lmm <- function(
  sobj,
  reduction = "pca",
  npcs = 30,
  match_id = NULL, # (옵션) 단지 결과에 key 붙이거나 sanity check용
  compare_var = "treatment",
  levels = c("pre", "post"),
  fixed_effects = c("sex", "response"),
  random_effects = c("emrid"),
  interactions = c("drug:response"),
  formula = NULL, # character or formula; 있으면 이걸 우선 사용
  use_lmerTest = TRUE,
  p_adjust_method = "BH",
  verbose = TRUE
) {
    stopifnot(inherits(sobj, "Seurat"))
    stopifnot(length(levels) == 2)

    # packages
    if (!requireNamespace("lme4", quietly = TRUE)) stop("Please install 'lme4'.")
    has_lmerTest <- requireNamespace("lmerTest", quietly = TRUE)

    # embeddings
    if (!reduction %in% names(sobj@reductions)) stop("Reduction not found: ", reduction)
    emb <- Seurat::Embeddings(sobj, reduction = reduction)
    if (ncol(emb) < npcs) npcs <- ncol(emb)

    # meta
    md <- sobj@meta.data
    # restrict to requested compare levels
    if (!compare_var %in% colnames(md)) stop("compare_var not in meta.data: ", compare_var)
    keep <- md[[compare_var]] %in% levels
    if (!any(keep)) stop("No cells/samples with compare_var in specified levels.")
    md <- md[keep, , drop = FALSE]
    emb <- emb[rownames(md), , drop = FALSE]

    # make compare_var a factor with desired order
    md[[compare_var]] <- factor(md[[compare_var]], levels = levels)

    # ensure covariates exist
    needed <- unique(c(compare_var, fixed_effects, random_effects))
    missing_cols <- setdiff(needed, colnames(md))
    if (length(missing_cols) > 0) stop("Missing meta.data columns: ", paste(missing_cols, collapse = ", "))

    # build formula
    build_formula <- function(resp_name) {
        if (!is.null(formula)) {
            if (inherits(formula, "formula")) {
                f <- formula
            } else if (is.character(formula)) {
                f <- stats::as.formula(formula)
            } else {
                stop("formula must be NULL, a formula, or a character string.")
            }
            # replace LHS with resp_name
            rhs <- as.character(f)[3]
            stats::as.formula(paste(resp_name, "~", rhs))
        } else {
            # base fixed: compare_var + fixed_effects
            fixed_terms <- unique(c(compare_var, fixed_effects))
            fixed_str <- paste(fixed_terms, collapse = " + ")

            # interactions
            inter_str <- ""
            if (!is.null(interactions) && length(interactions) > 0) {
                inter_str <- paste(interactions, collapse = " + ")
            }

            # random
            rand_str <- ""
            if (!is.null(random_effects) && length(random_effects) > 0) {
                rand_str <- paste0("(1|", random_effects, ")", collapse = " + ")
            }

            rhs_parts <- c(fixed_str, inter_str, rand_str)
            rhs_parts <- rhs_parts[rhs_parts != ""]
            stats::as.formula(paste(resp_name, "~", paste(rhs_parts, collapse = " + ")))
        }
    }

    # prepare data frame for modeling
    dat <- md
    # attach PCs as columns PC_1..PC_n
    for (i in seq_len(npcs)) {
        dat[[paste0("PC_", i)]] <- emb[, i]
    }
    dat$.row_id <- rownames(dat)

    # optional: basic sanity check for match_id uniqueness per level
    if (!is.null(match_id)) {
        miss_mid <- setdiff(match_id, colnames(dat))
        if (length(miss_mid) > 0) stop("match_id columns missing in meta.data: ", paste(miss_mid, collapse = ", "))
        # check duplicates within each level
        dup_tbl <- dat %>%
            dplyr::as_tibble() %>%
            dplyr::count(dplyr::across(dplyr::all_of(c(match_id, compare_var)))) %>%
            dplyr::filter(n > 1)
        if (nrow(dup_tbl) > 0 && verbose) {
            message("[warn] Duplicated match_id within levels detected. This is OK for centroid/aggregate analyses, but NOT 1:1 matching.")
        }
    }

    # fit models
    results <- vector("list", npcs)
    names(results) <- paste0("PC_", seq_len(npcs))

    tidy_one <- function(fit, pc_name) {
        # extract fixed effects table
        if (use_lmerTest && has_lmerTest) {
            sm <- summary(fit) # lmerTest provides df/t/p for fixed effects
            coefs <- as.data.frame(sm$coefficients)
            coefs$term <- rownames(coefs)
            # standardize column names
            # lmerTest: Estimate, Std. Error, df, t value, Pr(>|t|)
            out <- dplyr::tibble(
                PC = pc_name,
                term = coefs$term,
                estimate = coefs$Estimate,
                std_error = coefs$`Std. Error`,
                df = if ("df" %in% colnames(coefs)) coefs$df else NA_real_,
                statistic = coefs$`t value`,
                p_value = if ("Pr(>|t|)" %in% colnames(coefs)) coefs$`Pr(>|t|)` else NA_real_
            )
        } else {
            sm <- summary(fit) # lme4: no p-values
            coefs <- as.data.frame(sm$coefficients)
            coefs$term <- rownames(coefs)
            out <- dplyr::tibble(
                PC = pc_name,
                term = coefs$term,
                estimate = coefs$Estimate,
                std_error = coefs$`Std. Error`,
                df = NA_real_,
                statistic = coefs$`t value`,
                p_value = NA_real_
            )
        }
        out
    }

    all_tidy <- list()

    for (i in seq_len(npcs)) {
        pc <- paste0("PC_", i)
        f <- build_formula(pc)
        if (verbose) message("Fitting ", pc, ": ", deparse(f))

        fit <- tryCatch(
            {
                if (use_lmerTest && has_lmerTest) {
                    lmerTest::lmer(f, data = dat, REML = FALSE)
                } else {
                    lme4::lmer(f, data = dat, REML = FALSE)
                }
            },
            error = function(e) e
        )

        results[[pc]] <- list(formula = f, fit = fit)

        if (!inherits(fit, "error")) {
            all_tidy[[pc]] <- tidy_one(fit, pc)
        } else {
            all_tidy[[pc]] <- dplyr::tibble(
                PC = pc, term = NA_character_, estimate = NA_real_, std_error = NA_real_,
                df = NA_real_, statistic = NA_real_, p_value = NA_real_
            )
            if (verbose) message("[error] ", pc, ": ", fit$message)
        }
    }

    summary_tbl <- dplyr::bind_rows(all_tidy)

    # adjust p-values within term across PCs? or globally?
    # default: global adjustment across all PC×term rows with non-NA p
    summary_tbl <- summary_tbl %>%
        dplyr::mutate(
            adj_p_value = {
                pv <- p_value
                adj <- rep(NA_real_, length(pv))
                ok <- !is.na(pv)
                adj[ok] <- stats::p.adjust(pv[ok], method = p_adjust_method)
                adj
            }
        )

    # helpful: isolate compare_var terms quickly (incl. interaction)
    # you can filter later: grepl(compare_var, term)
    list(
        results = results,
        summary = summary_tbl,
        settings = list(
            reduction = reduction, npcs = npcs,
            compare_var = compare_var, levels = levels,
            fixed_effects = fixed_effects,
            random_effects = random_effects,
            interactions = interactions,
            formula = formula,
            use_lmerTest = use_lmerTest && has_lmerTest,
            p_adjust_method = p_adjust_method
        )
    )
}

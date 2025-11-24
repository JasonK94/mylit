#' Plot TML Model Comparison Metrics
#'
#' Visualizes the performance metrics of different L2 models across cross-validation folds.
#'
#' @param tml_object A result object from [TML7()].
#' @param metric The metric to plot (e.g., "ROC", "Sens", "Spec"). Default is the best metric used in training.
#' @return A ggplot object.
#' @export
plot_tml_metrics <- function(tml_object, metric = NULL) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")

    if (is.null(tml_object$model_comparison)) {
        stop("No model comparison data found in tml_object.")
    }

    df <- tml_object$model_comparison$values
    if (is.null(metric)) {
        metric <- tml_object$best_metric_name
    }

    cols <- colnames(df)
    metric_cols <- cols[grep(paste0("~", metric, "$"), cols)]

    if (length(metric_cols) == 0) {
        stop(sprintf("Metric '%s' not found in model comparison results.", metric))
    }

    plot_data <- data.frame()
    for (col in metric_cols) {
        model_name <- sub(paste0("~", metric), "", col)
        vals <- df[[col]]
        plot_data <- rbind(plot_data, data.frame(
            Model = model_name,
            Value = vals,
            Fold = df$Resample
        ))
    }

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Model, y = Value, fill = Model)) +
        ggplot2::geom_boxplot(alpha = 0.7, outlier.shape = NA) +
        ggplot2::geom_jitter(width = 0.2, height = 0, size = 2, alpha = 0.6) +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = paste("TML Model Comparison:", metric), y = metric) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

    return(p)
}

#' Analyze TML Outlier Folds
#'
#' Identifies cross-validation folds with poor performance and provides insights.
#' Supports standard CV, LOGO, and Repeated CV.
#'
#' @param tml_object A result object from [TML7()].
#' @param metric Metric to evaluate (default: "Sens").
#' @param threshold_method Method to define outliers: "absolute" (use threshold value) or "iqr" (Q1 - 1.5*IQR).
#' @param threshold Value for "absolute" method (default: 0.5).
#' @return A list containing outlier details and potentially problematic samples/groups.
#' @export
analyze_tml_outliers <- function(tml_object, metric = "Sens", threshold_method = "absolute", threshold = 0.5) {
    df <- tml_object$model_comparison$values
    cols <- colnames(df)
    metric_cols <- cols[grep(paste0("~", metric, "$"), cols)]

    if (length(metric_cols) == 0) {
        warning(sprintf("Metric '%s' not found.", metric))
        return(NULL)
    }

    outliers <- list()

    for (col in metric_cols) {
        model_name <- sub(paste0("~", metric), "", col)
        vals <- df[[col]]

        if (threshold_method == "iqr") {
            q1 <- stats::quantile(vals, 0.25, na.rm = TRUE)
            iqr <- stats::IQR(vals, na.rm = TRUE)
            cutoff <- q1 - 1.5 * iqr
            bad_idx <- which(vals < cutoff | is.na(vals))
        } else {
            bad_idx <- which(vals < threshold | is.na(vals))
        }

        if (length(bad_idx) > 0) {
            folds <- df$Resample[bad_idx]
            outliers[[model_name]] <- data.frame(
                Fold = folds,
                Value = vals[bad_idx]
            )
        }
    }

    if (length(outliers) == 0) {
        message("No outliers found.")
        return(NULL)
    }

    # Analyze problematic groups/samples if cv_folds info is available
    problematic_groups <- NULL
    if (!is.null(tml_object$cv_folds) && !is.null(tml_object$cv_folds$indexOut)) {
        fold_names <- names(tml_object$cv_folds$indexOut)
        all_bad_folds <- unique(unlist(lapply(outliers, function(x) x$Fold)))

        problematic_groups <- list()
        for (f in all_bad_folds) {
            if (f %in% fold_names) {
                idx <- tml_object$cv_folds$indexOut[[f]]
                problematic_groups[[f]] <- list(
                    indices = idx,
                    count = length(idx)
                )
            }
        }
    }

    return(list(
        outliers = outliers,
        problematic_folds = problematic_groups,
        threshold_used = if (threshold_method == "iqr") "Q1 - 1.5*IQR" else threshold
    ))
}

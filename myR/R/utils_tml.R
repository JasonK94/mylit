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
#' Identifies cross-validation folds with poor performance.
#'
#' @param tml_object A result object from [TML7()].
#' @param metric Metric to evaluate (default: "Sens").
#' @param threshold Value below which a fold is considered an outlier (default: 0.5).
#' @return A list containing outlier details.
#' @export
analyze_tml_outliers <- function(tml_object, metric = "Sens", threshold = 0.5) {
    df <- tml_object$model_comparison$values
    cols <- colnames(df)
    metric_cols <- cols[grep(paste0("~", metric, "$"), cols)]

    outliers <- list()

    for (col in metric_cols) {
        model_name <- sub(paste0("~", metric), "", col)
        vals <- df[[col]]
        bad_idx <- which(vals < threshold | is.na(vals))

        if (length(bad_idx) > 0) {
            folds <- df$Resample[bad_idx]
            outliers[[model_name]] <- data.frame(
                Fold = folds,
                Value = vals[bad_idx]
            )
        }
    }

    if (length(outliers) == 0) {
        message("No outliers found below threshold.")
        return(NULL)
    }

    return(list(outliers = outliers))
}

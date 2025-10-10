#' Composition Plotting Functions
#'
#' This module provides functions for visualizing cell/cluster composition
#' and proportions across samples or conditions.
#'
#' @name composition_plots
NULL

#' Create Proportional Stacked Bar Plot
#'
#' Creates a proportional stacked bar graph of clusters across samples.
#'
#' @param sobj Seurat object
#' @param cluster_col Cluster identity column (default: "seurat_clusters")
#' @param sample_col Sample grouping column
#' @param colors Color palette (optional)
#' @param title Plot title
#' @param xlab X-axis label
#' @param ylab Y-axis label (default: "Proportion")
#'
#' @return ggplot2 object
#'
#' @export
cmb <- function(sobj, 
                cluster_col = "seurat_clusters",
                sample_col,
                colors = NULL,
                title = "Cluster Composition by Sample",
                xlab = "Sample",
                ylab = "Proportion") {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }
  
  # Validate inputs
  if (!cluster_col %in% colnames(sobj@meta.data)) {
    stop("cluster_col '", cluster_col, "' not found in metadata")
  }
  if (!sample_col %in% colnames(sobj@meta.data)) {
    stop("sample_col '", sample_col, "' not found in metadata")
  }
  
  # Calculate proportions
  comp_data <- table(sobj@meta.data[[sample_col]], sobj@meta.data[[cluster_col]])
  comp_prop <- prop.table(comp_data, margin = 1)
  
  # Convert to data frame for plotting
  plot_data <- as.data.frame(comp_prop)
  colnames(plot_data) <- c("Sample", "Cluster", "Proportion")
  
  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Sample, y = Proportion, fill = Cluster)) +
    ggplot2::geom_bar(stat = "identity", position = "stack") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = title,
      x = xlab,
      y = ylab
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )
  
  # Add colors if provided
  if (!is.null(colors)) {
    p <- p + ggplot2::scale_fill_manual(values = colors)
  }
  
  return(p)
}

#' Create Absolute Count Stacked Bar Plot
#'
#' Creates an absolute count stacked bar graph of clusters across samples.
#'
#' @param sobj Seurat object
#' @param cluster_col Cluster identity column (default: "seurat_clusters")
#' @param sample_col Sample grouping column
#' @param colors Color palette (optional)
#' @param title Plot title
#' @param xlab X-axis label
#' @param ylab Y-axis label (default: "Cell Count")
#'
#' @return ggplot2 object
#'
#' @export
acmb <- function(sobj,
                 cluster_col = "seurat_clusters",
                 sample_col,
                 colors = NULL,
                 title = "Cluster Counts by Sample",
                 xlab = "Sample",
                 ylab = "Cell Count") {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }
  
  # Validate inputs
  if (!cluster_col %in% colnames(sobj@meta.data)) {
    stop("cluster_col '", cluster_col, "' not found in metadata")
  }
  if (!sample_col %in% colnames(sobj@meta.data)) {
    stop("sample_col '", sample_col, "' not found in metadata")
  }
  
  # Calculate counts
  comp_data <- table(sobj@meta.data[[sample_col]], sobj@meta.data[[cluster_col]])
  
  # Convert to data frame for plotting
  plot_data <- as.data.frame(comp_data)
  colnames(plot_data) <- c("Sample", "Cluster", "Count")
  
  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Sample, y = Count, fill = Cluster)) +
    ggplot2::geom_bar(stat = "identity", position = "stack") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = title,
      x = xlab,
      y = ylab
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )
  
  # Add colors if provided
  if (!is.null(colors)) {
    p <- p + ggplot2::scale_fill_manual(values = colors)
  }
  
  return(p)
}

#' Create Cumulative Line Plot of Cluster Proportions
#'
#' Produces a cumulative line plot of per-cluster cell proportions.
#'
#' @param sobj Seurat object
#' @param cluster_col Cluster identity column (default: "seurat_clusters")
#' @param sample_col Sample grouping column
#' @param title Plot title
#' @param xlab X-axis label (default: "Sample")
#' @param ylab Y-axis label (default: "Cumulative Proportion")
#'
#' @return ggplot2 object
#'
#' @export
cml <- function(sobj,
                cluster_col = "seurat_clusters",
                sample_col,
                title = "Cumulative Cluster Proportions",
                xlab = "Sample",
                ylab = "Cumulative Proportion") {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }
  
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr package required")
  }
  
  # Validate inputs
  if (!cluster_col %in% colnames(sobj@meta.data)) {
    stop("cluster_col '", cluster_col, "' not found in metadata")
  }
  if (!sample_col %in% colnames(sobj@meta.data)) {
    stop("sample_col '", sample_col, "' not found in metadata")
  }
  
  # Calculate proportions
  comp_data <- table(sobj@meta.data[[sample_col]], sobj@meta.data[[cluster_col]])
  comp_prop <- prop.table(comp_data, margin = 1)
  
  # Convert to data frame
  plot_data <- as.data.frame(comp_prop)
  colnames(plot_data) <- c("Sample", "Cluster", "Proportion")
  
  # Sort samples (numerically if possible)
  numeric_test <- suppressWarnings(as.numeric(as.character(plot_data$Sample)))
  if (!all(is.na(numeric_test))) {
    plot_data$Sample <- factor(plot_data$Sample,
                               levels = as.character(sort(unique(numeric_test))))
  }
  
  # Calculate cumulative proportions
  plot_data <- plot_data %>%
    dplyr::group_by(Cluster) %>%
    dplyr::arrange(Sample, .by_group = TRUE) %>%
    dplyr::mutate(Cumulative = cumsum(Proportion)) %>%
    dplyr::ungroup()
  
  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Sample, y = Cumulative, 
                                               color = Cluster, group = Cluster)) +
    ggplot2::geom_line(size = 1) +
    ggplot2::geom_point(size = 2) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = title,
      x = xlab,
      y = ylab
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  
  return(p)
}

#' Create UpSet Plot for Gene Lists
#'
#' Creates an UpSet plot for visualizing intersections of multiple gene lists.
#'
#' @param gene_lists Named list of character vectors (gene lists)
#' @param min_size Minimum intersection size to display (default: 0)
#' @param nsets Number of sets to display (default: all)
#' @param nintersects Number of intersections to display (default: 40)
#' @param order_by Order by: "freq" or "degree" (default: "freq")
#' @param ... Additional arguments passed to UpSetR::upset
#'
#' @return UpSet plot
#'
#' @export
upset_gene_lists <- function(gene_lists,
                             min_size = 0,
                             nsets = NULL,
                             nintersects = 40,
                             order_by = "freq",
                             ...) {
  
  if (!requireNamespace("UpSetR", quietly = TRUE)) {
    stop("UpSetR package required. Install with: install.packages('UpSetR')")
  }
  
  if (!is.list(gene_lists)) {
    stop("gene_lists must be a named list of character vectors")
  }
  
  if (is.null(names(gene_lists)) || any(names(gene_lists) == "")) {
    stop("gene_lists must have names for all elements")
  }
  
  # Get all unique genes
  all_genes <- unique(unlist(gene_lists))
  
  # Create binary matrix
  upset_matrix <- matrix(0, 
                        nrow = length(all_genes), 
                        ncol = length(gene_lists),
                        dimnames = list(all_genes, names(gene_lists)))
  
  for (i in seq_along(gene_lists)) {
    list_name <- names(gene_lists)[i]
    genes_in_list <- gene_lists[[i]]
    upset_matrix[genes_in_list, list_name] <- 1
  }
  
  # Convert to data frame
  upset_df <- as.data.frame(upset_matrix)
  
  # Determine nsets if not specified
  if (is.null(nsets)) {
    nsets <- length(gene_lists)
  }
  
  # Create UpSet plot
  UpSetR::upset(
    upset_df,
    nsets = nsets,
    nintersects = nintersects,
    order.by = order_by,
    decreasing = TRUE,
    mb.ratio = c(0.6, 0.4),
    mainbar.y.label = "Intersection Size",
    sets.x.label = "Gene Set Size",
    ...
  )
}

#' Create Violin Plot with Statistical Comparisons
#'
#' Generates violin plots with statistical testing using ggpubr.
#'
#' @param sobj Seurat object
#' @param features Features to plot
#' @param group_by Grouping variable
#' @param split_by Optional splitting variable
#' @param test_method Statistical test: "wilcox", "t.test", "anova", "kruskal" (default: "wilcox")
#' @param comparisons List of comparisons (default: all pairwise)
#' @param ncol Number of columns for faceting (default: NULL)
#'
#' @return ggplot2 object
#'
#' @export
vln_p <- function(sobj,
                  features,
                  group_by,
                  split_by = NULL,
                  test_method = "wilcox",
                  comparisons = NULL,
                  ncol = NULL) {
  
  if (!requireNamespace("ggpubr", quietly = TRUE)) {
    stop("ggpubr package required. Install with: install.packages('ggpubr')")
  }
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }
  
  # Validate inputs
  if (!group_by %in% colnames(sobj@meta.data)) {
    stop("group_by '", group_by, "' not found in metadata")
  }
  
  if (!is.null(split_by) && !split_by %in% colnames(sobj@meta.data)) {
    stop("split_by '", split_by, "' not found in metadata")
  }
  
  # Prepare data
  plot_list <- list()
  
  for (feat in features) {
    # Extract feature values
    if (feat %in% colnames(sobj@meta.data)) {
      values <- sobj@meta.data[[feat]]
    } else if (feat %in% rownames(sobj)) {
      values <- Seurat::FetchData(sobj, vars = feat)[[1]]
    } else {
      warning("Feature '", feat, "' not found. Skipping.")
      next
    }
    
    # Create plot data
    plot_data <- data.frame(
      value = values,
      group = sobj@meta.data[[group_by]]
    )
    
    if (!is.null(split_by)) {
      plot_data$split <- sobj@meta.data[[split_by]]
    }
    
    # Remove NAs
    plot_data <- plot_data[complete.cases(plot_data), ]
    
    # Create plot
    if (is.null(split_by)) {
      p <- ggpubr::ggviolin(
        plot_data,
        x = "group",
        y = "value",
        fill = "group",
        add = "boxplot",
        add.params = list(fill = "white")
      ) +
      ggpubr::stat_compare_means(
        comparisons = comparisons,
        method = test_method
      ) +
      ggplot2::labs(
        title = feat,
        x = group_by,
        y = "Expression"
      ) +
      ggplot2::theme(legend.position = "none")
    } else {
      p <- ggpubr::ggviolin(
        plot_data,
        x = "group",
        y = "value",
        fill = "split",
        add = "boxplot",
        add.params = list(fill = "white"),
        position = ggplot2::position_dodge(0.8)
      ) +
      ggplot2::labs(
        title = feat,
        x = group_by,
        y = "Expression",
        fill = split_by
      )
    }
    
    plot_list[[feat]] <- p
  }
  
  # Combine plots if multiple features
  if (length(plot_list) == 1) {
    return(plot_list[[1]])
  } else if (length(plot_list) > 1) {
    if (requireNamespace("patchwork", quietly = TRUE)) {
      if (is.null(ncol)) {
        ncol <- ceiling(sqrt(length(plot_list)))
      }
      return(patchwork::wrap_plots(plot_list, ncol = ncol))
    } else {
      warning("patchwork package not available. Returning list of plots.")
      return(plot_list)
    }
  } else {
    stop("No valid features to plot")
  }
}


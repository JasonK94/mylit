#' Generate boxplots of average feature expression from a Seurat object.
#'
#' This function calculates the average expression of specified features
#' (genes or metadata) per sample, grouped by a primary grouping variable
#' (e.g., cluster) and a secondary splitting variable. It then generates
#' boxplots showing the distribution of these sample averages.
#'
#' @param sobj A Seurat object.
#' @param features A character vector of features to plot (genes or metadata column names).
#' @param sample_col Character string: the name of the metadata column identifying samples/patients.
#' @param group.by Character string: the name of the metadata column for primary grouping on the x-axis (e.g., clusters, cell types). If NULL, a single group "Overall" is used.
#' @param split.by Character string: the name of the metadata column for secondary grouping, creating side-by-side boxplots.
#' @param idents Optional: A character vector of values from the `group.by` column to include. If NULL, all values are shown.
#' @param assay Optional: Character string: the assay to use for gene features (e.g., "RNA", "SCT"). Defaults to the current default assay in `sobj`.
#' @param layer Optional: Character string: the layer to use for gene features ("counts", "data", "scale.data"). Defaults to "data", which `VlnPlot` typically uses.
#' @param ncol Integer: number of columns for arranging plots if multiple features are provided. Default is 3.
#' @param pt.size Numeric: size of points for individual sample averages. If 0, points are not shown. Default is 0.
#' @param violin Logical: whether to overlay violin plots. Default is FALSE.
#'
#' @return A ggplot object (created by `patchwork`) combining plots for all features.
#'
#' @import Seurat
#' @import dplyr
#' @import ggplot2
#' @import patchwork
#'
#' @examples
#' \dontrun{
#' # Assuming 'pbmc_small' is a Seurat object with 'RNA_snn_res.0.7' for clusters,
#' # 'letter.idents' for sample IDs (dummy), and 'groups' for a split.by condition.
#' # pbmc_small$sample <- sample(paste0("Patient", 1:4), ncol(pbmc_small), replace = TRUE)
#' # pbmc_small$condition <- sample(c("Control", "Treatment"), ncol(pbmc_small), replace = TRUE)
#'
#' # Plot average CD3E expression per sample, grouped by cluster, split by condition
#' # mybox(pbmc_small,
#' #       features = "CD3E",
#' #       sample_col = "sample",
#' #       group.by = "seurat_clusters", # or another metadata column like cell type
#' #       split.by = "condition",
#' #       pt.size = 1,
#' #       violin = TRUE)
#'
#' # Plot average nCount_RNA per sample, group.by is NULL (overall), split by condition
#' # mybox(pbmc_small,
#' #       features = "nCount_RNA",
#' #       sample_col = "sample",
#' #       group.by = NULL,
#' #       split.by = "condition",
#' #       pt.size = 1)
#' }
#' @export
mybox <- function(sobj, features, sample_col, group.by, split.by,
                  idents = NULL, assay = NULL, layer = NULL,
                  ncol = 3, pt.size = 0, violin = FALSE) {
  
  # --- Argument Checks and Setup ---
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' needed for this function to work. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' needed for this function to work. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' needed for this function to work. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' needed for this function to work. Please install it.", call. = FALSE)
  }
  
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  assay_to_use <- assay %||% Seurat::DefaultAssay(sobj)
  layer_to_use <- layer %||% "data"
  
  if (!(sample_col %in% names(sobj@meta.data))) {
    stop(paste("sample_col '", sample_col, "' not found in Seurat object metadata.", sep=""))
  }
  if (!is.null(group.by) && !(group.by %in% names(sobj@meta.data))) {
    stop(paste("group.by '", group.by, "' not found in Seurat object metadata.", sep=""))
  }
  if (!(split.by %in% names(sobj@meta.data))) {
    stop(paste("split.by '", split.by, "' not found in Seurat object metadata.", sep=""))
  }
  
  all_plots <- list()
  
  for (feature_name in features) {
    metadata_cols_to_fetch <- c(sample_col, split.by)
    if (!is.null(group.by)) {
      metadata_cols_to_fetch <- c(metadata_cols_to_fetch, group.by)
    }
    metadata_cols_to_fetch <- unique(metadata_cols_to_fetch)
    
    if (feature_name %in% names(sobj@meta.data)) {
      data_for_feature <- Seurat::FetchData(object = sobj, vars = c(feature_name, metadata_cols_to_fetch))
    } else {
      data_for_feature <- Seurat::FetchData(object = sobj, vars = c(feature_name, metadata_cols_to_fetch),
                                            assay = assay_to_use, layer = layer_to_use)
    }
    
    group.by.internal <- group.by
    if (is.null(group.by)) {
      group.by.internal <- ".internal_placeholder_group"
      data_for_feature[[group.by.internal]] <- "Overall"
    }
    
    if (!is.null(group.by) && !is.null(idents)) {
      if (!(group.by.internal %in% names(data_for_feature))) {
        stop(paste("group.by column '", group.by.internal, "' not found in fetched data.", sep=""))
      }
      original_col_type <- class(data_for_feature[[group.by.internal]])
      data_for_feature[[group.by.internal]] <- as.character(data_for_feature[[group.by.internal]])
      data_for_feature <- data_for_feature[data_for_feature[[group.by.internal]] %in% idents, ]
      
      if ("factor" %in% original_col_type) {
        original_levels <- levels(sobj@meta.data[[group.by]])
        valid_subset_levels <- intersect(original_levels, idents)
        if(length(valid_subset_levels) > 0){
          data_for_feature[[group.by.internal]] <- factor(data_for_feature[[group.by.internal]], levels = valid_subset_levels)
        } else {
          data_for_feature[[group.by.internal]] <- factor(data_for_feature[[group.by.internal]]) 
        }
      } else {
        data_for_feature[[group.by.internal]] <- as.factor(data_for_feature[[group.by.internal]])
      }
      if (nrow(data_for_feature) == 0) {
        warning(paste("No data remaining for feature", feature_name, "after filtering by idents. Skipping this feature."))
        next
      }
    }
    
    data_for_feature[[group.by.internal]] <- as.factor(data_for_feature[[group.by.internal]])
    data_for_feature[[split.by]] <- as.factor(data_for_feature[[split.by]])
    data_for_feature[[sample_col]] <- as.factor(data_for_feature[[sample_col]])
    
    grouping_vars_for_agg <- unique(c(sample_col, group.by.internal, split.by))
    
    aggregated_df <- data_for_feature %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(grouping_vars_for_agg))) %>%
      dplyr::summarise(mean_value = mean(.data[[feature_name]], na.rm = TRUE), .groups = 'drop')
    
    if (nrow(aggregated_df) == 0) {
      warning(paste("No data available for plotting feature", feature_name, "after aggregation. Skipping this feature."))
      next
    }
    
    x_axis_label <- if(is.null(group.by)) "Overall" else group.by
    
    p <- ggplot2::ggplot(aggregated_df, ggplot2::aes(x = .data[[group.by.internal]], y = mean_value, fill = .data[[split.by]]))
    
    if (violin) {
      p <- p + ggplot2::geom_violin(position = ggplot2::position_dodge(preserve = "single", width = 0.9),
                                    alpha = 0.5, scale = "width", trim = FALSE, na.rm = TRUE)
    }
    
    # MODIFIED outlier.shape logic:
    current_outlier_shape <- if(pt.size > 0) NA else 19 
    p <- p + ggplot2::geom_boxplot(position = ggplot2::position_dodge(preserve = "single", width = 0.9),
                                   na.rm = TRUE, outlier.shape = current_outlier_shape) 
    
    if (pt.size > 0) {
      p <- p + ggplot2::geom_jitter(ggplot2::aes(color = .data[[split.by]]),
                                    position = ggplot2::position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9),
                                    size = pt.size, show.legend = FALSE, na.rm = TRUE)
    }
    
    p <- p + ggplot2::labs(title = feature_name,
                           x = x_axis_label,
                           y = paste("Average", feature_name),
                           fill = split.by, color = split.by) +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     plot.title = ggplot2::element_text(hjust = 0.5))
    
    all_plots[[feature_name]] <- p
  }
  
  if (length(all_plots) == 0) {
    message("No plots were generated. Check input data and parameters.")
    return(invisible(NULL))
  }
  
  return(patchwork::wrap_plots(all_plots, ncol = ncol))
}

#' Generate boxplots of average feature values from a data frame.
#'
#' This function takes a data frame where each row is a minimal data unit.
#' It calculates the average of specified features per sample, grouped by a
#' primary grouping variable and a secondary splitting variable. It then generates
#' boxplots showing the distribution of these sample averages.
#'
#' @param df A data frame. Each row should represent a minimal data unit (e.g., a cell).
#' @param features A character vector of feature column names to plot.
#' @param sample_col Character string: the name of the column identifying samples/patients.
#' @param group.by Character string: the name of the column for primary grouping on the x-axis. If NULL, a single group "Overall" is used.
#' @param split.by Character string: the name of the column for secondary grouping, creating side-by-side boxplots.
#' @param idents Optional: A character vector of values from the `group.by` column to include. If NULL, all values are shown.
#' @param ncol Integer: number of columns for arranging plots if multiple features are provided. Default is 3.
#' @param pt.size Numeric: size of points for individual sample averages. If 0, points are not shown. Default is 0.
#' @param violin Logical: whether to overlay violin plots. Default is FALSE.
#'
#' @return A ggplot object (created by `patchwork`) combining plots for all features.
#'
#' @import dplyr
#' @import ggplot2
#' @import patchwork
#'
#' @examples
#' \dontrun{
#' # Create a sample data frame
#' set.seed(123)
#' n_cells <- 500
#' example_df <- data.frame(
#'   feature1 = rnorm(n_cells, 10, 2),
#'   feature2 = rpois(n_cells, 5),
#'   sample_id = sample(paste0("Sample", 1:5), n_cells, replace = TRUE),
#'   cell_type = sample(c("TypeA", "TypeB", "TypeC"), n_cells, replace = TRUE),
#'   treatment = sample(c("DrugX", "Placebo"), n_cells, replace = TRUE)
#' )
#'
#' # Plot average feature1 per sample, grouped by cell_type, split by treatment
#' # mybox_df(example_df,
#' #          features = "feature1",
#' #          sample_col = "sample_id",
#' #          group.by = "cell_type",
#' #          split.by = "treatment",
#' #          pt.size = 1.5,
#' #          violin = TRUE)
#'
#' # Plot multiple features, no specific grouping (overall), split by treatment
#' # mybox_df(example_df,
#' #          features = c("feature1", "feature2"),
#' #          sample_col = "sample_id",
#' #          group.by = NULL,
#' #          split.by = "treatment",
#' #          ncol = 2,
#' #          pt.size = 1)
#' }
mybox_df <- function(df, features, sample_col, group.by, split.by,
                     idents = NULL, ncol = 3, pt.size = 0, violin = FALSE) {
  
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' needed for this function to work. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' needed for this function to work. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' needed for this function to work. Please install it.", call. = FALSE)
  }
  
  if (!is.data.frame(df)) {
    stop("Input 'df' must be a data frame.", call. = FALSE)
  }
  for (col_name in c(sample_col, if(!is.null(group.by)) group.by, split.by, features)) {
    if (!(col_name %in% names(df))) {
      stop(paste("Column '", col_name, "' not found in the data frame.", sep=""))
    }
  }
  
  all_plots <- list()
  
  for (feature_name in features) {
    current_df <- df 
    group.by.internal <- group.by
    if (is.null(group.by)) {
      group.by.internal <- ".internal_placeholder_group"
      current_df[[group.by.internal]] <- "Overall"
    }
    
    if (!is.null(group.by) && !is.null(idents)) {
      if (!(group.by.internal %in% names(current_df))) {
        stop(paste("group.by column '", group.by.internal, "' not found in data frame.", sep=""))
      }
      original_col_type <- class(current_df[[group.by.internal]])
      current_df[[group.by.internal]] <- as.character(current_df[[group.by.internal]])
      current_df <- current_df[current_df[[group.by.internal]] %in% idents, ]
      
      if ("factor" %in% original_col_type) {
        original_levels_df <- levels(df[[group.by]]) 
        valid_subset_levels <- intersect(original_levels_df, idents)
        if(length(valid_subset_levels) > 0) {
          current_df[[group.by.internal]] <- factor(current_df[[group.by.internal]], levels = valid_subset_levels)
        } else {
          current_df[[group.by.internal]] <- as.factor(current_df[[group.by.internal]]) 
        }
      } else {
        current_df[[group.by.internal]] <- as.factor(current_df[[group.by.internal]])
      }
      if (nrow(current_df) == 0) {
        warning(paste("No data remaining for feature", feature_name, "after filtering by idents. Skipping this feature."))
        next
      }
    }
    
    current_df[[group.by.internal]] <- as.factor(current_df[[group.by.internal]])
    current_df[[split.by]] <- as.factor(current_df[[split.by]])
    current_df[[sample_col]] <- as.factor(current_df[[sample_col]])
    
    grouping_vars_for_agg <- unique(c(sample_col, group.by.internal, split.by))
    
    aggregated_df <- current_df %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(grouping_vars_for_agg))) %>%
      dplyr::summarise(mean_value = mean(.data[[feature_name]], na.rm = TRUE), .groups = 'drop')
    
    if (nrow(aggregated_df) == 0) {
      warning(paste("No data available for plotting feature", feature_name, "after aggregation. Skipping this feature."))
      next
    }
    
    x_axis_label <- if(is.null(group.by)) "Overall" else group.by
    
    p <- ggplot2::ggplot(aggregated_df, ggplot2::aes(x = .data[[group.by.internal]], y = mean_value, fill = .data[[split.by]]))
    
    if (violin) {
      p <- p + ggplot2::geom_violin(position = ggplot2::position_dodge(preserve = "single", width = 0.9),
                                    alpha = 0.5, scale = "width", trim = FALSE, na.rm = TRUE)
    }
    
    # MODIFIED outlier.shape logic:
    current_outlier_shape <- if(pt.size > 0) NA else 19 
    p <- p + ggplot2::geom_boxplot(position = ggplot2::position_dodge(preserve = "single", width = 0.9),
                                   na.rm = TRUE, outlier.shape = current_outlier_shape)
    
    if (pt.size > 0) {
      p <- p + ggplot2::geom_jitter(ggplot2::aes(color = .data[[split.by]]), 
                                    position = ggplot2::position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9),
                                    size = pt.size, show.legend = FALSE, na.rm = TRUE)
    }
    
    p <- p + ggplot2::labs(title = feature_name,
                           x = x_axis_label,
                           y = paste("Average", feature_name),
                           fill = split.by, color = split.by) +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     plot.title = ggplot2::element_text(hjust = 0.5))
    
    all_plots[[feature_name]] <- p
  } 
  
  if (length(all_plots) == 0) {
    message("No plots were generated. Check input data and parameters.")
    return(invisible(NULL))
  }
  
  return(patchwork::wrap_plots(all_plots, ncol = ncol))
}

#' Create a Heatmap of Gene Set Expression per Cluster
#'
#' This function creates a heatmap showing the normalized expression of gene sets across clusters.
#' The expression values are z-score normalized for better visualization.
#'
#' @param sobj A Seurat object
#' @param group Character string specifying the identity to use for clustering.
#'             Default is "seurat_clusters"
#' @param value Character string specifying how to aggregate expression values.
#'             Options: "average" (default) or "sum"
#' @param assay Character string specifying which assay to use.
#'             Default is "SCT"
#' @param gene_sets A list of character vectors, where each vector contains gene names.
#'                 The list can be named, and names will be used as gene set labels.
#'                 If a single vector is provided, it will be converted to a list.
#' @param title Character string for the plot title.
#'             Default is "Normalized Gene Set Expression per Cluster"
#' @param x_label Character string for the x-axis label.
#'               Default is "Cluster"
#' @param y_label Character string for the y-axis label.
#'               Default is "Gene Set"
#'
#' @return A data frame containing the normalized expression values for each gene set
#'         across clusters, with an additional column indicating the assigned cell type
#'         based on highest expression
#'
#' @examples
#' \dontrun{
#' # Create heatmap with named gene sets
#' gene_sets <- list(
#'   T_cells = c("CD3D", "CD3E", "CD4", "CD8A"),
#'   B_cells = c("MS4A1", "CD19", "CD79A"),
#'   Myeloid = c("CD14", "FCGR3A", "LYZ")
#' )
#' result <- myhm_genesets2(sobj, gene_sets = gene_sets)
#' 
#' # Create heatmap with a single gene set
#' result <- myhm_genesets2(sobj, gene_sets = c("CD3D", "CD3E", "CD4"))
#' }
#' @export
myhm_genesets2 <- function(
    sobj,
    group = "seurat_clusters",
    value = "average",
    assay = "SCT",
    gene_sets = NULL,
    title="Normalized Gene Set Expression per Cluster",
    x_label="Cluster",
    y_label="Gene Set"
){
  library(Seurat)
  library(dplyr)
  library(reshape2)
  library(ggplot2)
  
  #-------------------------------'
  # (A) 유효성 체크
  #-------------------------------'
  if(is.null(gene_sets)){
    stop("gene_sets를 지정해 주세요. (ex: list(Immune=c('CD3D','CD3E'), Bcell=c('MS4A1','CD79A'))) ")
  }
  
  # 만약 gene_sets가 리스트가 아니라 벡터만 들어왔다면 리스트로 변환
  if(!is.list(gene_sets)){
    gene_sets <- list(GeneSet1 = gene_sets)
  }
  
  # 이름이 없는 리스트 원소가 있다면 자동으로 이름 부여
  if(is.null(names(gene_sets)) || any(names(gene_sets) == "")){
    for(i in seq_along(gene_sets)){
      if(is.null(names(gene_sets)[i]) || names(gene_sets)[i] == ""){
        names(gene_sets)[i] <- paste0("GeneSet", i)
      }
    }
  }
  
  #-------------------------------'
  # (B) Seurat 객체에 grouping 적용
  #-------------------------------'
  Idents(sobj) <- group
  
  
  
  #-------------------------------'
  # (C) 평균 발현량(또는 합계 등) 계산
  #-------------------------------'
  if(value == "average"){
    cluster_avg <- AverageExpression(sobj, assays = assay, slot = "data", group.by = group)[[assay]]
  } else {
    cluster_avg <- AggregateExpression(sobj, assays = assay, slot = "data", group.by = group)[[assay]]
  }
  
  # If it's a vector (one group), convert to matrix
  if(is.null(dim(cluster_avg)) || length(dim(cluster_avg)) < 2){
    cluster_avg <- as.matrix(cluster_avg)
    if(ncol(cluster_avg) == 1 && is.null(colnames(cluster_avg))){
      colnames(cluster_avg) <- unique(Idents(sobj))[1]
    }
  }
  
  #-------------------------------'
  # (D) Gene Set 별 발현량 계산
  #-------------------------------'
  cluster_names <- colnames(cluster_avg)
  
  # 결과를 담을 data.frame 생성
  gene_set_expression <- data.frame(Cluster = cluster_names, stringsAsFactors = FALSE)
  
  # gene_sets 각각에 대해 평균 발현량을 구함
  for(gset_name in names(gene_sets)){
    genes <- gene_sets[[gset_name]]
    genes_present <- genes[genes %in% rownames(cluster_avg)]
    
    if(length(genes_present) == 0){
      warning(paste("No genes from", gset_name, "found in the dataset."))
      gene_set_expression[[gset_name]] <- NA
      next
    }
    
    # Calculate means
    subset_data <- cluster_avg[genes_present, , drop = FALSE]
    if(nrow(subset_data) == 1){
      # If only one gene, the result is already the mean
      gene_set_expression[[gset_name]] <- as.numeric(subset_data)
    } else {
      gene_set_expression[[gset_name]] <- colMeans(subset_data)
    }
  }
  
  #-------------------------------'
  # (E) Z-score 정규화
  #-------------------------------'
  gene_set_expression_normalized <- gene_set_expression
  gene_set_expression_normalized[,-1] <- scale(gene_set_expression_normalized[,-1])
  
  # 각 Cluster에서 가장 높은 값을 가지는 gene set을 배정
  gene_set_expression_normalized$Assigned_CellType <- apply(
    gene_set_expression_normalized[,-1], 1, 
    function(x){
      if(all(is.na(x))) return(NA)
      names(x)[which.max(x)]
    }
  )
  
  #-------------------------------'
  # (F) 클러스터 순서 정렬
  #-------------------------------'
  numeric_test <- suppressWarnings(as.numeric(gene_set_expression_normalized$Cluster))
  
  if(!all(is.na(numeric_test))){
    if(sum(is.na(numeric_test)) == 0){
      sorted_levels <- sort(unique(numeric_test))
      gene_set_expression_normalized$Cluster <- factor(
        gene_set_expression_normalized$Cluster,
        levels = as.character(sorted_levels)
      )
    } else {
      sorted_levels <- sort(unique(gene_set_expression_normalized$Cluster))
      gene_set_expression_normalized$Cluster <- factor(
        gene_set_expression_normalized$Cluster,
        levels = sorted_levels
      )
    }
  } else {
    sorted_levels <- sort(unique(gene_set_expression_normalized$Cluster))
    gene_set_expression_normalized$Cluster <- factor(
      gene_set_expression_normalized$Cluster,
      levels = sorted_levels
    )
  }
  
  #-------------------------------'
  # (G) Heatmap용 long format 만들기
  #-------------------------------'
  melted_data <- melt(
    gene_set_expression_normalized,
    id.vars = c("Cluster","Assigned_CellType"),
    variable.name = "GeneSet",
    value.name = "Zscore"
  )
  
  #-------------------------------'
  # (H) Heatmap 그리기
  #-------------------------------'
  p <- ggplot(melted_data, aes(x = Cluster, y = GeneSet, fill = Zscore)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    theme_minimal() +
    labs(
      title = title,
      x = x_label,
      y = y_label
    ) +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=12),
          axis.text.y = element_text(hjust=1, size=12),
          axis.title.x = element_text(face="bold", size=14),
          axis.title.y = element_text(face="bold", size=14),
          plot.title = element_text(size=16, face="bold", hjust=0.5))
  
  print(p)
  
  #-------------------------------'
  # (I) 결과 반환
  #-------------------------------'
  return(gene_set_expression_normalized)
}

#' Create a Heatmap of Individual Gene Expression per Cluster
#'
#' This function creates a heatmap showing the normalized expression of individual genes
#' across clusters. The expression values are z-score normalized for better visualization.
#'
#' @param sobj A Seurat object
#' @param group Character string specifying the identity to use for clustering.
#'             Default is "seurat_clusters"
#' @param value Character string specifying how to aggregate expression values.
#'             Options: "average" (default) or "sum"
#' @param assay Character string specifying which assay to use.
#'             Default is "SCT"
#' @param genes Character vector containing gene names to plot
#' @param title Character string for the plot title.
#'             Default is "Normalized Gene Expression per Cluster"
#' @param x_label Character string for the x-axis label.
#'               Default is "Cluster"
#' @param y_label Character string for the y-axis label.
#'               Default is "Genes"
#'
#' @return A data frame containing the normalized expression values for each gene
#'         across clusters
#'
#' @examples
#' \dontrun{
#' # Create heatmap for a set of genes
#' genes <- c("CD3D", "CD3E", "CD4", "CD8A", "MS4A1", "CD19")
#' result <- myhm_genes2(sobj, genes = genes)
#' 
#' # Create heatmap with custom title and labels
#' result <- myhm_genes2(sobj, 
#'                      genes = genes,
#'                      title = "T and B Cell Markers",
#'                      x_label = "Cell Clusters",
#'                      y_label = "Marker Genes")
#' }
#' @export
myhm_genes2 <- function(
    sobj,
    group = "seurat_clusters",
    value = "average",
    assay = "SCT",
    genes = NULL,
    title = "Normalized Gene Expression per Cluster",
    x_label = "Cluster",
    y_label = "Genes"
){
  #-------------------------------'
  # (A) 유효성 체크
  #-------------------------------'
  if(is.null(genes)){
    stop("genes를 지정해 주세요. (예: genes=c('CD3D','CD3E','MS4A1'))")
  }
  
  #-------------------------------'
  # (B) Seurat 객체에 grouping 적용
  #-------------------------------'
  Seurat::Idents(sobj) <- group
  
  #-------------------------------'
  # (C) 평균 발현량(또는 합계 등) 계산
  #-------------------------------'
  if(value == "average"){
    cluster_avg <- Seurat::AverageExpression(sobj, assays = assay, slot = "data", group.by = group)[[assay]]
  } else {
    cluster_avg <- Seurat::AggregateExpression(sobj, assays = assay, slot = "data", group.by = group)[[assay]]
  }
  
  # IMPORTANT FIX: Handle case when cluster_avg is a vector (single group)
  if(is.null(dim(cluster_avg)) || length(dim(cluster_avg)) < 2){
    cluster_avg <- as.matrix(cluster_avg)
    if(ncol(cluster_avg) == 1 && is.null(colnames(cluster_avg))){
      colnames(cluster_avg) <- unique(Seurat::Idents(sobj))[1]
    }
  }
  
  #-------------------------------'
  # (D) 원하는 유전자만 필터
  #-------------------------------'
  genes_present <- genes[genes %in% rownames(cluster_avg)]
  if(length(genes_present) == 0){
    stop("지정하신 유전자 중 데이터셋에 존재하는 유전자가 없습니다.")
  }
  
  if(length(genes_present) < length(genes)){
    missing_genes <- genes[!genes %in% rownames(cluster_avg)]
    warning(paste("다음 유전자들이 데이터셋에 없습니다:", paste(missing_genes, collapse=", ")))
  }
  
  # Subset 
  gene_expression <- cluster_avg[genes_present, , drop=FALSE]
  
  # IMPORTANT FIX: Handle case when only one cluster exists
  if(ncol(gene_expression) == 1){
    # For single cluster, scale doesn't work properly
    # Create a data frame directly
    gene_expression_df <- as.data.frame(gene_expression)
    colnames(gene_expression_df) <- colnames(gene_expression)
    gene_expression_df$Cluster <- colnames(gene_expression)
    rownames(gene_expression_df) <- genes_present
    
    # Z-score는 단일 클러스터에서는 의미가 없으므로 0으로 설정
    gene_expression_scaled <- matrix(0, 
                                     nrow = length(genes_present), 
                                     ncol = 1,
                                     dimnames = list(genes_present, colnames(gene_expression)))
  } else {
    # Multiple clusters - proceed with normal scaling
    # Transpose for scaling (rows=clusters, cols=genes)
    gene_expression_t <- t(gene_expression)
    gene_expression_scaled <- scale(gene_expression_t)
    gene_expression_scaled <- t(gene_expression_scaled)  # Transpose back
  }
  
  # Convert to data.frame
  gene_expression_df <- as.data.frame(t(gene_expression_scaled))
  gene_expression_df$Cluster <- rownames(gene_expression_df)
  
  #-------------------------------'
  # (E) 클러스터 순서 정렬
  #-------------------------------'
  numeric_test <- suppressWarnings(as.numeric(gene_expression_df$Cluster))
  
  if(!all(is.na(numeric_test))){
    if(sum(is.na(numeric_test)) == 0){
      sorted_levels <- sort(unique(numeric_test))
      gene_expression_df$Cluster <- factor(
        gene_expression_df$Cluster,
        levels = as.character(sorted_levels)
      )
    } else {
      sorted_levels <- sort(unique(gene_expression_df$Cluster))
      gene_expression_df$Cluster <- factor(
        gene_expression_df$Cluster,
        levels = sorted_levels
      )
    }
  } else {
    sorted_levels <- sort(unique(gene_expression_df$Cluster))
    gene_expression_df$Cluster <- factor(
      gene_expression_df$Cluster,
      levels = sorted_levels
    )
  }
  
  #-------------------------------'
  # (F) long format으로 melt
  #-------------------------------'
  melted_data <- reshape2::melt(
    gene_expression_df,
    id.vars = "Cluster",
    variable.name = "Gene",
    value.name = "Zscore"
  )
  
  # 유전자 순서를 입력 순서대로 유지
  melted_data$Gene <- factor(melted_data$Gene, levels = genes_present)
  
  #-------------------------------'
  # (G) Heatmap 그리기
  #-------------------------------'
  p <- ggplot2::ggplot(melted_data, ggplot2::aes(x = Cluster, y = Gene, fill = Zscore)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = title,
      x = x_label,
      y = y_label
    ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1),
                   plot.title = ggplot2::element_text(size=14, face="bold", hjust=0.5))
  
  print(p)
  
  #-------------------------------'
  # (H) 결과 반환
  #-------------------------------'
  return(gene_expression_df)
}
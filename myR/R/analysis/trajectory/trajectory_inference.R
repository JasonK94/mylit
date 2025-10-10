#' Trajectory Inference and Analysis
#'
#' This module provides functions for trajectory inference using Slingshot
#' and gene dynamics analysis along pseudotime.
#'
#' @name trajectory_inference
NULL

#' Run Slingshot from Seurat Object
#'
#' Runs Slingshot trajectory inference from a Seurat object and returns a
#' SingleCellExperiment object with computed trajectories.
#'
#' @param seurat_obj Seurat object
#' @param cluster_col Cluster column in metadata (default: "seurat_clusters")
#' @param reduction Dimensionality reduction to use (default: "umap")
#' @param start_cluster Starting cluster for trajectory (optional)
#' @param end_cluster Ending cluster(s) for trajectory (optional)
#' @param assay Assay to use (default: "RNA")
#' @param slot Slot to use (default: "counts")
#'
#' @return SingleCellExperiment object with Slingshot results
#'
#' @examples
#' \dontrun{
#' sce <- run_slingshot_from_seurat(
#'   seurat_obj,
#'   cluster_col = "seurat_clusters",
#'   reduction = "umap",
#'   start_cluster = "0"
#' )
#' }
#'
#' @export
run_slingshot_from_seurat <- function(seurat_obj,
                                      cluster_col = "seurat_clusters",
                                      reduction = "umap",
                                      start_cluster = NULL,
                                      end_cluster = NULL,
                                      assay = "RNA",
                                      slot = "counts") {
  
  # Check required packages
  if (!requireNamespace("slingshot", quietly = TRUE)) {
    stop("slingshot package required. Install with: BiocManager::install('slingshot')")
  }
  
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("SingleCellExperiment package required. Install with: BiocManager::install('SingleCellExperiment')")
  }
  
  # Validate inputs
  if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
    stop("cluster_col '", cluster_col, "' not found in metadata")
  }
  
  if (!reduction %in% names(seurat_obj@reductions)) {
    stop("reduction '", reduction, "' not found in Seurat object")
  }
  
  # Get counts
  counts <- Seurat::GetAssayData(seurat_obj, assay = assay, slot = slot)
  
  # Get dimensionality reduction
  dim_red <- Seurat::Embeddings(seurat_obj, reduction = reduction)
  
  # Get cluster assignments
  clusters <- seurat_obj@meta.data[[cluster_col]]
  
  # Validate alignment
  if (!identical(colnames(counts), rownames(dim_red))) {
    warning("Cell names don't match. Attempting to align...")
    common_cells <- intersect(colnames(counts), rownames(dim_red))
    if (length(common_cells) == 0) {
      stop("No common cells found between counts and dimensionality reduction")
    }
    counts <- counts[, common_cells]
    dim_red <- dim_red[common_cells, ]
    clusters <- clusters[match(common_cells, colnames(seurat_obj))]
  }
  
  # Check for issues
  if (any(!is.finite(dim_red))) {
    stop("Non-finite values found in dimensionality reduction")
  }
  
  if (any(!is.finite(as.matrix(counts)))) {
    warning("Non-finite values found in counts matrix. These will be replaced with 0.")
    counts[!is.finite(as.matrix(counts))] <- 0
  }
  
  # Ensure counts are integers
  if (!all(counts == floor(counts))) {
    warning("Non-integer values in counts. Rounding to nearest integer.")
    counts <- round(counts)
  }
  
  # Create SingleCellExperiment
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts),
    reducedDims = list(reducedDim = dim_red),
    colData = data.frame(cluster = clusters)
  )
  
  # Run Slingshot
  message("Running Slingshot trajectory inference...")
  sce <- slingshot::slingshot(
    sce,
    clusterLabels = "cluster",
    reducedDim = "reducedDim",
    start.clus = start_cluster,
    end.clus = end_cluster
  )
  
  n_lineages <- length(slingshot::slingLineages(sce))
  message("Found ", n_lineages, " lineage(s)")
  
  return(sce)
}

#' Analyze Gene Dynamics Along Pseudotime
#'
#' Fits a Generalized Additive Model (GAM) to model gene expression changes
#' along pseudotime, accounting for conditions.
#'
#' @param sce SingleCellExperiment object from run_slingshot_from_seurat
#' @param gene Gene name to analyze
#' @param condition Optional condition variable for comparison
#' @param lineage Lineage number to analyze (default: 1)
#' @param plot Create plots (default: TRUE)
#' @param output_dir Directory for saving plots (optional)
#'
#' @return List containing:
#'   \item{model}{GAM model object}
#'   \item{summary}{Model summary statistics}
#'   \item{plot}{ggplot object (if plot=TRUE)}
#'
#' @export
analyze_gene_dynamics <- function(sce,
                                  gene,
                                  condition = NULL,
                                  lineage = 1,
                                  plot = TRUE,
                                  output_dir = NULL) {
  
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("mgcv package required. Install with: install.packages('mgcv')")
  }
  
  if (!requireNamespace("slingshot", quietly = TRUE)) {
    stop("slingshot package required")
  }
  
  # Get pseudotime
  pseudotime <- slingshot::slingPseudotime(sce)[, lineage]
  
  # Get expression
  if (gene %in% rownames(sce)) {
    expr <- SingleCellExperiment::counts(sce)[gene, ]
  } else {
    stop("Gene '", gene, "' not found in SCE object")
  }
  
  # Remove cells with NA pseudotime
  valid_cells <- !is.na(pseudotime)
  pseudotime <- pseudotime[valid_cells]
  expr <- expr[valid_cells]
  
  # Prepare model data
  model_data <- data.frame(
    pseudotime = pseudotime,
    expression = expr
  )
  
  # Add condition if provided
  if (!is.null(condition)) {
    if (length(condition) != nrow(sce)) {
      stop("condition length must match number of cells in SCE")
    }
    model_data$condition <- condition[valid_cells]
    
    # Fit GAM with condition interaction
    gam_model <- mgcv::gam(
      expression ~ s(pseudotime, by = condition) + condition,
      data = model_data,
      family = mgcv::nb()
    )
  } else {
    # Fit simple GAM
    gam_model <- mgcv::gam(
      expression ~ s(pseudotime),
      data = model_data,
      family = mgcv::nb()
    )
  }
  
  # Extract summary
  model_summary <- summary(gam_model)
  
  results <- list(
    model = gam_model,
    summary = model_summary,
    gene = gene,
    lineage = lineage
  )
  
  # Create plot if requested
  if (plot) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      warning("ggplot2 not available. Skipping plot.")
    } else {
      # Predictions
      pred_data <- data.frame(pseudotime = seq(min(pseudotime), max(pseudotime), length.out = 100))
      
      if (!is.null(condition)) {
        # Expand for each condition
        pred_list <- list()
        for (cond in unique(model_data$condition)) {
          pred_temp <- pred_data
          pred_temp$condition <- cond
          pred_temp$fit <- predict(gam_model, newdata = pred_temp, type = "response")
          pred_list[[cond]] <- pred_temp
        }
        pred_data <- do.call(rbind, pred_list)
        
        # Create plot with conditions
        p <- ggplot2::ggplot() +
          ggplot2::geom_point(
            data = model_data,
            ggplot2::aes(x = pseudotime, y = expression, color = condition),
            alpha = 0.3, size = 1
          ) +
          ggplot2::geom_line(
            data = pred_data,
            ggplot2::aes(x = pseudotime, y = fit, color = condition),
            size = 1.5
          )
      } else {
        pred_data$fit <- predict(gam_model, newdata = pred_data, type = "response")
        
        # Create simple plot
        p <- ggplot2::ggplot() +
          ggplot2::geom_point(
            data = model_data,
            ggplot2::aes(x = pseudotime, y = expression),
            alpha = 0.3, size = 1, color = "gray50"
          ) +
          ggplot2::geom_line(
            data = pred_data,
            ggplot2::aes(x = pseudotime, y = fit),
            color = "blue", size = 1.5
          )
      }
      
      p <- p +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = paste(gene, "- Lineage", lineage),
          x = "Pseudotime",
          y = "Expression",
          subtitle = paste0("Deviance explained: ", 
                           round(model_summary$dev.expl * 100, 2), "%")
        )
      
      results$plot <- p
      
      # Save if output_dir provided
      if (!is.null(output_dir)) {
        if (!dir.exists(output_dir)) {
          dir.create(output_dir, recursive = TRUE)
        }
        
        filename <- file.path(output_dir, paste0(gene, "_lineage", lineage, "_dynamics.png"))
        ggplot2::ggsave(filename, p, width = 7, height = 5, dpi = 300)
        message("Plot saved: ", filename)
      }
      
      print(p)
    }
  }
  
  return(results)
}

#' Process Multiple Genes for Dynamics Analysis
#'
#' Processes a list of genes using analyze_gene_dynamics, with optional
#' parallel processing.
#'
#' @param sce SingleCellExperiment object
#' @param genes Character vector of gene names
#' @param condition Optional condition variable
#' @param lineage Lineage number (default: 1)
#' @param plot Create plots (default: FALSE)
#' @param output_dir Directory for saving results
#' @param n_cores Number of cores for parallel processing (default: 1)
#'
#' @return List of results for each gene
#'
#' @export
process_gene_list_dynamics <- function(sce,
                                       genes,
                                       condition = NULL,
                                       lineage = 1,
                                       plot = FALSE,
                                       output_dir = NULL,
                                       n_cores = 1) {
  
  # Validate genes
  available_genes <- intersect(genes, rownames(sce))
  if (length(available_genes) == 0) {
    stop("None of the specified genes found in SCE object")
  }
  
  if (length(available_genes) < length(genes)) {
    missing <- setdiff(genes, available_genes)
    warning("Genes not found: ", paste(missing, collapse = ", "))
  }
  
  message("Processing ", length(available_genes), " genes...")
  
  # Function to process a single gene
  process_single_gene <- function(gene) {
    tryCatch({
      analyze_gene_dynamics(
        sce = sce,
        gene = gene,
        condition = condition,
        lineage = lineage,
        plot = plot,
        output_dir = output_dir
      )
    }, error = function(e) {
      warning("Failed for gene ", gene, ": ", e$message)
      return(NULL)
    })
  }
  
  # Process genes
  if (n_cores > 1 && requireNamespace("parallel", quietly = TRUE)) {
    message("Using ", n_cores, " cores for parallel processing")
    results <- parallel::mclapply(
      available_genes,
      process_single_gene,
      mc.cores = n_cores
    )
  } else {
    results <- lapply(available_genes, process_single_gene)
  }
  
  names(results) <- available_genes
  
  # Remove failed genes
  results <- results[!sapply(results, is.null)]
  
  message("Successfully processed ", length(results), " genes")
  
  return(results)
}

#' Analyze Gene Dynamics using tradeSeq
#'
#' Fits GAMs using tradeSeq for differential expression patterns along trajectories.
#'
#' @param sce SingleCellExperiment object from Slingshot
#' @param gene Gene name to analyze
#' @param n_knots Number of knots for spline (default: 6)
#' @param conditions Optional condition variable
#'
#' @return List with test results
#'
#' @export
analyze_gene_dynamics_tradeSeq <- function(sce,
                                           gene,
                                           n_knots = 6,
                                           conditions = NULL) {
  
  if (!requireNamespace("tradeSeq", quietly = TRUE)) {
    stop("tradeSeq package required. Install with: BiocManager::install('tradeSeq')")
  }
  
  # Get pseudotime and weights
  pseudotime <- slingshot::slingPseudotime(sce)
  weights <- slingshot::slingCurveWeights(sce)
  
  # Get expression
  if (gene %in% rownames(sce)) {
    expr <- SingleCellExperiment::counts(sce)[gene, , drop = FALSE]
  } else {
    stop("Gene '", gene, "' not found in SCE object")
  }
  
  # Fit GAM
  message("Fitting tradeSeq GAM for ", gene, "...")
  sce_fit <- tradeSeq::fitGAM(
    counts = expr,
    pseudotime = pseudotime,
    cellWeights = weights,
    nknots = n_knots,
    verbose = FALSE
  )
  
  # Run tests
  results <- list(
    fit = sce_fit,
    gene = gene
  )
  
  # Association test (is gene changing along pseudotime?)
  results$association_test <- tradeSeq::associationTest(sce_fit)
  
  # Start vs end test
  results$start_vs_end_test <- tradeSeq::startVsEndTest(sce_fit)
  
  # Pattern test (if multiple lineages)
  if (ncol(pseudotime) > 1) {
    results$pattern_test <- tradeSeq::patternTest(sce_fit)
  }
  
  # Condition test (if conditions provided)
  if (!is.null(conditions)) {
    results$condition_test <- tradeSeq::conditionTest(sce_fit, conditions = conditions)
  }
  
  return(results)
}


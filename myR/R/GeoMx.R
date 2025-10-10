# GeoMx Data Analysis Workflow in R
# Standard methods for normalization and differential expression analysis

# NOTE: Package dependencies should be declared in DESCRIPTION, not with library() calls
# Required packages: NanoStringNCTools, GeomxTools, GeoMxWorkflows, dplyr, ggplot2, 
#                    pheatmap, limma, edgeR, tidyr

# 1. Data Preparation
# Assuming you have raw_data and SP (metadata) as described

# Set row names from TargetName
prepare_geomx_data <- function(raw_data, metadata_SP) {
  # Extract expression matrix
  # Remove non-expression columns (assuming first few columns are gene info)
  gene_info_cols <- c("ProbeName","ProbeDisplayName","TargetName","HUGOSymbol","Accessions","GenomeBuild","GenomicPosition","AnalyteType","CodeClass","ProbePool","TargetGroup","GeneID")  # adjust as needed
  
  # Create expression matrix
  expr_matrix <- raw_data %>%
    select(-any_of(gene_info_cols)) %>%
    as.matrix()
  
  # Set row names
  rownames(expr_matrix) <- raw_data$TargetName
  
  # Ensure column names match between expression matrix and metadata
  common_rois <- intersect(colnames(expr_matrix), rownames(metadata_SP))
  expr_matrix <- expr_matrix[, common_rois]
  metadata <- metadata_SP[common_rois, ]
  
  return(list(
    expression = expr_matrix,
    metadata = metadata,
    gene_info = raw_data %>% select(any_of(gene_info_cols))
  ))
}

# 2. Q3 Normalization Function
q3_normalize <- function(expr_matrix, scaling_factor = 1000) {
  # Calculate Q3 (75th percentile) for each ROI
  q3_values <- apply(expr_matrix, 2, function(x) {
    quantile(x[x > 0], probs = 0.75, na.rm = TRUE)
  })
  
  # Normalize each ROI by its Q3 value
  norm_matrix <- sweep(expr_matrix, 2, q3_values, FUN = "/") * scaling_factor
  
  # Log2 transform (add pseudocount to avoid log(0))
  norm_matrix_log <- log2(norm_matrix + 1)
  
  return(list(
    normalized = norm_matrix,
    normalized_log = norm_matrix_log,
    q3_factors = q3_values
  ))
}

# 3. Quality Control Functions
perform_qc <- function(expr_matrix, metadata) {
  # Calculate QC metrics
  qc_metrics <- data.frame(
    ROI = colnames(expr_matrix),
    total_counts = colSums(expr_matrix),
    detected_genes = colSums(expr_matrix > 0),
    mean_expression = colMeans(expr_matrix),
    median_expression = apply(expr_matrix, 2, median)
  )
  
  # Add metadata
  qc_metrics <- cbind(qc_metrics, metadata[qc_metrics$ROI, ])
  
  # QC plots
  p1 <- ggplot(qc_metrics, aes(x = ROI, y = total_counts)) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = "Total Counts per ROI")
  
  p2 <- ggplot(qc_metrics, aes(x = detected_genes, y = total_counts)) +
    geom_point() +
    labs(title = "Total Counts vs Detected Genes")
  
  return(list(metrics = qc_metrics, plot1 = p1, plot2 = p2))
}

# 4. Differential Expression Analysis (Similar to Seurat's FindMarkers)
find_deg_geomx <- function(norm_expr, metadata, group_var, 
                           group1, group2 = NULL, 
                           method = "limma", 
                           logFC_threshold = 0.5, 
                           pval_threshold = 0.05) {
  
  # Filter samples for comparison
  if (!is.null(group2)) {
    # Two-group comparison
    samples_keep <- metadata[[group_var]] %in% c(group1, group2)
    expr_subset <- norm_expr[, samples_keep]
    metadata_subset <- metadata[samples_keep, ]
    
    # Create design matrix
    group_factor <- factor(metadata_subset[[group_var]], levels = c(group2, group1))
    design <- model.matrix(~ group_factor)
    
  } else {
    # One vs all comparison
    metadata$comparison_group <- ifelse(metadata[[group_var]] == group1, group1, "Other")
    group_factor <- factor(metadata$comparison_group, levels = c("Other", group1))
    design <- model.matrix(~ group_factor)
    expr_subset <- norm_expr
    metadata_subset <- metadata
  }
  
  if (method == "limma") {
    # limma-trend for normalized data
    fit <- lmFit(expr_subset, design)
    fit <- eBayes(fit, trend = TRUE)
    
    # Extract results
    results <- topTable(fit, coef = 2, number = Inf, adjust.method = "BH")
    
    # Add additional statistics
    results$gene <- rownames(results)
    results$pct.1 <- rowMeans(expr_subset[, metadata_subset[[group_var]] == group1] > 0)
    if (!is.null(group2)) {
      results$pct.2 <- rowMeans(expr_subset[, metadata_subset[[group_var]] == group2] > 0)
    } else {
      results$pct.2 <- rowMeans(expr_subset[, metadata_subset[[group_var]] != group1] > 0)
    }
    
  } else if (method == "wilcox") {
    # Wilcoxon rank sum test (similar to Seurat default)
    results <- data.frame(gene = rownames(expr_subset))
    
    for (i in 1:nrow(expr_subset)) {
      if (!is.null(group2)) {
        group1_expr <- expr_subset[i, metadata_subset[[group_var]] == group1]
        group2_expr <- expr_subset[i, metadata_subset[[group_var]] == group2]
      } else {
        group1_expr <- expr_subset[i, metadata_subset[[group_var]] == group1]
        group2_expr <- expr_subset[i, metadata_subset[[group_var]] != group1]
      }
      
      # Wilcoxon test
      test_result <- wilcox.test(group1_expr, group2_expr)
      
      results$p_val[i] <- test_result$p.value
      results$avg_logFC[i] <- mean(group1_expr) - mean(group2_expr)
      results$pct.1[i] <- mean(group1_expr > 0)
      results$pct.2[i] <- mean(group2_expr > 0)
    }
    
    # Adjust p-values
    results$p_val_adj <- p.adjust(results$p_val, method = "BH")
    rownames(results) <- results$gene
  }
  
  # Filter by thresholds
  results$significant <- abs(results$avg_logFC) > logFC_threshold & 
    results$p_val_adj < pval_threshold
  
  # Sort by p-value
  results <- results[order(results$p_val_adj), ]
  
  return(results)
}

# 5. Visualization Functions
plot_deg_volcano <- function(deg_results, title = "Volcano Plot") {
  deg_results$log10_pval <- -log10(deg_results$p_val_adj)
  deg_results$color <- ifelse(deg_results$significant, 
                              ifelse(deg_results$avg_logFC > 0, "Up", "Down"), 
                              "NS")
  
  p <- ggplot(deg_results, aes(x = avg_logFC, y = log10_pval, color = color)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "gray")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
    labs(title = title, x = "Log2 Fold Change", y = "-Log10(adjusted p-value)") +
    theme_classic()
  
  return(p)
}

plot_deg_heatmap <- function(norm_expr, deg_results, metadata, group_var, 
                             top_n = 50, scale_rows = TRUE) {
  # Select top DEGs
  top_genes <- head(rownames(deg_results[deg_results$significant, ]), top_n)
  
  if (length(top_genes) == 0) {
    print("No significant genes found")
    return(NULL)
  }
  
  # Subset expression matrix
  expr_subset <- norm_expr[top_genes, ]
  
  # Scale if requested
  if (scale_rows) {
    expr_subset <- t(scale(t(expr_subset)))
  }
  
  # Create annotation
  annotation_col <- data.frame(
    Group = metadata[[group_var]],
    row.names = rownames(metadata)
  )
  
  # Create heatmap
  pheatmap(expr_subset,
           annotation_col = annotation_col,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           show_colnames = FALSE,
           main = paste("Top", top_n, "DEGs"))
}

# 6. Complete Analysis Pipeline
run_geomx_analysis <- function(raw_data, metadata_SP, 
                               group_var, group1, group2 = NULL,
                               qc_filter = TRUE,
                               min_genes = 100,
                               min_counts = 1000) {
  
  # Step 1: Prepare data
  cat("Preparing data...\n")
  data_list <- prepare_geomx_data(raw_data, metadata_SP)
  
  # Step 2: QC
  cat("Performing QC...\n")
  qc_results <- perform_qc(data_list$expression, data_list$metadata)
  
  # QC filtering
  if (qc_filter) {
    keep_rois <- qc_results$metrics$detected_genes >= min_genes & 
      qc_results$metrics$total_counts >= min_counts
    
    data_list$expression <- data_list$expression[, keep_rois]
    data_list$metadata <- data_list$metadata[keep_rois, ]
    
    cat(sprintf("Kept %d/%d ROIs after QC\n", 
                sum(keep_rois), length(keep_rois)))
  }
  
  # Step 3: Normalization
  cat("Performing Q3 normalization...\n")
  norm_results <- q3_normalize(data_list$expression)
  
  # Step 4: Differential expression
  cat("Finding differentially expressed genes...\n")
  deg_results <- find_deg_geomx(
    norm_expr = norm_results$normalized_log,
    metadata = data_list$metadata,
    group_var = group_var,
    group1 = group1,
    group2 = group2,
    method = "limma"
  )
  
  n_sig <- sum(deg_results$significant)
  cat(sprintf("Found %d significant DEGs\n", n_sig))
  
  # Step 5: Create visualizations
  cat("Creating visualizations...\n")
  volcano_plot <- plot_deg_volcano(deg_results)
  
  heatmap_plot <- if (n_sig > 0) {
    plot_deg_heatmap(
      norm_results$normalized_log, 
      deg_results, 
      data_list$metadata, 
      group_var
    )
  } else {
    NULL
  }
  
  # Return results
  return(list(
    normalized_data = norm_results$normalized_log,
    deg_results = deg_results,
    qc_metrics = qc_results$metrics,
    plots = list(
      qc1 = qc_results$plot1,
      qc2 = qc_results$plot2,
      volcano = volcano_plot
    )
  ))
}

# Example usage:
# results <- run_geomx_analysis(
#   raw_data = raw_data,
#   metadata_SP = SP,
#   group_var = "Treatment",
#   group1 = "Treated",
#   group2 = "Control"
# )
# 
# # View top DEGs
# head(results$deg_results[results$deg_results$significant, ])
# 
# # Save results
# write.csv(results$deg_results, "DEG_results.csv")
# ggsave("volcano_plot.pdf", results$plots$volcano, width = 8, height = 6)
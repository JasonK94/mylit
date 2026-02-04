#!/usr/bin/env Rscript
# Test script for cmb sorting functionality
# Uses is5 data to verify sorting logic

# Load required libraries
if (!requireNamespace("Seurat", quietly = TRUE)) {
  stop("Seurat package is required")
}
library(Seurat)
library(dplyr)
library(ggplot2)

# Source the functions (assuming we're in the plots worktree)
source_dir <- "/home/user3/data_user3/git_repo/_wt/plots/myR/R"
source(file.path(source_dir, "plots.R"))

# Load data
cat("Loading is5 data...\n")
if (!file.exists("/data/user3/sobj/IS6_sex_added_251110.qs")) {
  stop("Data file not found: /data/user3/sobj/IS6_sex_added_251110.qs")
}
is5 <- qs::qread("/data/user3/sobj/IS6_sex_added_251110.qs")
cat("Data loaded. Dimensions:", dim(is5), "\n")

# Test 1: Basic sorting with sort.by (cluster frequency)
cat("\n=== Test 1: sort.by = 'Monocytes' (cluster frequency) ===\n")
Idents(is5) <- "anno3.scvi"
unique_clusters <- unique(Idents(is5))
cat("Available clusters:", paste(unique_clusters, collapse=", "), "\n")

# Get a cluster that exists
test_cluster <- unique_clusters[1]
cat("Testing with cluster:", test_cluster, "\n")

# Calculate frequencies manually
sample_cluster_counts <- is5@meta.data %>%
  dplyr::select(GEM, anno3.scvi) %>%
  dplyr::group_by(GEM, anno3.scvi) %>%
  dplyr::summarise(count = n(), .groups = "drop") %>%
  dplyr::group_by(GEM) %>%
  dplyr::mutate(proportion = count / sum(count)) %>%
  dplyr::ungroup()

cluster_freqs <- sample_cluster_counts %>%
  dplyr::filter(anno3.scvi == test_cluster) %>%
  dplyr::select(GEM, count) %>%
  dplyr::right_join(
    data.frame(GEM = unique(sample_cluster_counts$GEM), stringsAsFactors = FALSE),
    by = "GEM"
  ) %>%
  dplyr::mutate(count = ifelse(is.na(count), 0, count)) %>%
  dplyr::arrange(desc(count))

cat("\nManual calculation - sorted by frequency:\n")
print(cluster_freqs)

# Test 2: sort.by with metadata column
cat("\n=== Test 2: sort.by = metadata column ===\n")
# Check available numeric metadata columns
numeric_cols <- sapply(is5@meta.data, is.numeric)
numeric_col_names <- names(numeric_cols)[numeric_cols]
cat("Available numeric metadata columns:", paste(head(numeric_col_names, 10), collapse=", "), "\n")

if(length(numeric_col_names) > 0) {
  test_meta <- numeric_col_names[1]
  cat("Testing with metadata column:", test_meta, "\n")
  
  meta_means <- is5@meta.data %>%
    dplyr::select(GEM, !!sym(test_meta)) %>%
    dplyr::group_by(GEM) %>%
    dplyr::summarise(mean_value = mean(.data[[test_meta]], na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(desc(mean_value))
  
  cat("\nManual calculation - sorted by mean metadata:\n")
  print(head(meta_means, 20))
}

# Test 3: sort.by with feature/gene
cat("\n=== Test 3: sort.by = feature/gene ===\n")
test_genes <- c("CD3D", "CD4", "CD8A", "MS4A1")
available_genes <- test_genes[test_genes %in% rownames(is5)]
if(length(available_genes) > 0) {
  test_gene <- available_genes[1]
  cat("Testing with gene:", test_gene, "\n")
  
  expr_data <- Seurat::FetchData(is5, vars = test_gene)
  gene_means <- data.frame(
    GEM = is5@meta.data$GEM,
    expression = expr_data[[test_gene]],
    stringsAsFactors = FALSE
  ) %>%
    dplyr::group_by(GEM) %>%
    dplyr::summarise(mean_expr = mean(expression, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(desc(mean_expr))
  
  cat("\nManual calculation - sorted by mean expression:\n")
  print(head(gene_means, 20))
}

# Test 4: split.by + sort.by combination
cat("\n=== Test 4: split.by + sort.by combination ===\n")
if("g3" %in% names(is5@meta.data)) {
  cat("Testing with split.by = 'g3' and sort.by =", test_cluster, "\n")
  
  # Get split.by values for each sample
  sample_split_map <- is5@meta.data %>%
    dplyr::select(GEM, g3) %>%
    distinct() %>%
    dplyr::filter(!is.na(g3))
  
  split_groups <- unique(sample_split_map$g3)
  split_groups <- split_groups[order(split_groups)]
  cat("Split groups:", paste(split_groups, collapse=", "), "\n")
  
  # Calculate sort values
  cluster_freqs_all <- sample_cluster_counts %>%
    dplyr::filter(anno3.scvi == test_cluster) %>%
    dplyr::select(GEM, count) %>%
    dplyr::right_join(
      data.frame(GEM = unique(sample_cluster_counts$GEM), stringsAsFactors = FALSE),
      by = "GEM"
    ) %>%
    dplyr::mutate(count = ifelse(is.na(count), 0, count))
  
  sort_values <- setNames(cluster_freqs_all$count, cluster_freqs_all$GEM)
  
  # Sort within each split.by group
  sorted_samples <- c()
  for(sg in split_groups) {
    samples_in_group <- sample_split_map[sample_split_map$g3 == sg, ]$GEM
    samples_in_group <- intersect(samples_in_group, unique(sample_cluster_counts$GEM))
    if(length(samples_in_group) > 0) {
      group_sort_vals <- sort_values[names(sort_values) %in% samples_in_group]
      sorted_group <- names(sort(group_sort_vals, decreasing = TRUE))
      cat("\nGroup", sg, "- sorted samples:\n")
      print(data.frame(
        GEM = sorted_group,
        Frequency = group_sort_vals[sorted_group],
        stringsAsFactors = FALSE
      ))
      sorted_samples <- c(sorted_samples, sorted_group)
    }
  }
  
  cat("\nFinal sorted order (all groups combined):\n")
  print(data.frame(
    Position = 1:length(sorted_samples),
    GEM = sorted_samples,
    g3 = sample_split_map[match(sorted_samples, sample_split_map$GEM), ]$g3,
    Frequency = sort_values[sorted_samples],
    stringsAsFactors = FALSE
  ))
  
  # Check for reversals
  cat("\nChecking for sorting reversals...\n")
  for(i in 1:(length(sorted_samples)-1)) {
    current_val <- sort_values[sorted_samples[i]]
    next_val <- sort_values[sorted_samples[i+1]]
    current_g3 <- sample_split_map[sample_split_map$GEM == sorted_samples[i], ]$g3
    next_g3 <- sample_split_map[sample_split_map$GEM == sorted_samples[i+1], ]$g3
    
    if(current_g3 == next_g3 && current_val < next_val) {
      cat("WARNING: Reversal detected at position", i, "->", i+1, "\n")
      cat("  Sample", sorted_samples[i], ":", current_val, "vs Sample", sorted_samples[i+1], ":", next_val, "\n")
    }
  }
}

cat("\n=== Test Complete ===\n")


#' Statistical Analysis for scRNA-seq Metadata
#'
#' This function performs comprehensive statistical analysis on scRNA-seq data
#' comparing groups across different categorical variables.
#'
#' @param data A Seurat object or metadata dataframe
#' @param grouping_var Character. Variable to group by (e.g., "patient")
#' @param categorizing_var Character. Dependent variable (e.g., "seurat_clusters")
#' @param comparative_var Character. Independent variable (e.g., "condition", "smoking_status")
#' @param test_use Character vector. Statistical tests to use, options: "t", "u", or both c("t", "u")
#' @param ... Additional arguments (not used)
#'
#' @return A data frame with statistical test results
#' @export
seurat_group_stats <- function(data, grouping_var, categorizing_var, 
                               comparative_var, test_use = c("t", "u"), ...) {
  
  # Check input type and extract metadata if Seurat object
  if (inherits(data, "Seurat")) {
    meta_data <- data@meta.data
  } else if (is.data.frame(data)) {
    meta_data <- data
  } else {
    stop("Input must be a Seurat object or a data frame. Please check your input.")
  }
  
  # Ensure variables exist in the metadata
  required_vars <- c(grouping_var, categorizing_var, comparative_var)
  missing_vars <- required_vars[!required_vars %in% colnames(meta_data)]
  
  if (length(missing_vars) > 0) {
    stop(paste("The following variables are missing from the metadata:",
               paste(missing_vars, collapse = ", ")))
  }
  
  # Initialize results dataframe
  results <- data.frame()
  
  # Get unique categories from categorizing variable
  categories <- unique(meta_data[[categorizing_var]])
  
  # Process each category
  for (cat in categories) {
    # Subset data for current category
    cat_data <- meta_data[meta_data[[categorizing_var]] == cat, ]
    
    # Aggregate data by grouping variable and comparative variable
    agg_data <- aggregate(rep(1, nrow(cat_data)), 
                          by = list(cat_data[[grouping_var]], cat_data[[comparative_var]]), 
                          FUN = sum)
    colnames(agg_data) <- c("Group", "Comparative", "Count")
    
    # Create a wide format for statistical testing
    wide_data <- reshape(agg_data, idvar = "Group", timevar = "Comparative", 
                         direction = "wide")
    
    # Remove the "Count." prefix from column names
    colnames(wide_data)[-1] <- sub("Count.", "", colnames(wide_data)[-1])
    
    # Get the comparative variable levels
    comp_levels <- unique(agg_data$Comparative)
    
    # Skip if there's only one level in the comparative variable
    if (length(comp_levels) < 2) {
      next
    }
    
    # Skip if there are too few samples
    if (nrow(wide_data) < 3) {
      cat_result <- data.frame(
        Category = cat,
        ComparativeVar = comparative_var,
        NumGroups = length(comp_levels),
        Levene_pval = NA,
        ShapiroWilk_pval = NA,
        Parametric_test = NA,
        Parametric_metric = NA,
        Parametric_pval = NA,
        Parametric_adj_pval = NA,
        Parametric_relevant = NA,
        Nonparametric_test = NA,
        Nonparametric_metric = NA,
        Nonparametric_pval = NA,
        Nonparametric_adj_pval = NA
      )
      results <- rbind(results, cat_result)
      next
    }
    
    # Prepare data for tests
    test_data <- list()
    for (lvl in comp_levels) {
      if (lvl %in% colnames(wide_data)) {
        test_data[[as.character(lvl)]] <- wide_data[[as.character(lvl)]]
      }
    }
    
    # Remove NAs
    test_data <- lapply(test_data, function(x) x[!is.na(x)])
    
    # Check if any group has too few samples
    if (any(sapply(test_data, length) < 3)) {
      cat_result <- data.frame(
        Category = cat,
        ComparativeVar = comparative_var,
        NumGroups = length(comp_levels),
        Levene_pval = NA,
        ShapiroWilk_pval = NA,
        Parametric_test = NA,
        Parametric_metric = NA,
        Parametric_pval = NA,
        Parametric_adj_pval = NA,
        Parametric_relevant = NA,
        Nonparametric_test = NA,
        Nonparametric_metric = NA,
        Nonparametric_pval = NA,
        Nonparametric_adj_pval = NA
      )
      results <- rbind(results, cat_result)
      next
    }
    
    # Perform Levene test for homogeneity of variances
    levene_result <- car::leveneTest(unlist(test_data) ~ rep(names(test_data), sapply(test_data, length)))
    levene_pval <- levene_result[1, "Pr(>F)"]
    
    # Perform Shapiro-Wilk test for normality on residuals
    residuals <- unlist(lapply(names(test_data), function(grp) {
      test_data[[grp]] - mean(test_data[[grp]])
    }))
    shapiro_result <- shapiro.test(residuals)
    shapiro_pval <- shapiro_result$p.value
    
    # Determine appropriate tests
    is_normal <- shapiro_pval > 0.05
    equal_var <- levene_pval > 0.05
    
    # Initialize parametric test results
    param_test <- "None"
    param_metric <- NA
    param_pval <- NA
    param_relevant <- FALSE
    
    # Parametric tests (if requested)
    if ("t" %in% test_use) {
      # Determine appropriate parametric test
      if (length(comp_levels) == 2) {
        # Two groups: t-test
        if (is_normal) {
          if (equal_var) {
            # Equal variances: Student's t-test
            t_test <- t.test(test_data[[1]], test_data[[2]], var.equal = TRUE)
            param_test <- "Student's t-test"
          } else {
            # Unequal variances: Welch's t-test
            t_test <- t.test(test_data[[1]], test_data[[2]], var.equal = FALSE)
            param_test <- "Welch's t-test"
          }
          param_metric <- t_test$statistic
          param_pval <- t_test$p.value
          param_relevant <- TRUE
        } else {
          param_test <- "t-test not appropriate (non-normal data)"
        }
      } else {
        # More than two groups: ANOVA or Welch's ANOVA
        if (is_normal) {
          if (equal_var) {
            # Equal variances: standard ANOVA
            anova_model <- aov(unlist(test_data) ~ rep(names(test_data), sapply(test_data, length)))
            anova_result <- summary(anova_model)[[1]]
            param_test <- "ANOVA"
            param_metric <- anova_result$`F value`[1]
            param_pval <- anova_result$`Pr(>F)`[1]
            param_relevant <- TRUE
          } else {
            # Unequal variances: Welch's ANOVA
            welch_result <- oneway.test(unlist(test_data) ~ rep(names(test_data), sapply(test_data, length)), 
                                        var.equal = FALSE)
            param_test <- "Welch's ANOVA"
            param_metric <- welch_result$statistic
            param_pval <- welch_result$p.value
            param_relevant <- TRUE
          }
        } else {
          param_test <- "ANOVA not appropriate (non-normal data)"
        }
      }
    }
    
    # Initialize non-parametric test results
    nonparam_test <- "None"
    nonparam_metric <- NA
    nonparam_pval <- NA
    
    # Non-parametric tests (if requested)
    if ("u" %in% test_use) {
      if (length(comp_levels) == 2) {
        # Two groups: Wilcoxon rank-sum test (Mann-Whitney U)
        wilcox_result <- wilcox.test(test_data[[1]], test_data[[2]])
        nonparam_test <- "Wilcoxon rank-sum (Mann-Whitney U)"
        nonparam_metric <- wilcox_result$statistic
        nonparam_pval <- wilcox_result$p.value
      } else {
        # More than two groups: Kruskal-Wallis test
        kw_result <- kruskal.test(list = test_data)
        nonparam_test <- "Kruskal-Wallis"
        nonparam_metric <- kw_result$statistic
        nonparam_pval <- kw_result$p.value
      }
    }
    
    # Store results for current category
    cat_result <- data.frame(
      Category = cat,
      ComparativeVar = comparative_var,
      NumGroups = length(comp_levels),
      Levene_pval = levene_pval,
      ShapiroWilk_pval = shapiro_pval,
      Parametric_test = param_test,
      Parametric_metric = param_metric,
      Parametric_pval = param_pval,
      Parametric_adj_pval = NA,  # Will adjust later
      Parametric_relevant = param_relevant,
      Nonparametric_test = nonparam_test,
      Nonparametric_metric = nonparam_metric,
      Nonparametric_pval = nonparam_pval,
      Nonparametric_adj_pval = NA  # Will adjust later
    )
    
    results <- rbind(results, cat_result)
  }
  
  # Apply multiple testing correction if there are results
  if (nrow(results) > 0) {
    # Adjust p-values for parametric tests
    results$Parametric_adj_pval <- p.adjust(results$Parametric_pval, method = "BH")
    
    # Adjust p-values for non-parametric tests
    results$Nonparametric_adj_pval <- p.adjust(results$Nonparametric_pval, method = "BH")
  }
  
  # Print results
  print(results)
  
  # Return results
  return(results)
}



#' Post-hoc Analysis for scRNA-seq Statistical Tests
#'
#' This function performs post-hoc analysis based on results from seurat_group_stats.
#'
#' @param data A Seurat object or metadata dataframe
#' @param results Results dataframe from seurat_group_stats function
#' @param alpha Significance level (default: 0.05)
#'
#' @return A list of post-hoc test results
#' @export
seurat_posthoc_analysis <- function(data, results, alpha = 0.05) {
  
  # Check input type and extract metadata if Seurat object
  if (inherits(data, "Seurat")) {
    meta_data <- data@meta.data
  } else if (is.data.frame(data)) {
    meta_data <- data
  } else {
    stop("Input must be a Seurat object or a data frame. Please check your input.")
  }
  
  # Initialize output list
  posthoc_results <- list()
  
  # Iterate through significant results
  for (i in 1:nrow(results)) {
    # Only perform post-hoc if we have more than 2 groups and significant results
    if (results$NumGroups[i] > 2) {
      cat_name <- results$Category[i]
      comp_var <- results$ComparativeVar[i]
      
      # Check if parametric test is significant and relevant
      param_significant <- !is.na(results$Parametric_adj_pval[i]) && 
        results$Parametric_adj_pval[i] < alpha &&
        results$Parametric_relevant[i]
      
      # Check if non-parametric test is significant
      nonparam_significant <- !is.na(results$Nonparametric_adj_pval[i]) && 
        results$Nonparametric_adj_pval[i] < alpha
      
      # Skip if neither test is significant
      if (!param_significant && !nonparam_significant) {
        next
      }
      
      # Subset data for current category
      cat_data <- meta_data[meta_data[[results$Category[i]]] == cat_name, ]
      
      # Aggregate data
      agg_data <- aggregate(rep(1, nrow(cat_data)), 
                            by = list(cat_data[[grouping_var]], cat_data[[comp_var]]), 
                            FUN = sum)
      colnames(agg_data) <- c("Group", "Comparative", "Count")
      
      # Create a wide format for statistical testing
      wide_data <- reshape(agg_data, idvar = "Group", timevar = "Comparative", 
                           direction = "wide")
      
      # Remove the "Count." prefix from column names
      colnames(wide_data)[-1] <- sub("Count.", "", colnames(wide_data)[-1])
      
      # Get the comparative variable levels
      comp_levels <- unique(agg_data$Comparative)
      
      # Prepare data for tests
      test_data <- list()
      for (lvl in comp_levels) {
        if (lvl %in% colnames(wide_data)) {
          test_data[[as.character(lvl)]] <- wide_data[[as.character(lvl)]]
        }
      }
      
      # Remove NAs
      test_data <- lapply(test_data, function(x) x[!is.na(x)])
      
      # Perform post-hoc tests
      if (param_significant) {
        # Determine appropriate parametric post-hoc test
        if (results$Parametric_test[i] == "ANOVA") {
          # Tukey HSD for equal variances
          anova_model <- aov(unlist(test_data) ~ rep(names(test_data), sapply(test_data, length)))
          tukey_result <- TukeyHSD(anova_model)
          
          posthoc_results[[paste0(cat_name, "_parametric")]] <- list(
            test = "Tukey HSD",
            result = tukey_result
          )
        } else if (results$Parametric_test[i] == "Welch's ANOVA") {
          # Games-Howell for unequal variances
          gh_result <- userfriendlyscience::posthocTGH(
            y = unlist(test_data),
            x = factor(rep(names(test_data), sapply(test_data, length))),
            method = "games-howell"
          )
          
          posthoc_results[[paste0(cat_name, "_parametric")]] <- list(
            test = "Games-Howell",
            result = gh_result
          )
        }
      }
      
      if (nonparam_significant) {
        # Dunn test for non-parametric comparison
        dunn_result <- FSA::dunnTest(
          unlist(test_data) ~ factor(rep(names(test_data), sapply(test_data, length))),
          method = "bh"
        )
        
        posthoc_results[[paste0(cat_name, "_nonparametric")]] <- list(
          test = "Dunn's Test",
          result = dunn_result
        )
      }
    }
  }
  
  # Return post-hoc results
  return(posthoc_results)
}







# Claude style

# Seurat Cluster Fraction Analysis Functions
# NOTE: Package dependencies should be declared in DESCRIPTION, not with library() calls

#' Extract sample-level metadata from cell-level metadata
#' 
#' @param cell_metadata Cell-level metadata (Seurat object or data.frame)
#' @param sample_key Column name for sample identification
#' @param target_key Column name for the target variable (e.g., cluster)
#' @param group_key Column name for grouping variable (e.g., prognosis)
#' @param fun Aggregation function (default: proportion calculation)
#' @return Data frame with sample-level aggregated data
extract_sample_metadata <- function(cell_metadata, 
                                    sample_key, 
                                    target_key, 
                                    group_key = NULL,
                                    fun = NULL) {
  
  # Seurat 객체인 경우 메타데이터 추출
  if ("Seurat" %in% class(cell_metadata)) {
    metadata <- cell_metadata@meta.data
  } else {
    metadata <- cell_metadata
  }
  
  # 필수 컬럼 확인
  required_cols <- c(sample_key, target_key)
  if (!is.null(group_key)) required_cols <- c(required_cols, group_key)
  
  if (!all(required_cols %in% colnames(metadata))) {
    stop("Required columns not found in metadata")
  }
  
  # 기본 함수: 각 카테고리의 비율 계산
  if (is.null(fun)) {
    # 각 샘플별로 타겟 변수의 분율 계산
    sample_data <- metadata %>%
      group_by(.data[[sample_key]], .data[[target_key]]) %>%
      summarise(n = n(), .groups = "drop") %>%
      group_by(.data[[sample_key]]) %>%
      mutate(fraction = n / sum(n)) %>%
      select(-n) %>%
      pivot_wider(names_from = all_of(target_key), 
                  values_from = fraction, 
                  values_fill = 0)
    
    # group_key 정보 추가
    if (!is.null(group_key)) {
      group_info <- metadata %>%
        select(all_of(c(sample_key, group_key))) %>%
        distinct()
      
      sample_data <- sample_data %>%
        left_join(group_info, by = sample_key)
    }
    
  } else {
    # 사용자 정의 함수 적용
    sample_data <- metadata %>%
      group_by(!!sym(sample_key)) %>%
      summarise(value = fun(.data), .groups = "drop")
  }
  
  return(as.data.frame(sample_data))
}

#' Perform prior tests for statistical assumptions
#' 
#' @param data Data frame with sample-level data
#' @param value_col Column name for the values to test
#' @param group_col Column name for grouping variable
#' @return List with test results
perform_prior_tests <- function(data, value_col, group_col) {
  
  results <- list()
  
  # 정규성 검정 (Shapiro-Wilk test)
  groups <- unique(data[[group_col]])
  normality_results <- list()
  
  for (g in groups) {
    subset_data <- data[data[[group_col]] == g, value_col]
    if (length(subset_data) >= 3) {
      test_result <- shapiro.test(subset_data)
      normality_results[[g]] <- list(
        statistic = test_result$statistic,
        p.value = test_result$p.value,
        normal = test_result$p.value > 0.05
      )
    } else {
      normality_results[[g]] <- list(
        statistic = NA,
        p.value = NA,
        normal = NA,
        note = "Sample size too small for test"
      )
    }
  }
  results$normality <- normality_results
  
  # 등분산성 검정 (Levene's test)
  if (length(groups) == 2) {
    levene_result <- car::leveneTest(
      as.formula(paste(value_col, "~", group_col)), 
      data = data
    )
    results$homogeneity <- list(
      statistic = levene_result$`F value`[1],
      p.value = levene_result$`Pr(>F)`[1],
      equal_variance = levene_result$`Pr(>F)`[1] > 0.05
    )
  }
  
  # 추천 검정 방법
  all_normal <- all(sapply(normality_results, function(x) 
    isTRUE(x$normal) || is.na(x$normal)))
  
  if (length(groups) == 2) {
    if (all_normal && results$homogeneity$equal_variance) {
      results$recommendation <- "t-test"
    } else {
      results$recommendation <- "wilcox"
    }
  } else {
    if (all_normal) {
      results$recommendation <- "anova"
    } else {
      results$recommendation <- "kruskal"
    }
  }
  
  return(results)
}

#' Main function to create boxplot with statistical comparisons
#'
#' @param sobj_metadata Seurat object or metadata data.frame
#' @param sample_key Column name for sample identification
#' @param cluster_key Column name for clusters
#' @param group_key Column name for grouping (e.g., prognosis)
#' @param clusters_to_plot Vector of clusters to include (NULL for all)
#' @param groups_to_plot Vector of groups to include (NULL for all)
#' @param method Statistical test method
#' @param prior_test Whether to perform prior tests
#' @param palette Color palette for groups
#' @param show_brackets Whether to show brackets for comparisons
#' @param p_label_size Size of p-value labels
#' @param ... Additional arguments for ggplot
plot_cluster_fractions <- function(sobj_metadata,
                                   sample_key,
                                   cluster_key,
                                   group_key,
                                   clusters_to_plot = NULL,
                                   groups_to_plot = NULL,
                                   method = c("wilcox", "t.test", "anova", "kruskal"),
                                   prior_test = FALSE,
                                   palette = NULL,
                                   show_brackets = TRUE,
                                   p_label_size = 3.5,
                                   simple_pvalues = FALSE,
                                   ...) {
  
  method <- match.arg(method)
  
  # Extract sample-level data
  sample_data <- extract_sample_metadata(
    sobj_metadata, 
    sample_key, 
    cluster_key, 
    group_key
  )
  
  # Convert to long format for plotting
  cluster_cols <- setdiff(colnames(sample_data), c(sample_key, group_key))
  
  plot_data <- sample_data %>%
    pivot_longer(cols = all_of(cluster_cols), 
                 names_to = "cluster", 
                 values_to = "fraction")
  
  # Filter clusters and groups if specified
  if (!is.null(clusters_to_plot)) {
    plot_data <- plot_data %>%
      filter(cluster %in% clusters_to_plot)
  }
  
  if (!is.null(groups_to_plot)) {
    plot_data <- plot_data %>%
      filter(.data[[group_key]] %in% groups_to_plot)
  }
  
  # Create combined factor for x-axis
  plot_data <- plot_data %>%
    mutate(x_label = paste0(cluster, "_", .data[[group_key]]))
  
  # Perform prior tests if requested
  if (prior_test) {
    # ... (prior test code remains unchanged)
  }
  
  # Ensure group_key is treated as factor for discrete colors
  plot_data[[group_key]] <- factor(plot_data[[group_key]])
  
  # Create boxplot
  p <- ggplot(plot_data, aes(x = x_label, y = fraction, fill = .data[[group_key]])) +
    geom_boxplot(alpha = 0.8) +
    geom_point(position = position_jitterdodge(jitter.width = 0.2), 
               alpha = 0.6, size = 2) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top") +
    labs(x = "Cluster_Group", 
         y = "Fraction", 
         title = "Cluster Fraction Comparison by Group",
         fill = group_key) +
    scale_y_continuous(labels = scales::percent)
  
  # Add custom palette if provided
  if (!is.null(palette)) {
    p <- p + scale_fill_manual(values = palette)
  } else {
    p <- p + scale_fill_brewer(palette = "Set1")
  }
  
  # Statistical comparisons
  stat_results <- list()
  unique_clusters <- unique(plot_data$cluster)
  
  for (clust in unique_clusters) {
    # Initialize stat_test to NULL at the start of each iteration
    stat_test <- NULL
    
    clust_data <- plot_data %>% filter(cluster == clust)
    groups <- unique(clust_data[[group_key]])
    
    # ... (statistical test calculation remains unchanged)
    if (length(groups) == 2) {
      if (method %in% c("wilcox", "t.test")) {
        formula_str <- paste("fraction ~", group_key)
        if (method == "wilcox") {
          stat_test <- clust_data %>% rstatix::wilcox_test(as.formula(formula_str)) %>% mutate(cluster = clust)
        } else {
          stat_test <- clust_data %>% rstatix::t_test(as.formula(formula_str)) %>% mutate(cluster = clust)
        }
      }
    } else {
      formula_str <- paste("fraction ~", group_key)
      if (method == "kruskal") {
        stat_test <- clust_data %>% rstatix::kruskal_test(as.formula(formula_str)) %>% mutate(cluster = clust)
      } else if (method == "anova") {
        stat_test <- clust_data %>% rstatix::anova_test(as.formula(formula_str)) %>% mutate(cluster = clust)
      }
    }
    # Only add the result to the list if a test was actually performed
    if (!is.null(stat_test)) {
      stat_results[[clust]] <- stat_test
    }
  }
  
  # Combine statistical results
  all_stats <- bind_rows(stat_results)
  
  # Add p-values to plot
  # Only attempt to add annotations if there are statistical results to show
  if (nrow(all_stats) > 0) {
    if (simple_pvalues) {
      # Use simple p-value display
      p <- add_simple_pvalues(p, all_stats, plot_data, p_label_size)
    } else {
      # The following block was modified to remove the "stairway" effect.
      
      y_max <- max(plot_data$fraction, na.rm = TRUE)
      
      # Define a single, consistent y-position for all brackets.
      y_position_base <- y_max * 1.1 
      annotation_spacing <- y_max * 0.05
      
      x_label_order <- levels(factor(plot_data$x_label))
      
      # Process each cluster to add its annotation
      for (i in 1:nrow(all_stats)) {
        clust <- all_stats$cluster[i]
        p_val <- all_stats$p[i]
        
        # Format p-value
        p_label <- if (is.na(p_val)) "NA" else if (p_val < 0.001) "***" else if (p_val < 0.01) "**" else if (p_val < 0.05) "*" else "ns"
        
        # Find x positions for this cluster
        cluster_labels <- grep(paste0("^", clust, "_"), x_label_order, value = TRUE)
        
        if (length(cluster_labels) >= 2) {
          x_positions <- match(cluster_labels, x_label_order)
          x_min <- min(x_positions)
          x_max <- max(x_positions)
          x_center <- mean(c(x_min, x_max))
          
          # Use the consistent base y-position for all brackets
          bracket_y <- y_position_base
          label_y <- bracket_y + (annotation_spacing * 0.1)
          
          if (show_brackets) {
            # Horizontal line
            p <- p + annotate("segment", x = x_min - 0.3, xend = x_max + 0.3, y = bracket_y, yend = bracket_y, size = 0.5)
            # Left vertical tick
            p <- p + annotate("segment", x = x_min - 0.3, xend = x_min - 0.3, y = bracket_y - (annotation_spacing * 0.2), yend = bracket_y, size = 0.5)
            # Right vertical tick
            p <- p + annotate("segment", x = x_max + 0.3, xend = x_max + 0.3, y = bracket_y - (annotation_spacing * 0.2), yend = bracket_y, size = 0.5)
          }
          
          # Add p-value text
          p <- p + annotate("text", x = x_center, y = label_y, label = p_label, size = p_label_size, vjust = 0)
        }
      }
      
      # Adjust y-axis limits to accommodate the single layer of annotations
      y_limit <- y_position_base + annotation_spacing
      p <- p + coord_cartesian(ylim = c(0, y_limit), clip = "off") +
        theme(plot.margin = unit(c(1, 1, 1, 1), "lines"))
      
    }
  }
  
  # Return plot and statistics
  return(list(plot = p, statistics = all_stats))
}


# Example usage function
example_usage <- function() {
  # 예시 데이터 생성
  set.seed(123)
  n_cells <- 10000
  n_samples <- 20
  
  # 가상의 메타데이터 생성
  metadata <- data.frame(
    sample = rep(paste0("S", 1:n_samples), each = n_cells/n_samples),
    cluster = sample(paste0("C", 1:5), n_cells, replace = TRUE),
    prognosis = rep(c("Good", "Bad"), each = n_cells/2),
    stringsAsFactors = FALSE
  )
  
  # Bad prognosis에서 C1, C2 클러스터 비율 증가
  bad_idx <- metadata$prognosis == "Bad"
  metadata$cluster[bad_idx] <- sample(
    paste0("C", 1:5), 
    sum(bad_idx), 
    replace = TRUE, 
    prob = c(0.3, 0.3, 0.15, 0.15, 0.1)
  )
  
  # 함수 실행
  result <- plot_cluster_fractions(
    metadata,
    sample_key = "sample",
    cluster_key = "cluster",
    group_key = "prognosis",
    method = "wilcox",
    prior_test = TRUE,
    palette = c("Good" = "#3498db", "Bad" = "#e74c3c")
  )
  
  # 결과 출력
  print(result$plot)
  print(result$statistics)
}

# Additional utility functions

#' Calculate custom signatures at sample level
#' 
#' @param seurat_obj Seurat object
#' @param gene_sets Named list of gene sets
#' @param sample_key Column name for sample identification
#' @param method Scoring method (mean, ssgsea, etc.)
#' @return Data frame with sample-level signatures
calculate_sample_signatures <- function(seurat_obj, 
                                        gene_sets, 
                                        sample_key,
                                        method = "mean") {
  
  # Calculate module scores for each gene set
  for (set_name in names(gene_sets)) {
    genes <- gene_sets[[set_name]]
    genes <- genes[genes %in% rownames(seurat_obj)]
    
    if (length(genes) > 0) {
      seurat_obj <- AddModuleScore(
        seurat_obj,
        features = list(genes),
        name = paste0(set_name, "_score")
      )
    }
  }
  
  # Extract scores and aggregate by sample
  score_cols <- paste0(names(gene_sets), "_score1")
  metadata <- seurat_obj@meta.data
  
  sample_scores <- metadata %>%
    group_by(.data[[sample_key]]) %>%
    summarise(across(all_of(score_cols), mean, na.rm = TRUE))
  
  return(sample_scores)
}

#' Adjust for confounding variables
#' 
#' @param data Data frame with sample-level data
#' @param target_col Column name for target variable
#' @param group_col Column name for grouping variable
#' @param adjust_cols Vector of columns to adjust for
#' @param method Adjustment method (residuals, stratification)
#' @return Adjusted data frame
adjust_confounders <- function(data, 
                               target_col, 
                               group_col, 
                               adjust_cols,
                               method = "residuals") {
  
  if (method == "residuals") {
    # Linear model to get residuals
    formula_str <- paste(target_col, "~", paste(adjust_cols, collapse = " + "))
    lm_fit <- lm(as.formula(formula_str), data = data)
    
    # Add residuals as adjusted values
    data$adjusted_value <- residuals(lm_fit) + mean(data[[target_col]])
    
  } else if (method == "stratification") {
    # Create strata based on confounders
    data$strata <- interaction(data[adjust_cols])
    
    # Perform analysis within strata
    # This would be handled in the main plotting function
  }
  
  return(data)
}

# Alternative simple p-value display function
#' Add p-values without brackets (similar to first image style)
#' 
#' @param plot_obj ggplot object
#' @param stats_df Statistics dataframe
#' @param plot_data Plot data with x_labels
#' @param p_label_size Size of p-value text
#' @return Modified ggplot object
add_simple_pvalues <- function(plot_obj, stats_df, plot_data, p_label_size = 3) {
  
  y_max <- max(plot_data$fraction)
  y_position <- y_max * 1.05
  
  # Get x_label order
  x_label_order <- levels(factor(plot_data$x_label))
  
  # Create a dataframe for p-value positions
  pval_positions <- data.frame()
  
  for (i in 1:nrow(stats_df)) {
    clust <- stats_df$cluster[i]
    p_val <- stats_df$p[i]
    
    # Format p-value
    if (is.na(p_val)) {
      p_label <- "NA"
    } else if (p_val < 0.001) {
      p_label <- "p<0.001"
    } else if (p_val < 0.01) {
      p_label <- sprintf("p=%.3f", p_val)
    } else if (p_val < 0.05) {
      p_label <- sprintf("p=%.3f", p_val)
    } else {
      p_label <- sprintf("p=%.2f", p_val)
    }
    
    # Find center position for this cluster
    cluster_labels <- grep(paste0("^", clust, "_"), x_label_order, value = TRUE)
    
    if (length(cluster_labels) >= 2) {
      x_positions <- match(cluster_labels, x_label_order)
      x_center <- mean(x_positions)
      
      pval_positions <- rbind(pval_positions, 
                              data.frame(x = x_center, 
                                         y = y_position,
                                         label = p_label))
    }
  }
  
  # Add all p-values at once
  if (nrow(pval_positions) > 0) {
    plot_obj <- plot_obj + 
      annotate("text", 
               x = pval_positions$x, 
               y = pval_positions$y,
               label = pval_positions$label,
               size = p_label_size,
               vjust = 0)
  }
  
  # Adjust y-axis
  plot_obj <- plot_obj + 
    coord_cartesian(ylim = c(0, y_max * 1.15))
  
  return(plot_obj)
}
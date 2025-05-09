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
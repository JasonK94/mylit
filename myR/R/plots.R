#' Statistical Analysis for scRNA-seq Metadata
#'
#' This function performs comprehensive statistical analysis on scRNA-seq data
#' comparing groups across different categorical variables.
#' @import dplyr ggplot2
#' @param
#'
#' @return
#' @export
scatter_smooth = function(sobj, feature, clinical_variable = "nih_change", transpose = FALSE) {
  
  # 1) Make a per-cell data.frame with the expression for 'feature'
  df <- data.frame(
    sample_no = sobj@meta.data$sample_no,
    FEATURE = as.numeric(FetchData(sobj, vars = feature)[, feature])
  )
  
  # 2) Compute average expression per patient
  df_avg <- df %>%
    group_by(sample_no) %>%
    summarise(avg_FEATURE = mean(FEATURE, na.rm = TRUE))
  
  # 3) Merge with your clinical variable
  meta_patient <- sobj@meta.data %>%
    select(sample_no, all_of(clinical_variable)) %>%
    distinct(sample_no, .keep_all = TRUE)
  
  df_merged <- left_join(df_avg, meta_patient, by = "sample_no") %>%
    mutate(
      nih_change_numeric = as.numeric(as.character(.data[[clinical_variable]]))
    )
  
  # 4) Branch logic for transpose:
  #    transpose = FALSE -> x = clinical_variable, y = avg_FEATURE
  #    transpose = TRUE  -> x = avg_FEATURE,         y = clinical_variable
  
  if (!transpose) {
    # Model: avg_FEATURE ~ nih_change_numeric
    model <- lm(avg_FEATURE ~ nih_change_numeric, data = df_merged)
    summary(model)
    
    intercept <- round(coef(model)[1], 3)
    slope <- round(coef(model)[2], 3)
    pval <- signif(summary(model)$coefficients[2, 4], 3)
    label_text <- paste0("y = ", intercept, " + ", slope, " * x\np = ", pval)
    
    # Plot: x = nih_change_numeric, y = avg_FEATURE
    p = ggplot(df_merged, aes(x = nih_change_numeric, y = avg_FEATURE)) +
      geom_point() +
      geom_smooth(method = "lm", se = TRUE) +
      annotate(
        "text",
        x = min(df_merged$nih_change_numeric, na.rm = TRUE),
        y = max(df_merged$avg_FEATURE, na.rm = TRUE),
        label = label_text, hjust = 0, vjust = 1, size = 5
      ) +
      theme_bw() +
      xlab(clinical_variable) +
      ylab(paste0("Average ", feature, " Expression"))
    
  } else {
    # transpose = TRUE -> we flip the roles:
    # Model: nih_change_numeric ~ avg_FEATURE
    model <- lm(nih_change_numeric ~ avg_FEATURE, data = df_merged)
    summary(model)
    
    intercept <- round(coef(model)[1], 3)
    slope <- round(coef(model)[2], 3)
    pval <- signif(summary(model)$coefficients[2, 4], 3)
    label_text <- paste0("y = ", intercept, " + ", slope, " * x\np = ", pval)
    
    # Plot: x = avg_FEATURE, y = nih_change_numeric
    p = ggplot(df_merged, aes(x = avg_FEATURE, y = nih_change_numeric)) +
      geom_point() +
      geom_smooth(method = "lm", se = TRUE) +
      annotate(
        "text",
        x = min(df_merged$avg_FEATURE, na.rm = TRUE),
        y = max(df_merged$nih_change_numeric, na.rm = TRUE),
        label = label_text, hjust = 0, vjust = 1, size = 5
      ) +
      theme_bw() +
      xlab(paste0("Average ", feature, " Expression")) +
      ylab(clinical_variable)
  }
  
  return(p)
}
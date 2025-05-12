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


#' helper for cmb, acmb
sort_samples <- function(samples) {
  # Helper function to extract numbers from a string
  extract_numbers <- function(x) {
    nums <- as.numeric(strsplit(x, "\\+")[[1]])
    if (length(nums) == 0 || any(is.na(nums))) return(c(Inf, Inf))
    return(nums)
  }
  
  # Check if all samples are non-numeric (no '+' and can't be converted to number)
  all_non_numeric <- all(sapply(samples, function(x) {
    !grepl("\\+", x) && is.na(suppressWarnings(as.numeric(x)))
  }))
  
  if (all_non_numeric) {
    return(sort(samples))  # If all non-numeric, just return alphabetically sorted
  }
  
  # Separate samples into those with and without "+"
  single_samples <- samples[!grepl("\\+", samples)]
  doublet_samples <- samples[grepl("\\+", samples)]
  
  # Sort single samples
  single_sorted <- single_samples[order(sapply(single_samples, function(x) {
    num <- suppressWarnings(as.numeric(x))
    if (is.na(num)) Inf else num
  }))]
  
  # Sort doublet samples
  doublet_sorted <- doublet_samples[order(
    sapply(doublet_samples, function(x) extract_numbers(x)[1]),  # First number
    sapply(doublet_samples, function(x) extract_numbers(x)[2])   # Second number
  )]
  
  # Combine and return
  return(c(single_sorted, doublet_sorted))
}


#' This function shows relative frequency of each cluster per each sample(or group)
#' 
#' @import dplyr ggplot2
#' @param
#'
#' @return
#' @export
cmb <- function(sobj, identity = "seurat_clusters", group.by = "sample", idents = NULL, df=F, vlines=NULL, vline_color = "red") {
  Idents(sobj) <- identity
  cluster_ids <- Idents(sobj)
  sample_ids <- sobj@meta.data[[group.by]]
  
  data <- data.frame(cluster = cluster_ids, sample = sample_ids)
  
  # Filter identities if idents is provided
  if(!is.null(idents)) {
    data <- data[data$cluster %in% idents,]
  }
  
  summary_data <- data %>%
    group_by(sample, cluster) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(sample) %>%
    mutate(proportion = count / sum(count)) %>%
    ungroup()
  
  # Sort the samples using the new sort_samples function
  sorted_samples <- sort_samples(unique(as.character(summary_data$sample)))
  summary_data$sample <- factor(summary_data$sample, levels = sorted_samples)
  
  output=ggplot(summary_data, aes(x = sample, y = proportion, fill = cluster)) +
    geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.2) +
    labs(title = paste0("Proportional Bar Graph of Clusters for Each ", group.by),
         x = group.by,
         y = "Proportion of Cells") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right") +
    scale_fill_viridis_d() +
    scale_y_continuous(labels = scales::percent_format())
  
  if(!is.null(vlines)){
    for (line_pos in vlines) {
      output <- output + geom_vline(xintercept = line_pos + 0.5, color = vline_color, linetype = "dashed", linewidth = 0.5)
    }
  }
  
  if(df){
    return(summary_data)
  }
  return(output)
}


#' This function shows absolute counts of each cluster per each sample(or group)
#' 
#' @import dplyr ggplot2
#' @param
#'
#' @return
#' @export
acmb <- function(sobj, identity="seurat_clusters", group.by="sample", idents = NULL, df=F, vlines=NULL, vline_color = "red") {
  Idents(sobj) <- identity
  cluster_ids <- Idents(sobj)
  sample_ids <- sobj@meta.data[[group.by]]
  
  data <- data.frame(cluster = cluster_ids, sample = sample_ids)
  
  # Filter identities if idents is provided
  if(!is.null(idents)) {
    data <- data[data$cluster %in% idents,]
  }
  
  summary_data <- data %>%
    group_by(sample, cluster) %>%
    summarise(count = n(), .groups = "drop")
  
  # Sort the samples using the new sort_samples function
  sorted_samples <- sort_samples(unique(summary_data$sample))
  summary_data$sample <- factor(summary_data$sample, levels = sorted_samples)
  
  output=ggplot(summary_data, aes(x = sample, y = count, fill = cluster)) +
    geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.2) +
    labs(title = paste0("Cumulative Bar Graph of Clusters for Each ", group.by),
         x = group.by,
         y = "Number of Cells") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right") +
    scale_fill_viridis_d()
  
  if(!is.null(vlines)){
    for (line_pos in vlines) {
      output <- output + geom_vline(xintercept = line_pos + 0.5, color = vline_color, linetype = "dashed", linewidth = 0.5)
    }
  }
  
  if(df){
    return(summary_data)
  }
  return(output)
}

#' awesome function for making heatmap(tile plot) of gene modules
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
  
  #-------------------------------
  # (A) 유효성 체크
  #-------------------------------
  if(is.null(gene_sets)){
    stop("gene_sets를 지정해 주세요. (ex: list(Immune=c('CD3D','CD3E'), Bcell=c('MS4A1','CD79A'))) ")
  }
  
  # 만약 gene_sets가 리스트가 아니라 벡터만 들어왔다면 리스트로 변환
  # 예: c("CD3D","CD3E") -> list(GeneSet1 = c("CD3D","CD3E"))
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
  
  #-------------------------------
  # (B) Seurat 객체에 grouping 적용
  #-------------------------------
  Idents(sobj) <- group
  
  #-------------------------------
  # (C) 평균 발현량(또는 합계 등) 계산
  #-------------------------------
  if(value == "average"){
    # group.by = group 로 명시
    cluster_avg <- AverageExpression(sobj, assays = assay, slot = "data", group.by = group)[[assay]]
  } else {
    cluster_avg <- AggregateExpression(sobj, assays = assay, slot = "data", group.by = group)[[assay]]
  }
  
  #-------------------------------
  # (D) Gene Set 별 발현량 계산
  #-------------------------------
  # cluster_avg의 컬럼은 cluster 이름이 된다.
  cluster_names <- colnames(cluster_avg)
  
  # 결과를 담을 data.frame 생성
  gene_set_expression <- data.frame(Cluster = cluster_names, stringsAsFactors = FALSE)
  
  # gene_sets 각각에 대해 평균 발현량을 구함
  for(gset_name in names(gene_sets)){
    genes <- gene_sets[[gset_name]]
    genes_present <- genes[genes %in% rownames(cluster_avg)]
    
    if(length(genes_present) == 0){
      warning(paste("No genes from", gset_name, "found in the dataset."))
      # 데이터프레임에 NA 열을 넣고 다음으로 넘어감
      gene_set_expression[[gset_name]] <- NA
      next
    }
    
    # colMeans를 이용해, 해당 유전자들의 평균 발현량 계산
    gene_set_expression[[gset_name]] <- colMeans(cluster_avg[genes_present, , drop = FALSE])
  }
  
  #-------------------------------
  # (E) Z-score 정규화
  #-------------------------------
  # 첫 번째 열(Cluster)을 제외한 나머지를 scale()
  gene_set_expression_normalized <- gene_set_expression
  gene_set_expression_normalized[,-1] <- scale(gene_set_expression_normalized[,-1])
  
  # 각 Cluster에서 가장 높은 값을 가지는 gene set을 배정해보자(부가 기능)
  gene_set_expression_normalized$Assigned_CellType <- apply(
    gene_set_expression_normalized[,-1], 1, 
    function(x){
      names(x)[which.max(x)]
    }
  )
  
  #-------------------------------
  # (F) 클러스터 순서 정렬
  #-------------------------------
  # 사용자가 만든 cluster 이름이 꼭 숫자일 필요는 없으므로,
  # 1) 전부 숫자로 바꿀 수 있다면 numeric 정렬
  # 2) 아니면 문자 알파벳 순 정렬
  
  # 임시로 numeric 변환
  numeric_test <- suppressWarnings(as.numeric(gene_set_expression_normalized$Cluster))
  
  if(!all(is.na(numeric_test))){
    # NA가 아닌 값이 있다 => 전부 숫자로 파싱되는 경우
    # 실제로 모두 정상 변환인지 다시 확인 (NA가 하나라도 있으면 문자)
    if(sum(is.na(numeric_test)) == 0){
      # 전부 숫자면 해당 순서로 factor 설정
      sorted_levels <- (unique(numeric_test))
      gene_set_expression_normalized$Cluster <- factor(
        gene_set_expression_normalized$Cluster,
        levels = as.character(sorted_levels)
      )
    } else {
      # 일부만 숫자인 경우 => 그냥 문자 정렬
      sorted_levels <- (unique(gene_set_expression_normalized$Cluster))
      gene_set_expression_normalized$Cluster <- factor(
        gene_set_expression_normalized$Cluster,
        levels = sorted_levels
      )
    }
  } else {
    # 전부 NA => 아예 숫자로 파싱 불가 -> 문자 정렬
    sorted_levels <- (unique(gene_set_expression_normalized$Cluster))
    gene_set_expression_normalized$Cluster <- factor(
      gene_set_expression_normalized$Cluster,
      levels = sorted_levels
    )
  }
  
  #-------------------------------
  # (G) Heatmap용 long format 만들기
  #-------------------------------
  melted_data <- melt(
    gene_set_expression_normalized,
    id.vars = c("Cluster","Assigned_CellType"),
    variable.name = "GeneSet",
    value.name = "Zscore"
  )
  
  #-------------------------------
  # (H) Heatmap 그리기
  #-------------------------------
  p <- ggplot(melted_data, aes(x = Cluster, y = GeneSet, fill = Zscore)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    theme_minimal() +
    labs(
      title = title,
      x = x_label,
      y = y_label
    ) +
    theme(axis.text.x = element_text(hjust=0.5,size=16),
          axis.text.y = element_text(hjust=1,size=16),
          axis.title.x = element_text(hjust=0.5,face="bold",size=16),
          axis.title.y = element_text(hjust=0.5,face="bold",size=16),
          plot.title=element_text(size=16,face="bold",hjust=0.5))
  
  print(p)
  
  #-------------------------------
  # (I) 결과 반환
  #-------------------------------
  # Assigned_CellType 까지 붙어있는 최종 테이블 반환
  return(gene_set_expression_normalized)
}



#' title = "Normalized Gene Set Expression per Cluster",
#' x = "Cluster",
#' y = "Gene Set"
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
  library(Seurat)
  library(dplyr)
  library(reshape2)
  library(ggplot2)
  
  #-------------------------------
  # (A) 유효성 체크
  #-------------------------------
  if(is.null(genes)){
    stop("genes를 지정해 주세요. (예: genes=c('CD3D','CD3E','MS4A1'))")
  }
  
  #-------------------------------
  # (B) Seurat 객체에 grouping 적용
  #-------------------------------
  Idents(sobj) <- group
  
  #-------------------------------
  # (C) 평균 발현량(또는 합계 등) 계산
  #-------------------------------
  if(value == "average"){
    cluster_avg <- AverageExpression(sobj, assays = assay, slot = "data", group.by = group)[[assay]]
  } else {
    cluster_avg <- AggregateExpression(sobj, assays = assay, slot = "data", group.by = group)[[assay]]
  }
  
  #-------------------------------
  # (D) 원하는 유전자만 필터
  #-------------------------------
  genes_present <- genes[genes %in% rownames(cluster_avg)]
  if(length(genes_present) == 0){
    stop("지정하신 유전자 중 데이터셋에 존재하는 유전자가 없습니다.")
  }
  
  # Subset 후에 (행=유전자, 열=클러스터)
  # 행렬을 (열=클러스터, 행=유전자) 형태로 보고 싶으면 Transpose
  # Heatmap을 그릴 때는 보통 row=유전자, col=클러스터가 익숙하므로, 아래처럼 melt를 하려면
  # 먼저 t() 한 뒤 scale() 적용한 다음, 다시 melt 시 row는 cluster로, column은 gene이 되도록 했습니다.
  gene_expression <- cluster_avg[genes_present, , drop=FALSE]
  # gene_expression: rows=genes, cols=clusters
  
  # Z-score 정규화를 위해 t() (rows=clusters, cols=genes)
  gene_expression <- t(gene_expression)
  gene_expression <- scale(gene_expression)
  
  # data.frame으로 변환
  gene_expression <- as.data.frame(gene_expression)
  gene_expression$Cluster <- rownames(gene_expression)
  
  #-------------------------------
  # (E) 클러스터 순서 정렬
  #-------------------------------
  numeric_test <- suppressWarnings(as.numeric(gene_expression$Cluster))
  
  if(!all(is.na(numeric_test))){
    # 전부 숫자로 파싱되는 경우
    if(sum(is.na(numeric_test)) == 0){
      sorted_levels <- sort(unique(numeric_test))
      gene_expression$Cluster <- factor(
        gene_expression$Cluster,
        levels = as.character(sorted_levels)
      )
    } else {
      # 일부만 숫자인 경우 => 문자 정렬
      sorted_levels <- sort(unique(gene_expression$Cluster))
      gene_expression$Cluster <- factor(
        gene_expression$Cluster,
        levels = sorted_levels
      )
    }
  } else {
    # 전부 NA => 문자 정렬
    sorted_levels <- sort(unique(gene_expression$Cluster))
    gene_expression$Cluster <- factor(
      gene_expression$Cluster,
      levels = sorted_levels
    )
  }
  
  #-------------------------------
  # (F) long format으로 melt
  #-------------------------------
  melted_data <- melt(
    gene_expression,
    id.vars = "Cluster",
    variable.name = "Gene",
    value.name = "Zscore"
  )
  
  #-------------------------------
  # (G) Heatmap 그리기
  #-------------------------------
  p <- ggplot(melted_data, aes(x = Cluster, y = Gene, fill = Zscore)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    theme_minimal() +
    labs(
      title = title,
      x = x_label,
      y = y_label
    ) +
    theme(axis.text.x = element_text(angle=45, hjust=1),
          plot.title=element_text(size=14,face="bold",hjust=0.5))
  
  print(p)
  
  #-------------------------------
  # (H) 결과 반환
  #-------------------------------
  # 정규화된 수치를 담고 있는 wide-format data.frame 반환
  return(gene_expression)
}

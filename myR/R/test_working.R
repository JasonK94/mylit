# -------- 설명 --------
# https://www.notion.so/drgoze/205-myR-27f650f736bf80f08125e4b1f3896125?source=copy_link
# PIFG: patient by patient (pre-post)...
# PDH: drug by drug (delta)...

# MDMP: patient별로 delta matrix 구축. -> seurat화 -> clustering (엄청 오래걸림)

# DOGE: drug 별로 delta 반환 - one gene
# MDEM: drug 별로 delta 반환 - multiple genes -> ? -> PDRM (엄청 오래걸림)
# PDRM: delta matrix 시각화 (heatmap; 그 굉장히 빽뺵한거..)

# PIFG(sobj, gene, subse_...) res_IL1B_CK2=plot_interaction_for_gene(data.ck.adj, "IL1B",nested=FALSE,segment_col=NULL)
# PDH(pifg(sobj, gene_list(vector)))

# -------- PIFG: sobj -> spaghetti plot, emmeans plot --------
#' @export
plot_interaction_for_gene <- function(
  sobj, gene, assay = "RNA", layer = "data",
  drug_col = "drug", nested = FALSE, ref_drug = NULL, 
  time_col = "treatment", time_levels = c("pre","post"), time_labels = c("pre","post"),
  segment_col = "ck",subset_comp = NULL,
  patient_col = "emrid",
  optimizer = "bobyqa"
) {
  # ==== Data prep ====
  expr <- GetAssayData(sobj, assay = assay, slot = layer)
  meta <- sobj@meta.data %>%
    mutate(
      time    = factor(.data[[time_col]], levels = time_levels, labels = time_labels),
      drug    = factor(.data[[drug_col]]),
      patient = factor(.data[[patient_col]])
    )
  if (!is.null(segment_col)){
    meta=meta%>%mutate(
      comp=.data[[segment_col]]
    )
  }else{subset_comp=NULL} #현재 상태로는 comp가 있을지 여부가 갈리므로 후속 분석에 안 좋을 가능성이 크다. 현재는 문제없음.
  
  if (!is.null(subset_comp)) {
    keep <- meta$comp == subset_comp
  } else {
    keep <- rep(TRUE, nrow(meta))
  }
  
  if (!is.null(ref_drug)) meta$drug <- relevel(meta$drug, ref = ref_drug)
  
  y <- as.numeric(expr[gene, keep])
  df <- cbind.data.frame(y = y, meta[keep, , drop = FALSE])

  # ==== Checks ====
  if (all(is.na(y)) || sd(y, na.rm = TRUE) == 0) {
    stop(sprintf("Gene %s has no variance or all NA in this subset.", gene))
  }
  if (nlevels(df$time) < 2) stop("time has <2 levels in this subset.")
  if (nlevels(df$patient) < 2) stop("patient has <2 levels in this subset.")

  # ==== Model ====
  ctrl <- lmerControl(optimizer = optimizer, calc.derivs = TRUE,
                      check.conv.singular = "ignore")
  if (nested) {
    fit <- lmer(y ~ time * drug + (1 | drug/patient),
                data = df, REML = FALSE, control = ctrl)
  } else {
    fit <- lmer(y ~ time * drug + (1 | patient),
                data = df, REML = FALSE, control = ctrl)
  }

  # ==== emmeans ====
  emm <- emmeans(fit, ~ time * drug)
  emm_df <- as.data.frame(confint(emm)) %>%
    mutate(time = factor(time, levels = time_levels, labels = time_labels))

  p1 <- ggplot(emm_df, aes(x = time, y = emmean, group = drug, color = drug)) +
    geom_point(size = 3) +
    geom_line(size = 1) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1) +
    labs(
      title = sprintf("%s: time × drug (EMMeans ±95%% CI)", gene),
      x = "Time", y = "Estimated mean expression"
    ) +
    theme_bw()

  # ==== spaghetti ====
  df_plot <- df %>%
    dplyr::select(y, time, drug, patient) %>%
    mutate(time = factor(time, levels = time_levels, labels = time_labels))

  p2 <- ggplot(df_plot, aes(x = time, y = y, group = interaction(patient, drug), color = drug)) +
    geom_line(alpha = 0.35) +
    geom_point(size = 2, alpha = 0.7) +
    labs(
      title = sprintf("%s: individual pre→post (raw)", gene),
      x = "Time", y = "Expression (Q3/log scale)"
    ) +
    theme_bw()

  # ==== delta (post-pre) ====
  em_time_by_drug <- emmeans(fit, ~ time | drug)
  deltas <- contrast(em_time_by_drug, method = list(delta = c(-1, 1)))
  deltas_df <- as.data.frame(confint(deltas)) %>%
    mutate(gene = gene)

  p3 <- ggplot(deltas_df, aes(x = drug, y = estimate,
                              ymin = lower.CL, ymax = upper.CL, color = drug)) +
    geom_pointrange(position = position_dodge(width = 0.2)) +
    geom_hline(yintercept = 0, linetype = 2) +
    labs(
      title = sprintf("%s: Δ(post-pre) by drug (EMMeans ±95%% CI)", gene),
      x = "Drug", y = expression(Delta~"(post - pre)")
    ) +
    theme_bw() +
    theme(legend.position = "none")

  # Return
  list(
    df = df,
    fit = fit,
    emm = emm_df,
    deltas = deltas_df,
    plot_emmeans = p1,
    plot_spaghetti = p2,
    plot_delta = p3
  )
}

# -------- pifg: PIFG lapply  -> heatmap --------
pifg=function(sobj, gene_list){
  return(lapply(gene_list, function(g) plot_interaction_for_gene(sobj, g, segment_col=NULL)) %>% set_names(gene_list))
}

# -------- PDH: results_list from pifg  -> heatmap --------
#' @export
plot_delta_heatmap <- function(results_list) {
  # results_list: list of outputs from plot_interaction_for_gene (여러 gene)
  all_deltas <- purrr::map_dfr(results_list, "deltas")
  
  # drug x gene 매트릭스
  mat <- all_deltas %>%
    dplyr::select(gene, drug, estimate) %>%
    tidyr::pivot_wider(names_from = drug, values_from = estimate)
  
  df_long <- all_deltas
  
  p <- ggplot(df_long, aes(x = drug, y = gene, fill = estimate)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
    labs(title = "Δ(post-pre) heatmap by drug",
         x = "Drug", y = "Gene", fill = expression(Delta)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  list(matrix = mat, plot = p, raw = df_long)
}

# -------- MDMP --------
#' @export
make_delta_matrix_by_patient <- function(
  sobj,
  assay="RNA", layer="data",
  time_col="treatment", drug_col="drug", patient_col="emrid", segment_col="ck",
  time_levels=c("pre","post"), subset_comp=NULL,
  agg_fun = function(v) mean(v, na.rm=TRUE)
){
  X <- GetAssayData(sobj, assay=assay, layer=layer) # genes x AOIs
  meta <- sobj@meta.data %>%
    mutate(
      time    = factor(.data[[time_col]], levels=time_levels),
      drug    = factor(.data[[drug_col]]),
      patient = factor(.data[[patient_col]]),
      comp    = factor(.data[[segment_col]])
    )
  stopifnot(identical(colnames(X), rownames(meta)))

  if (!is.null(subset_comp)) {
    keep <- meta$comp == subset_comp
    X <- X[, keep, drop=FALSE]
    meta <- meta[keep, , drop=FALSE]
  }
  
  # meta column filtering
  gene_names <- rownames(X)
  
  run=0#
  while(length(which(gene_names%in%names(meta)))>0){
    num_meta_in_gene=which(names(meta)%in%gene_names)
    names(meta)[num_meta_in_gene]=paste0(names(meta)[num_meta_in_gene],"_re")
    run=run+1#
    print(run)#
    if(run>5)break#
  }
  # AOI별 표현을 long으로
  long <- as.data.frame(t(X)) %>%
    tibble::rownames_to_column("AOI") %>%
    cbind(meta, .)

  # pre/post를 같은 환자/약제/comp 안에서 집계 후 Δ = post − pre
  # 메모리 절약을 위해 gene-wise로 루프
  make_delta_one_gene <- function(g){
    dd <- long %>%
      dplyr::select(patient, drug, comp, time, !!sym(g)) %>%
      group_by(patient, drug, comp, time) %>%
      summarise(expr = agg_fun(.data[[g]]), .groups="drop") %>%
      tidyr::pivot_wider(names_from = time, values_from = expr) %>%
      mutate(delta = .data[[time_levels[2]]] - .data[[time_levels[1]]]) %>%
      dplyr::select(patient, drug, comp, delta)
    dd$gene <- g
    dd
  }

  delta_long <- purrr::map_dfr(gene_names, make_delta_one_gene)

  # 샘플 정의: 환자*약제*(comp)
  delta_long <- delta_long %>%
    mutate(sample_id = if (!is.null(subset_comp)) 
             paste(patient, drug, sep="|") else
             paste(patient, drug, comp, sep="|"))

  # gene x sample 행렬
  delta_mat <- delta_long %>%
    dplyr::select(gene, sample_id, delta) %>%
    tidyr::pivot_wider(names_from = sample_id, values_from = delta) %>%
    tibble::column_to_rownames("gene") %>%
    as.matrix()

  # 메타(샘플) 테이블
  sample_meta <- delta_long %>%
    distinct(sample_id, patient, drug, comp) %>%
    tibble::column_to_rownames("sample_id")

  list(delta_mat = delta_mat, sample_meta = sample_meta)
}

# -------- DOGE -> MDEM --------
# 유전자 1개 -> deltas_df (drug별 Δ와 CI, p) 반환
#' @export
delta_one_gene_emm <- function(
  sobj, gene,
  nested = FALSE, subset_comp = NULL,
  assay = "RNA", layer = "data",
  time_col = "treatment", time_levels = c("pre","post"),
  drug_col = "drug", patient_col = "emrid", segment_col = "ck",
  ref_drug = NULL, optimizer = "bobyqa"
){
  X <- GetAssayData(sobj, assay=assay, slot=layer)
  meta <- sobj@meta.data %>%
    mutate(
      time    = factor(.data[[time_col]], levels=time_levels),
      drug    = factor(.data[[drug_col]]),
      patient = factor(.data[[patient_col]]),
      comp    = factor(.data[[segment_col]])
    )
  if (!is.null(ref_drug)) meta$drug <- relevel(meta$drug, ref=ref_drug)
  stopifnot(identical(colnames(X), rownames(meta)))

  keep <- if (is.null(subset_comp)) rep(TRUE, nrow(meta)) else meta$comp == subset_comp
  y <- as.numeric(X[gene, keep])
  df <- cbind.data.frame(y = y, meta[keep, , drop=FALSE])

  if (all(is.na(y)) || sd(y, na.rm=TRUE) == 0) return(NULL)
  if (nlevels(df$time) < 2 || nlevels(df$patient) < 2) return(NULL)

  ctrl <- lmerControl(optimizer=optimizer, calc.derivs=TRUE, check.conv.singular="ignore")
  form <- if (nested) y ~ time * drug + (1|drug/patient) else y ~ time * drug + (1|patient)
  fit <- try(lmer(form, data=df, REML=FALSE, control=ctrl), silent=TRUE)
  if (inherits(fit, "try-error")) return(NULL)

  em_time_by_drug <- emmeans(fit, ~ time | drug)
  deltas <- contrast(em_time_by_drug, method=list(delta=c(-1,1)))
  deltas_df <- as.data.frame(confint(deltas))
  deltas_df$gene <- gene
  deltas_df
}

# 여러 유전자 병렬 처리하여 Δ(post-pre) 행렬 (genes x drugs) 구축
#' @export
make_delta_emm_matrix <- function(
  sobj, genes,
  subset_comp = NULL, nested=FALSE,
  assay="RNA", layer="data",
  time_col="treatment", time_levels=c("pre","post"),
  drug_col="drug", patient_col="emrid", segment_col="ck",
  ref_drug=NULL, optimizer="bobyqa",
  mc.cores = 1
){
  fun <- function(g)
    delta_one_gene_emm(sobj, g, nested, subset_comp, assay, layer,
                       time_col, time_levels, drug_col, patient_col, segment_col,
                       ref_drug, optimizer)
  res_list <- if (mc.cores > 1 && .Platform$OS.type != "windows") {
    parallel::mclapply(genes, fun, mc.cores=mc.cores)
  } else {
    lapply(genes, fun)
  }
  res <- dplyr::bind_rows(res_list)
  if (nrow(res) == 0) stop("No gene produced a delta result.")

  # 행렬화
  mat <- res %>%
    dplyr::select(gene, drug, estimate) %>%
    tidyr::pivot_wider(names_from = drug, values_from = estimate) %>%
    tibble::column_to_rownames("gene") %>%
    as.matrix()

  list(delta_df = res, delta_mat = mat)
}

# -------- helper: pathway membership matrix (gene x pathway, 0/1) --------
.make_pathway_mat <- function(genes, gene_sets=list()){
  if (length(gene_sets)==0) return(NULL)
  pws <- names(gene_sets)
  m <- sapply(gene_sets, function(gv) as.integer(genes %in% gv))
  colnames(m) <- pws
  rownames(m) <- genes
  m
}

# -------- helper: p-value 별표 --------
.stars <- function(p){ ifelse(p < 0.001, "***",
                       ifelse(p < 0.01,  "**",
                       ifelse(p < 0.05,  "*", ""))) }

# -------- PDRM(MDEM$delta_df); options: top variable genes? or gene list. main plotting function --------
#' @export
plot_drug_response_map <- function(
  df_long,                      # gene, drug, estimate
  drug_annot = NULL,            # optional: drug, class, route ...
  gene_sets = list(),           # optional: list(pathway = c(genes))
  sig_df = NULL,                # optional: gene, drug, p_adj
  clip_q = c(0.01, 0.99),       # quantile clipping to tame outliers
  scale_row = TRUE,             # gene-wise z
  show_rownames = FALSE,
  show_colnames = TRUE,
  top_n_genes = NULL,           # 선택: 상위 변이 gene만 (예: 2000)
  row_split_k = NULL,           # 선택: 행 클러스터 k-split
  col_split_k = NULL            # 선택: 열 클러스터 k-split
){
  # 1) wide matrix
  mat <- df_long %>%
    select(gene, drug, estimate) %>%
    pivot_wider(names_from = drug, values_from = estimate) %>%
    tibble::column_to_rownames("gene") %>%
    as.matrix()

  # 2) optional gene filter by variance
  if (!is.null(top_n_genes)) {
    v <- apply(mat, 1, function(x) var(x, na.rm=TRUE))
    keep <- names(sort(v, decreasing=TRUE))[seq_len(min(top_n_genes, length(v)))]
    mat <- mat[keep, , drop=FALSE]
  }

  # 3) clipping
  qlims <- quantile(mat, probs = clip_q, na.rm = TRUE)
  mat_clip <- pmin(pmax(mat, qlims[1]), qlims[2])

  # 4) row scaling (gene-wise z)
  if (scale_row) {
    mat_use <- t(scale(t(mat_clip)))
  } else {
    mat_use <- mat_clip
  }
  mat_use[is.na(mat_use)] <- 0  # 안전장치

  # 5) color scale (대칭)
  col_fun <- colorRamp2(c(min(mat_use), 0, max(mat_use)), c("blue","white","red"))

  # 6) column annotation (drug info)
  ha_col <- NULL
  if (!is.null(drug_annot)) {
    drug_annot <- as.data.frame(drug_annot)
    rownames(drug_annot) <- drug_annot$drug
    drug_annot <- drug_annot[colnames(mat_use), setdiff(colnames(drug_annot),"drug"), drop=FALSE]
    # 범주형에 대한 색 자동 생성
    ann_cols <- lapply(drug_annot, function(v){
      if (is.numeric(v)) NULL else
        structure(scales::hue_pal()(length(unique(v))), names=unique(v))
    })
    ha_col <- HeatmapAnnotation(df = drug_annot, col = ann_cols, which = "column")
  }

  # 7) row annotation: pathway membership (bar/dot)
  ha_row <- NULL
  pw_mat <- .make_pathway_mat(rownames(mat_use), gene_sets)
  if (!is.null(pw_mat)) {
    # pathway count bar
    pw_count <- rowSums(pw_mat)
    ha_row <- rowAnnotation(
      pw_count = anno_barplot(pw_count, gp = gpar(fill="#888888", col=NA), width = unit(0.8,"cm"))
    )
  }

  # 8) 유의성 별표 overlay (선택)
  layer_fun <- NULL
  if (!is.null(sig_df)) {
    sig_df <- sig_df %>% mutate(star = .stars(p_adj))
    # matrix of stars aligned to mat_use
    star_mat <- matrix("", nrow=nrow(mat_use), ncol=ncol(mat_use),
                       dimnames = dimnames(mat_use))
    common_genes <- intersect(rownames(mat_use), sig_df$gene)
    common_drugs <- intersect(colnames(mat_use), sig_df$drug)
    sig_sub <- sig_df %>% filter(gene %in% common_genes, drug %in% common_drugs)
    idx_r <- match(sig_sub$gene, rownames(mat_use))
    idx_c <- match(sig_sub$drug, colnames(mat_use))
    star_mat[cbind(idx_r, idx_c)] <- sig_sub$star

    layer_fun <- function(j, i, x, y, w, h, fill){
      grid.text(star_mat[i, j], x, y, gp = gpar(col = "black", fontsize = 8))
    }
  }

  # 9) Heatmap 생성
  ht <- Heatmap(
    mat_use,
    name = "Δ (scaled)",
    col = col_fun,
    cluster_rows = TRUE, cluster_columns = TRUE,
    show_row_names = show_rownames,
    show_column_names = show_colnames,
    top_annotation = ha_col,
    left_annotation = ha_row,
    row_split = row_split_k,
    column_split = col_split_k,
    border = FALSE,
    rect_gp = gpar(col = NA),
    cell_fun = layer_fun,
    column_names_rot = 45
  )

  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  invisible(list(matrix = mat_use, heatmap = ht))
}



# -------- PCP: Plot cluster proportions, and other functions --------
#' 클러스터 비율 비교 박스플롯 생성 (Pseudobulk 방식)
#'
#' @param sobj Seurat 객체
#' @param patient_col 메타데이터 내 환자 또는 샘플 ID 컬럼명
#' @param group.by 비교할 클러스터 컬럼명
#' @param split.by 비교할 그룹 컬럼명
#' @param test.use 사용할 통계 검정 ("wilcox.test", "t.test", "kruskal.test", "anova")
#' @param p.adj p-value 조정 방법 ("BH", "bonferroni", "none" 등)
#' @param split.na.use 'split.by' 컬럼의 NA 값 처리 방법 (FALSE: 제외, TRUE: 별도 그룹으로 취급)
#' @param colors 그룹별로 사용할 색상 팔레트 (선택 사항)
#'
#' @return ggplot 객체와 통계 결과가 포함된 리스트
#'
#' @export
PlotClusterProportions <- function(sobj, 
                                   patient_col, 
                                   group.by, 
                                   split.by, 
                                   test.use = "wilcox.test",
                                   p.adj = "BH",
                                   split.na.use = FALSE,
                                   colors = NULL) {
  
  # 1. 메타데이터 추출
  meta_data <- sobj@meta.data
  
  # 2. 필수 컬럼 확인
  required_cols <- c(patient_col, group.by, split.by)
  if (!all(required_cols %in% colnames(meta_data))) {
    stop("제공된 컬럼명 중 일부가 메타데이터에 없습니다: ", 
         paste(required_cols[!required_cols %in% colnames(meta_data)], collapse = ", "))
  }
  
  # 3. 환자(샘플)별 그룹 정보 추출 (dplyr:: 명시적 호출)
  patient_group_info <- meta_data %>%
    dplyr::select(all_of(c(patient_col, split.by))) %>%
    dplyr::distinct()
  
  # 4. Pseudobulk 비율 계산 (dplyr:: 명시적 호출)
  prop_data <- meta_data %>%
    dplyr::group_by(!!sym(patient_col), !!sym(group.by)) %>%
    dplyr::summarise(n = n(), .groups = "drop") %>%
    dplyr::group_by(!!sym(patient_col)) %>%
    dplyr::mutate(proportion = n / sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(patient_group_info, by = patient_col)
  
  # 5. NA 처리 로직 (dplyr:: 명시적 호출)
  any_na <- any(is.na(prop_data[[split.by]]))
  if (any_na) {
    if (split.na.use == FALSE) {
      message("NA values found in '", split.by, "'. Ignoring samples with NA.")
      prop_data <- prop_data %>% dplyr::filter(!is.na(!!sym(split.by)))
    } else { 
      message("NA values found in '", split.by, "'. Treating as a separate group: 'NA_as_Group'.")
      prop_data <- prop_data %>%
        dplyr::mutate(!!sym(split.by) := as.character(!!sym(split.by))) %>%
        dplyr::mutate(!!sym(split.by) := dplyr::if_else(is.na(!!sym(split.by)), "NA_as_Group", !!sym(split.by))) %>%
        dplyr::mutate(!!sym(split.by) := factor(!!sym(split.by)))
    }
  }
  
  # 5.5. 팩터 레벨 정리
  prop_data <- droplevels(prop_data)
  
  # 5.6. split.by 컬럼을 factor로 변환 (V8 수정)
  prop_data[[split.by]] <- as.factor(prop_data[[split.by]])
  
  # *** 5.7. group.by 컬럼을 factor로 변환 (V10 수정) ***
  # 이 부분이 .subset2 오류의 최종 원인일 가능성이 높습니다.
  prop_data[[group.by]] <- as.factor(prop_data[[group.by]])
  
  # 6. 통계 검정을 위한 데이터 필터링
  n_groups <- length(unique(prop_data[[split.by]]))
  stat_test <- data.frame() # 통계 결과 초기화
  
  if (n_groups < 2) {
    warning("NA 제거 후 비교할 그룹이 2개 미만입니다. 통계 검정을 생략합니다.")
  } else {
    # 7. 비교 가능한 클러스터 탐색 (dplyr:: 명시적 호출)
    groups_per_cluster <- prop_data %>%
      dplyr::group_by(!!sym(group.by)) %>%
      dplyr::summarise(n_groups_in_cluster = dplyr::n_distinct(!!sym(split.by)))
    
    clusters_to_test <- groups_per_cluster %>%
      dplyr::filter(n_groups_in_cluster == n_groups) %>%
      dplyr::pull(!!sym(group.by))
    
    # 8. 통계 검정 수행 (dplyr:: 명시적 호출)
    if (length(clusters_to_test) > 0) {
      
      prop_data_for_stat_test <- prop_data %>%
        dplyr::filter(!!sym(group.by) %in% clusters_to_test)
      
      if (n_groups > 2 && test.use == "t.test") {
        message("그룹이 3개 이상이므로 'anova'를 사용합니다.")
        test.use <- "anova"
      }
      if (n_groups > 2 && test.use == "wilcox.test") {
        message("그룹이 3개 이상이므로 'kruskal.test'를 사용합니다.")
        test.use <- "kruskal.test"
      }
      
      stat_test <- prop_data_for_stat_test %>%
        dplyr::group_by(!!sym(group.by)) %>%
        dplyr::do(
          if(test.use == "wilcox.test") {
            rstatix::wilcox_test(data = ., as.formula(paste("proportion ~", split.by)))
          } else if (test.use == "t.test") {
            rstatix::t_test(data = ., as.formula(paste("proportion ~", split.by)))
          } else if (test.use == "kruskal.test") {
            rstatix::kruskal_test(data = ., as.formula(paste("proportion ~", split.by)))
          } else if (test.use == "anova") {
            rstatix::anova_test(data = ., as.formula(paste("proportion ~", split.by)))
          } else {
            stop("지원하지 않는 test.use입니다.")
          }
        ) %>%
        dplyr::ungroup()
      
      # 9. P-value 조정
      if (nrow(stat_test) > 0 && p.adj != "none") {
        stat_test <- stat_test %>%
          rstatix::adjust_pvalue(method = p.adj) %>%
          rstatix::add_significance(p.col = "p.adj")
        stat_test$p_label_col <- "p.adj.signif"
      } else if (nrow(stat_test) > 0) {
        stat_test <- stat_test %>%
          rstatix::add_significance(p.col = "p")
        stat_test$p_label_col <- "p.signif"
      }
      
      # 10. P-value 브라켓 위치 계산
      if (nrow(stat_test) > 0 && n_groups == 2) {
        stat_test <- stat_test %>%
          rstatix::add_xy_position(
            data = prop_data, 
            formula = as.formula(paste("proportion ~", split.by)), 
            x = group.by, 
            dodge = 0.8
          )
      }
      
    } else {
      warning("NA 제거 후, 모든 그룹 간 비교가 가능한 클러스터가 없습니다. 통계 검정을 생략합니다.")
    }
  }
  
  # 11. 박스플롯 생성
  p <- ggplot2::ggplot(prop_data, ggplot2::aes(x = !!sym(group.by), y = proportion, fill = !!sym(split.by))) +
    ggplot2::geom_boxplot(position = ggplot2::position_dodge(0.8), outlier.shape = NA) +
    ggplot2::geom_point(position = ggplot2::position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), 
                        alpha = 0.7, size = 1.5) +
    ggplot2::labs(
      title = "Cluster Proportion by Group (Pseudobulk)",
      x = group.by,
      y = "Proportion per Sample",
      fill = split.by
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    ) +
    ggplot2::scale_y_continuous(labels = scales::percent)
  
  if (!is.null(colors)) {
    p <- p + ggplot2::scale_fill_manual(values = colors)
  }
  
  # 12. P-value 브라켓 추가
  if (nrow(stat_test) > 0) {
    if (n_groups > 2) {
      p <- p + ggpubr::stat_pvalue_manual(
        data = stat_test,
        aes(x = !!sym(group.by)), 
        label = "{p_label_col}",
        y.position = max(prop_data$proportion, na.rm = TRUE) * 1.05,
        hide.ns = TRUE,
        inherit.aes = FALSE 
      )
    } else {
      p <- p + ggpubr::stat_pvalue_manual(
        stat_test, 
        label = "{p_label_col}",
        tip.length = 0.01,
        hide.ns = TRUE
      )
    }
  }
  
  # 13. 반환
  return(list(plot = p, statistics = stat_test))
}

# -------- PTML: utility. --------
#' FindAllMarkers 결과에서 Top N 유전자를 LLM용으로 출력
#'
#' @param markers_df FindAllMarkers()의 결과인 데이터프레임
#' @param top_n 각 클러스터별로 추출할 상위 유전자 개수
#' @param cluster_col 'cluster' 정보가 있는 컬럼명 (기본값: "cluster")
#' @param gene_col 'gene' 정보가 있는 컬럼명 (기본값: "gene")
#' @param p_val_col p-value 컬럼명 (기본값: "p_val")
#' @param p_val_adj_col p-value (adjusted) 컬럼명 (기본값: "p_val_adj")
#' @param avg_log2FC_col log-fold change 컬럼명 (기본값: "avg_log2FC")
#'
#' @return 콘솔에 클러스터별 Top N 유전자 목록을 출력합니다.
#'
#' @export
PrintTopMarkersForLLM <- function(markers_df,
                                  top_n = 50,
                                  cluster_col = "cluster",
                                  gene_col = "gene",
                                  p_val_col = "p_val",
                                  p_val_adj_col = "p_val_adj",
                                  avg_log2FC_col = "avg_log2FC") {
  
  # 1. 필수 컬럼 존재 여부 확인
  required_cols <- c(cluster_col, gene_col, avg_log2FC_col)
  
  if (!all(required_cols %in% colnames(markers_df))) {
    missing_cols <- required_cols[!required_cols %in% colnames(markers_df)]
    stop("필수 컬럼이 markers_df에 없습니다: ", paste(missing_cols, collapse = ", "),
         "\n'cluster_col', 'gene_col', 'avg_log2FC_col' 인자를 확인하세요.")
  }
  
  # 2. 정렬을 위한 데이터 전처리 (p-value 0 처리)
  markers_processed <- markers_df
  
  use_pval <- p_val_col %in% colnames(markers_processed)
  use_padj <- p_val_adj_col %in% colnames(markers_processed)
  
  # p_val_mod 컬럼을 안전하게 추가
  if (use_pval) {
    markers_processed$p_val_mod <- ifelse(markers_processed[[p_val_col]] == 0, 
                                          1e-300, # 0 대신 매우 작은 값
                                          markers_processed[[p_val_col]])
  }
  
  # 3. 정렬 기준에 따라 Top N 마커 필터링
  markers_grouped <- markers_processed %>%
    group_by(!!sym(cluster_col))
  
  # 정렬 로직
  if (use_padj && use_pval) {
    message(paste("Sorting by", p_val_adj_col, ",", p_val_col, ", and", avg_log2FC_col))
    sorted_markers <- markers_grouped %>%
      arrange(!!sym(p_val_adj_col), p_val_mod, desc(!!sym(avg_log2FC_col)), .by_group = TRUE)
    
  } else if (use_pval) {
    message(paste("Sorting by", p_val_col, "and", avg_log2FC_col))
    sorted_markers <- markers_grouped %>%
      arrange(p_val_mod, desc(!!sym(avg_log2FC_col)), .by_group = TRUE)
    
  } else {
    message(paste("Sorting by", avg_log2FC_col, "only."))
    sorted_markers <- markers_grouped %>%
      arrange(desc(!!sym(avg_log2FC_col)), .by_group = TRUE)
  }
  
  final_markers <- sorted_markers %>%
    slice_head(n = top_n) %>%
    ungroup()
  
  # 4. LLM에 복사/붙여넣기 편한 형태로 콘솔에 출력
  
  # 클러스터 목록 (숫자형/문자형 모두 정렬되도록)
  clusters <- unique(final_markers[[cluster_col]])
  tryCatch({
    clusters_sorted <- sort(as.numeric(as.character(clusters)))
  }, warning = function(w) {
    clusters_sorted <<- sort(as.character(clusters))
  })
  
  
  cat("--- Top", top_n, "Markers per Cluster ---\n\n")
  
  for (cl in clusters_sorted) {
    cat("--- Cluster:", cl, "---\n")
    
    genes <- final_markers %>%
      filter(!!sym(cluster_col) == cl) %>%
      pull(!!sym(gene_col))
    
    # 쉼표로 구분된 한 줄의 문자열로 출력
    cat(paste(genes, collapse = ", "))
    cat("\n\n") # 클러스터 간 구분을 위한 공백
  }
  
  message("--- End of Marker List ---")
}

# --------- analyses: LISI, PERMANOVA ---------

#' LISI 점수에 대해 의사 반복을 피해 통계 검정 수행
#'
#' @param lisi_scores LISI 점수 (compute_lisi의 결과. Rowname이 cell barcode)
#' @param metadata Seurat 객체의 메타데이터 (is@meta.data)
#' @param lisi_col_name LISI 점수 컬럼명 (e.g., "GEM", "batch")
#' @param patient_var 환자/샘플 ID가 있는 메타데이터 컬럼명 (e.g., "patient_ID")
#' @param celltype_var 세포 유형 정보가 있는 메타데이터 컬럼명 (e.g., "cell_type")
#' @param biological_group_var 비교하려는 그룹 (e.g., "disease_status", "set")
#'
#' @return cell_type별 통계 검정 결과 (adj_p_value로 정렬됨)
#' @export
TestLISI = function(lisi_scores, 
                    metadata, 
                    lisi_col_name, 
                    patient_var, 
                    celltype_var, 
                    biological_group_var) {
  
  # 1. 데이터 준비 (메타데이터와 LISI 점수 결합)
  
  # LISI 점수 컬럼이 있는지 확인
  if (!lisi_col_name %in% colnames(lisi_scores)) {
    stop(paste0("Error: '", lisi_col_name, "' 컬럼이 lisi_scores에 없습니다."))
  }
  
  # LISI 점수 컬럼을 'LISI_score'라는 공통 이름으로 변경
  lisi_scores_renamed <- lisi_scores %>%
    select(LISI_score = !!sym(lisi_col_name))
  
  # 메타데이터와 LISI 점수 결합 (rownames 기준)
  meta_data_with_lisi <- metadata %>%
    rownames_to_column("cell_barcode_temp") %>%
    left_join(
      lisi_scores_renamed %>% rownames_to_column("cell_barcode_temp"),
      by = "cell_barcode_temp"
    ) %>%
    select(-cell_barcode_temp) # 임시 바코드 컬럼 제거

  # 2. 세포 유형별 루프 시작
  all_cell_types <- unique(meta_data_with_lisi[[celltype_var]])
  
  message(paste("Starting LISI statistical test for", length(all_cell_types), "cell types..."))
  
  stats_results <- purrr::map_dfr(all_cell_types, function(current_cell_type) {
    
    # 3. 환자 레벨로 LISI 점수 집계 (Aggregation)
    # (핵심) 각 환자별 LISI 점수의 중앙값을 계산
    median_lisi_per_patient <- meta_data_with_lisi %>%
      filter(!!sym(celltype_var) == current_cell_type) %>%
      filter(!is.na(!!sym(biological_group_var))) %>% # 비교 그룹이 NA인 셀 제외
      group_by(!!sym(patient_var), !!sym(biological_group_var)) %>%
      summarise(
        median_LISI = median(LISI_score, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      # 한 환자가 여러 메타데이터 행을 가질 경우(e.g. 여러 배치) 
      # group_by/summarise가 이미 환자당 1줄로 만들지만,
      # 만일을 대비해 환자 ID로 중복 제거
      distinct(!!sym(patient_var), .keep_all = TRUE) 
      
    # 4. 통계 검정을 위한 데이터 유효성 검사
    if (nrow(median_lisi_per_patient) < 3 || 
        length(unique(median_lisi_per_patient[[biological_group_var]])) < 2) {
      
      message(paste("  [Skipping]", current_cell_type, "| Not enough samples or < 2 groups"))
      return(tibble(
        cell_type = current_cell_type,
        statistic = NA,
        p_value = NA,
        error = "Not enough samples or only one group"
      ))
    }
    
    # 5. Wilcoxon Rank-Sum Test 수행
    test_formula <- as.formula(paste("median_LISI ~", biological_group_var))
    
    tryCatch({
      wilcox_res <- wilcox.test(test_formula, data = median_lisi_per_patient)
      
      return(tibble(
        cell_type = current_cell_type,
        statistic_W = wilcox_res$statistic,
        p_value = wilcox_res$p.value,
        error = NA
      ))
    }, error = function(e) {
      # (e.g., 한 그룹의 모든 값이 동일할 때 오류 발생 가능)
      message(paste("  [Error]", current_cell_type, ":", e$message))
      return(tibble(
        cell_type = current_cell_type,
        statistic_W = NA,
        p_value = NA,
        error = e$message
      ))
    })
  })
  
  # 6. P-value 보정 (BH) 및 정렬
  stats_results <- stats_results %>%
    mutate(adj_p_value = p.adjust(p_value, method = "BH")) %>%
    arrange(adj_p_value)
    
  message("LISI statistical test finished.")
  print(stats_results)
  return(stats_results)
}

#' @export
PERMANOVA=function(sobj, patient_var, celltype_var, group_var, assay="RNA", na.drop=TRUE){
  # 환자 ID가 들어있는 컬럼명 (예: "p1", "p2"...)
  # patient_var
  # 세포 유형 정보가 있는 컬럼명
  # celltype_var <- "seurat_clusters"
  # 테스트하려는 생물학적 그룹 컬럼명 (예: "set1", "set2", "set3")
  # group_var <- "set"
  # 사용할 Assay
  # assay <- "SCT"
  # na.drop: group_var에 NA가 있는 셀을 제거할지 여부
  
  # --- 분석 시작 ---
  
  # 1. 분석할 모든 세포 유형 목록 가져오기
  all_cell_types <- unique(sobj@meta.data[[celltype_var]])
  
  # 2. 각 세포 유형에 대해 반복적으로 PERMANOVA 수행
  # map_dfr은 결과를 data.frame으로 깔끔하게 묶어줍니다.
  pb_permanova_results <- purrr::map_dfr(all_cell_types, function(current_cell_type) {
    message(paste("Processing:", current_cell_type))
    # 3. 현재 세포 유형으로만 Seurat 객체 필터링
    seurat_subset <- subset(sobj, !!sym(celltype_var) == current_cell_type)
    
    # --- [수정] group_var에 NA가 있으면 해당 셀 제거 (na.drop=TRUE일 경우) ---
    if (na.drop) {
      seurat_subset <- subset(seurat_subset, !is.na(!!sym(group_var)))
    }
    
    # --- [추가] 필터링 후 셀이 남아있는지 확인 ---
    if (nrow(seurat_subset@meta.data) == 0) {
      message(paste("  [Skipping]", current_cell_type, "| No cells remaining after NA drop"))
      return(data.frame(
        cell_type = current_cell_type,
        R2 = NA,
        p_value = NA,
        error = "No cells remaining after NA drop"
      ))
    }
    # --- [추가 끝] ---
    
    # 4. 환자(patient_var)별로 Pseudo-bulk 수행
    pb_subset <- AggregateExpression(
      seurat_subset,
      group.by = patient_var,
      assays = assay,
      slot = "data" # SCT/RNA의 정규화된 데이터를 사용
    )
    
    # 5. 매트릭스 준비 (샘플 x 유전자)
    pb_matrix_t <- t(pb_subset[[assay]])
    
    # 6. Pseudo-bulk 샘플(환자)에 대한 메타데이터 생성
    # 환자 ID와 그 환자가 속한 'set' 그룹을 매칭
    meta_pb_subset <- seurat_subset@meta.data %>%
      select(!!sym(patient_var), !!sym(group_var)) %>%
      filter(!is.na(!!sym(patient_var))) %>% #NA 값을 가진 행을 미리 제거합니다.
      distinct() # 환자-그룹 매핑 정보만 남김

    # 매트릭스 row 순서와 메타데이터 순서 맞추기
    # AggregateExpression은 숫자로 시작하는 ID 앞에 'g'를 붙일 수 있으므로,
    # join_key를 만들어 두 경우 모두 처리합니다. (숫자이거나, 캐릭터이거나.)
    
    # --- [수정] dplyr의 NSE(Non-Standard Evaluation) 오류를 피하기 위해 base R로 변경 ---
    patient_ids <- as.character(meta_pb_subset[[patient_var]])
    needs_g_prefix <- grepl("^[0-9]", patient_ids)
    join_keys <- ifelse(needs_g_prefix, paste0("g", patient_ids), patient_ids)
    
    meta_pb_subset$join_key <- join_keys
    # column_to_rownames를 사용하기 전에 기존의 행 이름(cell barcodes)을 제거합니다.
    rownames(meta_pb_subset) <- NULL
    meta_pb_subset <- tibble::column_to_rownames(meta_pb_subset, "join_key")
    
    # AggregateExpression이 생성한 실제 샘플 이름으로 메타데이터 순서를 맞춥니다.
    meta_pb_subset <- meta_pb_subset[rownames(pb_matrix_t), , drop = FALSE]


    # # 매트릭스 row 순서와 메타데이터 순서 맞추기
    # rownames(meta_pb_subset) <- paste0("g",meta_pb_subset[[patient_var]])
    # meta_pb_subset <- meta_pb_subset[rownames(pb_matrix_t), ]
    
    # (오류 방지) 해당 세포 유형에 샘플이 너무 적거나 그룹이 1개면 스킵
    if (nrow(meta_pb_subset) < 3 || length(unique(meta_pb_subset[[group_var]])) < 2) {
      # --- [확인용 코드 추가] ---
      message(paste(" [Skipping]", current_cell_type,
                    "| Samples:", nrow(meta_pb_subset),
                    "| Unique Groups:", length(unique(meta_pb_subset[[group_var]]))))
      # --- [추가 끝] ---
      return(data.frame(
        cell_type = current_cell_type,
        R2 = NA,
        p_value = NA,
        error = "Not enough samples or only one group"
      ))
    }
    
    # 7. 유클리드 거리 계산
    dist_pb <- dist(pb_matrix_t, method = "euclidean")
    
    # 8. PERMANOVA 수행 (formula: 거리 ~ 그룹)
    set_formula <- as.formula(paste("dist_pb ~", group_var))
    
    res <- adonis2(
      set_formula,
      data = meta_pb_subset,
      permutations = 999)
    
    # 9. 결과 추출
    data.frame(
      cell_type = current_cell_type,
      R2 = res$R2[1], # 'set' 그룹이 설명하는 분산의 양
      p_value = res$`Pr(>F)`[1], # 'set' 그룹 간 차이의 p-value
      error = NA
    )}
  )
  
  # --- 10. P-value 보정 (Benjamini-Hochberg) ---
  pb_permanova_results$adj_p_value <- p.adjust(pb_permanova_results$p_value, method = "BH")
  
  # --- 11. adj_p_value 기준으로 정렬 ---
  pb_permanova_results <- pb_permanova_results %>%
    dplyr::arrange(adj_p_value)
  
  # --- 최종 결과 확인 ---
  print(pb_permanova_results)
  
  return(pb_permanova_results)
}


# --------- analyses: NEBULA, MUSCAT, CLMM, Spearman ---------

#' Seurat 객체에서 유전자별 CLMM 분석 수행 (v2: 기능 추가)
#'
#' @param sobj Seurat 객체
#' @param ordinal_var 메타데이터의 순서형 반응 변수 (예: "heal")
#' @param patient_col 메타데이터의 환자 ID(random effect) (예: "emrid")
#' @param assay 사용할 Seurat assay (기본값: "RNA")
#' @param layer 사용할 Seurat layer (v4/v3의 'slot') (기본값: "data")
#' @param top_n 분석할 상위 유전자 개수 (기본값: NULL = 모든 유전자)
#' @param sort_by 'top_n'을 사용할 경우 정렬 기준 ("variance" 또는 "mean")
#' @param decreasing 'top_n' 정렬 시 내림차순 (기본값: TRUE)
#' @param ... clmm() 함수에 전달할 추가 인자 (예: control = clmm.control(maxIter = 200))
#'
#' @return 'gene' 및 CLMM 결과 (Estimate, p-value 등) data.frame
#'
run_clmm_analysis <- function(sobj, ordinal_var, patient_col, assay = "RNA", layer = "data", 
                              top_n = NULL, sort_by = "variance", decreasing = TRUE, ...) {
  
  # 0. 필수 패키지 확인
  if (!requireNamespace("ordinal", quietly = TRUE)) {
    stop("'ordinal' 패키지가 필요합니다: install.packages('ordinal')")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("'dplyr' 패키지가 필요합니다: install.packages('dplyr')")
  }
  
  # 1. 메타데이터 준비
  meta <- sobj@meta.data
  unique_levels <- sort(unique(meta[[ordinal_var]]))
  meta[[ordinal_var]] <- factor(meta[[ordinal_var]], levels = unique_levels, ordered = TRUE)
  cat("Note: '", ordinal_var, "'를 ordered factor로 변환했습니다: ", paste(unique_levels, collapse = "<"), "\n")
  meta[[patient_col]] <- as.factor(meta[[patient_col]])
  
  # 2. 유전자 발현 데이터 및 유전자 목록 준비
  expr_matrix <- GetAssayData(sobj, assay = assay, slot = layer)
  
  if (!is.null(top_n)) {
    cat(paste("'top_n'이 설정되었습니다.", top_n, "개 유전자를", sort_by, "기준으로 정렬합니다...\n"))
    
    if (sort_by == "variance") {
      metrics <- apply(expr_matrix, 1, var)
    } else if (sort_by == "mean") {
      metrics <- rowMeans(expr_matrix)
    } else {
      stop("'sort_by'는 'variance' 또는 'mean'이어야 합니다.")
    }
    
    ordered_indices <- order(metrics, decreasing = decreasing)
    genes_to_test <- rownames(expr_matrix)[ordered_indices][1:min(top_n, length(metrics))]
    
  } else {
    genes_to_test <- rownames(expr_matrix)
  }
  
  total_genes <- length(genes_to_test)
  cat("총", total_genes, "개의 유전자에 대해 CLMM 모델을 실행합니다...\n")
  
  # 3. clmm 모델 실행 (for loop로 변경: 진행 상황 및 시간 측정)
  clmm_results_list <- list()
  overall_start_time <- Sys.time()
  
  for (i in seq_along(genes_to_test)) {
    gene <- genes_to_test[i]
    
    model_data <- data.frame(
      ord_var = meta[[ordinal_var]],
      pat_var = meta[[patient_col]],
      gene_expr = as.numeric(expr_matrix[gene, ])
    )
    
    # 모델 피팅 (오류 발생 시 NULL 반환)
    model_fit <- tryCatch({
      # '...' 인자를 통해 clmm.control(maxIter=200) 등 전달 가능
      ordinal::clmm(ord_var ~ gene_expr + (1|pat_var), data = model_data, ...)
    }, error = function(e) {
      NULL
    })
    
    # 결과 추출
    if (!is.null(model_fit)) {
      coef_summary <- summary(model_fit)$coefficients
      if ("gene_expr" %in% rownames(coef_summary)) {
        clmm_results_list[[gene]] <- coef_summary["gene_expr", ]
      } else {
        clmm_results_list[[gene]] <- rep(NA, 4)
      }
    } else {
      clmm_results_list[[gene]] <- rep(NA, 4)
    }
    
    # 4. 진행 상황 및 예상 시간 리포트
    if (i == 100) {
      first_100_time <- Sys.time()
      time_per_100 <- as.numeric(difftime(first_100_time, overall_start_time, units = "secs"))
      total_estimated_secs <- (time_per_100 / 100) * total_genes
      total_estimated_mins <- total_estimated_secs / 60
      
      cat(paste0("[진행 상황] 첫 100개 유전자 완료 (", round(time_per_100, 1), "초 소요)\n"))
      cat(paste0("   > 총 예상 소요 시간: 약 ", round(total_estimated_mins, 1), " 분\n"))
      
    } else if (i > 100 && i %% 100 == 0) {
      current_time <- Sys.time()
      elapsed_mins <- as.numeric(difftime(current_time, overall_start_time, units = "mins"))
      cat(paste0("[진행 상황] ", i, " / ", total_genes, " (", round(i/total_genes*100), "%) 완료. (", round(elapsed_mins, 1), "분 경과)\n"))
    }
  } # end for loop
  
  # 5. 결과 data.frame으로 변환
  clmm_results_df <- do.call(rbind, clmm_results_list)
  rownames(clmm_results_df) <- names(clmm_results_list) # genes_to_test 순서
  colnames(clmm_results_df) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  
  clmm_results_df <- as.data.frame(clmm_results_df)
  clmm_results_df$gene <- rownames(clmm_results_df)
  
  clmm_results_df <- clmm_results_df[!is.na(clmm_results_df$`Pr(>|z|)`), ]
  clmm_results_df <- clmm_results_df %>% 
    dplyr::arrange(`Pr(>|z|)`) %>%
    dplyr::select(gene, dplyr::everything())
  
  cat("CLMM 분석 완료.\n")
  return(clmm_results_df)
}

#' Seurat 객체에서 유전자별 스피어만 상관관계 분석 (Pseudobulk) 수행 (v2: 버그 수정 및 기능 추가)
#'
#' @param sobj Seurat 객체
#' @param ordinal_var 메타데이터의 순서형 변수 (예: "heal")
#' @param patient_col 메타데이터의 환자 ID (예: "emrid")
#' @param assay 사용할 Seurat assay (기본값: "RNA")
#' @param layer 사용할 Seurat layer (v4/v3의 'slot') (기본값: "data")
#' @param top_n 분석할 상위 유전자 개수 (기본값: NULL = 모든 유전자)
#' @param sort_by 'top_n'을 사용할 경우 정렬 기준 ("variance" 또는 "mean")
#' @param decreasing 'top_n' 정렬 시 내림차순 (기본값: TRUE)
#'
#' @return 'gene', 'rho', 'p.value'를 포함하는 data.frame
#'
run_spearman_pseudobulk <- function(sobj, ordinal_var, patient_col, assay = "RNA", layer = "data",
                                    top_n = NULL, sort_by = "variance", decreasing = TRUE) {
  
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("'dplyr' 패키지가 필요합니다: install.packages('dplyr')")
  }

  # 1. (오류 수정) 숫자로 시작하는 ID에 'g'를 덧붙인 임시 컬럼 생성
  #    AggregateExpression이 내부적으로 하는 작업을 미리 수행하여 ID 불일치 방지
  temp_group_by_col <- "temp_spearman_group_by"
  sobj@meta.data[[temp_group_by_col]] <- ifelse(
    grepl("^[0-9]", sobj@meta.data[[patient_col]]),
    paste0("g", sobj@meta.data[[patient_col]]),
    as.character(sobj@meta.data[[patient_col]])
  )

  # 2. Pseudobulk 발현 매트릭스 생성
  cat("환자별 Pseudobulk 프로필을 생성합니다 (임시 ID 사용)...\n")
  avg_expr <- AggregateExpression(
    sobj,
    group.by = temp_group_by_col, # 수정된 임시 컬럼 사용
    assays = assay,
    slot = layer, 
    return.seurat = FALSE
  )
  
  pseudobulk_matrix <- avg_expr[[assay]]
  
  # (참고) AggregateExpression의 'Inf' 경고 메시지:
  # 'layer = "data"' (보통 log-normalized)를 사용할 때 발생할 수 있습니다.
  # 이는 'AggregateExpression'이 내부적으로 'exp(x)-1'을 시도할 때 발생할 수 있으나,
  # 'slot = "data"'의 경우 평균(mean)을 계산하므로 결과(평균 발현 값)는 유효합니다.

  # 3. 환자별 메타데이터(ordinal_var) 준비
  meta <- sobj@meta.data
  
  # 임시 ID를 기준으로 유니크한 메타데이터 추출
  patient_meta <- meta %>%
    dplyr::select(dplyr::all_of(c(temp_group_by_col, ordinal_var))) %>%
    dplyr::distinct(!!dplyr::sym(temp_group_by_col), .keep_all = TRUE)
  
  # pseudobulk_matrix의 컬럼 이름(환자 ID) 순서대로 patient_meta 정렬
  pseudobulk_patients <- colnames(pseudobulk_matrix)
  match_indices <- match(pseudobulk_patients, patient_meta[[temp_group_by_col]])
  
  # 오류 검사 (NA가 없는지)
  if (any(is.na(match_indices))) {
      stop(paste("Patient ID 매칭 실패. 임시 ID 생성 후에도 NA가 발생:",
                 paste(pseudobulk_patients[is.na(match_indices)], collapse = ",")))
  }
  
  ordered_patient_meta <- patient_meta[match_indices, ]
  
  # (오류 수정) 원본 코드의 if문 검사를 더 강력하게 수정
  if (!all(pseudobulk_patients == ordered_patient_meta[[temp_group_by_col]])) {
    stop("환자 ID 정렬 오류: Pseudobulk 순서와 메타데이터 순서가 일치하지 않습니다.")
  }
  
  heal_values <- ordered_patient_meta[[ordinal_var]]
  
  # 4. 유전자 필터링 (Pseudobulk Matrix 기준)
  if (!is.null(top_n)) {
    cat(paste("'top_n'이 설정되었습니다.", top_n, "개 유전자를", sort_by, "기준으로 정렬합니다 (Pseudobulk 기준)...\n"))
    
    if (sort_by == "variance") {
      metrics <- apply(pseudobulk_matrix, 1, var)
    } else if (sort_by == "mean") {
      metrics <- rowMeans(pseudobulk_matrix)
    } else {
      stop("'sort_by'는 'variance' 또는 'mean'이어야 합니다.")
    }
    
    ordered_indices <- order(metrics, decreasing = decreasing)
    genes_to_test <- rownames(pseudobulk_matrix)[ordered_indices][1:min(top_n, length(metrics))]
    
    # 분석할 매트릭스 서브셋
    pseudobulk_matrix_subset <- pseudobulk_matrix[genes_to_test, ]
    
  } else {
    genes_to_test <- rownames(pseudobulk_matrix)
    pseudobulk_matrix_subset <- pseudobulk_matrix
  }

  # 5. 스피어만 상관관계 계산 (for loop 사용)
  total_genes <- length(genes_to_test)
  cat("총", total_genes, "개의 유전자에 대해 스피어만 상관관계를 계산합니다...\n")
  
  spearman_results_list <- list()
  overall_start_time <- Sys.time()
  
  for (i in seq_along(genes_to_test)) {
    gene <- genes_to_test[i]
    gene_avg_expr <- pseudobulk_matrix_subset[gene, ]
    
    # 발현에 분산이 있는지 확인
    if (var(gene_avg_expr, na.rm=TRUE) == 0) {
      spearman_results_list[[gene]] <- c(rho = NA, p.value = NA)
      next # 다음 유전자로
    }
    
    cor_test_result <- tryCatch({
      stats::cor.test(gene_avg_expr, heal_values, method = "spearman")
    }, error = function(e) { NULL })
    
    if (!is.null(cor_test_result)) {
      spearman_results_list[[gene]] <- c(
        rho = cor_test_result$estimate,
        p.value = cor_test_result$p.value
      )
    } else {
      spearman_results_list[[gene]] <- c(rho = NA, p.value = NA)
    }
    
    # 진행 상황 리포트
    if (i == 100) {
      first_100_time <- Sys.time()
      time_per_100 <- as.numeric(difftime(first_100_time, overall_start_time, units = "secs"))
      total_estimated_secs <- (time_per_100 / 100) * total_genes
      total_estimated_mins <- total_estimated_secs / 60
      
      cat(paste0("[진행 상황] 첫 100개 유전자 완료 (", round(time_per_100, 1), "초 소요)\n"))
      cat(paste0("   > 총 예상 소요 시간: 약 ", round(total_estimated_mins, 1), " 분\n"))
      
    } else if (i > 100 && i %% 100 == 0) {
      current_time <- Sys.time()
      elapsed_mins <- as.numeric(difftime(current_time, overall_start_time, units = "mins"))
      cat(paste0("[진행 상황] ", i, " / ", total_genes, " (", round(i/total_genes*100), "%) 완료. (", round(elapsed_mins, 1), "분 경과)\n"))
    }
  } # end for loop

  # 6. 결과 data.frame으로 변환
  spearman_results_df <- do.call(rbind, spearman_results_list)
  spearman_results_df <- as.data.frame(spearman_results_df)
  spearman_results_df$gene <- rownames(spearman_results_df)
  
  spearman_results_df <- spearman_results_df[!is.na(spearman_results_df$p.value), ]
  spearman_results_df <- spearman_results_df %>%
    dplyr::arrange(p.value) %>%
    dplyr::select(gene, rho, p.value)
  
  cat("스피어만 상관관계 분석 완료.\n")
  return(spearman_results_df)
}

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
#' @title NEBULA 파이프라인 함수
#' @description Seurat 객체를 받아 NA 처리, 데이터 정렬, NEBULA GLMM을 실행.
#'              Modeling에서, y~NB(mu,phi)이며, log(mu)=log(o)+BX+random effect 에서 "o"는 offset, random effect는 patient_col, X는 fixed effects와 covar_effects.
#'
#' @param sobj (Seurat) Seurat 객체
#' @param layer (character) Raw count가 저장된 layer (예: "counts")
#' @param fixed_effects (character vector) 고정 효과로 사용할 meta.data 컬럼명
#'                    (예: c("target_col", "celltype_col"))
#' @param covar_effects (character vector) 임의 효과처럼 보정하고 싶은 공변량.
#'                    NEBULA의 한계로 인해 '고정 효과'로 처리됩니다.
#'                    (예: c("batch_col", "percent.mt"))
#' @param patient_col (character) 환자 ID 컬럼. NEBULA의 'id' (유일한 임의 효과)로 사용됨.
#' @param offset (character) Offset으로 사용할 'numeric' 컬럼 (예: "nCount_RNA")
#' @param min_count (numeric) 유전자 필터링 기준 (최소 발현 세포 수)
#'
#' @return (nebula) NEBULA 실행 결과 객체
#'
#' @export 
runNEBULA <- function(sobj,
                                layer = "counts",
                                fixed_effects = c("target_col", "celltype_col"),
                                covar_effects = c("batch_col"),
                                patient_col = "patient_col",
                                offset = "nCount_RNA",
                                min_count = 10) {
  
  # --- 0. 데이터 추출 ---
  meta <- sobj@meta.data
  counts <- GetAssayData(sobj, layer = layer) # dgCMatrix (희소 행렬)
  
  # --- 1. 유전자 필터링 ---
  message(sprintf("1/6: 유전자 필터링 (min %d cells)...", min_count))
  keep_genes <- rowSums(counts > 0) >= min_count
  counts_filtered <- counts[keep_genes, ]
  message(sprintf("... %d / %d 유전자 통과", sum(keep_genes), nrow(counts)))
  
  # --- 2. NA 값 확인 및 제거 ---
  # [수정] NEBULA는 'covar_effects'도 고정 효과로 처리해야 함
  all_fixed_vars <- c(fixed_effects, covar_effects)
  
  # [수정] NA 체크할 모든 변수 목록
  vars_to_check <- c(all_fixed_vars, patient_col, offset)
  
  message("2/6: 모델 변수에서 NA 값 확인 중...")
  message(paste("... 확인 대상:", paste(vars_to_check, collapse = ", ")))
  
  keep_cells_idx <- complete.cases(meta[, vars_to_check])
  n_removed <- sum(!keep_cells_idx)
  
  if (n_removed > 0) {
    message(sprintf("... NA 값으로 인해 %d 개의 세포를 제거합니다.", n_removed))
  }
  
  # NA가 없는 '깨끗한' 데이터 생성
  meta_clean <- meta[keep_cells_idx, ]
  counts_clean <- counts_filtered[, keep_cells_idx]
  
  message(sprintf("... 최종 분석 대상 세포: %d 개", nrow(meta_clean)))
  
  # --- 3. 디자인 행렬 및 벡터 생성 ---
  message("3/6: 디자인 행렬 및 벡터 생성 중...")
  
  # [수정] 범주형 변수와 숫자형 변수 분리
  # 'offset'은 숫자로 유지해야 함!
  factor_vars <- c(all_fixed_vars, patient_col)
  
  # [수정] lapply + as.factor 적용 (오타 수정 as.factors -> as.factor)
  meta_clean[factor_vars] <- lapply(meta_clean[factor_vars], as.factor)
  
  # [수정] 'offset'은 numeric으로 강제 변환 (안전장치)
  meta_clean[[offset]] <- as.numeric(meta_clean[[offset]])
  if (any(is.na(meta_clean[[offset]]))) {
      stop(sprintf("'%s' 컬럼에 NA가 아닌 숫자형 값만 있어야 합니다.", offset))
  }

  # [수정] 'formula' 인수를 없애고, 변수명으로부터 동적 생성
  formula_str <- paste("~", paste(all_fixed_vars, collapse = " + "))
  message(sprintf("... 사용할 고정 효과 포뮬러: %s", formula_str))
  
  design_matrix <- model.matrix(as.formula(formula_str), data = meta_clean)
  
  # id 및 offset 벡터 추출
  id_vector <- meta_clean[[patient_col]]
  offset_vector <- meta_clean[[offset]]
  
  # --- 4. group_cell()로 데이터 정렬 ---
  message(sprintf("4/6: NEBULA를 위해 id (%s) 기준으로 데이터 정렬 중...", patient_col))
  data_grouped <- nebula::group_cell(
    count = counts_clean,
    id = id_vector,
    pred = design_matrix,
    offset = offset_vector
  )
  
  # --- 5. NEBULA 실행 ---
  message("5/6: NEBULA 실행 중 (NBLMM)...")
  re_nebula <- nebula::nebula(
    count = data_grouped$count,
    id = data_grouped$id,
    pred = data_grouped$pred,
    offset = data_grouped$offset,
    model = "NBLMM" # lme4::glmer.nb와 유사한 Lognormal random effect 권장
  )
  
  # --- 6. 결과 반환 ---
  message("6/6: 분석 완료.")
  
  # [수정] 함수는 결과를 'return' 해야 함
  return(re_nebula)
}

#' @title Muscat (Pseudo-bulking) 파이프라인 함수 (V5 - prepData 사용)
#'
#' @description prepData를 사용하여 SCE 객체를 포맷하고, pbDS로 GLM(edgeR/DESeq2) 실행
#' @note 이 함수는 'pbDS'의 한계로 인해 혼합 효과(Random Effects/block)를 지원하지 않습니다.
#'
#' @param sobj (Seurat) Seurat 객체
#' @param cluster_id (character) 세포 유형 컬럼명
#' @param sample_id (character) 환자/샘플 ID 컬럼명
#' @param group_id (character) 비교할 그룹 컬럼명
#' @param formula_str (character) 포뮬러 (예: "~ 0 + group + set")
#'                 'group' 키워드는 'group_id'로 자동 변환됩니다.
#' @param contrast (character) Contrast (예: "groupA-groupB")
#' @param method (character) "edgeR", "DESeq2", "limma-trend", "limma-voom" (block 없음)
#'
#' @export
runMUSCAT <- function(...) {
  .Deprecated("runMUSCAT_v5", package = "myR", msg = "runMUSCAT in test_working.R is deprecated. Use runMUSCAT_v5 or runMUSCAT2_v1 from test_analysis.R instead.")
  if (exists("runMUSCAT_v5", envir = asNamespace("myR"), inherits = FALSE)) {
    fun <- get("runMUSCAT_v5", envir = asNamespace("myR"))
    return(fun(...))
  } else {
    stop("runMUSCAT_v5 from test_analysis.R not found. Please ensure test_analysis.R is loaded.")
  }
}

# Original runMUSCAT function body removed - use runMUSCAT_v5 from test_analysis.R instead
# The following code block is kept for reference but will not execute:
if (FALSE) {
  runMUSCAT_original <- function(
  sobj,
  cluster_id = "seurat_clusters",
  sample_id  = "hos_no",
  group_id   = "type",
  batch_id   = NULL,                 # ex) "exp_batch"
  contrast   = NULL,                 # ex) "IS - SAH"
  method     = "edgeR",
  pb_min_cells = 3,
  filter_genes = c("none","genes","both","edgeR"),
  keep_clusters = NULL,
  cluster_label_map = NULL
){
  if (is.null(contrast)) stop("'contrast'를 지정하세요. 예: 'IS - SAH'")
  filter_genes <- match.arg(filter_genes)

  # deps
  req <- c("Seurat","muscat","SingleCellExperiment","SummarizedExperiment","S4Vectors","limma","dplyr")
  miss <- req[!vapply(req, requireNamespace, logical(1), quietly=TRUE)]
  if (length(miss)) stop("필요 패키지 설치: ", paste(miss, collapse=", "))

  # 1) Seurat -> SCE, prepSCE
  sce <- Seurat::as.SingleCellExperiment(sobj)
  sce <- muscat::prepSCE(sce, kid = cluster_id, sid = sample_id, gid = group_id)

  # factor 보장
  sce$cluster_id <- droplevels(factor(SummarizedExperiment::colData(sce)$cluster_id))
  sce$sample_id  <- droplevels(factor(SummarizedExperiment::colData(sce)$sample_id))
  sce$group_id   <- droplevels(factor(SummarizedExperiment::colData(sce)$group_id))
  if (!is.null(batch_id) && batch_id %in% colnames(SummarizedExperiment::colData(sce))) {
    sce[[batch_id]] <- droplevels(factor(SummarizedExperiment::colData(sce)[[batch_id]]))
  }

  # 2) Pseudobulk
  pb <- muscat::aggregateData(sce, assay = "counts", by = c("cluster_id","sample_id"))

  # (선택) 특정 클러스터만
  if (!is.null(keep_clusters)) {
    keep_clusters <- as.character(keep_clusters)
    pb <- pb[names(SummarizedExperiment::assays(pb)) %in% keep_clusters]
    if (length(SummarizedExperiment::assays(pb)) == 0L) stop("keep_clusters에 해당하는 클러스터가 없습니다.")
  }

  # 2-1) pb 메타 보강 (sample_id / group_id / batch)
  pb_meta <- as.data.frame(SummarizedExperiment::colData(pb))

  # sample_id 없으면 assay의 colnames로 복구
  if (!"sample_id" %in% names(pb_meta)) {
    first_assay <- names(SummarizedExperiment::assays(pb))[1]
    sid_guess <- colnames(SummarizedExperiment::assays(pb)[[first_assay]])
    if (is.null(sid_guess)) stop("pb에 sample_id가 없습니다.")
    pb_meta$sample_id <- sid_guess
    rownames(pb_meta) <- pb_meta$sample_id
    SummarizedExperiment::colData(pb) <- S4Vectors::DataFrame(pb_meta)
  }

  # sce에서 (sample_id -> group_id / batch) map
  sce_meta <- as.data.frame(SummarizedExperiment::colData(sce))
  map_cols <- c("sample_id","group_id")
  if (!is.null(batch_id) && batch_id %in% names(sce_meta)) map_cols <- c(map_cols, batch_id)
  sce_map <- unique(sce_meta[, map_cols, drop=FALSE])

  # pb에 group_id / batch 보강
  pb_meta <- as.data.frame(SummarizedExperiment::colData(pb))
  need_fix <- (!"group_id" %in% names(pb_meta)) ||
              (length(unique(pb_meta$group_id)) < 2) ||
              (all(unique(pb_meta$group_id) %in% c("type","group","group_id", NA, "")))
  if (need_fix || (!is.null(batch_id) && !batch_id %in% names(pb_meta))) {
    pb_meta2 <- dplyr::left_join(pb_meta, sce_map, by = "sample_id")
    if ("group_id.x" %in% names(pb_meta2) && "group_id.y" %in% names(pb_meta2)) {
      pb_meta2$group_id <- ifelse(is.na(pb_meta2$group_id.y), pb_meta2$group_id.x, pb_meta2$group_id.y)
      pb_meta2$group_id.x <- NULL; pb_meta2$group_id.y <- NULL
    }
    rownames(pb_meta2) <- rownames(pb_meta)
    SummarizedExperiment::colData(pb) <- S4Vectors::DataFrame(pb_meta2)
  }

  # factor화
  pb$sample_id <- droplevels(factor(SummarizedExperiment::colData(pb)$sample_id))
  pb$group_id  <- droplevels(factor(SummarizedExperiment::colData(pb)$group_id))
  if (!is.null(batch_id) && batch_id %in% colnames(SummarizedExperiment::colData(pb))) {
    pb[[batch_id]] <- droplevels(factor(SummarizedExperiment::colData(pb)[[batch_id]]))
  }

  # 3) contrast 그룹만 자동 subset
  extract_groups <- function(contrast_str, levels_available){
    z <- gsub("\\s+", "", contrast_str)
    toks <- unique(gsub("^group(_id)?", "", unlist(strsplit(z, "[^A-Za-z0-9_]+"))))
    toks <- toks[nchar(toks) > 0]
    keep <- intersect(toks, levels_available)
    if (length(keep) < 1) {
      g2 <- levels_available[vapply(levels_available, function(g) grepl(g, z), logical(1))]
      keep <- unique(g2)
    }
    keep
  }
  grp_lvls <- levels(SummarizedExperiment::colData(pb)$group_id)
  tg <- extract_groups(contrast, grp_lvls)
  if (length(tg) < 2) stop(sprintf("contrast에서 추출한 그룹이 부족합니다. contrast='%s', 사용가능레벨=%s",
                                   contrast, paste(grp_lvls, collapse=", ")))

  keep_idx <- SummarizedExperiment::colData(pb)$group_id %in% tg
  pb_sub <- pb[, keep_idx]
  pb_sub$group_id <- droplevels(factor(SummarizedExperiment::colData(pb_sub)$group_id))

  # **sce도 동일 기준으로 subset (resDS용 필수)**
  sce_sub <- sce[, sce$sample_id %in% SummarizedExperiment::colData(pb_sub)$sample_id &
                    sce$group_id  %in% tg]
  sce_sub$cluster_id <- droplevels(factor(sce_sub$cluster_id))
  sce_sub$sample_id  <- droplevels(factor(sce_sub$sample_id))
  sce_sub$group_id   <- droplevels(factor(sce_sub$group_id))
  if (!is.null(batch_id) && batch_id %in% colnames(SummarizedExperiment::colData(sce_sub))) {
    sce_sub[[batch_id]] <- droplevels(factor(sce_sub[[batch_id]]))
  }

  # 4) design/contrast (batch는 'batch'로 복사해서 사용)
  pb_sub$group <- pb_sub$group_id
  if (!is.null(batch_id) && batch_id %in% colnames(SummarizedExperiment::colData(pb_sub))) {
    pb_sub$batch <- droplevels(factor(SummarizedExperiment::colData(pb_sub)[[batch_id]]))
    design <- stats::model.matrix(~ 0 + group + batch,
                                  data = as.data.frame(SummarizedExperiment::colData(pb_sub)))
  } else {
    design <- stats::model.matrix(~ 0 + group,
                                  data = as.data.frame(SummarizedExperiment::colData(pb_sub)))
  }

  fix_contrast <- function(contrast_str, design_cols){
    z <- gsub("\\s+", "", contrast_str)
    toks <- unlist(strsplit(z, "([+\\-])", perl=TRUE))
    ops  <- unlist(regmatches(z, gregexpr("([+\\-])", z, perl=TRUE)))
    rebuild <- function(tok){
      tok <- gsub("^group(_id)?", "group", tok)
      if (!grepl("^group", tok)) tok <- paste0("group", tok)
      tok
    }
    toks2 <- vapply(toks, rebuild, character(1))
    out <- toks2[1]; if (length(ops)) for (i in seq_along(ops)) out <- paste0(out, ops[i], toks2[i+1])
    out
  }
  contrast_fixed <- fix_contrast(contrast, colnames(design))
  contrast_matrix <- limma::makeContrasts(contrasts = contrast_fixed, levels = design)

  # 5) pbDS
  res <- muscat::pbDS(
    pb_sub,
    design    = design,
    method    = method,
    contrast  = contrast_matrix,
    min_cells = pb_min_cells,
    filter    = filter_genes,
    verbose   = TRUE
  )

  # 6) 결과 평탄화: **sce_sub가 먼저, res가 다음**
  combined <- muscat::resDS(sce_sub, res)

  # cluster_id 정리 + 라벨
  if ("cluster" %in% names(combined) && !"cluster_id" %in% names(combined)) {
    combined$cluster_id <- combined$cluster
  }
  if (!"cluster_id" %in% names(combined)) stop("resDS 결과에 'cluster_id'가 없습니다.")
  if (!is.null(cluster_label_map)) {
    combined$cluster_label <- cluster_label_map[as.character(combined$cluster_id)]
    combined$cluster_label[is.na(combined$cluster_label)] <- as.character(combined$cluster_id)
  } else {
    combined$cluster_label <- as.character(combined$cluster_id)
  }

  # 7) 요약 헬퍼
  .pick_cols_for <- function(tab, contrast_fixed) {
    cstr <- gsub("\\s+", "", contrast_fixed)
    patt <- paste0("(^|\\.)", gsub("([+\\-])", "\\\\\\1", cstr), "$")

    # 1) 접미사 있는 형태 먼저 시도
    logfc <- grep("^logFC(\\.|_)?", names(tab), value=TRUE); logfc <- logfc[grep(patt, logfc)]
    padj  <- grep("^(p_adj(\\.loc|\\.glb)?|FDR)(\\.|_)?", names(tab), value=TRUE); padj <- padj[grep(patt, padj)]
    pcol  <- grep("^(p_val|PValue)(\\.|_)?", names(tab), value=TRUE); pcol <- pcol[grep(patt, pcol)]

    # 2) 접미사 매칭 실패 시 기본 컬럼으로 폴백
    if (!length(logfc) && "logFC" %in% names(tab)) logfc <- "logFC"
    if (!length(pcol)  && "p_val" %in% names(tab)) pcol  <- "p_val"
    if (!length(padj)) {
      if ("p_adj.loc" %in% names(tab))      padj <- "p_adj.loc"
      else if ("p_adj.glb" %in% names(tab)) padj <- "p_adj.glb"
      else if ("FDR" %in% names(tab))       padj <- "FDR"
    }

    if (!length(logfc) || !length(padj)) {
      stop("결과 테이블에서 logFC/adj.p 컬럼을 찾지 못했습니다. resDS 출력 컬럼명을 확인하세요.")
    }
    list(logfc = logfc[1], padj = padj[1], p = if (length(pcol)) pcol[1] else NULL)
  }

  # --- [REPLACE] 클러스터별 Top-N
  top_by_cluster <- function(n=25){
    tab_use <- combined
    # contrast 컬럼이 있으면 현재 contrast만 사용
    if ("contrast" %in% names(tab_use)) {
      target <- gsub("\\s+", "", contrast_fixed)
      tab_use <- tab_use[gsub("\\s+", "", tab_use$contrast) == target, , drop=FALSE]
    }
    cols <- .pick_cols_for(tab_use, contrast_fixed)

    tab_use |>
      dplyr::mutate(
        logFC_view = .data[[cols$logfc]],
        padj_view  = .data[[cols$padj]],
        p_view     = if (!is.null(cols$p)) .data[[cols$p]] else NA_real_
      ) |>
      dplyr::arrange(padj_view) |>
      dplyr::group_by(cluster_label) |>
      dplyr::slice_head(n = n) |>
      dplyr::ungroup()
  }

  # --- [REPLACE] 전체 Top-N
  top_overall <- function(n=100){
    tab_use <- combined
    if ("contrast" %in% names(tab_use)) {
      target <- gsub("\\s+", "", contrast_fixed)
      tab_use <- tab_use[gsub("\\s+", "", tab_use$contrast) == target, , drop=FALSE]
    }
    cols <- .pick_cols_for(tab_use, contrast_fixed)

    tab_use |>
      dplyr::mutate(
        logFC_view = .data[[cols$logfc]],
        padj_view  = .data[[cols$padj]],
        p_view     = if (!is.null(cols$p)) .data[[cols$p]] else NA_real_
      ) |>
      dplyr::arrange(padj_view) |>
      dplyr::slice_head(n = n)
  }

  list(
    sce        = sce,
    sce_sub    = sce_sub,
    pb         = pb,
    pb_sub     = pb_sub,
    design     = design,
    contrast_fixed = contrast_fixed,
    res_raw    = res,
    combined   = combined,
    top_by_cluster = top_by_cluster,
    top_overall    = top_overall
  )
}
} # End of if (FALSE) block

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

# --------- Analyses for GeoMX Seurat ----------

#' @title Limma-Dream-SVA (LDS) 파이프라인
#' @description SVA로 숨겨진 공변량을 찾고, limma-dream으로 다중 임의 효과
#'              (예: 환자, 배치)를 포함한 LMM을 피팅합니다. 기본적으로 GeoMx data와 같이 샘플 수가 많지 않은 경우에 사용됩니다.
#'
#' @param sobj Seurat 객체, DGEList 객체, 또는 Raw Count Matrix.
#' @param formula (formula) 모델 포뮬러. LMM 수식(lme4)을 따릅니다.
#'                (예: ~ response + (1|emrid) + (1|set))
#' @param meta.data (data.frame) 메타데이터. `sobj`가 Matrix일 경우 필수.
#' @param layer (character) Seurat 객체 사용 시, Raw Count가 저장된 layer.
#' @param n_sv (numeric) 사용할 SV의 개수. `NULL` (기본값)일 경우,
#'             `sv_var_cutoff`에 따라 자동으로 결정됩니다.
#' @param sv_var_cutoff (numeric) `n_sv=NULL`일 때, SV가 설명해야 할
#'                      *잔차 분산(residual variance)*의 누적 비율 (0 ~ 1).
#' @param n_cores (numeric) 병렬 처리에 사용할 CPU 코어 수.
#'
#' @return (list) 'fit' (MArrayLM 객체), 'voom' (EList 객체),
#'         'sva' (SVA 객체), 'final_formula' (사용된 최종 포뮬러)
#'
LDS <- function(sobj,
                formula,
                meta.data = NULL,
                layer = "counts",
                n_sv = NULL,
                sv_var_cutoff = 0.5,
                n_cores = max(1, parallel::detectCores() - 2)) {

  # --- 1. 입력값 검증 및 데이터 추출 ---
  message("1/7: 입력 데이터 처리 중...")
  if (inherits(sobj, "Seurat")) {
    counts_matrix <- GetAssayData(sobj, layer = layer)
    if (is.null(meta.data)) {
      meta.data <- sobj@meta.data
    }
  } else if (inherits(sobj, "DGEList")) {
    counts_matrix <- sobj$counts
    if (is.null(meta.data)) {
      meta.data <- sobj$samples
    }
  } else if (is.matrix(sobj) || inherits(sobj, "dgCMatrix")) {
    counts_matrix <- sobj
    if (is.null(meta.data)) {
      stop("`sobj`가 Matrix일 경우, `meta.data`를 반드시 제공해야 합니다.")
    }
  } else {
    stop("`sobj`는 Seurat, DGEList, 또는 Matrix여야 합니다.")
  }
  
  if (ncol(counts_matrix) != nrow(meta.data)) {
    stop("Count matrix의 열 수(샘플)와 meta.data의 행 수(샘플)가 일치하지 않습니다.")
  }
  
  # 병렬 백엔드 설정
  BPPARAM_SETUP <- MulticoreParam(n_cores)

  # --- 2. SVA를 위한 포뮬러 파싱 ---
  message("2/7: 포뮬러 파싱 중...")
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("LMM 포뮬러 파싱을 위해 'lme4' 패키지가 필요합니다.")
  }
  
  # `formula`에서 고정 효과(Fixed Effects)만 추출
  # SVA의 'mod' 인자 및 `filterByExpr`에 사용됨
  fixed_effects_formula <- lme4::nobars(formula)
  
  if (is.null(fixed_effects_formula)) {
     # 예: ~ (1|emrid) + (1|set) 처럼 고정 효과가 아예 없는 경우
     fixed_effects_formula <- ~ 1 
  }

  # --- 3. DGEList 생성, 필터링, 정규화 ---
  message("3/7: DGEList 생성 및 필터링 중...")
  
  # `filterByExpr`를 위한 디자인 행렬 (고정 효과 기반)
  design_for_filter <- model.matrix(fixed_effects_formula, data = meta.data)
  
  dge <- DGEList(counts_matrix, samples = meta.data)
  
  keep_genes <- filterByExpr(dge, design = design_for_filter)
  dge <- dge[keep_genes, , keep.lib.sizes = FALSE]
  
  if (sum(keep_genes) == 0) {
    stop("모든 유전자가 필터링되었습니다. `filterByExpr` 조건을 확인하십시오.")
  }
  
  dge <- calcNormFactors(dge)
  
  message(sprintf("... 유전자 필터링 완료: %d / %d 개 통과", 
                  sum(keep_genes), nrow(counts_matrix)))

  # --- 4. SVA 실행 (임시 voom 기반) ---
  message("4/7: SVA 실행 (숨겨진 변동성 탐색)...")
  
  mod_sva <- model.matrix(fixed_effects_formula, data = meta.data)
  mod0_sva <- model.matrix(~ 1, data = meta.data)
  
  # SVA는 voom 변환된 데이터에 실행
  v_sva <- voom(dge, mod_sva, plot = FALSE)
  
  # 1) n.sv=NULL로 SVA를 실행하여 *최대* SV 개수 및 SV *모두* 찾기
  # (n_sv=NULL일 때만 SVD를 수행하여 sva_obj$svd를 반환함. 
  #  -> 이 부분이 sva 버전마다 다를 수 있어, 수동 잔차 계산이 더 안정적임)
  sva_obj <- sva(v_sva$E, mod = mod_sva, mod0 = mod0_sva, n.sv = NULL)
  
  n_sv_max <- sva_obj$n.sv
  
  if (n_sv_max == 0) {
    message("... SVA가 유의미한 대리 변수(SV)를 찾지 못했습니다.")
    svs_final <- NULL
    n_sv_final <- 0
  } else {
    message(sprintf("... SVA가 최대 %d개의 SV를 탐지했습니다.", n_sv_max))
    
    # 2) 사용할 SV 개수 결정
    if (!is.null(n_sv)) {
      # 사용자가 SV 개수 명시
      n_sv_final <- min(n_sv, n_sv_max)
      
    } else {
      # (사용자 요청) 잔차 분산 기반으로 SV 개수 자동 결정
      message(sprintf("... 잔차 분산의 %.0f%%를 설명하는 SV 개수 자동 탐색 중...", 
                      sv_var_cutoff * 100))
                      
      # SVA가 사용한 것과 동일한 잔차(Residuals)를 수동 계산
      # lm.fit이 (샘플 x 유전자) 형태의 t(v_sva$E)를 입력받음
      res_matrix <- t(resid(lm.fit(mod_sva, t(v_sva$E))))
      
      # 잔차 행렬의 SVD (Singular Value Decomposition)
      svd_res <- svd(res_matrix)
      
      # 각 SV가 설명하는 분산 비율
      percent_var_explained <- (svd_res$d^2) / sum(svd_res$d^2)
      cumulative_var <- cumsum(percent_var_explained)
      
      # Cutoff를 만족하는 최소 SV 개수 찾기
      n_sv_auto <- which(cumulative_var >= sv_var_cutoff)[1]
      
      if (is.na(n_sv_auto)) { # 모든 SV를 합쳐도 cutoff 미만일 경우
        n_sv_final <- n_sv_max
      } else {
        n_sv_final <- n_sv_auto
      }
      
      message(sprintf("... SV %d개가 잔차 분산의 %.1f%%를 설명합니다.", 
                      n_sv_final, cumulative_var[n_sv_final] * 100))
    }
    
    # 3) 최종 SV 추출 및 메타데이터에 추가
    if (n_sv_final > 0) {
      svs_final <- sva_obj$sv[, 1:n_sv_final, drop = FALSE]
      colnames(svs_final) <- paste0("SV", 1:n_sv_final)
      meta.data <- cbind(meta.data, svs_final)
    }
  }

  # --- 5. 최종 포뮬러 생성 ---
  message("5/7: 최종 모델 포뮬러 생성 중...")
  
  original_formula_str <- paste(deparse(formula), collapse = "")
  
  if (n_sv_final > 0) {
    sv_str <- paste(colnames(svs_final), collapse = " + ")
    # 포뮬러에 SV 추가
    final_formula_str <- paste(original_formula_str, sv_str, sep = " + ")
  } else {
    final_formula_str <- original_formula_str
  }
  
  final_formula <- as.formula(final_formula_str)
  message(sprintf("... 최종 포뮬러: %s", final_formula_str))

  # --- 6. limma-dream 파이프라인 실행 ---
  message(sprintf("6/7: limma-dream 실행 (Core: %d개)...", n_cores))
  
  # 1) 유효 가중치 계산 (voomWithDreamWeights)
  v_dream <- voomWithDreamWeights(
    dge, 
    final_formula, 
    meta.data, 
    BPPARAM = BPPARAM_SETUP
  )
  
  # 2) LMM 피팅 (dream)
  fit_dream <- dream(
    v_dream, 
    final_formula, 
    meta.data, 
    BPPARAM = BPPARAM_SETUP
  )
  
  # 3) Empirical Bayes 조정
  fit_ebayes <- eBayes(fit_dream)
  
  message("7/7: 분석 완료.")

  # --- 7. 결과 반환 ---
  return(list(
    fit = fit_ebayes,       # 최종 MArrayLM 객체 (topTable 사용 가능)
    voom = v_dream,         # 가중치가 포함된 EList 객체
    sva_obj = sva_obj,      # 원본 SVA 객체
    svs_used = svs_final,   # 모델에 실제 사용된 SV 매트릭스
    final_formula = final_formula, # 최종 사용된 포뮬러
    dge = dge               # 필터링/정규화된 DGEList
  ))
}
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
  time_col = "treatment", time_levels = c("pre", "post"), time_labels = c("pre", "post"),
  segment_col = "ck", subset_comp = NULL,
  patient_col = "emrid",
  optimizer = "bobyqa"
) {
  # ==== Data prep ====
  expr <- GetAssayData(sobj, assay = assay, layer = layer)
  meta <- sobj@meta.data %>%
    mutate(
      time    = factor(.data[[time_col]], levels = time_levels, labels = time_labels),
      drug    = factor(.data[[drug_col]]),
      patient = factor(.data[[patient_col]])
    )
  if (!is.null(segment_col)) {
    meta <- meta %>% mutate(
      comp = .data[[segment_col]]
    )
  } else {
    subset_comp <- NULL
  } # 현재 상태로는 comp가 있을지 여부가 갈리므로 후속 분석에 안 좋을 가능성이 크다. 현재는 문제없음.

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
  ctrl <- lmerControl(
    optimizer = optimizer, calc.derivs = TRUE,
    check.conv.singular = "ignore"
  )
  if (nested) {
    fit <- lmer(y ~ time * drug + (1 | drug / patient),
      data = df, REML = FALSE, control = ctrl
    )
  } else {
    fit <- lmer(y ~ time * drug + (1 | patient),
      data = df, REML = FALSE, control = ctrl
    )
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

  p3 <- ggplot(deltas_df, aes(
    x = drug, y = estimate,
    ymin = lower.CL, ymax = upper.CL, color = drug
  )) +
    geom_pointrange(position = position_dodge(width = 0.2)) +
    geom_hline(yintercept = 0, linetype = 2) +
    labs(
      title = sprintf("%s: Δ(post-pre) by drug (EMMeans ±95%% CI)", gene),
      x = "Drug", y = expression(Delta ~ "(post - pre)")
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
pifg <- function(sobj, gene_list) {
  return(lapply(gene_list, function(g) plot_interaction_for_gene(sobj, g, segment_col = NULL)) %>% set_names(gene_list))
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
    labs(
      title = "Δ(post-pre) heatmap by drug",
      x = "Drug", y = "Gene", fill = expression(Delta)
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  list(matrix = mat, plot = p, raw = df_long)
}

# -------- MDMP --------
#' @export
make_delta_matrix_by_patient <- function(
  sobj,
  assay = "RNA", layer = "data",
  time_col = "treatment", drug_col = "drug", patient_col = "emrid", segment_col = "ck",
  time_levels = c("pre", "post"), subset_comp = NULL,
  agg_fun = function(v) mean(v, na.rm = TRUE)
) {
  X <- GetAssayData(sobj, assay = assay, layer = layer) # genes x AOIs
  meta <- sobj@meta.data %>%
    mutate(
      time    = factor(.data[[time_col]], levels = time_levels),
      drug    = factor(.data[[drug_col]]),
      patient = factor(.data[[patient_col]]),
      comp    = factor(.data[[segment_col]])
    )
  stopifnot(identical(colnames(X), rownames(meta)))

  if (!is.null(subset_comp)) {
    keep <- meta$comp == subset_comp
    X <- X[, keep, drop = FALSE]
    meta <- meta[keep, , drop = FALSE]
  }

  # meta column filtering
  gene_names <- rownames(X)

  run <- 0 #
  while (length(which(gene_names %in% names(meta))) > 0) {
    num_meta_in_gene <- which(names(meta) %in% gene_names)
    names(meta)[num_meta_in_gene] <- paste0(names(meta)[num_meta_in_gene], "_re")
    run <- run + 1 #
    print(run) #
    if (run > 5) break #
  }
  # AOI별 표현을 long으로
  long <- as.data.frame(t(X)) %>%
    tibble::rownames_to_column("AOI") %>%
    cbind(meta, .)

  # pre/post를 같은 환자/약제/comp 안에서 집계 후 Δ = post − pre
  # 메모리 절약을 위해 gene-wise로 루프
  make_delta_one_gene <- function(g) {
    dd <- long %>%
      dplyr::select(patient, drug, comp, time, !!sym(g)) %>%
      group_by(patient, drug, comp, time) %>%
      summarise(expr = agg_fun(.data[[g]]), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = time, values_from = expr) %>%
      mutate(delta = .data[[time_levels[2]]] - .data[[time_levels[1]]]) %>%
      dplyr::select(patient, drug, comp, delta)
    dd$gene <- g
    dd
  }

  delta_long <- purrr::map_dfr(gene_names, make_delta_one_gene)

  # 샘플 정의: 환자*약제*(comp)
  delta_long <- delta_long %>%
    mutate(sample_id = if (!is.null(subset_comp)) {
      paste(patient, drug, sep = "|")
    } else {
      paste(patient, drug, comp, sep = "|")
    })

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
  time_col = "treatment", time_levels = c("pre", "post"),
  drug_col = "drug", patient_col = "emrid", segment_col = "ck",
  ref_drug = NULL, optimizer = "bobyqa"
) {
  X <- GetAssayData(sobj, assay = assay, layer = layer)
  meta <- sobj@meta.data %>%
    mutate(
      time    = factor(.data[[time_col]], levels = time_levels),
      drug    = factor(.data[[drug_col]]),
      patient = factor(.data[[patient_col]]),
      comp    = factor(.data[[segment_col]])
    )
  if (!is.null(ref_drug)) meta$drug <- relevel(meta$drug, ref = ref_drug)
  stopifnot(identical(colnames(X), rownames(meta)))

  keep <- if (is.null(subset_comp)) rep(TRUE, nrow(meta)) else meta$comp == subset_comp
  y <- as.numeric(X[gene, keep])
  df <- cbind.data.frame(y = y, meta[keep, , drop = FALSE])

  if (all(is.na(y)) || sd(y, na.rm = TRUE) == 0) {
    return(NULL)
  }
  if (nlevels(df$time) < 2 || nlevels(df$patient) < 2) {
    return(NULL)
  }

  ctrl <- lmerControl(optimizer = optimizer, calc.derivs = TRUE, check.conv.singular = "ignore")
  form <- if (nested) y ~ time * drug + (1 | drug / patient) else y ~ time * drug + (1 | patient)
  fit <- try(lmer(form, data = df, REML = FALSE, control = ctrl), silent = TRUE)
  if (inherits(fit, "try-error")) {
    return(NULL)
  }

  em_time_by_drug <- emmeans(fit, ~ time | drug)
  deltas <- contrast(em_time_by_drug, method = list(delta = c(-1, 1)))
  deltas_df <- as.data.frame(confint(deltas))
  deltas_df$gene <- gene
  deltas_df
}

# 여러 유전자 병렬 처리하여 Δ(post-pre) 행렬 (genes x drugs) 구축
#' @export
make_delta_emm_matrix <- function(
  sobj, genes,
  subset_comp = NULL, nested = FALSE,
  assay = "RNA", layer = "data",
  time_col = "treatment", time_levels = c("pre", "post"),
  drug_col = "drug", patient_col = "emrid", segment_col = "ck",
  ref_drug = NULL, optimizer = "bobyqa",
  mc.cores = 1
) {
  fun <- function(g) {
    delta_one_gene_emm(
      sobj, g, nested, subset_comp, assay, layer,
      time_col, time_levels, drug_col, patient_col, segment_col,
      ref_drug, optimizer
    )
  }
  res_list <- if (mc.cores > 1 && .Platform$OS.type != "windows") {
    parallel::mclapply(genes, fun, mc.cores = mc.cores)
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
.make_pathway_mat <- function(genes, gene_sets = list()) {
  if (length(gene_sets) == 0) {
    return(NULL)
  }
  pws <- names(gene_sets)
  m <- sapply(gene_sets, function(gv) as.integer(genes %in% gv))
  colnames(m) <- pws
  rownames(m) <- genes
  m
}

# -------- helper: p-value 별표 --------
.stars <- function(p) {
  ifelse(p < 0.001, "***",
    ifelse(p < 0.01, "**",
      ifelse(p < 0.05, "*", "")
    )
  )
}

# -------- PDRM(MDEM$delta_df); options: top variable genes? or gene list. main plotting function --------
#' @export
plot_drug_response_map <- function(
  df_long, # gene, drug, estimate
  drug_annot = NULL, # optional: drug, class, route ...
  gene_sets = list(), # optional: list(pathway = c(genes))
  sig_df = NULL, # optional: gene, drug, p_adj
  clip_q = c(0.01, 0.99), # quantile clipping to tame outliers
  scale_row = TRUE, # gene-wise z
  show_rownames = FALSE,
  show_colnames = TRUE,
  top_n_genes = NULL, # 선택: 상위 변이 gene만 (예: 2000)
  row_split_k = NULL, # 선택: 행 클러스터 k-split
  col_split_k = NULL # 선택: 열 클러스터 k-split
) {
  # 1) wide matrix
  mat <- df_long %>%
    select(gene, drug, estimate) %>%
    pivot_wider(names_from = drug, values_from = estimate) %>%
    tibble::column_to_rownames("gene") %>%
    as.matrix()

  # 2) optional gene filter by variance
  if (!is.null(top_n_genes)) {
    v <- apply(mat, 1, function(x) var(x, na.rm = TRUE))
    keep <- names(sort(v, decreasing = TRUE))[seq_len(min(top_n_genes, length(v)))]
    mat <- mat[keep, , drop = FALSE]
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
  mat_use[is.na(mat_use)] <- 0 # 안전장치

  # 5) color scale (대칭)
  col_fun <- colorRamp2(c(min(mat_use), 0, max(mat_use)), c("blue", "white", "red"))

  # 6) column annotation (drug info)
  ha_col <- NULL
  if (!is.null(drug_annot)) {
    drug_annot <- as.data.frame(drug_annot)
    rownames(drug_annot) <- drug_annot$drug
    drug_annot <- drug_annot[colnames(mat_use), setdiff(colnames(drug_annot), "drug"), drop = FALSE]
    # 범주형에 대한 색 자동 생성
    ann_cols <- lapply(drug_annot, function(v) {
      if (is.numeric(v)) {
        NULL
      } else {
        structure(scales::hue_pal()(length(unique(v))), names = unique(v))
      }
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
      pw_count = anno_barplot(pw_count, gp = gpar(fill = "#888888", col = NA), width = unit(0.8, "cm"))
    )
  }

  # 8) 유의성 별표 overlay (선택)
  layer_fun <- NULL
  if (!is.null(sig_df)) {
    sig_df <- sig_df %>% mutate(star = .stars(p_adj))
    # matrix of stars aligned to mat_use
    star_mat <- matrix("",
      nrow = nrow(mat_use), ncol = ncol(mat_use),
      dimnames = dimnames(mat_use)
    )
    common_genes <- intersect(rownames(mat_use), sig_df$gene)
    common_drugs <- intersect(colnames(mat_use), sig_df$drug)
    sig_sub <- sig_df %>% filter(gene %in% common_genes, drug %in% common_drugs)
    idx_r <- match(sig_sub$gene, rownames(mat_use))
    idx_c <- match(sig_sub$drug, colnames(mat_use))
    star_mat[cbind(idx_r, idx_c)] <- sig_sub$star

    layer_fun <- function(j, i, x, y, w, h, fill) {
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
    stop(
      "제공된 컬럼명 중 일부가 메타데이터에 없습니다: ",
      paste(required_cols[!required_cols %in% colnames(meta_data)], collapse = ", ")
    )
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
          if (test.use == "wilcox.test") {
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
    ggplot2::geom_point(
      position = ggplot2::position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
      alpha = 0.7, size = 1.5
    ) +
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

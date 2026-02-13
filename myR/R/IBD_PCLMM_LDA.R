PCLMM <- function(
  sobj,
  npcs = 10,
  fixed = c("ck_plot", "response", "treatment"),
  random = c("emrid")
) {
    library(lme4)
    library(lmerTest)
    library(dplyr)
    library(broom.mixed)

    pcs <- paste0("PC_", 1:npcs)

    meta <- sobj@meta.data

    # PCA metadata 추가
    if (!all(pcs %in% colnames(meta))) {
        sobj <- AddMetaData(sobj, Embeddings(sobj, "pca"))
        meta <- sobj@meta.data
    }

    fixed <- intersect(fixed, colnames(meta))
    random <- intersect(random, colnames(meta))

    if (length(fixed) == 0) {
        stop("No fixed factors found")
    }

    fixed_str <- paste(fixed, collapse = " + ")
    random_str <- if (length(random) > 0) {
        paste(paste0("(1|", random, ")"), collapse = " + ")
    } else {
        NULL
    }

    rhs <- paste(c(fixed_str, random_str), collapse = " + ")

    res_list <- lapply(pcs, function(pc) {
        form <- as.formula(paste(pc, "~", rhs))

        fit <- tryCatch(
            lmer(form, data = meta),
            error = function(e) NULL
        )
        if (is.null(fit)) {
            return(NULL)
        }

        # ---------- ANOVA (factor-level p) ----------
        aov_tab <- as.data.frame(anova(fit))
        aov_tab$term <- rownames(aov_tab)

        aov_keep <- aov_tab %>%
            dplyr::filter(term %in% fixed) %>%
            dplyr::select(term, anova_p = `Pr(>F)`)

        # ---------- coefficient table ----------
        coef_tab <- broom.mixed::tidy(fit, effects = "fixed") %>%
            dplyr::rename(coef_term = term)

        # coef_term 예:
        # responseResponder
        # treatmentDrugA
        # ck_plotTRUE

        # factor 이름 매칭 (startsWith)
        coef_tab$term <- sapply(coef_tab$coef_term, function(x) {
            hit <- fixed[startsWith(x, fixed)]
            if (length(hit) == 0) NA else hit[1]
        })

        coef_tab <- coef_tab %>%
            dplyr::filter(!is.na(term)) %>%
            dplyr::select(
                term,
                coef_term,
                estimate,
                std.error,
                statistic,
                p.value
            )

        # ---------- merge ----------
        out <- left_join(coef_tab, aov_keep, by = "term") %>%
            dplyr::mutate(
                PC = pc,
                singular = isSingular(fit)
            )

        out
    })

    res <- bind_rows(res_list)

    if (nrow(res) == 0) {
        stop("No model results produced")
    }

    # FDR — factor별
    res <- res %>%
        dplyr::group_by(term) %>%
        dplyr::mutate(padj = p.adjust(p.value, "BH")) %>%
        dplyr::ungroup()

    return(res)
}

LDAbyFeature <- function(
  sobj,
  features = c("sex", "ck_plot"),
  npcs = 10,
  reduction = "pca",
  prefix = "LD1_",
  drop_na = TRUE,
  verbose = TRUE
) {
    stopifnot(inherits(sobj, "Seurat"))
    if (!requireNamespace("MASS", quietly = TRUE)) stop("Please install MASS")

    # PCA embeddings (cells x PCs)
    E <- Embeddings(sobj, reduction = reduction)
    if (is.null(E) || nrow(E) == 0) stop("No PCA embeddings found. RunPCA first.")

    npcs <- min(npcs, ncol(E))
    pcs_use <- colnames(E)[seq_len(npcs)]

    meta <- sobj@meta.data

    out_info <- list()

    for (feat in features) {
        if (!feat %in% colnames(meta)) {
            if (verbose) message("Skipping '", feat, "': not found in metadata")
            next
        }

        y0 <- meta[[feat]]

        # NA 처리
        keep <- rep(TRUE, length(y0))
        if (drop_na) keep <- keep & !is.na(y0)

        X <- E[keep, pcs_use, drop = FALSE]
        y <- y0[keep]

        # factor로
        y <- as.factor(y)
        y <- droplevels(y)

        # 클래스 수 체크
        k <- nlevels(y)
        if (k < 2) {
            if (verbose) message("Skipping '", feat, "': <2 levels after filtering")
            next
        }

        # 소표본/고차원 방지: PC 수를 (min_class_n - 1)로 제한
        min_class_n <- min(table(y))
        max_p <- max(1, min(npcs, min_class_n - 1))
        if (max_p < 1) {
            if (verbose) message("Skipping '", feat, "': too few samples per class")
            next
        }
        X <- X[, seq_len(max_p), drop = FALSE]

        # LDA (PC space)
        lda_fit <- MASS::lda(x = X, grouping = y)

        # LD1 scores
        ld <- predict(lda_fit)$x
        ld1 <- as.numeric(ld[, 1])

        # 전체 cell 길이에 맞춰 넣기 (NA 포함)
        ld1_full <- rep(NA_real_, nrow(meta))
        ld1_full[which(keep)] <- ld1

        new_col <- paste0(prefix, feat)
        sobj[[new_col]] <- ld1_full

        out_info[[feat]] <- list(
            feature = feat,
            ld1_col = new_col,
            n_used = sum(keep),
            n_levels = k,
            levels = levels(y),
            n_pcs_used = ncol(X)
        )

        if (verbose) {
            message(
                "OK: ", feat, " → ", new_col,
                " (n=", sum(keep), ", levels=", k, ", PCs=", ncol(X), ")"
            )
        }
    }

    return(list(sobj = sobj, info = out_info))
}
LD1TopGenes <- function(
  sobj,
  feature,
  npcs = 10,
  reduction = "pca",
  top_n = 30,
  abs_rank = TRUE,
  drop_na = TRUE
) {
    stopifnot(inherits(sobj, "Seurat"))
    if (!requireNamespace("MASS", quietly = TRUE)) stop("Please install MASS")

    meta <- sobj@meta.data
    if (!feature %in% colnames(meta)) stop("feature not found in metadata: ", feature)

    # PCA embedding + loading
    E <- Embeddings(sobj, reduction = reduction)
    L <- Loadings(sobj, reduction = reduction) # (genes x PCs)
    if (is.null(E) || nrow(E) == 0) stop("No embeddings found. RunPCA first.")
    if (is.null(L) || nrow(L) == 0) stop("No loadings found. RunPCA first.")

    npcs <- min(npcs, ncol(E), ncol(L))
    pcs_use <- colnames(E)[seq_len(npcs)]

    y0 <- meta[[feature]]
    keep <- rep(TRUE, length(y0))
    if (drop_na) keep <- keep & !is.na(y0)

    y <- droplevels(as.factor(y0[keep]))
    if (nlevels(y) < 2) stop("feature has <2 levels after filtering.")

    X <- E[keep, pcs_use, drop = FALSE]

    # 소표본 방지: PC 수를 (min_class_n - 1)로 제한
    min_class_n <- min(table(y))
    max_p <- max(1, min(npcs, min_class_n - 1))
    X <- X[, seq_len(max_p), drop = FALSE]
    pcs_use2 <- colnames(X)

    # LDA fit
    lda_fit <- MASS::lda(x = X, grouping = y)

    # LDA scaling (PC weights): (p x (k-1)) ; 우리는 LD1만
    w <- lda_fit$scaling[, 1]
    w <- as.numeric(w)
    names(w) <- rownames(lda_fit$scaling)

    # PCA gene loadings subset: (genes x p)
    Lsub <- L[, pcs_use2, drop = FALSE]

    # gene weight for LD1: (genes x p) %*% (p)
    gene_w <- as.numeric(Lsub %*% w)
    names(gene_w) <- rownames(Lsub)

    df <- data.frame(
        gene = names(gene_w),
        ld1_gene_weight = gene_w,
        abs_weight = abs(gene_w),
        stringsAsFactors = FALSE
    )

    if (abs_rank) {
        df <- df[order(df$abs_weight, decreasing = TRUE), ]
    } else {
        df <- df[order(df$ld1_gene_weight, decreasing = TRUE), ]
    }

    df$rank <- seq_len(nrow(df))
    df <- head(df, top_n)

    # 방향 구분(참고용)
    df$direction <- ifelse(df$ld1_gene_weight >= 0, "pos", "neg")
    df$feature <- feature
    df$n_pcs_used <- length(pcs_use2)
    df$levels <- paste(levels(y), collapse = " vs ")

    rownames(df) <- NULL
    return(df)
}
GeneAUC <- function(
  sobj,
  genes,
  label = "ck_plot", # 이진 factor/0-1 변수
  assay = NULL, # NULL이면 DefaultAssay(sobj)
  layer = "data", # "data"(log-normalized), "counts", "scale.data"
  positive = NULL # 양성 클래스 레벨 지정(예: "TRUE" 또는 "CK+")
) {
    if (!requireNamespace("pROC", quietly = TRUE)) {
        stop("Please install pROC: install.packages('pROC')")
    }

    if (is.null(assay)) assay <- Seurat::DefaultAssay(sobj)
    Seurat::DefaultAssay(sobj) <- assay

    genes <- intersect(genes, rownames(sobj))
    if (length(genes) == 0) stop("No genes found in sobj rownames.")

    if (!label %in% colnames(sobj@meta.data)) stop("label not in metadata: ", label)

    # expression + label 가져오기
    df <- Seurat::FetchData(sobj, vars = c(label, genes), layer = layer)

    y0 <- df[[label]]
    # 이진화
    if (!is.factor(y0)) y0 <- as.factor(y0)
    y0 <- droplevels(y0)

    if (nlevels(y0) != 2) stop("label must be binary (2 levels). Current levels: ", paste(levels(y0), collapse = ", "))

    # positive 클래스 지정
    if (is.null(positive)) positive <- levels(y0)[2]
    if (!positive %in% levels(y0)) stop("positive not in label levels.")

    # pROC은 response를 factor로 주면 levels 순서를 참고하므로, levels를 (neg, pos)로 맞춤
    neg <- setdiff(levels(y0), positive)
    y <- factor(y0, levels = c(neg, positive))

    res <- lapply(genes, function(g) {
        x <- df[[g]]

        keep <- is.finite(x) & !is.na(y)
        x <- x[keep]
        yy <- y[keep]

        # 변이가 없으면 AUC 불가
        if (length(unique(x)) < 2) {
            return(data.frame(gene = g, auc = NA_real_, n = length(x), n_pos = sum(yy == positive), n_neg = sum(yy == neg)))
        }

        roc_obj <- pROC::roc(response = yy, predictor = x, quiet = TRUE, direction = "auto")
        data.frame(
            gene = g,
            auc = as.numeric(pROC::auc(roc_obj)),
            n = length(x),
            n_pos = sum(yy == positive),
            n_neg = sum(yy == neg)
        )
    })

    out <- do.call(rbind, res)
    out <- out[order(out$auc, decreasing = TRUE), ]
    rownames(out) <- NULL
    out
}

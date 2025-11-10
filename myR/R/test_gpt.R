#' Set parallel & threading profile for RStudio vs Terminal
#' @export
#'
#' @param mode One of c("rstudio","terminal","auto").
#'   - "rstudio": use PSOCK / multisession (no forking).
#'   - "terminal": use Multicore where possible.
#'   - "auto": detect RStudio and choose automatically.
#' @param workers Integer number of workers. Default: parallel::detectCores()-1.
#' @param blas_threads Integer BLAS threads (via RhpcBLASctl if available).
#' @return Invisibly returns a list with chosen backends.
#' @examples
#' set_parallel_profile(mode="auto", workers=8, blas_threads=1)
#' @export
set_parallel_profile <- function(mode = c("auto","rstudio","terminal"),
                                 workers = max(1, parallel::detectCores()-1),
                                 blas_threads = 1) {
  mode <- match.arg(mode)
  # Detect RStudio
  in_rstudio <- FALSE
  if (mode == "auto") {
    in_rstudio <- isTRUE(Sys.getenv("RSTUDIO") == "1")
  } else if (mode == "rstudio") {
    in_rstudio <- TRUE
  }

  # BLAS threads
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(blas_threads)
    RhpcBLASctl::omp_set_num_threads(blas_threads)
  }

  # future backend
  if (requireNamespace("future", quietly = TRUE)) {
    if (in_rstudio) {
      future::plan(future::multisession, workers = workers)
    } else {
      # Terminal: prefer multicore on Unix; fallback to multisession
      if (.Platform$OS.type == "unix") {
        future::plan(future::multicore, workers = workers)
      } else {
        future::plan(future::multisession, workers = workers)
      }
    }
  }

  # BiocParallel backend
  if (requireNamespace("BiocParallel", quietly = TRUE)) {
    if (in_rstudio || .Platform$OS.type != "unix") {
      BiocParallel::register(BiocParallel::SnowParam(workers = workers, type = "SOCK"))
    } else {
      BiocParallel::register(BiocParallel::MulticoreParam(workers = workers))
    }
  }

  invisible(list(rstudio = in_rstudio, workers = workers, blas_threads = blas_threads))
}

#' Set parallel and threading profile safely (RStudio / Terminal)
#'
#' Uses only \code{future} by default. Optionally registers a BiocParallel backend
#' with correct types (\code{type="SOCK"}). Keeps BLAS/OMP to 1 thread to avoid oversubscription.
#'
#' @param mode One of c("auto","rstudio","terminal").
#' @param workers Number of workers for \code{future} (default: detectCores()-1).
#' @param blas_threads BLAS/OMP threads (default 1).
#' @param use_biocparallel Logical; if TRUE, also register a BiocParallel backend.
#' @return Invisibly returns a list with chosen backends.
#' @examples
#' set_parallel_profile_v2("auto", workers = 16, use_biocparallel = FALSE)
set_parallel_profile_v2 <- function(mode = c("auto","rstudio","terminal"),
                                    workers = max(1, parallel::detectCores()-1),
                                    blas_threads = 1,
                                    use_biocparallel = FALSE) {
  mode <- match.arg(mode)

  # --- Threads to 1 ---
  Sys.setenv(OMP_NUM_THREADS = as.character(blas_threads),
             OPENBLAS_NUM_THREADS = as.character(blas_threads),
             MKL_NUM_THREADS = as.character(blas_threads),
             VECLIB_MAXIMUM_THREADS = as.character(blas_threads),
             NUMEXPR_NUM_THREADS = as.character(blas_threads))
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(blas_threads)
    RhpcBLASctl::omp_set_num_threads(blas_threads)
  }

  # --- Detect RStudio ---
  in_rstudio <- if (mode == "auto") isTRUE(Sys.getenv("RSTUDIO") == "1") else (mode == "rstudio")

  # --- future backend ---
  if (requireNamespace("future", quietly = TRUE)) {
    if (in_rstudio) {
      future::plan(future::multisession, workers = workers)
    } else {
      # Unix 터미널이면 multicore가 약간 더 가볍지만, 호환성 위해 multisession로 통일해도 OK
      if (.Platform$OS.type == "unix") {
        tryCatch(
          future::plan(future::multicore, workers = workers),
          error = function(e) future::plan(future::multisession, workers = workers)
        )
      } else {
        future::plan(future::multisession, workers = workers)
      }
    }
  }

  # --- (옵션) BiocParallel 등록: type 이름은 "SOCK"이 맞습니다 ---
  if (use_biocparallel && requireNamespace("BiocParallel", quietly = TRUE)) {
    bp <- try({
      if (in_rstudio || .Platform$OS.type != "unix") {
        BiocParallel::SnowParam(workers = workers, type = "SOCK")
      } else {
        # 터미널/유닉스에서는 MulticoreParam 사용 가능
        BiocParallel::MulticoreParam(workers = workers)
      }
    }, silent = TRUE)
    if (!inherits(bp, "try-error")) BiocParallel::register(bp)
  }

  invisible(list(
    rstudio = in_rstudio,
    workers = workers,
    blas_threads = blas_threads,
    future_plan = if (requireNamespace("future", quietly = TRUE)) future::plan() else NULL,
    biocparallel = use_biocparallel
  ))
}

set_parallel_profile_v3 <- function(mode = c("auto","rstudio","terminal"),
                                    workers = max(1, parallel::detectCores()-1),
                                    blas_threads = 1,
                                    use_biocparallel = FALSE) {
  mode <- match.arg(mode)

  # --- Threads to 1 ---
  Sys.setenv(OMP_NUM_THREADS = as.character(blas_threads),
             OPENBLAS_NUM_THREADS = as.character(blas_threads),
             MKL_NUM_THREADS = as.character(blas_threads),
             VECLIB_MAXIMUM_THREADS = as.character(blas_threads),
             NUMEXPR_NUM_THREADS = as.character(blas_threads))
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(blas_threads)
    RhpcBLASctl::omp_set_num_threads(blas_threads)
  }

  # --- Detect RStudio ---
  in_rstudio <- if (mode == "auto") isTRUE(Sys.getenv("RSTUDIO") == "1") else (mode == "rstudio")

  # --- future backend ---
  if (requireNamespace("future", quietly = TRUE)) {
    if (in_rstudio) {
      future::plan(future::multisession, workers = workers)
    } else {
      # Unix 터미널이면 multicore가 약간 더 가볍지만, 호환성 위해 multisession로 통일해도 OK
      if (.Platform$OS.type == "unix") {
        tryCatch(
          future::plan(future::multicore, workers = workers),
          error = function(e) future::plan(future::multisession, workers = workers)
        )
      } else {
        future::plan(future::multisession, workers = workers)
      }
    }
  }

  # --- (옵션) BiocParallel 등록: type 이름은 "SOCK"이 맞습니다 ---
  if (use_biocparallel && requireNamespace("BiocParallel", quietly = TRUE)) {
    bp <- try({
      if (in_rstudio || .Platform$OS.type != "unix") {
        BiocParallel::SnowParam(workers = workers, type = "SOCK")
      } else {
        # 터미널/유닉스에서는 MulticoreParam 사용 가능
        BiocParallel::MulticoreParam(workers = workers)
      }
    }, silent = TRUE)
    if (!inherits(bp, "try-error")) BiocParallel::register(bp)
  }

  invisible(list(
    rstudio = in_rstudio,
    workers = workers,
    blas_threads = blas_threads,
    future_plan = if (requireNamespace("future", quietly = TRUE)) future::plan() else NULL,
    biocparallel = use_biocparallel
  ))
}

#' Minimal, safe parallel setup (future-only; no BiocParallel)
#'
#' @param workers integer, number of workers for \code{future}. Default: max(1, parallel::detectCores()-1).
#' @param blas_threads integer, set BLAS/OMP threads. Default 1.
#' @return Invisibly returns list(workers, plan).
#' @examples
#' set_parallel_minimal(workers = 16, blas_threads = 1)
set_parallel_minimal <- function(workers = max(1, parallel::detectCores()-1), blas_threads = 1) {
  Sys.setenv(
    OMP_NUM_THREADS = as.character(blas_threads),
    OPENBLAS_NUM_THREADS = as.character(blas_threads),
    MKL_NUM_THREADS = as.character(blas_threads),
    VECLIB_MAXIMUM_THREADS = as.character(blas_threads),
    NUMEXPR_NUM_THREADS = as.character(blas_threads)
  )
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(blas_threads)
    RhpcBLASctl::omp_set_num_threads(blas_threads)
  }
  if (requireNamespace("future", quietly = TRUE)) {
    future::plan(future::multisession, workers = workers)
  }
  invisible(list(workers = workers, plan = if (requireNamespace("future", quietly = TRUE)) future::plan() else NULL))
}



#' Fast neighborhood-by-sample counts using sparse matrix multiplication
#' @export
#'
#' @description
#' Replaces miloR::countCells() with a vectorized sparse multiplication:
#' counts(H × S) = t(nhoods(cells × H)) %*% model.matrix(~0 + samples) (cells × S).
#'
#' @param milo A \code{Milo} object with neighborhoods already built.
#' @param samples Factor vector of length = ncol(milo) giving sample IDs
#'   (e.g., patient IDs like \code{hos_no}).
#' @return A matrix of neighborhood-by-sample cell counts (H × S).
#' @examples
#' counts <- fast_count_cells(milo, samples = factor(colData(milo)$hos_no))
#' @export
fast_count_cells <- function(milo, samples) {
  stopifnot(inherits(milo, "Milo"))
  X <- miloR::nhoods(milo)              # cells x nhoods (dgCMatrix)
  stopifnot(nrow(X) == ncol(milo))      # rows=cells, columns=nhoods
  if (!is.factor(samples)) samples <- factor(samples)
  # cells x S one-hot sparse matrix
  S <- Matrix::sparse.model.matrix(~ 0 + samples)
  # (H x S) = t(cells x H) %*% (cells x S)
  counts <- Matrix::t(X) %*% S
  as.matrix(counts)
}

#' Fast neighborhood-by-sample counts (internal)
#' @keywords internal
.fast_count_cells <- function(milo, samples) {
  stopifnot(inherits(milo, "Milo"))
  X <- miloR::nhoods(milo)                   # cells x H (dgCMatrix)
  stopifnot(nrow(X) == ncol(milo))
  if (!is.factor(samples)) samples <- factor(samples)
  S <- Matrix::sparse.model.matrix(~ 0 + samples)  # cells x S
  Matrix::t(X) %*% S |> as.matrix()               # H x S
}

milo_neighborhood_da <- function(
  sobj,
  target_var   = "g3",
  batch_var    = "GEM",
  patient_var  = "hos_no",
  cluster_var  = "anno3.scvi",
  reduced_dim  = "PCA",
  k            = 30,
  d            = 50,
  prop         = 0.1,
  fixed_effects  = NULL,            # 예: c("g3","GEM")
  random_effects = NULL,            # 예: c("hos_no")
  glmm_family    = c("nb","poisson"),
  p_adj_method   = "BH",
  fdr_thresh     = 0.1,
  logfc_thresh   = 0.5,
  seed           = 1
){
  suppressPackageStartupMessages({
    library(Seurat)
    library(SingleCellExperiment)
    library(miloR)
    library(Matrix)
    library(dplyr)
    library(tibble)
    library(purrr)
  })
  set.seed(seed)

  stopifnot(all(c(target_var,batch_var,patient_var,cluster_var) %in% colnames(sobj@meta.data)))

  # 0) 차원축소 준비
  if (!reduced_dim %in% c("PCA","UMAP","TSNE")) reduced_dim <- "PCA"
  if (!"pca" %in% Reductions(sobj)) {
    DefaultAssay(sobj) <- "RNA"
    sobj <- NormalizeData(sobj, verbose=FALSE)
    sobj <- FindVariableFeatures(sobj, nfeatures=3000, verbose=FALSE)
    sobj <- ScaleData(sobj, verbose=FALSE)
    sobj <- RunPCA(sobj, npcs=max(50,d), verbose=FALSE)
  }
  if (!"umap" %in% Reductions(sobj)) {
    sobj <- RunUMAP(sobj, dims=1:min(30,ncol(Embeddings(sobj,"pca"))), reduction="pca", verbose=FALSE)
  }

  # 1) Seurat -> SCE
  sce <- as.SingleCellExperiment(sobj, assay="RNA")
  reducedDim(sce, "PCA")  <- Embeddings(sobj, "pca")
  if ("umap" %in% Reductions(sobj)) reducedDim(sce,"UMAP") <- Embeddings(sobj,"umap")

  # 2) Milo 객체 & 이웃 구성
  milo <- Milo(sce)
  milo <- buildGraph(milo, k = k, d = d, reduced.dim = "PCA")
  milo <- makeNhoods(milo, prop = prop, k = k, refined = TRUE, reduced_dims = "PCA")

  # 3) 이웃별 샘플 카운트
  md <- as.data.frame(colData(milo))
  md[[patient_var]] <- factor(md[[patient_var]])
  md[[target_var]]  <- factor(md[[target_var]])
  md[[batch_var]]   <- factor(md[[batch_var]])
  md[[cluster_var]] <- factor(md[[cluster_var]])

  milo <- countCells(milo, meta.data = md, samples = patient_var)

  # 4) 기본 설계식 구성
  if (is.null(fixed_effects)) fixed_effects <- c(target_var, batch_var)

  # ===== 분기: 랜덤효과 유무 =====
  use_glmm <- length(random_effects) > 0

  if (!use_glmm) {
    # ---- (A) miloR 표준: edgeR GLM (고정효과만) ----
    design.df <- md[, unique(c(patient_var, fixed_effects, random_effects)), drop=FALSE]
    rownames(design.df) <- colnames(milo)  # samples
    # 기준레벨 설정(첫 레벨이 reference)
    for (v in fixed_effects) design.df[[v]] <- droplevels(factor(design.df[[v]]))
    mm <- model.matrix(as.formula(paste("~", paste(fixed_effects, collapse = " + "))), data = design.df)
    da.res <- testNhoods(milo, design = mm, reduced.dim = "PCA")
    da.res <- graphSpatialFDR(milo, da.res)
    # 정리
    nh_mat <- nhoods(milo)  # cells x nhoods
    cell_clusters <- md[[cluster_var]]
    major_cluster <- apply(nh_mat, 2, function(idx){
      labs <- cell_clusters[as.logical(idx)]
      if (!length(labs)) NA else names(sort(table(labs), decreasing=TRUE))[1]
    })
    da.res$major_cluster <- major_cluster[match(da.res$nhood, as.integer(colnames(nh_mat)))]
    da.res <- da.res %>% mutate(is_sig = spatialFDR < fdr_thresh & abs(logFC) > logfc_thresh)

    cluster_summary <- da.res %>% filter(is_sig) %>%
      group_by(major_cluster) %>%
      summarise(n_sig_nhoods = n(),
                mean_logFC = mean(logFC),
                min_spatialFDR = min(spatialFDR),
                .groups="drop") %>%
      arrange(desc(n_sig_nhoods))

    return(list(
      mode                  = "edgeR_fixed_effects",
      milo                  = milo,
      da_nhood_with_cluster = da.res,
      da_by_cluster         = cluster_summary
    ))
  }

  # ---- (B) GLMM 모드: glmmTMB로 랜덤효과 지원 ----
  # 주의: miloR의 공간 FDR 보정은 edgeR 경로에 붙습니다.
  # GLMM 경로에선 per-nhood 회귀 후 p값 보정만 수행(보수적 해석 권장).

  suppressPackageStartupMessages({
    library(glmmTMB)
    library(broom.mixed)
    library(future.apply)
  })

  # nhood x sample 카운트 행렬 추출
  nh_counts <- nhoodCounts(milo)    # nhoods x samples
  samples   <- colnames(nh_counts)
  stopifnot(identical(samples, rownames(md)))

  # 샘플별 라이브러리(셀) 크기 → offset로 사용
  sample_total_cells <- Matrix::colSums(nhoods(milo))  # 각 sample의 총 셀수
  # 안정화를 위해 size-factor: 총셀수 / 평균총셀수
  size_factor <- sample_total_cells / mean(sample_total_cells)
  log_offset  <- log(size_factor + 1e-8)

  # 긴형 데이터 생성(메모리 절약 위해 필요 열만)
  sample_df <- tibble(
    sample = samples,
    !!target_var  := droplevels(factor(md[[target_var]])),
    !!batch_var   := droplevels(factor(md[[batch_var]])),
    !!patient_var := droplevels(factor(md[[patient_var]])),
    log_offset    = as.numeric(log_offset)
  )

  fam <- match.arg(glmm_family)
  glmm_family_obj <- if (fam == "nb") glmmTMB::nbinom2(link="log") else poisson(link="log")

  # formula 구성: fixed + random
  # 예) count ~ g3 + GEM + (1|hos_no) + offset(log_offset)
  if (is.null(fixed_effects)) fixed_effects <- c(target_var, batch_var)
  fixed_part  <- paste(fixed_effects, collapse = " + ")
  random_part <- paste(sprintf("(1|%s)", random_effects), collapse = " + ")
  rhs <- paste(c(fixed_part, random_part), collapse = " + ")
  full_formula <- as.formula(paste("count ~", rhs, "+ offset(log_offset)"))

  # nhood별 GLMM 피팅
  nh_ids <- rownames(nh_counts)
  plan_old <- future::plan()
  on.exit(try(future::plan(plan_old), silent=TRUE), add=TRUE)
  future::plan(future::sequential)

  fit_one <- function(i){
    y <- as.numeric(nh_counts[i, ])
    dat <- dplyr::mutate(sample_df, count = y)
    # rare case: 모두 0이면 스킵
    if (sum(y) == 0) return(tibble(nhood = nh_ids[i], term = NA, estimate = NA, std.error = NA, statistic = NA, p.value = NA))
    fit <- try(glmmTMB::glmmTMB(formula = full_formula, data = dat, family = glmm_family_obj), silent = TRUE)
    if (inherits(fit, "try-error")) {
      return(tibble(nhood = nh_ids[i], term = NA, estimate = NA, std.error = NA, statistic = NA, p.value = NA))
    }
    co <- try(broom.mixed::tidy(fit, effects = "fixed"), silent = TRUE)
    if (inherits(co, "try-error")) {
      return(tibble(nhood = nh_ids[i], term = NA, estimate = NA, std.error = NA, statistic = NA, p.value = NA))
    }
    co$nhood <- nh_ids[i]
    as_tibble(co)
  }

  glmm_tab <- future.apply::future_lapply(seq_along(nh_ids), fit_one)
  glmm_tab <- bind_rows(glmm_tab)

  # 관심 계수: target_var의 각 대비(참조수준 대비). term 필터
  target_terms <- grep(paste0("^", target_var), unique(glmm_tab$term), value=TRUE)
  res <- glmm_tab %>%
    filter(term %in% target_terms) %>%
    group_by(nhood) %>%
    slice_head(n = 1) %>%  # 대표 대비 1개(필요시 다 보존 가능)
    ungroup() %>%
    mutate(p_adj = p.adjust(p.value, method = p_adj_method),
           logFC = estimate)  # 로그링크 → 계수=log rate ratio

  # 이웃→클러스터 다수표 매핑
  nh_mat <- nhoods(milo)  # cells x nhoods
  cell_clusters <- md[[cluster_var]]
  major_cluster <- apply(nh_mat, 2, function(idx){
    labs <- cell_clusters[as.logical(idx)]
    if (!length(labs)) NA else names(sort(table(labs), decreasing=TRUE))[1]
  })
  # nh_counts rownames = nhood ids (character). nh_mat colnames = nhood ids
  res$major_cluster <- major_cluster[match(res$nhood, colnames(nh_mat))]

  res <- res %>%
    mutate(spatialFDR = p_adj,                        # 공간보정 없음 → adj.p를 보수적으로 사용
           is_sig = (p_adj < fdr_thresh) & (abs(logFC) > logfc_thresh)) %>%
    arrange(p_adj)

  cluster_summary <- res %>% filter(is_sig) %>%
    group_by(major_cluster) %>%
    summarise(n_sig_nhoods = n(),
              mean_logFC = mean(logFC, na.rm=TRUE),
              min_p_adj  = min(p_adj, na.rm=TRUE),
              .groups="drop") %>%
    arrange(desc(n_sig_nhoods))

  list(
    mode                  = "glmm_random_effects",
    formula               = deparse(full_formula),
    milo                  = milo,
    da_nhood_with_cluster = res,
    da_by_cluster         = cluster_summary
  )
}


#' MILO neighborhood DA with fixed/random effects and fast counting
#'
#' @description
#' Runs neighborhood-level differential abundance using MiloR. Provides two modes:
#' (A) edgeR GLM with spatial FDR (fixed effects only), or
#' (B) GLMM via glmmTMB with random effects (no spatial FDR; uses adj.p conservatively).
#' Also replaces \code{countCells()} by fast sparse multiplication for large data.
#'
#' @param sobj Seurat object (uses \code{RNA} assay by default).
#' @param target_var Character. Group variable to compare (e.g., "g3").
#' @param batch_var Character. Batch factor (e.g., "GEM").
#' @param patient_var Character. Sample/pseudobulk unit (e.g., "hos_no").
#' @param cluster_var Character. Cluster labels (e.g., "anno3.scvi").
#' @param reduced_dim One of "PCA","UMAP","TSNE" to build graph on. Default "PCA".
#' @param k,d,prop Milo graph/neighborhood params.
#' @param fixed_effects Character vector. Defaults to c(target_var, batch_var).
#' @param random_effects Character vector. If length > 0, GLMM mode is used.
#' @param glmm_family "nb" or "poisson". Default "nb".
#' @param p_adj_method P value adjustment method. Default "BH".
#' @param fdr_thresh,logfc_thresh Thresholds for significance flags.
#' @param seed RNG seed.
#' @param parallel_profile One of c("auto","rstudio","terminal"). Passed to \code{set_parallel_profile()}.
#' @param workers Integer parallel workers for GLMM path.
#' @return A list with components:
#'   \item{mode}{ "edgeR_fixed_effects" or "glmm_random_effects" }
#'   \item{milo}{ Milo object }
#'   \item{da_nhood_with_cluster}{ per-neighborhood results }
#'   \item{da_by_cluster}{ cluster-level summary }
#' @examples
#' set_parallel_profile("auto", workers=8, blas_threads=1)
#' out <- milo_neighborhood_da_v2(
#'   sobj = data_seurat,
#'   target_var="g3", batch_var="GEM", patient_var="hos_no", cluster_var="anno3.scvi",
#'   random_effects=c("hos_no")
#' )
#' @export
milo_neighborhood_da_v2 <- function(
  sobj,
  target_var   = "g3",
  batch_var    = "GEM",
  patient_var  = "hos_no",
  cluster_var  = "anno3.scvi",
  reduced_dim  = "PCA",
  k            = 30,
  d            = 50,
  prop         = 0.1,
  fixed_effects  = NULL,
  random_effects = NULL,
  glmm_family    = c("nb","poisson"),
  p_adj_method   = "BH",
  fdr_thresh     = 0.1,
  logfc_thresh   = 0.5,
  seed           = 1,
  parallel_profile = c("auto","rstudio","terminal"),
  workers          = max(1, parallel::detectCores()-1)
){
  suppressPackageStartupMessages({
    library(Seurat)
    library(SingleCellExperiment)
    library(miloR)
    library(Matrix)
    library(dplyr)
    library(tibble)
    library(purrr)
  })
  set.seed(seed)
  parallel_profile <- match.arg(parallel_profile)
  set_parallel_profile(mode = parallel_profile, workers = workers, blas_threads = 1)

  stopifnot(all(c(target_var,batch_var,patient_var,cluster_var) %in% colnames(sobj@meta.data)))

  # ---------- DimRed ----------
  if (!reduced_dim %in% c("PCA","UMAP","TSNE")) reduced_dim <- "PCA"
  if (!"pca" %in% Reductions(sobj)) {
    DefaultAssay(sobj) <- "RNA"
    sobj <- NormalizeData(sobj, verbose=FALSE)
    sobj <- FindVariableFeatures(sobj, nfeatures=3000, verbose=FALSE)
    sobj <- ScaleData(sobj, verbose=FALSE)
    sobj <- RunPCA(sobj, npcs=max(50,d), verbose=FALSE)
  }
  if (!"umap" %in% Reductions(sobj)) {
    sobj <- RunUMAP(sobj, dims=1:min(30,ncol(Embeddings(sobj,"pca"))), reduction="pca", verbose=FALSE)
  }

  # ---------- Seurat -> SCE ----------
  sce <- as.SingleCellExperiment(sobj, assay="RNA")
  reducedDim(sce, "PCA")  <- Embeddings(sobj, "pca")
  if ("umap" %in% Reductions(sobj)) reducedDim(sce,"UMAP") <- Embeddings(sobj,"umap")

  # ---------- Milo Graph & Nhoods ----------
  milo <- Milo(sce)
  milo <- buildGraph(milo, k = k, d = d, reduced.dim = "PCA")
  milo <- makeNhoods(milo, prop = prop, k = k, refined = TRUE, reduced_dims = "PCA")

  # ---------- Meta prep ----------
  md <- as.data.frame(colData(milo))
  md[[patient_var]] <- factor(md[[patient_var]])
  md[[target_var]]  <- factor(md[[target_var]])
  md[[batch_var]]   <- factor(md[[batch_var]])
  md[[cluster_var]] <- factor(md[[cluster_var]])

  # ---------- FAST count (replace countCells) ----------
  counts <- fast_count_cells(milo, samples = md[[patient_var]]) # H x S

  # attach counts into milo object (so downstream milo funcs can find them)
  # nhoodCounts<- setter is not exported; most downstream steps do not REQUIRE setter
  # we'll keep 'counts' locally for GLM/GLMM; edgeR path can still use testNhoods if we rebuild DGEList

  if (is.null(fixed_effects)) fixed_effects <- c(target_var, batch_var)
  use_glmm <- length(random_effects) > 0
  glmm_family <- match.arg(glmm_family)

  if (!use_glmm) {
    # ---- edgeR GLM with spatial FDR (fixed effects only) ----
    # Reconstruct edgeR design at SAMPLE-level:
    design.df <- md[, unique(c(patient_var, fixed_effects)), drop=FALSE]
    rownames(design.df) <- colnames(milo)  # samples

    for (v in fixed_effects) design.df[[v]] <- droplevels(factor(design.df[[v]]))
    mm <- model.matrix(as.formula(paste("~", paste(fixed_effects, collapse = " + "))), data = design.df)

    # Build DGEList from our fast counts
    suppressPackageStartupMessages(library(edgeR))
    dge <- edgeR::DGEList(counts = counts)
    dge <- edgeR::calcNormFactors(dge)  # TMM

    fit <- edgeR::glmQLFit(dge, design = mm)
    # contrast: 첫 레벨 대비(필요시 makeContrasts로 확장)
    qlf <- edgeR::glmQLFTest(fit)
    da.res <- tibble::tibble(
      nhood = seq_len(nrow(counts)),
      logFC = qlf$table$logFC,
      pvalue = qlf$table$PValue
    )
    da.res$p.adjust <- p.adjust(da.res$pvalue, method = p_adj_method)

    # 공간 FDR 근사: miloR의 graphSpatialFDR는 내부 그래프를 활용.
    # 여기서는 milo 객체의 그래프를 이용해 동일 함수 호출 (pvalue 컬럼명 기대치를 맞춰줌)
    names(da.res)[names(da.res)=="pvalue"] <- "pval"
    da.res <- miloR::graphSpatialFDR(milo, da.res)

    # 이웃→클러스터 매핑
    nh_mat <- miloR::nhoods(milo)  # cells x H
    cell_clusters <- md[[cluster_var]]
    major_cluster <- apply(nh_mat, 2, function(idx){
      labs <- cell_clusters[as.logical(idx)]
      if (!length(labs)) NA else names(sort(table(labs), decreasing=TRUE))[1]
    })
    da.res$major_cluster <- major_cluster[da.res$nhood]
    da.res <- da.res %>% mutate(is_sig = spatialFDR < fdr_thresh & abs(logFC) > logfc_thresh)

    cluster_summary <- da.res %>% filter(is_sig) %>%
      group_by(major_cluster) %>%
      summarise(n_sig_nhoods = n(),
                mean_logFC = mean(logFC),
                min_spatialFDR = min(spatialFDR),
                .groups="drop") %>%
      arrange(desc(n_sig_nhoods))

    return(list(
      mode                  = "edgeR_fixed_effects_fastcount",
      milo                  = milo,
      da_nhood_with_cluster = da.res,
      da_by_cluster         = cluster_summary
    ))
  }

  # ---- GLMM with random effects ----
  suppressPackageStartupMessages({
    library(glmmTMB); library(broom.mixed); library(future.apply)
  })

  samples   <- colnames(counts)             # sample IDs (levels of patient_var)
  sample_df <- tibble::tibble(
    sample = samples,
    !!target_var  := droplevels(factor(md[match(samples, rownames(md)), target_var][[1]])),
    !!batch_var   := droplevels(factor(md[match(samples, rownames(md)), batch_var][[1]])),
    !!patient_var := droplevels(factor(md[match(samples, rownames(md)), patient_var][[1]]))
  )

  # offset: sample total cells (for rate)
  sample_total_cells <- Matrix::colSums(miloR::nhoods(milo)) # per-sample total cells
  size_factor <- sample_total_cells / mean(sample_total_cells)
  sample_df$log_offset <- as.numeric(log(size_factor + 1e-8))

  fam <- if (glmm_family == "nb") glmmTMB::nbinom2(link="log") else poisson(link="log")
  fixed_part  <- paste(fixed_effects, collapse = " + ")
  random_part <- paste(sprintf("(1|%s)", random_effects), collapse = " + ")
  rhs <- paste(c(fixed_part, random_part), collapse = " + ")
  full_formula <- as.formula(paste("count ~", rhs, "+ offset(log_offset)"))

  nh_ids <- seq_len(nrow(counts))
  fit_one <- function(i){
    y <- as.numeric(counts[i, ])
    dat <- dplyr::mutate(sample_df, count = y)
    if (sum(y) == 0) return(tibble(nhood = i, term = NA, estimate = NA, std.error = NA, statistic = NA, p.value = NA))
    fit <- try(glmmTMB::glmmTMB(full_formula, data = dat, family = fam), silent = TRUE)
    if (inherits(fit, "try-error")) return(tibble(nhood = i, term = NA, estimate = NA, std.error = NA, statistic = NA, p.value = NA))
    co <- try(broom.mixed::tidy(fit, effects = "fixed"), silent = TRUE)
    if (inherits(co, "try-error")) return(tibble(nhood = i, term = NA, estimate = NA, std.error = NA, statistic = NA, p.value = NA))
    co$nhood <- i
    as_tibble(co)
  }

  # 병렬
  oplan <- future::plan()
  on.exit(try(future::plan(oplan), silent=TRUE), add=TRUE)
  glmm_tab <- future.apply::future_lapply(nh_ids, fit_one)
  glmm_tab <- dplyr::bind_rows(glmm_tab)

  target_terms <- grep(paste0("^", target_var), unique(glmm_tab$term), value=TRUE)
  res <- glmm_tab %>%
    dplyr::filter(term %in% target_terms) %>%
    dplyr::group_by(nhood) %>%
    dplyr::slice_head(n=1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(p_adj = p.adjust(p.value, method = p_adj_method),
                  logFC = estimate)

  # 이웃→클러스터
  nh_mat <- miloR::nhoods(milo)
  cell_clusters <- md[[cluster_var]]
  major_cluster <- apply(nh_mat, 2, function(idx){
    labs <- cell_clusters[as.logical(idx)]
    if (!length(labs)) NA else names(sort(table(labs), decreasing=TRUE))[1]
  })
  res$major_cluster <- major_cluster[res$nhood]

  res <- res %>%
    mutate(spatialFDR = p_adj,  # 보수적으로 adj.p 사용
           is_sig = (p_adj < fdr_thresh) & (abs(logFC) > logfc_thresh)) %>%
    arrange(p_adj)

  cluster_summary <- res %>% filter(is_sig) %>%
    group_by(major_cluster) %>%
    summarise(n_sig_nhoods = n(),
              mean_logFC = mean(logFC, na.rm=TRUE),
              min_p_adj  = min(p_adj, na.rm=TRUE),
              .groups="drop") %>%
    arrange(desc(n_sig_nhoods))

  list(
    mode                  = "glmm_random_effects_fastcount",
    formula               = deparse(full_formula),
    milo                  = milo,
    da_nhood_with_cluster = res,
    da_by_cluster         = cluster_summary
  )
}

#' MILO neighborhood DA (self-contained v3)
#'
#' Neighborhood-level differential abundance with two modes:
#' (A) edgeR GLM + spatial FDR (fixed effects only),
#' (B) GLMM via glmmTMB with random effects (adj.p used; no spatial FDR).
#' Internally replaces countCells() by a fast sparse-matrix method.
#'
#' @param sobj Seurat object (uses RNA assay by default).
#' @param target_var Character: group variable to compare (e.g., "g3").
#' @param batch_var Character: batch factor (e.g., "GEM").
#' @param patient_var Character: pseudobulk unit / sample ID (e.g., "hos_no").
#' @param cluster_var Character: cluster labels (e.g., "anno3.scvi").
#' @param reduced_dim "PCA","UMAP","TSNE" (graph built on PCA internally).
#' @param k,d,prop Milo graph/neighborhood params (defaults: 30,50,0.1).
#' @param fixed_effects Character vector; defaults to c(target_var, batch_var).
#' @param random_effects Character vector; if length>0, GLMM mode is used.
#' @param glmm_family "nb" (default) or "poisson".
#' @param p_adj_method P-value adjustment method (default "BH").
#' @param fdr_thresh,logfc_thresh Thresholds for significance flags.
#' @param seed RNG seed.
#' @param use_future Logical; if TRUE and GLMM mode, use future.apply for parallel fits.
#' @param workers Integer workers for future (only if use_future=TRUE).
#' @return list(mode, milo, da_nhood_with_cluster, da_by_cluster, formula(optional))
#' @examples
#' out <- milo_neighborhood_da_v3(
#'   sobj = data_seurat,
#'   target_var="g3", batch_var="GEM", patient_var="hos_no", cluster_var="anno3.scvi",
#'   fixed_effects=c("g3","GEM")                 # edgeR 고정효과 경로
#' )
#' out2 <- milo_neighborhood_da_v3(
#'   sobj = data_seurat,
#'   target_var="g3", batch_var="GEM", patient_var="hos_no", cluster_var="anno3.scvi",
#'   fixed_effects=c("g3","GEM"), random_effects=c("hos_no"), # GLMM 경로
#'   use_future=TRUE, workers=16
#' )
milo_neighborhood_da_v3 <- function(
  sobj,
  target_var   = "g3",
  batch_var    = "GEM",
  patient_var  = "hos_no",
  cluster_var  = "anno3.scvi",
  reduced_dim  = "PCA",
  k            = 30,
  d            = 50,
  prop         = 0.1,
  fixed_effects  = NULL,
  random_effects = NULL,
  glmm_family    = c("nb","poisson"),
  p_adj_method   = "BH",
  fdr_thresh     = 0.1,
  logfc_thresh   = 0.5,
  seed           = 1,
  use_future     = FALSE,
  workers        = max(1, parallel::detectCores()-1)
){
  suppressPackageStartupMessages({
    library(Seurat)
    library(SingleCellExperiment)
    library(miloR)
    library(Matrix)
    library(dplyr)
    library(tibble)
  })
  set.seed(seed)

  # ---- basic checks ----
  req <- c(target_var, batch_var, patient_var, cluster_var)
  stopifnot(all(req %in% colnames(sobj@meta.data)))
  if (!reduced_dim %in% c("PCA","UMAP","TSNE")) reduced_dim <- "PCA"
  if (is.null(fixed_effects)) fixed_effects <- c(target_var, batch_var)
  glmm_family <- match.arg(glmm_family)
  use_glmm <- length(random_effects) > 0

  # ---- dimred (PCA mandatory for Milo graph) ----
  if (!"pca" %in% Reductions(sobj)) {
    DefaultAssay(sobj) <- "RNA"
    sobj <- NormalizeData(sobj, verbose=FALSE)
    sobj <- FindVariableFeatures(sobj, nfeatures=3000, verbose=FALSE)
    sobj <- ScaleData(sobj, verbose=FALSE)
    sobj <- RunPCA(sobj, npcs = max(50, d), verbose=FALSE)
  }
  if (!"umap" %in% Reductions(sobj)) {
    sobj <- RunUMAP(sobj, dims=1:min(30, ncol(Embeddings(sobj,"pca"))), reduction="pca", verbose=FALSE)
  }

  # ---- Seurat -> SCE -> Milo ----
  sce <- as.SingleCellExperiment(sobj, assay="RNA")
  reducedDim(sce, "PCA")  <- Embeddings(sobj, "pca")
  if ("umap" %in% Reductions(sobj)) reducedDim(sce, "UMAP") <- Embeddings(sobj, "umap")

  milo <- Milo(sce)
  milo <- buildGraph(milo, k = k, d = d, reduced.dim = "PCA")
  milo <- makeNhoods(milo, prop = prop, k = k, refined = TRUE, reduced_dims = "PCA")

  # ---- metadata ----
  md <- as.data.frame(colData(milo))
  md[[patient_var]] <- factor(md[[patient_var]])
  md[[target_var]]  <- factor(md[[target_var]])
  md[[batch_var]]   <- factor(md[[batch_var]])
  md[[cluster_var]] <- factor(md[[cluster_var]])

  # ---- super-fast counts HxS ----
  counts <- .fast_count_cells(milo, samples = md[[patient_var]])  # H x S
  H <- nrow(counts)

  # ---- helper: majority cluster per nhood ----
  nh_mat <- miloR::nhoods(milo)   # cells x H
  cell_clusters <- md[[cluster_var]]
  majority_cluster <- apply(nh_mat, 2, function(idx){
    labs <- cell_clusters[as.logical(idx)]
    if (!length(labs)) NA else names(sort(table(labs), decreasing = TRUE))[1]
  })

  if (!use_glmm) {
    # =========================
    # (A) edgeR GLM + spatial FDR
    # =========================
    suppressPackageStartupMessages(library(edgeR))
    design.df <- md[, unique(c(patient_var, fixed_effects)), drop=FALSE]
    rownames(design.df) <- colnames(milo)  # samples
    for (v in fixed_effects) design.df[[v]] <- droplevels(factor(design.df[[v]]))

    mm <- model.matrix(as.formula(paste("~", paste(fixed_effects, collapse = " + "))), data = design.df)

    dge <- edgeR::DGEList(counts = counts)
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    fit <- edgeR::glmQLFit(dge, design = mm)
    qlf <- edgeR::glmQLFTest(fit)  # default: overall test on last coef; 필요시 makeContrasts 사용

    da.res <- tibble(
      nhood = seq_len(H),
      logFC = qlf$table$logFC,
      pval  = qlf$table$PValue
    )

    # spatial FDR using Milo graph
    da.res <- miloR::graphSpatialFDR(milo, da.res)

    da.res$major_cluster <- majority_cluster[da.res$nhood]
    da.res <- da.res %>%
      mutate(is_sig = spatialFDR < fdr_thresh & abs(logFC) > logfc_thresh) %>%
      arrange(spatialFDR, desc(abs(logFC)))

    cluster_summary <- da.res %>% filter(is_sig) %>%
      group_by(major_cluster) %>%
      summarise(
        n_sig_nhoods = n(),
        mean_logFC = mean(logFC),
        min_spatialFDR = min(spatialFDR),
        .groups="drop"
      ) %>% arrange(desc(n_sig_nhoods))

    return(list(
      mode                  = "edgeR_fixed_effects_fastcount",
      milo                  = milo,
      da_nhood_with_cluster = da.res,
      da_by_cluster         = cluster_summary
    ))
  }

  # =========================
  # (B) GLMM with random effects
  # =========================
  suppressPackageStartupMessages({
    library(glmmTMB)
    library(broom.mixed)
  })

  # sample-level frame aligned to counts columns (samples)
  samples <- colnames(counts)
  sample_df <- tibble(
    sample = samples,
    !!target_var  := droplevels(factor(md[match(samples, rownames(md)), target_var][[1]])),
    !!batch_var   := droplevels(factor(md[match(samples, rownames(md)), batch_var][[1]])),
    !!patient_var := droplevels(factor(md[match(samples, rownames(md)), patient_var][[1]]))
  )

  # offset: per-sample total cells across all nhoods
  sample_total_cells <- Matrix::colSums(miloR::nhoods(milo))
  size_factor <- sample_total_cells / mean(sample_total_cells)
  sample_df$log_offset <- as.numeric(log(size_factor + 1e-8))

  fam <- if (glmm_family == "nb") glmmTMB::nbinom2(link="log") else poisson(link="log")
  fixed_part  <- paste(fixed_effects, collapse = " + ")
  random_part <- paste(sprintf("(1|%s)", random_effects), collapse = " + ")
  rhs <- paste(c(fixed_part, random_part), collapse = " + ")
  full_formula <- as.formula(paste("count ~", rhs, "+ offset(log_offset)"))

  # optional parallel with future.apply
  if (use_future) {
    oplan <- NULL
    if (requireNamespace("future", quietly = TRUE)) {
      oplan <- future::plan()
      # multisession은 RStudio/터미널 어디서나 안전
      future::plan(future::multisession, workers = workers)
    }
    on.exit({
      if (!is.null(oplan)) try(future::plan(oplan), silent = TRUE)
    }, add = TRUE)
    suppressPackageStartupMessages(library(future.apply))
    map_fun <- future.apply::future_lapply
  } else {
    map_fun <- lapply
  }

  fit_one <- function(i){
    y <- as.numeric(counts[i, ])
    dat <- dplyr::mutate(sample_df, count = y)
    if (sum(y) == 0) {
      return(tibble(nhood = i, term = NA_character_, estimate = NA_real_, std.error = NA_real_,
                    statistic = NA_real_, p.value = NA_real_))
    }
    fit <- try(glmmTMB::glmmTMB(full_formula, data = dat, family = fam), silent = TRUE)
    if (inherits(fit, "try-error")) {
      return(tibble(nhood = i, term = NA_character_, estimate = NA_real_, std.error = NA_real_,
                    statistic = NA_real_, p.value = NA_real_))
    }
    co <- try(broom.mixed::tidy(fit, effects = "fixed"), silent = TRUE)
    if (inherits(co, "try-error")) {
      return(tibble(nhood = i, term = NA_character_, estimate = NA_real_, std.error = NA_real_,
                    statistic = NA_real_, p.value = NA_real_))
    }
    co$nhood <- i
    as_tibble(co)
  }

  glmm_tab <- map_fun(seq_len(H), fit_one)
  glmm_tab <- dplyr::bind_rows(glmm_tab)

  target_terms <- grep(paste0("^", target_var), unique(glmm_tab$term), value=TRUE)
  res <- glmm_tab %>%
    dplyr::filter(term %in% target_terms) %>%
    dplyr::group_by(nhood) %>%
    dplyr::slice_head(n=1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      p_adj = p.adjust(p.value, method = p_adj_method),
      logFC = estimate,
      major_cluster = majority_cluster[nhood],
      spatialFDR = p_adj,  # 보수적으로 adj.p 사용
      is_sig = (p_adj < fdr_thresh) & (abs(logFC) > logfc_thresh)
    ) %>% arrange(p_adj)

  cluster_summary <- res %>% filter(is_sig) %>%
    group_by(major_cluster) %>%
    summarise(
      n_sig_nhoods = n(),
      mean_logFC   = mean(logFC, na.rm=TRUE),
      min_p_adj    = min(p_adj, na.rm=TRUE),
      .groups="drop"
    ) %>% arrange(desc(n_sig_nhoods))

  list(
    mode                  = "glmm_random_effects_fastcount",
    formula               = deparse(full_formula),
    milo                  = milo,
    da_nhood_with_cluster = res,
    da_by_cluster         = cluster_summary
  )
}


#' MILO neighborhood DA (flexible reductions; v4)
#'
#' Neighborhood DA with:
#' - (A) edgeR GLM + spatial FDR (fixed effects only), or
#' - (B) GLMM via glmmTMB (random effects; adj.p as conservative spatialFDR).
#' Uses arbitrary Seurat reductions (e.g., "integrated.scvi", "umap.scvi").
#'
#' @param sobj Seurat object.
#' @param target_var,batch_var,patient_var,cluster_var Character; column names in meta.data.
#' @param graph_reduction Character; Seurat reduction name used for Milo graph (e.g., "pca","integrated.scvi"). Required.
#' @param layout_reduction Character or NULL; reduction for plotting/positions (e.g., "umap","umap.scvi"). Optional.
#' @param k,d,prop Milo graph/neighborhood params. If \code{d=NULL}, uses all columns of \code{graph_reduction}.
#' @param fixed_effects Character vector; defaults to c(target_var, batch_var).
#' @param random_effects Character vector; if length>0, GLMM mode.
#' @param glmm_family "nb" (default) or "poisson".
#' @param p_adj_method "BH" etc.
#' @param fdr_thresh,logfc_thresh Thresholds for significance.
#' @param seed RNG seed.
#' @param use_future Logical; GLMM 병렬 시 \code{future.apply} 사용.
#' @param workers Integer workers for future (if \code{use_future=TRUE}).
#' @return list(mode, milo, da_nhood_with_cluster, da_by_cluster, formula(optional))
#' @examples
#' # scVI 임베딩 사용
#' out <- milo_neighborhood_da_v4(
#'   sobj = data_seurat,
#'   target_var="g3", batch_var="GEM", patient_var="hos_no", cluster_var="anno3.scvi",
#'   graph_reduction="integrated.scvi", layout_reduction="umap.scvi",
#'   fixed_effects=c("g3","GEM")
#' )
milo_neighborhood_da_v4 <- function(
  sobj,
  target_var,
  batch_var,
  patient_var,
  cluster_var,
  graph_reduction,                 # e.g., "integrated.scvi" (required)
  layout_reduction = NULL,         # e.g., "umap.scvi"
  k = 30,
  d = NULL,                        # if NULL -> ncol(Embeddings(graph_reduction))
  prop = 0.1,
  fixed_effects  = NULL,
  random_effects = NULL,
  glmm_family    = c("nb","poisson"),
  p_adj_method   = "BH",
  fdr_thresh     = 0.1,
  logfc_thresh   = 0.5,
  seed           = 1,
  use_future     = FALSE,
  workers        = max(1, parallel::detectCores()-1)
){
  suppressPackageStartupMessages({
    library(Seurat)
    library(SingleCellExperiment)
    library(miloR)
    library(Matrix)
    library(dplyr)
    library(tibble)
  })
  set.seed(seed)

  # ---- checks ----
  req_cols <- c(target_var, batch_var, patient_var, cluster_var)
  stopifnot(all(req_cols %in% colnames(sobj@meta.data)))
  stopifnot(graph_reduction %in% Reductions(sobj))
  if (!is.null(layout_reduction)) stopifnot(layout_reduction %in% Reductions(sobj))
  if (is.null(fixed_effects)) fixed_effects <- c(target_var, batch_var)
  glmm_family <- match.arg(glmm_family)
  use_glmm <- length(random_effects) > 0

  # ---- pull embeddings flexibly ----
  G_emb <- Embeddings(sobj, reduction = graph_reduction)   # cells x dims_G
  if (is.null(d)) d <- ncol(G_emb)
  if (d > ncol(G_emb)) d <- ncol(G_emb)
  G_use <- G_emb[, seq_len(d), drop = FALSE]

  L_use <- NULL
  if (!is.null(layout_reduction)) {
    L_emb <- Embeddings(sobj, reduction = layout_reduction)
    L_use <- L_emb
  }

  # ---- SCE & reducedDims ----
  sce <- as.SingleCellExperiment(sobj, assay = DefaultAssay(sobj))
  reducedDim(sce, "GRAPH") <- as.matrix(G_use)
  if (!is.null(L_use)) reducedDim(sce, "LAYOUT") <- as.matrix(L_use)

  # ---- Milo graph/neighborhoods ----
  milo <- Milo(sce)
  milo <- buildGraph(milo, k = k, d = ncol(G_use), reduced.dim = "GRAPH")
  milo <- makeNhoods(milo, prop = prop, k = k, refined = TRUE, reduced_dims = "GRAPH")

  # ---- metadata ----
  md <- as.data.frame(colData(milo))
  md[[patient_var]] <- factor(md[[patient_var]])
  md[[target_var]]  <- factor(md[[target_var]])
  md[[batch_var]]   <- factor(md[[batch_var]])
  md[[cluster_var]] <- factor(md[[cluster_var]])

  # ---- fast counts (H x S) ----
  counts <- .fast_count_cells(milo, samples = md[[patient_var]])
  H <- nrow(counts)

  # ---- majority cluster per nhood ----
  nh_mat <- miloR::nhoods(milo)  # cells x H
  cell_clusters <- md[[cluster_var]]
  majority_cluster <- apply(nh_mat, 2, function(idx){
    labs <- cell_clusters[as.logical(idx)]
    if (!length(labs)) NA else names(sort(table(labs), decreasing = TRUE))[1]
  })

  if (!use_glmm) {
    # ===== (A) edgeR GLM + spatial FDR =====
    suppressPackageStartupMessages(library(edgeR))
    design.df <- md[, unique(c(patient_var, fixed_effects)), drop = FALSE]
    rownames(design.df) <- colnames(milo)
    for (v in fixed_effects) design.df[[v]] <- droplevels(factor(design.df[[v]]))

    mm <- model.matrix(as.formula(paste("~", paste(fixed_effects, collapse = " + "))), data = design.df)
    dge <- edgeR::DGEList(counts = counts)
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    fit <- edgeR::glmQLFit(dge, design = mm)
    qlf <- edgeR::glmQLFTest(fit)

    da.res <- tibble(
      nhood = seq_len(H),
      logFC = qlf$table$logFC,
      pval  = qlf$table$PValue
    )
    da.res <- miloR::graphSpatialFDR(milo, da.res)  # uses Milo graph (GRAPH)

    da.res$major_cluster <- majority_cluster[da.res$nhood]
    da.res <- da.res %>%
      mutate(is_sig = spatialFDR < fdr_thresh & abs(logFC) > logfc_thresh) %>%
      arrange(spatialFDR, desc(abs(logFC)))

    cluster_summary <- da.res %>% filter(is_sig) %>%
      group_by(major_cluster) %>%
      summarise(n_sig_nhoods = n(),
                mean_logFC = mean(logFC),
                min_spatialFDR = min(spatialFDR),
                .groups="drop") %>%
      arrange(desc(n_sig_nhoods))

    return(list(
      mode                  = "edgeR_fixed_effects_fastcount",
      milo                  = milo,
      da_nhood_with_cluster = da.res,
      da_by_cluster         = cluster_summary
    ))
  }

  # ===== (B) GLMM with random effects =====
  suppressPackageStartupMessages({
    library(glmmTMB)
    library(broom.mixed)
  })
  samples <- colnames(counts)
  # align sample-level covariates to 'counts' columns
  sample_df <- tibble(
    sample = samples,
    !!target_var  := droplevels(factor(md[match(samples, rownames(md)), target_var][[1]])),
    !!batch_var   := droplevels(factor(md[match(samples, rownames(md)), batch_var][[1]])),
    !!patient_var := droplevels(factor(md[match(samples, rownames(md)), patient_var][[1]]))
  )
  # offset: per-sample total cells (across all nhoods)
  sample_total_cells <- Matrix::colSums(miloR::nhoods(milo))
  size_factor <- sample_total_cells / mean(sample_total_cells)
  sample_df$log_offset <- as.numeric(log(size_factor + 1e-8))

  fam <- if (glmm_family == "nb") glmmTMB::nbinom2(link="log") else poisson(link="log")
  fixed_part  <- paste(fixed_effects, collapse = " + ")
  random_part <- paste(sprintf("(1|%s)", random_effects), collapse = " + ")
  rhs <- paste(c(fixed_part, random_part), collapse = " + ")
  full_formula <- as.formula(paste("count ~", rhs, "+ offset(log_offset)"))

  # optional parallel via future.apply
  if (use_future && requireNamespace("future", quietly = TRUE)) {
    oplan <- future::plan()
    on.exit(try(future::plan(oplan), silent = TRUE), add = TRUE)
    future::plan(future::multisession, workers = workers)
    suppressPackageStartupMessages(library(future.apply))
    map_fun <- future.apply::future_lapply
  } else {
    map_fun <- lapply
  }

  fit_one <- function(i){
    y <- as.numeric(counts[i, ])
    dat <- dplyr::mutate(sample_df, count = y)
    if (sum(y) == 0) {
      return(tibble(nhood = i, term = NA_character_, estimate = NA_real_, std.error = NA_real_,
                    statistic = NA_real_, p.value = NA_real_))
    }
    fit <- try(glmmTMB::glmmTMB(full_formula, data = dat, family = fam), silent = TRUE)
    if (inherits(fit, "try-error")) {
      return(tibble(nhood = i, term = NA_character_, estimate = NA_real_, std.error = NA_real_,
                    statistic = NA_real_, p.value = NA_real_))
    }
    co <- try(broom.mixed::tidy(fit, effects = "fixed"), silent = TRUE)
    if (inherits(co, "try-error")) {
      return(tibble(nhood = i, term = NA_character_, estimate = NA_real_, std.error = NA_real_,
                    statistic = NA_real_, p.value = NA_real_))
    }
    co$nhood <- i
    as_tibble(co)
  }

  glmm_tab <- map_fun(seq_len(H), fit_one) |> dplyr::bind_rows()

  target_terms <- grep(paste0("^", target_var), unique(glmm_tab$term), value = TRUE)
  res <- glmm_tab %>%
    dplyr::filter(term %in% target_terms) %>%
    dplyr::group_by(nhood) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      p_adj = p.adjust(p.value, method = p_adj_method),
      logFC = estimate,
      major_cluster = majority_cluster[nhood],
      spatialFDR = p_adj,                                # 보수적으로 adj.p 사용
      is_sig = (p_adj < fdr_thresh) & (abs(logFC) > logfc_thresh)
    ) %>% arrange(p_adj)

  cluster_summary <- res %>% filter(is_sig) %>%
    group_by(major_cluster) %>%
    summarise(n_sig_nhoods = n(),
              mean_logFC   = mean(logFC, na.rm = TRUE),
              min_p_adj    = min(p_adj, na.rm = TRUE),
              .groups="drop") %>%
    arrange(desc(n_sig_nhoods))

  list(
    mode                  = "glmm_random_effects_fastcount",
    formula               = deparse(full_formula),
    milo                  = milo,
    da_nhood_with_cluster = res,
    da_by_cluster         = cluster_summary
  )
}


milo_opus1 <- function(
        sobj,
        target_var,
        batch_var,
        patient_var,
        cluster_var,
        graph_reduction,                 # e.g., "integrated.scvi" (required)
        layout_reduction = NULL,         # e.g., "umap.scvi"
        k = 30,
        d = NULL,                        # if NULL -> ncol(Embeddings(graph_reduction))
        prop = 0.1,
        fixed_effects  = NULL,
        random_effects = NULL,
        glmm_family    = c("nb","poisson"),
        p_adj_method   = "BH",
        fdr_thresh     = 0.1,
        logfc_thresh   = 0.5,
        seed           = 1,
        use_future     = FALSE,
        workers        = max(1, parallel::detectCores()-1),
        verbose        = TRUE,
        min_cells      = 3                # minimum cells per neighborhood
){
    suppressPackageStartupMessages({
        library(Seurat)
        library(SingleCellExperiment)
        library(miloR)
        library(Matrix)
        library(dplyr)
        library(tibble)
    })
    set.seed(seed)
    
    # ---- Verbose helper ----
    vcat <- function(...) if(verbose) cat(..., "\n")
    
    # ---- checks ----
    req_cols <- c(target_var, batch_var, patient_var, cluster_var)
    stopifnot("Required columns missing from metadata" = all(req_cols %in% colnames(sobj@meta.data)))
    stopifnot("Graph reduction not found" = graph_reduction %in% Reductions(sobj))
    if (!is.null(layout_reduction)) {
        stopifnot("Layout reduction not found" = layout_reduction %in% Reductions(sobj))
    }
    if (is.null(fixed_effects)) fixed_effects <- c(target_var, batch_var)
    glmm_family <- match.arg(glmm_family)
    use_glmm <- length(random_effects) > 0
    
    vcat("Starting Milo DA analysis...")
    vcat("Mode:", ifelse(use_glmm, "GLMM with random effects", "edgeR fixed effects"))
    
    # ---- pull embeddings flexibly ----
    G_emb <- Embeddings(sobj, reduction = graph_reduction)   # cells x dims_G
    if (is.null(d)) d <- ncol(G_emb)
    if (d > ncol(G_emb)) d <- ncol(G_emb)
    G_use <- G_emb[, seq_len(d), drop = FALSE]
    
    L_use <- NULL
    if (!is.null(layout_reduction)) {
        L_emb <- Embeddings(sobj, reduction = layout_reduction)
        L_use <- L_emb
    }
    
    vcat(sprintf("Using %d dimensions from %s reduction", d, graph_reduction))
    
    # ---- SCE & reducedDims ----
    sce <- as.SingleCellExperiment(sobj, assay = DefaultAssay(sobj))
    reducedDim(sce, "GRAPH") <- as.matrix(G_use)
    if (!is.null(L_use)) reducedDim(sce, "LAYOUT") <- as.matrix(L_use)
    
    # ---- Milo graph/neighborhoods ----
    vcat("Building kNN graph and neighborhoods...")
    milo <- Milo(sce)
    milo <- buildGraph(milo, k = k, d = ncol(G_use), reduced.dim = "GRAPH")
    milo <- makeNhoods(milo, prop = prop, k = k, refined = TRUE, reduced_dims = "GRAPH")
    
    vcat(sprintf("Created %d neighborhoods", ncol(nhoods(milo))))
    
    # ---- metadata ----
    md <- as.data.frame(colData(milo))
    md[[patient_var]] <- factor(md[[patient_var]])
    md[[target_var]]  <- factor(md[[target_var]])
    md[[batch_var]]   <- factor(md[[batch_var]])
    md[[cluster_var]] <- factor(md[[cluster_var]])
    
    # ---- fast counts (H x S) ----
    vcat("Computing neighborhood-by-sample counts...")
    counts <- .fast_count_cells(milo, samples = md[[patient_var]])
    H <- nrow(counts)
    
    # Filter out neighborhoods with too few cells
    nh_sizes <- rowSums(counts)
    keep_nh <- nh_sizes >= min_cells
    if (sum(!keep_nh) > 0) {
        vcat(sprintf("Filtering out %d neighborhoods with < %d cells", sum(!keep_nh), min_cells))
        counts <- counts[keep_nh, , drop = FALSE]
    }
    
    # ---- majority cluster per nhood ----
    nh_mat <- miloR::nhoods(milo)  # cells x H
    cell_clusters <- md[[cluster_var]]
    majority_cluster <- apply(nh_mat, 2, function(idx){
        labs <- cell_clusters[as.logical(idx)]
        if (!length(labs)) NA else names(sort(table(labs), decreasing = TRUE))[1]
    })
    
    if (!use_glmm) {
        # ===== (A) edgeR GLM + spatial FDR =====
        vcat("Running edgeR analysis...")
        suppressPackageStartupMessages(library(edgeR))
        
        # Create design data frame - one row per sample
        samples <- colnames(counts)
        design.df <- md[!duplicated(md[[patient_var]]), , drop = FALSE]
        rownames(design.df) <- design.df[[patient_var]]
        design.df <- design.df[samples, unique(c(patient_var, fixed_effects)), drop = FALSE]
        
        # Ensure factors are clean
        for (v in fixed_effects) {
            design.df[[v]] <- droplevels(factor(design.df[[v]]))
        }
        
        # Create design matrix
        mm <- tryCatch({
            model.matrix(as.formula(paste("~", paste(fixed_effects, collapse = " + "))), 
                         data = design.df)
        }, error = function(e) {
            stop("Failed to create design matrix: ", e$message, 
                 "\nCheck for singular design or missing levels")
        })
        
        # Check for full rank
        if (qr(mm)$rank < ncol(mm)) {
            warning("Design matrix is not full rank. Some coefficients may not be estimable.")
        }
        
        # EdgeR workflow with proper dispersion estimation
        dge <- edgeR::DGEList(counts = counts[keep_nh, , drop = FALSE])
        dge <- edgeR::calcNormFactors(dge, method = "TMM")
        
        # CRITICAL: Estimate dispersions before fitting
        vcat("Estimating dispersions...")
        dge <- edgeR::estimateDisp(dge, design = mm, robust = TRUE)
        
        # Now fit the model
        vcat("Fitting GLM...")
        fit <- edgeR::glmQLFit(dge, design = mm, robust = TRUE)
        
        # Test for differential abundance
        # Find the coefficient for target_var
        target_coef <- grep(paste0("^", target_var), colnames(mm), value = TRUE)[1]
        if (is.na(target_coef)) {
            stop("Could not find target variable coefficient in design matrix")
        }
        
        qlf <- edgeR::glmQLFTest(fit, coef = target_coef)
        
        # Build results
        nh_idx <- which(keep_nh)
        da.res <- tibble(
            nhood = nh_idx,
            logFC = qlf$table$logFC,
            pval  = qlf$table$PValue
        )
        
        vcat("Computing spatial FDR...")
        da.res <- miloR::graphSpatialFDR(milo, da.res)
        
        da.res$major_cluster <- majority_cluster[da.res$nhood]
        da.res <- da.res %>%
            mutate(is_sig = spatialFDR < fdr_thresh & abs(logFC) > logfc_thresh) %>%
            arrange(spatialFDR, desc(abs(logFC)))
        
        cluster_summary <- da.res %>% 
            filter(is_sig) %>%
            group_by(major_cluster) %>%
            summarise(
                n_sig_nhoods = n(),
                mean_logFC = mean(logFC),
                min_spatialFDR = min(spatialFDR),
                .groups = "drop"
            ) %>%
            arrange(desc(n_sig_nhoods))
        
        vcat(sprintf("Found %d significant neighborhoods", sum(da.res$is_sig)))
        
        return(list(
            mode                  = "edgeR_fixed_effects",
            milo                  = milo,
            da_nhood_with_cluster = da.res,
            da_by_cluster         = cluster_summary,
            design_matrix         = mm,
            dge                   = dge
        ))
    }
    
    # ===== (B) GLMM with random effects =====
    vcat("Running GLMM analysis...")
    suppressPackageStartupMessages({
        library(glmmTMB)
        library(broom.mixed)
    })
    
    samples <- colnames(counts)
    # align sample-level covariates to 'counts' columns
    sample_info <- md[!duplicated(md[[patient_var]]), , drop = FALSE]
    rownames(sample_info) <- sample_info[[patient_var]]
    
    sample_df <- tibble(
        sample = samples,
        !!target_var  := droplevels(factor(sample_info[samples, target_var])),
        !!batch_var   := droplevels(factor(sample_info[samples, batch_var])),
        !!patient_var := droplevels(factor(samples))
    )
    
    # Add any additional random effects variables
    for (re in random_effects) {
        if (re != patient_var) {
            sample_df[[re]] <- droplevels(factor(sample_info[samples, re]))
        }
    }
    
    # Compute offset
    sample_total_cells <- colSums(counts)
    sample_df$log_offset <- log(sample_total_cells + 1)
    
    fam <- if (glmm_family == "nb") glmmTMB::nbinom2(link = "log") else poisson(link = "log")
    fixed_part  <- paste(fixed_effects, collapse = " + ")
    random_part <- paste(sprintf("(1|%s)", random_effects), collapse = " + ")
    rhs <- paste(c(fixed_part, random_part), collapse = " + ")
    full_formula <- as.formula(paste("count ~", rhs, "+ offset(log_offset)"))
    
    vcat(sprintf("GLMM formula: %s", deparse(full_formula)))
    
    # Set up parallel processing
    if (use_future && requireNamespace("future", quietly = TRUE)) {
        oplan <- future::plan()
        on.exit(try(future::plan(oplan), silent = TRUE), add = TRUE)
        future::plan(future::multisession, workers = workers)
        suppressPackageStartupMessages(library(future.apply))
        map_fun <- future.apply::future_lapply
        vcat(sprintf("Using parallel processing with %d workers", workers))
    } else {
        map_fun <- lapply
    }
    
    # Fit models
    nh_to_fit <- which(keep_nh)
    n_to_fit <- length(nh_to_fit)
    
    fit_one <- function(idx) {
        i <- nh_to_fit[idx]
        y <- as.numeric(counts[idx, ])
        dat <- dplyr::mutate(sample_df, count = y)
        
        if (sum(y) == 0) {
            return(tibble(
                nhood = i, term = NA_character_, estimate = NA_real_, 
                std.error = NA_real_, statistic = NA_real_, p.value = NA_real_
            ))
        }
        
        fit <- tryCatch({
            glmmTMB::glmmTMB(full_formula, data = dat, family = fam,
                             control = glmmTMBControl(optimizer = nlminb))
        }, error = function(e) NULL, warning = function(w) NULL)
        
        if (is.null(fit)) {
            return(tibble(
                nhood = i, term = NA_character_, estimate = NA_real_,
                std.error = NA_real_, statistic = NA_real_, p.value = NA_real_
            ))
        }
        
        co <- tryCatch({
            broom.mixed::tidy(fit, effects = "fixed")
        }, error = function(e) NULL)
        
        if (is.null(co)) {
            return(tibble(
                nhood = i, term = NA_character_, estimate = NA_real_,
                std.error = NA_real_, statistic = NA_real_, p.value = NA_real_
            ))
        }
        
        co$nhood <- i
        as_tibble(co)
    }
    
    vcat(sprintf("Fitting %d GLMM models...", n_to_fit))
    glmm_tab <- map_fun(seq_len(n_to_fit), fit_one) %>% dplyr::bind_rows()
    
    # Extract results for target variable
    target_terms <- grep(paste0("^", target_var), unique(glmm_tab$term), value = TRUE)
    
    res <- glmm_tab %>%
        dplyr::filter(term %in% target_terms) %>%
        dplyr::group_by(nhood) %>%
        dplyr::slice_head(n = 1) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
            p_adj = p.adjust(p.value, method = p_adj_method),
            logFC = estimate,
            major_cluster = majority_cluster[nhood],
            spatialFDR = p_adj,
            is_sig = (p_adj < fdr_thresh) & (abs(logFC) > logfc_thresh)
        ) %>% 
        arrange(p_adj)
    
    cluster_summary <- res %>% 
        filter(is_sig) %>%
        group_by(major_cluster) %>%
        summarise(
            n_sig_nhoods = n(),
            mean_logFC = mean(logFC, na.rm = TRUE),
            min_p_adj = min(p_adj, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        arrange(desc(n_sig_nhoods))
    
    vcat(sprintf("Found %d significant neighborhoods", sum(res$is_sig, na.rm = TRUE)))
    
    list(
        mode                  = "glmm_random_effects",
        formula               = deparse(full_formula),
        milo                  = milo,
        da_nhood_with_cluster = res,
        da_by_cluster         = cluster_summary
    )
}

milo_opus2 <- function(
    sobj,
    target_var,
    batch_var,
    patient_var,
    cluster_var,
    graph_reduction,                 
    layout_reduction = NULL,         
    k = 30,
    d = NULL,                        
    prop = 0.1,
    fixed_effects  = NULL,
    random_effects = NULL,
    glmm_family    = c("nb","poisson"),
    p_adj_method   = "BH",
    fdr_thresh     = 0.1,
    logfc_thresh   = 0.5,
    seed           = 1,
    use_future     = FALSE,
    workers        = max(1, parallel::detectCores()-1),
    verbose        = TRUE,
    min_cells      = 3,
    force_simple_design = FALSE      # New option for nested designs
){
    suppressPackageStartupMessages({
        library(Seurat)
        library(SingleCellExperiment)
        library(miloR)
        library(Matrix)
        library(dplyr)
        library(tibble)
    })
    set.seed(seed)
    
    # ---- Verbose helper ----
    vcat <- function(...) if(verbose) cat(..., "\n")
    
    # ---- checks ----
    req_cols <- c(target_var, batch_var, patient_var, cluster_var)
    stopifnot("Required columns missing from metadata" = all(req_cols %in% colnames(sobj@meta.data)))
    stopifnot("Graph reduction not found" = graph_reduction %in% Reductions(sobj))
    if (!is.null(layout_reduction)) {
        stopifnot("Layout reduction not found" = layout_reduction %in% Reductions(sobj))
    }
    if (is.null(fixed_effects)) fixed_effects <- c(target_var, batch_var)
    glmm_family <- match.arg(glmm_family)
    use_glmm <- length(random_effects) > 0
    
    vcat("Starting Milo DA analysis...")
    vcat("Mode:", ifelse(use_glmm, "GLMM with random effects", "edgeR fixed effects"))
    
    # ---- Check for nested design ----
    md_check <- sobj@meta.data[, c(patient_var, target_var, batch_var), drop = FALSE]
    md_check <- unique(md_check)
    
    # Check if each patient has only one value for target and batch
    patients_per_target <- md_check %>% 
        group_by(!!sym(patient_var)) %>% 
        summarise(n_target = n_distinct(!!sym(target_var)), 
                  n_batch = n_distinct(!!sym(batch_var))) %>%
        ungroup()
    
    is_nested <- all(patients_per_target$n_target == 1) & all(patients_per_target$n_batch == 1)
    
    if (is_nested && !force_simple_design) {
        vcat("WARNING: Detected nested design - each patient has only one value for target and batch variables")
        vcat("This creates a confounded design. Adjusting analysis strategy...")
        
        # Check if batch is nested within target
        batch_target_table <- unique(md_check[, c(target_var, batch_var)])
        batch_nested_in_target <- nrow(batch_target_table) == length(unique(md_check[[batch_var]]))
        
        if (batch_nested_in_target) {
            vcat("Batch is nested within target variable. Removing batch from fixed effects.")
            fixed_effects <- setdiff(fixed_effects, batch_var)
            if (length(fixed_effects) == 0) fixed_effects <- target_var
        }
    }
    
    # ---- pull embeddings ----
    G_emb <- Embeddings(sobj, reduction = graph_reduction)
    if (is.null(d)) d <- ncol(G_emb)
    if (d > ncol(G_emb)) d <- ncol(G_emb)
    G_use <- G_emb[, seq_len(d), drop = FALSE]
    
    L_use <- NULL
    if (!is.null(layout_reduction)) {
        L_emb <- Embeddings(sobj, reduction = layout_reduction)
        L_use <- L_emb
    }
    
    vcat(sprintf("Using %d dimensions from %s reduction", d, graph_reduction))
    
    # ---- SCE & reducedDims ----
    sce <- as.SingleCellExperiment(sobj, assay = DefaultAssay(sobj))
    reducedDim(sce, "GRAPH") <- as.matrix(G_use)
    if (!is.null(L_use)) reducedDim(sce, "LAYOUT") <- as.matrix(L_use)
    
    # ---- Milo graph/neighborhoods ----
    vcat("Building kNN graph and neighborhoods...")
    milo <- Milo(sce)
    milo <- buildGraph(milo, k = k, d = ncol(G_use), reduced.dim = "GRAPH")
    milo <- makeNhoods(milo, prop = prop, k = k, refined = TRUE, reduced_dims = "GRAPH")
    
    vcat(sprintf("Created %d neighborhoods", ncol(nhoods(milo))))
    
    # ---- metadata ----
    md <- as.data.frame(colData(milo))
    md[[patient_var]] <- factor(md[[patient_var]])
    md[[target_var]]  <- factor(md[[target_var]])
    md[[batch_var]]   <- factor(md[[batch_var]])
    md[[cluster_var]] <- factor(md[[cluster_var]])
    
    # ---- fast counts (H x S) ----
    vcat("Computing neighborhood-by-sample counts...")
    counts <- .fast_count_cells(milo, samples = md[[patient_var]])
    H <- nrow(counts)
    
    # Filter neighborhoods
    nh_sizes <- rowSums(counts)
    keep_nh <- nh_sizes >= min_cells
    if (sum(!keep_nh) > 0) {
        vcat(sprintf("Filtering out %d neighborhoods with < %d cells", sum(!keep_nh), min_cells))
    }
    
    # ---- majority cluster per nhood ----
    nh_mat <- miloR::nhoods(milo)
    cell_clusters <- md[[cluster_var]]
    majority_cluster <- apply(nh_mat, 2, function(idx){
        labs <- cell_clusters[as.logical(idx)]
        if (!length(labs)) NA else names(sort(table(labs), decreasing = TRUE))[1]
    })
    
    if (!use_glmm) {
        # ===== edgeR with proper sample-level design =====
        vcat("Running edgeR analysis...")
        vcat(sprintf("Fixed effects in model: %s", paste(fixed_effects, collapse = ", ")))
        suppressPackageStartupMessages(library(edgeR))
        
        # Get unique sample-level metadata
        samples <- colnames(counts)
        sample_md <- md[!duplicated(md[[patient_var]]), , drop = FALSE]
        rownames(sample_md) <- as.character(sample_md[[patient_var]])
        
        # Ensure samples are in correct order
        sample_md <- sample_md[samples, , drop = FALSE]
        
        # Create design data frame with only the variables we need
        design.df <- sample_md[, unique(c(patient_var, fixed_effects)), drop = FALSE]
        
        # Clean up factors
        for (v in fixed_effects) {
            if (v %in% colnames(design.df)) {
                design.df[[v]] <- droplevels(factor(design.df[[v]]))
            }
        }
        
        # Create design matrix
        if (length(fixed_effects) == 1 && fixed_effects[1] == target_var) {
            # Simple design with just target variable
            formula_str <- paste("~ 0 +", target_var)
        } else {
            # Standard design
            formula_str <- paste("~", paste(fixed_effects, collapse = " + "))
        }
        
        vcat(sprintf("Design formula: %s", formula_str))
        
        mm <- tryCatch({
            model.matrix(as.formula(formula_str), data = design.df)
        }, error = function(e) {
            vcat("Error creating design matrix. Trying alternative approach...")
            # Fallback to simpler design
            model.matrix(as.formula(paste("~ 0 +", target_var)), data = design.df)
        })
        
        vcat(sprintf("Design matrix dimensions: %d samples x %d coefficients", nrow(mm), ncol(mm)))
        vcat(sprintf("Design matrix rank: %d", qr(mm)$rank))
        
        if (qr(mm)$rank < ncol(mm)) {
            warning("Design matrix is not full rank. Results may be unreliable.")
        }
        
        # EdgeR workflow
        dge <- edgeR::DGEList(counts = counts[keep_nh, , drop = FALSE])
        dge <- edgeR::calcNormFactors(dge, method = "TMM")
        
        # Estimate dispersions
        vcat("Estimating dispersions...")
        dge <- tryCatch({
            edgeR::estimateDisp(dge, design = mm, robust = TRUE)
        }, error = function(e) {
            vcat("Robust dispersion estimation failed, trying standard...")
            edgeR::estimateDisp(dge, design = mm, robust = FALSE)
        })
        
        # Fit model
        vcat("Fitting GLM...")
        fit <- edgeR::glmQLFit(dge, design = mm, robust = TRUE)
        
        # Determine which coefficient(s) to test
        if (grepl("~ 0 \\+", formula_str)) {
            # No intercept model - need to set up contrasts
            target_levels <- levels(design.df[[target_var]])
            if (length(target_levels) == 2) {
                # Create contrast for two-group comparison
                contrast <- makeContrasts(
                    contrasts = paste0(target_var, target_levels[2], " - ", 
                                     target_var, target_levels[1]),
                    levels = mm
                )
                qlf <- edgeR::glmQLFTest(fit, contrast = contrast[,1])
                vcat(sprintf("Testing contrast: %s vs %s", target_levels[2], target_levels[1]))
            } else {
                # Multi-level factor - test all coefficients
                qlf <- edgeR::glmQLFTest(fit, coef = 2:ncol(mm))
                vcat("Testing all non-reference levels")
            }
        } else {
            # Standard model with intercept
            target_coefs <- grep(paste0("^", target_var), colnames(mm), value = TRUE)
            if (length(target_coefs) > 0) {
                qlf <- edgeR::glmQLFTest(fit, coef = target_coefs[1])
                vcat(sprintf("Testing coefficient: %s", target_coefs[1]))
            } else {
                # Fallback to testing second coefficient
                qlf <- edgeR::glmQLFTest(fit, coef = 2)
                vcat(sprintf("Testing coefficient: %s", colnames(mm)[2]))
            }
        }
        
        # Build results
        nh_idx <- which(keep_nh)
        da.res <- tibble(
            nhood = nh_idx,
            logFC = qlf$table$logFC,
            pval  = qlf$table$PValue
        )
        
        # Spatial FDR
        vcat("Computing spatial FDR...")
        da.res <- miloR::graphSpatialFDR(milo, da.res)
        
        da.res$major_cluster <- majority_cluster[da.res$nhood]
        da.res <- da.res %>%
            mutate(is_sig = spatialFDR < fdr_thresh & abs(logFC) > logfc_thresh) %>%
            arrange(spatialFDR, desc(abs(logFC)))
        
        cluster_summary <- da.res %>% 
            filter(is_sig) %>%
            group_by(major_cluster) %>%
            summarise(
                n_sig_nhoods = n(),
                mean_logFC = mean(logFC),
                min_spatialFDR = min(spatialFDR),
                .groups = "drop"
            ) %>%
            arrange(desc(n_sig_nhoods))
        
        vcat(sprintf("Found %d significant neighborhoods", sum(da.res$is_sig)))
        
        return(list(
            mode                  = "edgeR_fixed_effects",
            formula               = formula_str,
            milo                  = milo,
            da_nhood_with_cluster = da.res,
            da_by_cluster         = cluster_summary,
            design_matrix         = mm,
            design_df             = design.df,
            is_nested_design      = is_nested
        ))
    }
    
    # GLMM section remains similar but with nested design awareness...
    # [Rest of GLMM code would go here]
}

milo_opus2_cursor <- function(
    sobj,
    target_var,
    batch_var,
    patient_var,
    cluster_var,
    graph_reduction,                 
    layout_reduction = NULL,         
    k = 30,
    d = NULL,                        
    prop = 0.1,
    fixed_effects  = NULL,
    random_effects = NULL,
    glmm_family    = c("nb","poisson"),
    p_adj_method   = "BH",
    fdr_thresh     = 0.1,
    logfc_thresh   = 0.5,
    seed           = 1,
    use_future     = FALSE,
    workers        = max(1, parallel::detectCores()-1),
    verbose        = TRUE,
    min_cells      = 3,
    force_simple_design = FALSE      # New option for nested designs
){
    suppressPackageStartupMessages({
        library(Seurat)
        library(SingleCellExperiment)
        library(miloR)
        library(Matrix)
        library(dplyr)
        library(tibble)
    })
    set.seed(seed)
    
    # ---- Verbose helper ----
    vcat <- function(...) if(verbose) cat(..., "\n")
    
    # ---- checks ----
    req_cols <- c(target_var, batch_var, patient_var, cluster_var)
    stopifnot("Required columns missing from metadata" = all(req_cols %in% colnames(sobj@meta.data)))
    stopifnot("Graph reduction not found" = graph_reduction %in% Reductions(sobj))
    if (!is.null(layout_reduction)) {
        stopifnot("Layout reduction not found" = layout_reduction %in% Reductions(sobj))
    }
    if (is.null(fixed_effects)) fixed_effects <- c(target_var, batch_var)
    glmm_family <- match.arg(glmm_family)
    use_glmm <- length(random_effects) > 0
    
    vcat("Starting Milo DA analysis...")
    vcat("Mode:", ifelse(use_glmm, "GLMM with random effects", "edgeR fixed effects"))
    
    # ---- Check for nested design ----
    md_check <- sobj@meta.data[, c(patient_var, target_var, batch_var), drop = FALSE]
    md_check <- unique(md_check)
    colnames(md_check) <- c("patient_var_col", "target_var_col", "batch_var_col")
    
    # Check if each patient has only one value for target and batch
    patients_per_target <- md_check %>% 
        group_by(patient_var_col) %>% 
        summarise(n_target = n_distinct(target_var_col), 
                  n_batch = n_distinct(batch_var_col),
                  .groups = "drop")
    
    is_nested <- all(patients_per_target$n_target == 1) & all(patients_per_target$n_batch == 1)
    
    if (is_nested && !force_simple_design) {
        vcat("WARNING: Detected nested design - each patient has only one value for target and batch variables")
        vcat("This creates a confounded design. Adjusting analysis strategy...")
        
        # Check if batch is nested within target
        batch_target_table <- unique(md_check[, c("target_var_col", "batch_var_col")])
        batch_nested_in_target <- nrow(batch_target_table) == length(unique(md_check[["batch_var_col"]]))
        
        if (batch_nested_in_target) {
            vcat("Batch is nested within target variable. Removing batch from fixed effects.")
            fixed_effects <- setdiff(fixed_effects, batch_var)
            if (length(fixed_effects) == 0) fixed_effects <- target_var
        }
    }
    
    # ---- pull embeddings ----
    G_emb <- Embeddings(sobj, reduction = graph_reduction)
    if (is.null(d)) d <- ncol(G_emb)
    if (d > ncol(G_emb)) d <- ncol(G_emb)
    G_use <- G_emb[, seq_len(d), drop = FALSE]
    
    L_use <- NULL
    if (!is.null(layout_reduction)) {
        L_emb <- Embeddings(sobj, reduction = layout_reduction)
        L_use <- L_emb
    }
    
    vcat(sprintf("Using %d dimensions from %s reduction", d, graph_reduction))
    
    # ---- SCE & reducedDims ----
    sce <- as.SingleCellExperiment(sobj, assay = DefaultAssay(sobj))
    reducedDim(sce, "GRAPH") <- as.matrix(G_use)
    if (!is.null(L_use)) reducedDim(sce, "LAYOUT") <- as.matrix(L_use)
    
    # ---- Milo graph/neighborhoods ----
    vcat("Building kNN graph and neighborhoods...")
    milo <- Milo(sce)
    milo <- buildGraph(milo, k = k, d = ncol(G_use), reduced.dim = "GRAPH")
    milo <- makeNhoods(milo, prop = prop, k = k, refined = TRUE, reduced_dims = "GRAPH")
    
    vcat(sprintf("Created %d neighborhoods", ncol(nhoods(milo))))
    
    # ---- metadata ----
    md <- as.data.frame(colData(milo))
    md[[patient_var]] <- factor(md[[patient_var]])
    md[[target_var]]  <- factor(md[[target_var]])
    md[[batch_var]]   <- factor(md[[batch_var]])
    md[[cluster_var]] <- factor(md[[cluster_var]])
    
    # ---- fast counts (H x S) ----
    vcat("Computing neighborhood-by-sample counts...")
    counts <- .fast_count_cells(milo, samples = md[[patient_var]])
    H <- nrow(counts)
    
    # Filter neighborhoods
    nh_sizes <- rowSums(counts)
    keep_nh <- nh_sizes >= min_cells
    if (sum(!keep_nh) > 0) {
        vcat(sprintf("Filtering out %d neighborhoods with < %d cells", sum(!keep_nh), min_cells))
    }
    
    # ---- majority cluster per nhood ----
    nh_mat <- miloR::nhoods(milo)
    cell_clusters <- md[[cluster_var]]
    majority_cluster <- apply(nh_mat, 2, function(idx){
        labs <- cell_clusters[as.logical(idx)]
        if (!length(labs)) NA else names(sort(table(labs), decreasing = TRUE))[1]
    })
    
    if (!use_glmm) {
        # ===== edgeR with proper sample-level design =====
        vcat("Running edgeR analysis...")
        vcat(sprintf("Fixed effects in model: %s", paste(fixed_effects, collapse = ", ")))
        suppressPackageStartupMessages(library(edgeR))
        
        # Get unique sample-level metadata
        samples <- colnames(counts)
        sample_md <- md[!duplicated(md[[patient_var]]), , drop = FALSE]
        rownames(sample_md) <- as.character(sample_md[[patient_var]])
        
        # Ensure samples are in correct order
        sample_md <- sample_md[samples, , drop = FALSE]
        
        # Create design data frame with only the variables we need
        design.df <- sample_md[, unique(c(patient_var, fixed_effects)), drop = FALSE]
        
        # Clean up factors
        for (v in fixed_effects) {
            if (v %in% colnames(design.df)) {
                design.df[[v]] <- droplevels(factor(design.df[[v]]))
            }
        }
        
        # Create design matrix
        if (length(fixed_effects) == 1 && fixed_effects[1] == target_var) {
            # Simple design with just target variable
            formula_str <- paste("~ 0 +", target_var)
        } else {
            # Standard design
            formula_str <- paste("~", paste(fixed_effects, collapse = " + "))
        }
        
        vcat(sprintf("Design formula: %s", formula_str))
        
        mm <- tryCatch({
            model.matrix(as.formula(formula_str), data = design.df)
        }, error = function(e) {
            vcat("Error creating design matrix. Trying alternative approach...")
            # Fallback to simpler design
            model.matrix(as.formula(paste("~ 0 +", target_var)), data = design.df)
        })
        
        vcat(sprintf("Design matrix dimensions: %d samples x %d coefficients", nrow(mm), ncol(mm)))
        vcat(sprintf("Design matrix rank: %d", qr(mm)$rank))
        
        if (qr(mm)$rank < ncol(mm)) {
            warning("Design matrix is not full rank. Results may be unreliable.")
        }
        
        # EdgeR workflow
        dge <- edgeR::DGEList(counts = counts[keep_nh, , drop = FALSE])
        dge <- edgeR::calcNormFactors(dge, method = "TMM")
        
        # Estimate dispersions
        vcat("Estimating dispersions...")
        dge <- tryCatch({
            edgeR::estimateDisp(dge, design = mm, robust = TRUE)
        }, error = function(e) {
            vcat("Robust dispersion estimation failed, trying standard...")
            edgeR::estimateDisp(dge, design = mm, robust = FALSE)
        })
        
        # Fit model
        vcat("Fitting GLM...")
        fit <- edgeR::glmQLFit(dge, design = mm, robust = TRUE)
        
        # Determine which coefficient(s) to test
        if (grepl("~ 0 \\+", formula_str)) {
            # No intercept model - need to set up contrasts
            target_levels <- levels(design.df[[target_var]])
            if (length(target_levels) == 2) {
                # Create contrast for two-group comparison
                contrast <- makeContrasts(
                    contrasts = paste0(target_var, target_levels[2], " - ", 
                                     target_var, target_levels[1]),
                    levels = mm
                )
                qlf <- edgeR::glmQLFTest(fit, contrast = contrast[,1])
                vcat(sprintf("Testing contrast: %s vs %s", target_levels[2], target_levels[1]))
            } else {
                # Multi-level factor - test all coefficients
                qlf <- edgeR::glmQLFTest(fit, coef = 2:ncol(mm))
                vcat("Testing all non-reference levels")
            }
        } else {
            # Standard model with intercept
            target_coefs <- grep(paste0("^", target_var), colnames(mm), value = TRUE)
            if (length(target_coefs) > 0) {
                qlf <- edgeR::glmQLFTest(fit, coef = target_coefs[1])
                vcat(sprintf("Testing coefficient: %s", target_coefs[1]))
            } else {
                # Fallback to testing second coefficient
                qlf <- edgeR::glmQLFTest(fit, coef = 2)
                vcat(sprintf("Testing coefficient: %s", colnames(mm)[2]))
            }
        }
        
        # Build results
        nh_idx <- which(keep_nh)
        da.res <- tibble(
            nhood = nh_idx,
            logFC = qlf$table$logFC,
            pval  = qlf$table$PValue
        )
        
        # Spatial FDR
        vcat("Computing spatial FDR...")
        da.res <- miloR::graphSpatialFDR(milo, da.res)
        
        da.res$major_cluster <- majority_cluster[da.res$nhood]
        da.res <- da.res %>%
            mutate(is_sig = spatialFDR < fdr_thresh & abs(logFC) > logfc_thresh) %>%
            arrange(spatialFDR, desc(abs(logFC)))
        
        cluster_summary <- da.res %>% 
            filter(is_sig) %>%
            group_by(major_cluster) %>%
            summarise(
                n_sig_nhoods = n(),
                mean_logFC = mean(logFC),
                min_spatialFDR = min(spatialFDR),
                .groups = "drop"
            ) %>%
            arrange(desc(n_sig_nhoods))
        
        vcat(sprintf("Found %d significant neighborhoods", sum(da.res$is_sig)))
        
        return(list(
            mode                  = "edgeR_fixed_effects",
            formula               = formula_str,
            milo                  = milo,
            da_nhood_with_cluster = da.res,
            da_by_cluster         = cluster_summary,
            design_matrix         = mm,
            design_df             = design.df,
            is_nested_design      = is_nested
        ))
    }
    
    # GLMM section remains similar but with nested design awareness...
    # [Rest of GLMM code would go here]
}

milo_opus3 <- function(
    sobj,
    target_var,
    batch_var,
    patient_var,
    cluster_var,
    graph_reduction,                 
    layout_reduction = NULL,         
    k = 30,
    d = NULL,                        
    prop = 0.1,
    fixed_effects  = NULL,
    random_effects = NULL,
    glmm_family    = c("nb","poisson"),
    p_adj_method   = "BH",
    fdr_thresh     = 0.1,
    logfc_thresh   = 0.5,
    seed           = 1,
    use_future     = FALSE,
    workers        = max(1, parallel::detectCores()-1),
    verbose        = TRUE,
    min_cells      = 3,
    auto_adjust_design = TRUE        # Automatically adjust for nested designs
){
    suppressPackageStartupMessages({
        library(Seurat)
        library(SingleCellExperiment)
        library(miloR)
        library(Matrix)
        library(dplyr)
        library(tibble)
    })
    set.seed(seed)
    
    # ---- Verbose helper ----
    vcat <- function(...) if(verbose) cat(..., "\n")
    
    # ---- checks ----
    req_cols <- c(target_var, batch_var, patient_var, cluster_var)
    stopifnot("Required columns missing from metadata" = all(req_cols %in% colnames(sobj@meta.data)))
    stopifnot("Graph reduction not found" = graph_reduction %in% Reductions(sobj))
    if (!is.null(layout_reduction)) {
        stopifnot("Layout reduction not found" = layout_reduction %in% Reductions(sobj))
    }
    if (is.null(fixed_effects)) fixed_effects <- c(target_var, batch_var)
    glmm_family <- match.arg(glmm_family)
    use_glmm <- length(random_effects) > 0
    
    vcat("Starting Milo DA analysis...")
    vcat("Mode:", ifelse(use_glmm, "GLMM with random effects", "edgeR fixed effects"))
    
    # ---- Analyze experimental design ----
    md_temp <- sobj@meta.data[, c(patient_var, target_var, batch_var), drop = FALSE]
    md_unique <- unique(md_temp)
    
    # Check confounding
    patient_summary <- md_unique %>%
        group_by(across(all_of(patient_var))) %>%
        summarise(
            n_target = n_distinct(!!sym(target_var)),
            n_batch = n_distinct(!!sym(batch_var)),
            target_val = first(!!sym(target_var)),
            batch_val = first(!!sym(batch_var)),
            .groups = "drop"
        )
    
    # Check if design is nested
    is_nested <- all(patient_summary$n_target == 1) & all(patient_summary$n_batch == 1)
    
    if (is_nested && auto_adjust_design) {
        vcat("WARNING: Detected nested/confounded design:")
        vcat("  - Each patient belongs to exactly one target group and one batch")
        
        # Check relationship between batch and target
        batch_target_df <- unique(patient_summary[, c("target_val", "batch_val")])
        n_batches <- n_distinct(batch_target_df$batch_val)
        n_targets <- n_distinct(batch_target_df$target_val)
        
        # Check if batches are nested within targets
        batches_per_target <- batch_target_df %>%
            group_by(target_val) %>%
            summarise(n_batches = n_distinct(batch_val), .groups = "drop")
        
        if (nrow(batch_target_df) == n_batches) {
            # Each batch appears in only one target group
            vcat("  - Batches are nested within target groups (each batch belongs to one target group)")
            vcat("  - Cannot separate batch effects from target effects at the patient level")
            vcat("  - Adjusting model to use only target variable as fixed effect")
            fixed_effects <- target_var
        } else {
            # Batches span multiple targets but patients are still confounded
            vcat("  - Batches span multiple target groups but are confounded with patients")
            vcat("  - Will attempt to include batch as a covariate, but interpretation should be cautious")
        }
    }
    
    # ---- pull embeddings ----
    G_emb <- Embeddings(sobj, reduction = graph_reduction)
    if (is.null(d)) d <- ncol(G_emb)
    if (d > ncol(G_emb)) d <- ncol(G_emb)
    G_use <- G_emb[, seq_len(d), drop = FALSE]
    
    L_use <- NULL
    if (!is.null(layout_reduction)) {
        L_emb <- Embeddings(sobj, reduction = layout_reduction)
        L_use <- L_emb
    }
    
    vcat(sprintf("Using %d dimensions from %s reduction", d, graph_reduction))
    
    # ---- SCE & reducedDims ----
    sce <- as.SingleCellExperiment(sobj, assay = DefaultAssay(sobj))
    reducedDim(sce, "GRAPH") <- as.matrix(G_use)
    if (!is.null(L_use)) reducedDim(sce, "LAYOUT") <- as.matrix(L_use)
    
    # ---- Milo graph/neighborhoods ----
    vcat("Building kNN graph and neighborhoods...")
    milo <- Milo(sce)
    milo <- buildGraph(milo, k = k, d = ncol(G_use), reduced.dim = "GRAPH")
    milo <- makeNhoods(milo, prop = prop, k = k, refined = TRUE, reduced_dims = "GRAPH")
    
    vcat(sprintf("Created %d neighborhoods", ncol(nhoods(milo))))
    
    # ---- metadata ----
    md <- as.data.frame(colData(milo))
    md[[patient_var]] <- factor(md[[patient_var]])
    md[[target_var]]  <- factor(md[[target_var]])
    md[[batch_var]]   <- factor(md[[batch_var]])
    md[[cluster_var]] <- factor(md[[cluster_var]])
    
    # ---- fast counts (H x S) ----
    vcat("Computing neighborhood-by-sample counts...")
    counts <- .fast_count_cells(milo, samples = md[[patient_var]])
    H <- nrow(counts)
    S <- ncol(counts)
    vcat(sprintf("Count matrix: %d neighborhoods x %d samples", H, S))
    
    # Filter neighborhoods
    nh_sizes <- rowSums(counts)
    keep_nh <- nh_sizes >= min_cells
    if (sum(!keep_nh) > 0) {
        vcat(sprintf("Filtering out %d neighborhoods with < %d cells", sum(!keep_nh), min_cells))
    }
    
    # ---- majority cluster per nhood ----
    nh_mat <- miloR::nhoods(milo)
    cell_clusters <- md[[cluster_var]]
    majority_cluster <- apply(nh_mat, 2, function(idx){
        labs <- cell_clusters[as.logical(idx)]
        if (!length(labs)) NA else names(sort(table(labs), decreasing = TRUE))[1]
    })
    
    if (!use_glmm) {
        # ===== edgeR with proper handling of nested design =====
        vcat("Running edgeR analysis...")
        vcat(sprintf("Fixed effects in model: %s", paste(fixed_effects, collapse = ", ")))
        suppressPackageStartupMessages(library(edgeR))
        
        # Create sample-level metadata
        samples <- colnames(counts)
        sample_df <- patient_summary[match(samples, patient_summary[[patient_var]]), ]
        rownames(sample_df) <- samples
        
        # Rename columns for design matrix
        colnames(sample_df)[colnames(sample_df) == "target_val"] <- target_var
        colnames(sample_df)[colnames(sample_df) == "batch_val"] <- batch_var
        
        # Ensure factors
        for (v in c(target_var, batch_var)) {
            if (v %in% colnames(sample_df)) {
                sample_df[[v]] <- factor(sample_df[[v]])
            }
        }
        
        # Build design matrix based on what's actually possible
        design_built <- FALSE
        mm <- NULL
        
        # Try the requested design first
        if (!design_built && length(fixed_effects) > 1) {
            # For multiple fixed effects, check if they're not perfectly confounded
            formula_str <- paste("~", paste(fixed_effects, collapse = " + "))
            vcat(sprintf("Attempting design: %s", formula_str))
            
            mm <- tryCatch({
                # Only include variables that are in fixed_effects
                design_vars <- intersect(fixed_effects, colnames(sample_df))
                if (length(design_vars) == 0) stop("No valid design variables")
                
                temp_df <- sample_df[, design_vars, drop = FALSE]
                # Ensure all are factors with >1 level
                for (v in design_vars) {
                    temp_df[[v]] <- droplevels(factor(temp_df[[v]]))
                    if (nlevels(temp_df[[v]]) < 2) {
                        stop(sprintf("Variable %s has only one level", v))
                    }
                }
                model.matrix(as.formula(paste("~", paste(design_vars, collapse = " + "))), 
                           data = temp_df)
            }, error = function(e) {
                vcat(sprintf("  Failed: %s", e$message))
                NULL
            })
            
            if (!is.null(mm) && qr(mm)$rank == ncol(mm)) {
                design_built <- TRUE
                vcat("  Success: Design matrix is full rank")
            } else if (!is.null(mm)) {
                vcat("  Design matrix is rank deficient, trying simpler model")
                mm <- NULL
            }
        }
        
        # Fallback to target variable only
        if (!design_built) {
            formula_str <- paste("~ 0 +", target_var)
            vcat(sprintf("Using simplified design: %s", formula_str))
            
            mm <- model.matrix(as.formula(formula_str), 
                             data = sample_df[, target_var, drop = FALSE])
            design_built <- TRUE
        }
        
        vcat(sprintf("Final design matrix: %d samples x %d coefficients", nrow(mm), ncol(mm)))
        vcat(sprintf("Coefficients: %s", paste(colnames(mm), collapse = ", ")))
        
        # EdgeR workflow
        dge <- edgeR::DGEList(counts = counts[keep_nh, , drop = FALSE])
        dge <- edgeR::calcNormFactors(dge, method = "TMM")
        
        # Estimate dispersions
        vcat("Estimating dispersions...")
        dge <- tryCatch({
            edgeR::estimateDisp(dge, design = mm, robust = TRUE)
        }, error = function(e) {
            vcat("Robust dispersion estimation failed, trying standard...")
            edgeR::estimateDisp(dge, design = mm, robust = FALSE)
        })
        
        vcat(sprintf("Common dispersion: %.4f", dge$common.dispersion))
        
        # Fit model
        vcat("Fitting GLM...")
        fit <- edgeR::glmQLFit(dge, design = mm, robust = TRUE)
        
        # Determine contrast/coefficient to test
        if (grepl("~ 0 \\+", formula_str)) {
            # No intercept model - create contrast
            target_levels <- levels(sample_df[[target_var]])
            if (length(target_levels) == 2) {
                contrast_str <- paste0(target_var, target_levels[2], " - ", 
                                     target_var, target_levels[1])
                vcat(sprintf("Testing contrast: %s", contrast_str))
                
                contrast <- makeContrasts(
                    contrasts = contrast_str,
                    levels = mm
                )
                qlf <- edgeR::glmQLFTest(fit, contrast = contrast[,1])
            } else {
                vcat("Testing all non-reference levels")
                qlf <- edgeR::glmQLFTest(fit, coef = 2:ncol(mm))
            }
        } else {
            # Model with intercept
            target_coefs <- grep(paste0("^", target_var), colnames(mm), value = TRUE)
            if (length(target_coefs) > 0) {
                vcat(sprintf("Testing coefficient: %s", target_coefs[1]))
                qlf <- edgeR::glmQLFTest(fit, coef = target_coefs[1])
            } else {
                vcat(sprintf("Testing coefficient: %s", colnames(mm)[2]))
                qlf <- edgeR::glmQLFTest(fit, coef = 2)
            }
        }
        
        # Build results
        nh_idx <- which(keep_nh)
        da.res <- tibble(
            nhood = nh_idx,
            logFC = qlf$table$logFC,
            logCPM = qlf$table$logCPM,
            F = qlf$table$F,
            pval = qlf$table$PValue
        )
        
        # Spatial FDR
        vcat("Computing spatial FDR...")
        da.res <- miloR::graphSpatialFDR(milo, da.res, 
                                         pvalues = "pval",
                                         indices = "nhood")
        
        da.res$major_cluster <- majority_cluster[da.res$nhood]
        da.res <- da.res %>%
            mutate(is_sig = spatialFDR < fdr_thresh & abs(logFC) > logfc_thresh) %>%
            arrange(spatialFDR, desc(abs(logFC)))
        
        cluster_summary <- da.res %>% 
            filter(is_sig) %>%
            group_by(major_cluster) %>%
            summarise(
                n_sig_nhoods = n(),
                mean_logFC = mean(logFC),
                median_logFC = median(logFC),
                min_spatialFDR = min(spatialFDR),
                .groups = "drop"
            ) %>%
            arrange(desc(n_sig_nhoods))
        
        vcat(sprintf("Found %d significant neighborhoods", sum(da.res$is_sig, na.rm = TRUE)))
        if (nrow(cluster_summary) > 0) {
            vcat(sprintf("Affecting %d clusters", nrow(cluster_summary)))
        }
        
        return(list(
            mode = "edgeR_fixed_effects",
            formula = formula_str,
            milo = milo,
            da_nhood_with_cluster = da.res,
            da_by_cluster = cluster_summary,
            design_matrix = mm,
            design_info = list(
                is_nested = is_nested,
                n_samples = S,
                n_neighborhoods = H,
                n_filtered = sum(!keep_nh)
            )
        ))
    }
    
    # [GLMM section would go here if needed...]
}

milo_opus4 <- function(
    sobj,
    target_var,
    batch_var,
    patient_var,
    cluster_var,
    graph_reduction,                 
    layout_reduction = NULL,         
    k = 30,
    d = NULL,                        
    prop = 0.1,
    fixed_effects  = NULL,
    random_effects = NULL,
    glmm_family    = c("nb","poisson"),
    p_adj_method   = "BH",
    fdr_thresh     = 0.1,
    logfc_thresh   = 0.5,
    seed           = 1,
    use_future     = FALSE,
    workers        = max(1, parallel::detectCores()-1),
    verbose        = TRUE,
    min_cells      = 3
){
    suppressPackageStartupMessages({
        library(Seurat)
        library(SingleCellExperiment)
        library(miloR)
        library(Matrix)
        library(dplyr)
        library(tibble)
    })
    set.seed(seed)
    
    # ---- Verbose helper ----
    vcat <- function(...) if(verbose) cat(..., "\n")
    
    # ---- checks ----
    req_cols <- c(target_var, batch_var, patient_var, cluster_var)
    stopifnot("Required columns missing from metadata" = all(req_cols %in% colnames(sobj@meta.data)))
    stopifnot("Graph reduction not found" = graph_reduction %in% Reductions(sobj))
    if (!is.null(layout_reduction)) {
        stopifnot("Layout reduction not found" = layout_reduction %in% Reductions(sobj))
    }
    if (is.null(fixed_effects)) fixed_effects <- c(target_var, batch_var)
    glmm_family <- match.arg(glmm_family)
    use_glmm <- length(random_effects) > 0
    
    vcat("Starting Milo DA analysis...")
    vcat("Mode:", ifelse(use_glmm, "GLMM with random effects", "edgeR fixed effects"))
    
    # ---- Analyze experimental design WITHOUT tidyeval ----
    md_temp <- sobj@meta.data[, c(patient_var, target_var, batch_var), drop = FALSE]
    md_unique <- unique(md_temp)
    
    # Create patient summary using base R approach
    patient_ids <- unique(md_unique[[patient_var]])
    patient_summary <- data.frame(
        patient = patient_ids,
        stringsAsFactors = FALSE
    )
    colnames(patient_summary)[1] <- patient_var
    
    # Calculate summaries for each patient
    for (i in seq_along(patient_ids)) {
        pid <- patient_ids[i]
        patient_rows <- md_unique[[patient_var]] == pid
        patient_data <- md_unique[patient_rows, , drop = FALSE]
        
        patient_summary$n_target[i] <- length(unique(patient_data[[target_var]]))
        patient_summary$n_batch[i] <- length(unique(patient_data[[batch_var]]))
        patient_summary$target_val[i] <- patient_data[[target_var]][1]
        patient_summary$batch_val[i] <- patient_data[[batch_var]][1]
    }
    
    # Check if design is nested
    is_nested <- all(patient_summary$n_target == 1) & all(patient_summary$n_batch == 1)
    
    if (is_nested) {
        vcat("WARNING: Detected nested/confounded design:")
        vcat("  - Each patient belongs to exactly one target group and one batch")
        
        # Simplified: just use target variable if nested
        vcat("  - Simplifying model to use only target variable")
        fixed_effects <- target_var
    }
    
    # ---- pull embeddings ----
    G_emb <- Embeddings(sobj, reduction = graph_reduction)
    if (is.null(d)) d <- ncol(G_emb)
    if (d > ncol(G_emb)) d <- ncol(G_emb)
    G_use <- G_emb[, seq_len(d), drop = FALSE]
    
    L_use <- NULL
    if (!is.null(layout_reduction)) {
        L_emb <- Embeddings(sobj, reduction = layout_reduction)
        L_use <- L_emb
    }
    
    vcat(sprintf("Using %d dimensions from %s reduction", d, graph_reduction))
    
    # ---- SCE & reducedDims ----
    sce <- as.SingleCellExperiment(sobj, assay = DefaultAssay(sobj))
    reducedDim(sce, "GRAPH") <- as.matrix(G_use)
    if (!is.null(L_use)) reducedDim(sce, "LAYOUT") <- as.matrix(L_use)
    
    # ---- Milo graph/neighborhoods ----
    vcat("Building kNN graph and neighborhoods...")
    milo <- Milo(sce)
    milo <- buildGraph(milo, k = k, d = ncol(G_use), reduced.dim = "GRAPH")
    milo <- makeNhoods(milo, prop = prop, k = k, refined = TRUE, reduced_dims = "GRAPH")
    
    vcat(sprintf("Created %d neighborhoods", ncol(nhoods(milo))))
    
    # ---- metadata ----
    md <- as.data.frame(colData(milo))
    md[[patient_var]] <- factor(md[[patient_var]])
    md[[target_var]]  <- factor(md[[target_var]])
    md[[batch_var]]   <- factor(md[[batch_var]])
    md[[cluster_var]] <- factor(md[[cluster_var]])
    
    # ---- fast counts (H x S) ----
    vcat("Computing neighborhood-by-sample counts...")
    counts <- .fast_count_cells(milo, samples = md[[patient_var]])
    H <- nrow(counts)
    S <- ncol(counts)
    vcat(sprintf("Count matrix: %d neighborhoods x %d samples", H, S))
    
    # Filter neighborhoods
    nh_sizes <- rowSums(counts)
    keep_nh <- nh_sizes >= min_cells
    if (sum(!keep_nh) > 0) {
        vcat(sprintf("Filtering out %d neighborhoods with < %d cells", sum(!keep_nh), min_cells))
    }
    
    # ---- majority cluster per nhood ----
    nh_mat <- miloR::nhoods(milo)
    cell_clusters <- md[[cluster_var]]
    majority_cluster <- vapply(seq_len(ncol(nh_mat)), function(i) {
        idx <- as.logical(nh_mat[, i])
        labs <- cell_clusters[idx]
        if (length(labs) == 0) return(NA_character_)
        tbl <- table(labs)
        names(tbl)[which.max(tbl)]
    }, character(1))
    
    if (!use_glmm) {
        # ===== edgeR analysis =====
        vcat("Running edgeR analysis...")
        suppressPackageStartupMessages(library(edgeR))
        
        # Create sample-level design matrix
        samples <- colnames(counts)
        
        # Match samples to patient summary
        sample_idx <- match(samples, patient_summary[[patient_var]])
        sample_df <- patient_summary[sample_idx, ]
        rownames(sample_df) <- samples
        
        # Rename columns
        colnames(sample_df)[colnames(sample_df) == "target_val"] <- target_var
        colnames(sample_df)[colnames(sample_df) == "batch_val"] <- batch_var
        
        # Create design based on what's possible
        vcat(sprintf("Fixed effects in model: %s", paste(fixed_effects, collapse = ", ")))
        
        # Build design matrix - try simple approach first
        if (target_var %in% fixed_effects) {
            sample_df[[target_var]] <- factor(sample_df[[target_var]])
            
            if (length(fixed_effects) == 1) {
                # Simple one-factor design
                formula_str <- paste("~ 0 +", target_var)
                vcat(sprintf("Using design: %s", formula_str))
                mm <- model.matrix(as.formula(formula_str), data = sample_df)
                
            } else if (batch_var %in% fixed_effects && !is_nested) {
                # Try to include batch if not nested
                sample_df[[batch_var]] <- factor(sample_df[[batch_var]])
                formula_str <- paste("~", target_var, "+", batch_var)
                vcat(sprintf("Attempting design: %s", formula_str))
                
                mm <- tryCatch({
                    model.matrix(as.formula(formula_str), data = sample_df)
                }, error = function(e) {
                    vcat("  Batch adjustment failed, using target only")
                    formula_str <<- paste("~ 0 +", target_var)
                    model.matrix(as.formula(formula_str), data = sample_df)
                })
            } else {
                # Default to target only
                formula_str <- paste("~ 0 +", target_var)
                vcat(sprintf("Using design: %s", formula_str))
                mm <- model.matrix(as.formula(formula_str), data = sample_df)
            }
        } else {
            stop("Target variable must be in fixed effects")
        }
        
        vcat(sprintf("Design matrix: %d samples x %d coefficients", nrow(mm), ncol(mm)))
        vcat(sprintf("Coefficients: %s", paste(colnames(mm), collapse = ", ")))
        
        # EdgeR workflow
        dge <- edgeR::DGEList(counts = counts[keep_nh, , drop = FALSE])
        dge <- edgeR::calcNormFactors(dge, method = "TMM")
        
        # Estimate dispersions
        vcat("Estimating dispersions...")
        dge <- edgeR::estimateDisp(dge, design = mm, robust = TRUE)
        vcat(sprintf("Common dispersion: %.4f", dge$common.dispersion))
        
        # Fit model
        vcat("Fitting GLM...")
        fit <- edgeR::glmQLFit(dge, design = mm, robust = TRUE)
        
        # Test for differential abundance
        if (grepl("~ 0 \\+", formula_str)) {
            # No intercept model
            target_levels <- levels(sample_df[[target_var]])
            if (length(target_levels) == 2) {
                # Two-group comparison
                contrast_formula <- paste0(target_var, target_levels[2], " - ", 
                                         target_var, target_levels[1])
                vcat(sprintf("Testing contrast: %s", contrast_formula))
                
                contrast <- edgeR::makeContrasts(
                    contrasts = contrast_formula,
                    levels = mm
                )
                qlf <- edgeR::glmQLFTest(fit, contrast = contrast[,1])
            } else {
                # Multi-level comparison
                vcat("Testing all non-reference levels")
                qlf <- edgeR::glmQLFTest(fit, coef = 2:ncol(mm))
            }
        } else {
            # Model with intercept - test target coefficient
            target_coef <- grep(target_var, colnames(mm), value = TRUE)[1]
            if (!is.na(target_coef)) {
                vcat(sprintf("Testing coefficient: %s", target_coef))
                qlf <- edgeR::glmQLFTest(fit, coef = target_coef)
            } else {
                qlf <- edgeR::glmQLFTest(fit, coef = 2)
            }
        }
        
        # Build results
        nh_idx <- which(keep_nh)
        da.res <- tibble(
            nhood = nh_idx,
            logFC = qlf$table$logFC,
            logCPM = qlf$table$logCPM,
            F = qlf$table$F,
            pval = qlf$table$PValue
        )
        
        # Spatial FDR
        vcat("Computing spatial FDR...")
        da.res <- miloR::graphSpatialFDR(milo, da.res, 
                                         pvalues = "pval",
                                         indices = "nhood")
        
        da.res$major_cluster <- majority_cluster[da.res$nhood]
        da.res <- da.res %>%
            mutate(is_sig = spatialFDR < fdr_thresh & abs(logFC) > logfc_thresh) %>%
            arrange(spatialFDR, desc(abs(logFC)))
        
        # Cluster summary
        if (sum(da.res$is_sig, na.rm = TRUE) > 0) {
            cluster_summary <- da.res %>% 
                filter(is_sig) %>%
                group_by(major_cluster) %>%
                summarise(
                    n_sig_nhoods = n(),
                    mean_logFC = mean(logFC),
                    median_logFC = median(logFC),
                    min_spatialFDR = min(spatialFDR),
                    .groups = "drop"
                ) %>%
                arrange(desc(n_sig_nhoods))
        } else {
            cluster_summary <- tibble()
        }
        
        vcat(sprintf("Found %d significant neighborhoods", sum(da.res$is_sig, na.rm = TRUE)))
        
        return(list(
            mode = "edgeR_fixed_effects",
            formula = formula_str,
            milo = milo,
            da_nhood_with_cluster = da.res,
            da_by_cluster = cluster_summary,
            design_matrix = mm,
            design_info = list(
                is_nested = is_nested,
                n_samples = S,
                n_neighborhoods = H,
                n_filtered = sum(!keep_nh)
            )
        ))
    }
    
    # [GLMM section would continue here...]
}

milo_opus5<- function(
    sobj,
    target_var,
    batch_var,
    patient_var,
    cluster_var,
    graph_reduction,                 
    layout_reduction = NULL,         
    k = 30,
    d = NULL,                        
    prop = 0.1,
    fixed_effects  = NULL,
    random_effects = NULL,
    glmm_family    = c("nb","poisson"),
    p_adj_method   = "BH",
    fdr_thresh     = 0.1,
    logfc_thresh   = 0.5,
    seed           = 1,
    use_future     = FALSE,
    workers        = max(1, parallel::detectCores()-1),
    verbose        = TRUE,
    min_cells      = 3
){
    suppressPackageStartupMessages({
        library(Seurat)
        library(SingleCellExperiment)
        library(miloR)
        library(Matrix)
        library(dplyr)
        library(tibble)
        library(edgeR)
    })
    set.seed(seed)
    
    vcat <- function(...) if(verbose) cat(..., "\n")
    
    # ---- Initial checks ----
    req_cols <- c(target_var, batch_var, patient_var, cluster_var)
    stopifnot("Required columns missing" = all(req_cols %in% colnames(sobj@meta.data)))
    stopifnot("Graph reduction not found" = graph_reduction %in% Reductions(sobj))
    
    if (is.null(fixed_effects)) fixed_effects <- c(target_var, batch_var)
    glmm_family <- match.arg(glmm_family)
    use_glmm <- length(random_effects) > 0
    
    vcat("Starting Milo DA analysis...")
    vcat("Mode:", ifelse(use_glmm, "GLMM with random effects", "edgeR fixed effects"))
    
    # ---- Pull embeddings ----
    G_emb <- Embeddings(sobj, reduction = graph_reduction)
    if (is.null(d)) d <- ncol(G_emb)
    d <- min(d, ncol(G_emb))
    G_use <- G_emb[, seq_len(d), drop = FALSE]
    
    L_use <- NULL
    if (!is.null(layout_reduction)) {
        L_emb <- Embeddings(sobj, reduction = layout_reduction)
        L_use <- L_emb
    }
    
    vcat(sprintf("Using %d dimensions from %s reduction", d, graph_reduction))
    
    # ---- Create SCE ----
    sce <- as.SingleCellExperiment(sobj, assay = DefaultAssay(sobj))
    reducedDim(sce, "GRAPH") <- as.matrix(G_use)
    if (!is.null(L_use)) reducedDim(sce, "LAYOUT") <- as.matrix(L_use)
    
    # ---- Build Milo object ----
    vcat("Building kNN graph and neighborhoods...")
    milo <- Milo(sce)
    milo <- buildGraph(milo, k = k, d = ncol(G_use), reduced.dim = "GRAPH")
    milo <- makeNhoods(milo, prop = prop, k = k, refined = TRUE, reduced_dims = "GRAPH")
    
    vcat(sprintf("Created %d neighborhoods", ncol(nhoods(milo))))
    
    # ---- Get metadata ----
    md <- as.data.frame(colData(milo))
    
    # ---- Compute counts ----
    vcat("Computing neighborhood-by-sample counts...")
    counts <- .fast_count_cells(milo, samples = md[[patient_var]])
    H <- nrow(counts)
    S <- ncol(counts)
    vcat(sprintf("Count matrix: %d neighborhoods x %d samples", H, S))
    
    # Filter small neighborhoods
    nh_sizes <- rowSums(counts)
    keep_nh <- nh_sizes >= min_cells
    if (sum(!keep_nh) > 0) {
        vcat(sprintf("Filtering out %d neighborhoods with < %d cells", sum(!keep_nh), min_cells))
    }
    
    # ---- Get majority cluster ----
    nh_mat <- miloR::nhoods(milo)
    cell_clusters <- md[[cluster_var]]
    majority_cluster <- vapply(seq_len(ncol(nh_mat)), function(i) {
        idx <- as.logical(nh_mat[, i])
        labs <- cell_clusters[idx]
        if (length(labs) == 0) return(NA_character_)
        tbl <- table(labs)
        names(tbl)[which.max(tbl)]
    }, character(1))
    
    # ===== CRITICAL FIX: Create proper sample-level design =====
    vcat("Creating sample-level design matrix...")
    
    # Get sample IDs from count matrix
    sample_ids <- colnames(counts)
    
    # Create sample-level metadata directly from original data
    sample_design <- data.frame(
        sample_id = sample_ids,
        stringsAsFactors = FALSE
    )
    
    # For each sample, get its target and batch values
    for (i in seq_along(sample_ids)) {
        sid <- sample_ids[i]
        # Find first cell from this patient
        cell_idx <- which(md[[patient_var]] == sid)[1]
        sample_design[[target_var]][i] <- as.character(md[[target_var]][cell_idx])
        sample_design[[batch_var]][i] <- as.character(md[[batch_var]][cell_idx])
    }
    
    # Convert to factors and check levels
    sample_design[[target_var]] <- factor(sample_design[[target_var]])
    sample_design[[batch_var]] <- factor(sample_design[[batch_var]])
    
    vcat(sprintf("Target variable levels: %s", paste(levels(sample_design[[target_var]]), collapse=", ")))
    vcat(sprintf("Batch variable levels: %s", paste(levels(sample_design[[batch_var]]), collapse=", ")))
    vcat(sprintf("Number of samples: %d", nrow(sample_design)))
    
    # Check if design is nested
    design_check <- table(sample_design[[target_var]], sample_design[[batch_var]])
    vcat("Design structure (target x batch):")
    print(design_check)
    
    is_nested <- all(apply(design_check > 0, 2, sum) == 1)  # Each batch in only one target
    
    if (is_nested) {
        vcat("WARNING: Batches are completely nested within target groups")
        vcat("This means batch effects are confounded with target effects")
        vcat("Consider using batch as a random effect or patient as blocking factor")
    }
    
    # ===== EdgeR analysis =====
    vcat("Running edgeR analysis...")
    
    # Strategy for nested design:
    # Option 1: Just target (loses batch info)
    # Option 2: Use duplicateCorrelation to account for batch
    # Option 3: Include batch anyway (will be partially confounded)
    
    rownames(sample_design) <- sample_ids
    
    # Try to build design matrix with both factors
    design_success <- FALSE
    mm <- NULL
    formula_str <- ""
    
    # First attempt: both factors
    if (length(fixed_effects) > 1 && all(c(target_var, batch_var) %in% fixed_effects)) {
        vcat("Attempting to include both target and batch...")
        formula_str <- paste("~", target_var, "+", batch_var)
        
        mm <- tryCatch({
            # Ensure we have >1 level for each factor
            if (nlevels(sample_design[[target_var]]) < 2) {
                stop(sprintf("%s has only one level", target_var))
            }
            if (nlevels(sample_design[[batch_var]]) < 2) {
                stop(sprintf("%s has only one level", batch_var))  
            }
            
            m <- model.matrix(as.formula(formula_str), data = sample_design)
            
            # Check rank
            if (qr(m)$rank < ncol(m)) {
                vcat("  Design is rank-deficient, but proceeding...")
            }
            m
        }, error = function(e) {
            vcat(sprintf("  Failed: %s", e$message))
            NULL
        })
        
        if (!is.null(mm)) design_success <- TRUE
    }
    
    # Fallback: target only with no intercept
    if (!design_success) {
        vcat("Using target-only design (batch effects will be ignored)")
        formula_str <- paste("~ 0 +", target_var)
        
        # Make sure target has 2+ levels
        if (nlevels(sample_design[[target_var]]) < 2) {
            stop(sprintf("Target variable %s must have at least 2 levels, found: %s",
                        target_var, paste(levels(sample_design[[target_var]]), collapse=", ")))
        }
        
        mm <- model.matrix(as.formula(formula_str), data = sample_design)
    }
    
    vcat(sprintf("Final design: %s", formula_str))
    vcat(sprintf("Design matrix: %d samples x %d coefficients", nrow(mm), ncol(mm)))
    vcat(sprintf("Coefficients: %s", paste(colnames(mm), collapse=", ")))
    
    # EdgeR workflow
    dge <- DGEList(counts = counts[keep_nh, , drop = FALSE])
    dge <- calcNormFactors(dge, method = "TMM")
    
    # Estimate dispersions
    vcat("Estimating dispersions...")
    dge <- estimateDisp(dge, design = mm, robust = TRUE)
    vcat(sprintf("Common dispersion: %.4f", dge$common.dispersion))
    
    # Fit model
    vcat("Fitting GLM...")
    fit <- glmQLFit(dge, design = mm, robust = TRUE)
    
    # Test for differential abundance
    if (grepl("~ 0 \\+", formula_str)) {
        # No intercept model - need contrast
        target_levels <- levels(sample_design[[target_var]])
        
        if (length(target_levels) == 2) {
            contrast_str <- paste0(colnames(mm)[2], " - ", colnames(mm)[1])
            vcat(sprintf("Testing contrast: %s", contrast_str))
            
            contrast_vec <- numeric(ncol(mm))
            contrast_vec[1] <- -1
            contrast_vec[2] <- 1
            
            qlf <- glmQLFTest(fit, contrast = contrast_vec)
        } else {
            vcat("Multiple target levels - testing all pairwise")
            qlf <- glmQLFTest(fit, coef = 2:ncol(mm))
        }
    } else {
        # Model with intercept
        target_coef <- grep(target_var, colnames(mm), value = TRUE)[1]
        vcat(sprintf("Testing coefficient: %s", target_coef))
        qlf <- glmQLFTest(fit, coef = target_coef)
    }
    
    # Results
    nh_idx <- which(keep_nh)
    da.res <- tibble(
        nhood = nh_idx,
        logFC = qlf$table$logFC,
        logCPM = qlf$table$logCPM,
        F = qlf$table$F,
        pval = qlf$table$PValue
    )
    
    # Spatial FDR
    vcat("Computing spatial FDR...")
    da.res <- graphSpatialFDR(milo, da.res, pvalues = "pval", indices = "nhood")
    
    da.res$major_cluster <- majority_cluster[da.res$nhood]
    da.res <- da.res %>%
        mutate(is_sig = spatialFDR < fdr_thresh & abs(logFC) > logfc_thresh) %>%
        arrange(spatialFDR, desc(abs(logFC)))
    
    # Summary
    cluster_summary <- if (sum(da.res$is_sig, na.rm = TRUE) > 0) {
        da.res %>%
            filter(is_sig) %>%
            group_by(major_cluster) %>%
            summarise(
                n_sig_nhoods = n(),
                mean_logFC = mean(logFC),
                median_logFC = median(logFC),
                min_spatialFDR = min(spatialFDR),
                .groups = "drop"
            ) %>%
            arrange(desc(n_sig_nhoods))
    } else {
        tibble()
    }
    
    vcat(sprintf("\nResults: %d significant neighborhoods out of %d tested",
                sum(da.res$is_sig, na.rm = TRUE), nrow(da.res)))
    
    if (is_nested) {
        vcat("\nNOTE: Due to nested design, batch effects could not be separated from target effects.")
        vcat("Consider alternative approaches:")
        vcat("  1. Use GLMM with batch as random effect")
        vcat("  2. Use limma's duplicateCorrelation")
        vcat("  3. Perform sensitivity analysis by stratifying by batch")
    }
    
    return(list(
        mode = "edgeR_fixed_effects",
        formula = formula_str,
        milo = milo,
        da_nhood_with_cluster = da.res,
        da_by_cluster = cluster_summary,
        design_matrix = mm,
        sample_design = sample_design,
        design_info = list(
            is_nested = is_nested,
            n_samples = S,
            n_neighborhoods = H,
            n_filtered = sum(!keep_nh)
        )
    ))
}

milo_opus6 <- function(
    sobj,
    target_var,
    batch_var,
    patient_var,
    cluster_var,
    graph_reduction,                 
    layout_reduction = NULL,         
    k = 30,
    d = NULL,                        
    prop = 0.1,
    fixed_effects  = NULL,
    random_effects = NULL,
    glmm_family    = c("nb","poisson"),
    p_adj_method   = "BH",
    fdr_thresh     = 0.1,
    logfc_thresh   = 0.5,
    seed           = 1,
    use_future     = FALSE,
    workers        = max(1, parallel::detectCores()-1),
    verbose        = TRUE,
    min_cells      = 3
){
    suppressPackageStartupMessages({
        library(Seurat)
        library(SingleCellExperiment)
        library(miloR)
        library(Matrix)
        library(dplyr)
        library(tibble)
        library(edgeR)
    })
    set.seed(seed)
    
    vcat <- function(...) if(verbose) cat(..., "\n")
    
    # ---- Initial checks ----
    req_cols <- c(target_var, batch_var, patient_var, cluster_var)
    stopifnot("Required columns missing" = all(req_cols %in% colnames(sobj@meta.data)))
    stopifnot("Graph reduction not found" = graph_reduction %in% Reductions(sobj))
    
    if (is.null(fixed_effects)) fixed_effects <- c(target_var, batch_var)
    glmm_family <- match.arg(glmm_family)
    use_glmm <- length(random_effects) > 0
    
    vcat("Starting Milo DA analysis...")
    vcat("Mode:", ifelse(use_glmm, "GLMM with random effects", "edgeR fixed effects"))
    
    # ---- Convert metadata to proper types ----
    md <- sobj@meta.data
    
    # Convert to character first to preserve values, then to factor
    md[[target_var]] <- as.character(md[[target_var]])
    md[[batch_var]] <- as.character(md[[batch_var]])
    md[[patient_var]] <- as.character(md[[patient_var]])
    md[[cluster_var]] <- as.character(md[[cluster_var]])
    
    vcat(sprintf("Target variable (%s) unique values: %s", 
                target_var, paste(unique(md[[target_var]]), collapse=", ")))
    vcat(sprintf("Batch variable (%s) unique values: %s", 
                batch_var, paste(unique(md[[batch_var]]), collapse=", ")))
    
    # ---- Pull embeddings ----
    G_emb <- Embeddings(sobj, reduction = graph_reduction)
    if (is.null(d)) d <- ncol(G_emb)
    d <- min(d, ncol(G_emb))
    G_use <- G_emb[, seq_len(d), drop = FALSE]
    
    L_use <- NULL
    if (!is.null(layout_reduction)) {
        L_emb <- Embeddings(sobj, reduction = layout_reduction)
        L_use <- L_emb
    }
    
    vcat(sprintf("Using %d dimensions from %s reduction", d, graph_reduction))
    
    # ---- Create SCE with updated metadata ----
    sce <- as.SingleCellExperiment(sobj, assay = DefaultAssay(sobj))
    
    # Update colData with converted metadata
    colData(sce)[[target_var]] <- md[[target_var]]
    colData(sce)[[batch_var]] <- md[[batch_var]]
    colData(sce)[[patient_var]] <- md[[patient_var]]
    colData(sce)[[cluster_var]] <- md[[cluster_var]]
    
    reducedDim(sce, "GRAPH") <- as.matrix(G_use)
    if (!is.null(L_use)) reducedDim(sce, "LAYOUT") <- as.matrix(L_use)
    
    # ---- Build Milo object ----
    vcat("Building kNN graph and neighborhoods...")
    milo <- Milo(sce)
    milo <- buildGraph(milo, k = k, d = ncol(G_use), reduced.dim = "GRAPH")
    milo <- makeNhoods(milo, prop = prop, k = k, refined = TRUE, reduced_dims = "GRAPH")
    
    vcat(sprintf("Created %d neighborhoods", ncol(nhoods(milo))))
    
    # ---- Get metadata from milo ----
    milo_md <- as.data.frame(colData(milo))
    
    # ---- Compute counts ----
    vcat("Computing neighborhood-by-sample counts...")
    counts <- .fast_count_cells(milo, samples = milo_md[[patient_var]])
    H <- nrow(counts)
    S <- ncol(counts)
    vcat(sprintf("Count matrix: %d neighborhoods x %d samples", H, S))
    
    # Filter small neighborhoods
    nh_sizes <- rowSums(counts)
    keep_nh <- nh_sizes >= min_cells
    if (sum(!keep_nh) > 0) {
        vcat(sprintf("Filtering out %d neighborhoods with < %d cells", sum(!keep_nh), min_cells))
    }
    
    # ---- Get majority cluster ----
    nh_mat <- miloR::nhoods(milo)
    cell_clusters <- milo_md[[cluster_var]]
    majority_cluster <- vapply(seq_len(ncol(nh_mat)), function(i) {
        idx <- as.logical(nh_mat[, i])
        labs <- cell_clusters[idx]
        if (length(labs) == 0) return(NA_character_)
        tbl <- table(labs)
        names(tbl)[which.max(tbl)]
    }, character(1))
    
    # ===== Create sample-level design matrix =====
    vcat("Creating sample-level design matrix...")
    
    # Get sample IDs from count matrix
    sample_ids <- colnames(counts)
    
    # Create sample metadata by aggregating cell-level data
    sample_design <- data.frame(
        sample_id = character(length(sample_ids)),
        stringsAsFactors = FALSE
    )
    
    # Initialize columns
    sample_design[[target_var]] <- character(length(sample_ids))
    sample_design[[batch_var]] <- character(length(sample_ids))
    
    # For each sample (patient), get their metadata
    for (i in seq_along(sample_ids)) {
        sid <- sample_ids[i]
        sample_design$sample_id[i] <- sid
        
        # Get all cells from this patient
        patient_cells <- which(milo_md[[patient_var]] == sid)
        
        if (length(patient_cells) > 0) {
            # Get unique values for this patient
            target_vals <- unique(milo_md[[target_var]][patient_cells])
            batch_vals <- unique(milo_md[[batch_var]][patient_cells])
            
            # Should be only one value per patient (nested design)
            sample_design[[target_var]][i] <- target_vals[1]
            sample_design[[batch_var]][i] <- batch_vals[1]
            
            if (length(target_vals) > 1) {
                warning(sprintf("Patient %s has multiple target values: %s", 
                              sid, paste(target_vals, collapse=", ")))
            }
            if (length(batch_vals) > 1) {
                warning(sprintf("Patient %s has multiple batch values: %s",
                              sid, paste(batch_vals, collapse=", ")))
            }
        }
    }
    
    # Convert to factors
    sample_design[[target_var]] <- factor(sample_design[[target_var]])
    sample_design[[batch_var]] <- factor(sample_design[[batch_var]])
    
    vcat(sprintf("Target variable levels: %s", 
                paste(levels(sample_design[[target_var]]), collapse=", ")))
    vcat(sprintf("Batch variable levels: %s", 
                paste(levels(sample_design[[batch_var]]), collapse=", ")))
    vcat(sprintf("Number of samples: %d", nrow(sample_design)))
    
    # Check design structure
    design_check <- table(sample_design[[target_var]], sample_design[[batch_var]])
    vcat("Design structure (target x batch):")
    print(design_check)
    
    is_nested <- all(apply(design_check > 0, 2, sum) == 1)  # Each batch in only one target
    
    if (is_nested) {
        vcat("WARNING: Batches are completely nested within target groups")
        vcat("Batch effects are confounded with target effects")
    }
    
    # ===== EdgeR analysis =====
    vcat("Running edgeR analysis...")
    
    rownames(sample_design) <- sample_ids
    
    # Check that we have at least 2 levels for target
    if (nlevels(sample_design[[target_var]]) < 2) {
        stop(sprintf("Target variable '%s' must have at least 2 levels, found: %s",
                    target_var, paste(levels(sample_design[[target_var]]), collapse=", ")))
    }
    
    # Build design matrix
    design_success <- FALSE
    mm <- NULL
    formula_str <- ""
    
    # Try both factors if not nested
    if (!is_nested && length(fixed_effects) > 1) {
        vcat("Attempting to include both target and batch...")
        formula_str <- paste("~", target_var, "+", batch_var)
        
        mm <- tryCatch({
            model.matrix(as.formula(formula_str), data = sample_design)
        }, error = function(e) {
            vcat(sprintf("  Failed: %s", e$message))
            NULL
        })
        
        if (!is.null(mm) && qr(mm)$rank == ncol(mm)) {
            design_success <- TRUE
            vcat("  Success: Including both factors")
        } else if (!is.null(mm)) {
            vcat("  Design is rank-deficient")
        }
    }
    
    # Fallback to target only
    if (!design_success) {
        if (is_nested) {
            vcat("Using target-only design due to nested structure")
        } else {
            vcat("Using target-only design")
        }
        
        # Use no-intercept model for cleaner contrasts
        formula_str <- paste("~ 0 +", target_var)
        mm <- model.matrix(as.formula(formula_str), data = sample_design)
        design_success <- TRUE
    }
    
    vcat(sprintf("Final design: %s", formula_str))
    vcat(sprintf("Design matrix: %d samples x %d coefficients", nrow(mm), ncol(mm)))
    vcat(sprintf("Coefficients: %s", paste(colnames(mm), collapse=", ")))
    
    # EdgeR workflow
    dge <- DGEList(counts = counts[keep_nh, , drop = FALSE])
    dge <- calcNormFactors(dge, method = "TMM")
    
    # Estimate dispersions
    vcat("Estimating dispersions...")
    dge <- estimateDisp(dge, design = mm, robust = TRUE)
    vcat(sprintf("Common dispersion: %.4f", dge$common.dispersion))
    
    # Fit model
    vcat("Fitting GLM...")
    fit <- glmQLFit(dge, design = mm, robust = TRUE)
    
    # Test for differential abundance
    if (grepl("~ 0 \\+", formula_str)) {
        # No intercept model - test contrast
        target_levels <- levels(sample_design[[target_var]])
        
        if (length(target_levels) == 2) {
            # Create contrast: group2 - group1
            contrast_str <- paste0(colnames(mm)[2], " - ", colnames(mm)[1])
            vcat(sprintf("Testing contrast: %s", contrast_str))
            
            contrast_vec <- numeric(ncol(mm))
            contrast_vec[1] <- -1
            contrast_vec[2] <- 1
            
            qlf <- glmQLFTest(fit, contrast = contrast_vec)
        } else {
            vcat("Testing all groups against baseline")
            qlf <- glmQLFTest(fit, coef = 2:ncol(mm))
        }
    } else {
        # Model with intercept
        target_coef <- grep(target_var, colnames(mm), value = TRUE)[1]
        vcat(sprintf("Testing coefficient: %s", target_coef))
        qlf <- glmQLFTest(fit, coef = target_coef)
    }
    
    # Build results
    nh_idx <- which(keep_nh)
    da.res <- tibble(
        nhood = nh_idx,
        logFC = qlf$table$logFC,
        logCPM = qlf$table$logCPM,
        F = qlf$table$F,
        pval = qlf$table$PValue
    )
    
    # Spatial FDR
    vcat("Computing spatial FDR...")
    da.res <- graphSpatialFDR(milo, da.res, pvalues = "pval", indices = "nhood")
    
    da.res$major_cluster <- majority_cluster[da.res$nhood]
    da.res <- da.res %>%
        mutate(is_sig = spatialFDR < fdr_thresh & abs(logFC) > logfc_thresh) %>%
        arrange(spatialFDR, desc(abs(logFC)))
    
    # Cluster summary
    cluster_summary <- if (sum(da.res$is_sig, na.rm = TRUE) > 0) {
        da.res %>%
            filter(is_sig) %>%
            group_by(major_cluster) %>%
            summarise(
                n_sig_nhoods = n(),
                mean_logFC = mean(logFC),
                median_logFC = median(logFC),
                min_spatialFDR = min(spatialFDR),
                .groups = "drop"
            ) %>%
            arrange(desc(n_sig_nhoods))
    } else {
        tibble()
    }
    
    vcat(sprintf("\nResults: %d significant neighborhoods out of %d tested",
                sum(da.res$is_sig, na.rm = TRUE), nrow(da.res)))
    
    if (is_nested) {
        vcat("\nNOTE: Batch effects could not be separated due to nested design")
    }
    
    return(list(
        mode = "edgeR_fixed_effects",
        formula = formula_str,
        milo = milo,
        da_nhood_with_cluster = da.res,
        da_by_cluster = cluster_summary,
        design_matrix = mm,
        sample_design = sample_design,
        design_info = list(
            is_nested = is_nested,
            n_samples = S,
            n_neighborhoods = H,
            n_filtered = sum(!keep_nh)
        )
    ))
}

milo_opus7 <- function(
    sobj,
    target_var,
    batch_var, 
    patient_var,
    cluster_var,
    graph_reduction,                 
    layout_reduction = NULL,         
    k = 30,
    d = NULL,                        
    prop = 0.1,
    fixed_effects  = NULL,
    random_effects = NULL,
    glmm_family    = c("nb","poisson"),
    p_adj_method   = "BH",
    fdr_thresh     = 0.1,
    logfc_thresh   = 0.5,
    seed           = 1,
    use_future     = FALSE,
    workers        = max(1, parallel::detectCores()-1),
    verbose        = TRUE,
    min_cells      = 3
){
    suppressPackageStartupMessages({
        library(Seurat)
        library(SingleCellExperiment)
        library(miloR)
        library(Matrix)
        library(dplyr)
        library(tibble)
        library(edgeR)
    })
    set.seed(seed)
    
    vcat <- function(...) if(verbose) cat(..., "\n")
    
    # ---- Checks ----
    req_cols <- c(target_var, batch_var, patient_var, cluster_var)
    stopifnot("Required columns missing" = all(req_cols %in% colnames(sobj@meta.data)))
    stopifnot("Graph reduction not found" = graph_reduction %in% Reductions(sobj))
    
    if (is.null(fixed_effects)) fixed_effects <- c(target_var, batch_var)
    glmm_family <- match.arg(glmm_family)
    use_glmm <- length(random_effects) > 0
    
    vcat("Starting Milo DA analysis...")
    vcat("Mode:", ifelse(use_glmm, "GLMM with random effects", "edgeR fixed effects"))
    
    # ---- First, create the sample-level design DIRECTLY from Seurat object ----
    vcat("Extracting sample-level information...")
    
    # Get unique patient-level information
    patient_info <- sobj@meta.data[, c(patient_var, target_var, batch_var)] %>%
        unique() %>%
        as.data.frame()
    
    # Convert to appropriate types
    patient_info[[patient_var]] <- as.character(patient_info[[patient_var]])
    patient_info[[target_var]] <- as.character(patient_info[[target_var]])
    patient_info[[batch_var]] <- as.character(patient_info[[batch_var]])
    
    vcat(sprintf("Found %d unique patients", nrow(patient_info)))
    vcat(sprintf("Target values in data: %s", 
                paste(sort(unique(patient_info[[target_var]])), collapse=", ")))
    vcat(sprintf("Batch values in data: %s", 
                paste(sort(unique(patient_info[[batch_var]])), collapse=", ")))
    
    # ---- Embeddings ----
    G_emb <- Embeddings(sobj, reduction = graph_reduction)
    if (is.null(d)) d <- ncol(G_emb)
    d <- min(d, ncol(G_emb))
    G_use <- G_emb[, seq_len(d), drop = FALSE]
    
    L_use <- NULL
    if (!is.null(layout_reduction)) {
        L_emb <- Embeddings(sobj, reduction = layout_reduction)
        L_use <- L_emb
    }
    
    vcat(sprintf("Using %d dimensions from %s reduction", d, graph_reduction))
    
    # ---- Create SCE ----
    sce <- as.SingleCellExperiment(sobj, assay = DefaultAssay(sobj))
    reducedDim(sce, "GRAPH") <- as.matrix(G_use)
    if (!is.null(L_use)) reducedDim(sce, "LAYOUT") <- as.matrix(L_use)
    
    # ---- Build Milo ----
    vcat("Building kNN graph and neighborhoods...")
    milo <- Milo(sce)
    milo <- buildGraph(milo, k = k, d = ncol(G_use), reduced.dim = "GRAPH")
    milo <- makeNhoods(milo, prop = prop, k = k, refined = TRUE, reduced_dims = "GRAPH")
    
    vcat(sprintf("Created %d neighborhoods", ncol(nhoods(milo))))
    
    # ---- Get metadata ----
    milo_md <- as.data.frame(colData(milo))
    
    # ---- Compute counts ----
    vcat("Computing neighborhood-by-sample counts...")
    counts <- .fast_count_cells(milo, samples = milo_md[[patient_var]])
    H <- nrow(counts)
    S <- ncol(counts)
    vcat(sprintf("Count matrix: %d neighborhoods x %d samples", H, S))
    
    # ---- Create design matrix using patient_info ----
    sample_ids <- colnames(counts)
    
    # Match the sample IDs to patient_info
    sample_idx <- match(sample_ids, patient_info[[patient_var]])
    
    if (any(is.na(sample_idx))) {
        missing <- sample_ids[is.na(sample_idx)]
        stop(sprintf("Could not find patient info for samples: %s", 
                    paste(missing, collapse=", ")))
    }
    
    sample_design <- patient_info[sample_idx, ]
    rownames(sample_design) <- sample_ids
    
    # Convert to factors
    sample_design[[target_var]] <- factor(sample_design[[target_var]])
    sample_design[[batch_var]] <- factor(sample_design[[batch_var]])
    
    vcat(sprintf("\nSample design summary:"))
    vcat(sprintf("  Target (%s) levels: %s", target_var,
                paste(levels(sample_design[[target_var]]), collapse=", ")))
    vcat(sprintf("  Batch (%s) levels: %s", batch_var,
                paste(levels(sample_design[[batch_var]]), collapse=", ")))
    
    # Check design
    design_table <- table(sample_design[[target_var]], sample_design[[batch_var]])
    vcat("\nDesign structure (rows=target, cols=batch):")
    print(design_table)
    
    # ---- Filter neighborhoods ----
    nh_sizes <- rowSums(counts)
    keep_nh <- nh_sizes >= min_cells
    if (sum(!keep_nh) > 0) {
        vcat(sprintf("Filtering %d neighborhoods with < %d cells", sum(!keep_nh), min_cells))
    }
    
    # ---- Majority cluster ----
    nh_mat <- miloR::nhoods(milo)
    cell_clusters <- milo_md[[cluster_var]]
    majority_cluster <- vapply(seq_len(ncol(nh_mat)), function(i) {
        idx <- as.logical(nh_mat[, i])
        labs <- cell_clusters[idx]
        if (length(labs) == 0) return(NA_character_)
        tbl <- table(labs)
        names(tbl)[which.max(tbl)]
    }, character(1))
    
    # ---- EdgeR analysis ----
    vcat("\nRunning edgeR analysis...")
    
    # Check levels
    if (nlevels(sample_design[[target_var]]) < 2) {
        stop(sprintf("Target variable needs at least 2 levels, found: %s",
                    paste(levels(sample_design[[target_var]]), collapse=", ")))
    }
    
    # Simple design with target only (since batch is nested)
    formula_str <- paste("~ 0 +", target_var)
    vcat(sprintf("Using formula: %s", formula_str))
    
    mm <- model.matrix(as.formula(formula_str), data = sample_design)
    vcat(sprintf("Design matrix: %d x %d", nrow(mm), ncol(mm)))
    vcat(sprintf("Column names: %s", paste(colnames(mm), collapse=", ")))
    
    # EdgeR workflow
    dge <- DGEList(counts = counts[keep_nh, , drop = FALSE])
    dge <- calcNormFactors(dge, method = "TMM")
    
    vcat("Estimating dispersions...")
    dge <- estimateDisp(dge, design = mm, robust = TRUE)
    vcat(sprintf("Common dispersion: %.4f", dge$common.dispersion))
    
    vcat("Fitting GLM...")
    fit <- glmQLFit(dge, design = mm, robust = TRUE)
    
    # Test contrast
    if (ncol(mm) == 2) {
        # Two groups - simple contrast
        contrast_name <- paste(colnames(mm)[2], "vs", colnames(mm)[1])
        vcat(sprintf("Testing: %s", contrast_name))
        
        contrast <- numeric(ncol(mm))
        contrast[1] <- -1
        contrast[2] <- 1
        qlf <- glmQLFTest(fit, contrast = contrast)
    } else {
        # Multiple groups
        vcat("Testing all groups")
        qlf <- glmQLFTest(fit, coef = 2:ncol(mm))
    }
    
    # Results
    nh_idx <- which(keep_nh)
    da.res <- tibble(
        nhood = nh_idx,
        logFC = qlf$table$logFC,
        logCPM = qlf$table$logCPM,
        F = qlf$table$F,
        pval = qlf$table$PValue
    )
    
    vcat("Computing spatial FDR...")
    da.res <- graphSpatialFDR(milo, da.res, pvalues = "pval", indices = "nhood")
    
    da.res$major_cluster <- majority_cluster[da.res$nhood]
    da.res <- da.res %>%
        mutate(is_sig = spatialFDR < fdr_thresh & abs(logFC) > logfc_thresh) %>%
        arrange(spatialFDR, desc(abs(logFC)))
    
    # Summary
    cluster_summary <- if (sum(da.res$is_sig, na.rm = TRUE) > 0) {
        da.res %>%
            filter(is_sig) %>%
            group_by(major_cluster) %>%
            summarise(
                n_sig_nhoods = n(),
                mean_logFC = mean(logFC),
                median_logFC = median(logFC),
                min_spatialFDR = min(spatialFDR),
                .groups = "drop"
            ) %>%
            arrange(desc(n_sig_nhoods))
    } else {
        tibble()
    }
    
    vcat(sprintf("\nFound %d significant neighborhoods", sum(da.res$is_sig, na.rm = TRUE)))
    if (nrow(cluster_summary) > 0) {
        vcat(sprintf("Affecting %d cell types", nrow(cluster_summary)))
    }
    
    return(list(
        mode = "edgeR_fixed_effects",
        formula = formula_str,
        milo = milo,
        da_nhood_with_cluster = da.res,
        da_by_cluster = cluster_summary,
        design_matrix = mm,
        sample_design = sample_design
    ))
}

.fast_count_cells <- function(milo, samples) {
  stopifnot(inherits(milo, "Milo"))
  X <- miloR::nhoods(milo)           # cells x H (dgCMatrix)
  stopifnot(nrow(X) == ncol(milo))
  
  if (!is.factor(samples)) samples <- factor(samples)
  
  S <- Matrix::sparse.model.matrix(~ 0 + samples)   # cells x S
  
  # H x S 매트릭스 계산
  counts_matrix <- Matrix::t(X) %*% S |> as.matrix() 
  
  # !! 컬럼명 오류 수정 !!
  if (ncol(counts_matrix) == length(levels(samples))) {
      colnames(counts_matrix) <- levels(samples)
  } else {
      warning("Count matrix columns do not match sample levels. Using default names.")
  }

  return(counts_matrix)
}

milo_opus_final <- function(
    sobj,
    target_var,
    batch_var, 
    patient_var,
    cluster_var,
    graph_reduction,              
    layout_reduction = NULL,      
    k = 30,
    d = NULL,                     
    prop = 0.1,
    fixed_effects  = NULL,
    random_effects = NULL,
    glmm_family    = c("nb","poisson"),
    p_adj_method   = "BH",
    fdr_thresh     = 0.1,
    logfc_thresh   = 0.5,
    seed           = 1,
    use_future     = FALSE,
    workers        = max(1, parallel::detectCores()-1),
    verbose        = TRUE,
    min_cells      = 3
){
    suppressPackageStartupMessages({
        library(Seurat)
        library(SingleCellExperiment)
        library(miloR)
        library(Matrix)
        library(dplyr)
        library(tibble)
        library(edgeR)
    })
    set.seed(seed)
    
    vcat <- function(...) if(verbose) cat(..., "\n")
    
    # ---- Checks ----
    req_cols <- c(target_var, batch_var, patient_var, cluster_var)
    stopifnot("Required columns missing" = all(req_cols %in% colnames(sobj@meta.data)))
    stopifnot("Graph reduction not found" = graph_reduction %in% Reductions(sobj))
    
    if (is.null(fixed_effects)) fixed_effects <- c(target_var, batch_var)
    glmm_family <- match.arg(glmm_family)
    use_glmm <- length(random_effects) > 0
    
    vcat("Starting Milo DA analysis...")
    vcat("Mode:", ifelse(use_glmm, "GLMM with random effects", "edgeR fixed effects"))
    
    # -------------------------------------------------------------------
    # !! 핵심 수정 사항 (BUGFIX) !!
    # - sce와 patient_info가 생성되기 *전*에 원본 sobj 메타데이터 타입을 통일합니다.
    # - 이것이 ver7의 'match' 오류와 ver6의 'nlevels' 오류를 모두 해결합니다.
    # -------------------------------------------------------------------
    vcat("Coercing metadata columns to character in source object...")
    sobj@meta.data[[patient_var]] <- as.character(sobj@meta.data[[patient_var]])
    sobj@meta.data[[target_var]]  <- as.character(sobj@meta.data[[target_var]])
    sobj@meta.data[[batch_var]]   <- as.character(sobj@meta.data[[batch_var]])
    sobj@meta.data[[cluster_var]] <- as.character(sobj@meta.data[[cluster_var]])
    
    # ---- First, create the sample-level design DIRECTLY from Seurat object ----
    vcat("Extracting sample-level information...")
    
    # 이제 sobj@meta.data가 이미 character이므로 patient_info는 안전합니다.
    patient_info <- sobj@meta.data[, c(patient_var, target_var, batch_var)] %>%
        unique() %>%
        as.data.frame()
    
    vcat(sprintf("Found %d unique patients", nrow(patient_info)))
    vcat(sprintf("Target values in data: %s", 
                 paste(sort(unique(patient_info[[target_var]])), collapse=", ")))
    vcat(sprintf("Batch values in data: %s", 
                 paste(sort(unique(patient_info[[batch_var]])), collapse=", ")))
    
    # ---- Embeddings ----
    G_emb <- Embeddings(sobj, reduction = graph_reduction)
    if (is.null(d)) d <- ncol(G_emb)
    d <- min(d, ncol(G_emb))
    G_use <- G_emb[, seq_len(d), drop = FALSE]
    
    L_use <- NULL
    if (!is.null(layout_reduction)) {
        L_emb <- Embeddings(sobj, reduction = layout_reduction)
        L_use <- L_emb
    }
    
    vcat(sprintf("Using %d dimensions from %s reduction", d, graph_reduction))
    
    # ---- Create SCE ----
    # sobj의 메타데이터가 이미 수정되었으므로, sce의 colData도 올바른 타입을 갖게 됩니다.
    sce <- as.SingleCellExperiment(sobj, assay = DefaultAssay(sobj))
    reducedDim(sce, "GRAPH") <- as.matrix(G_use)
    if (!is.null(L_use)) reducedDim(sce, "LAYOUT") <- as.matrix(L_use)
    
    # ---- Build Milo ----
    vcat("Building kNN graph and neighborhoods...")
    milo <- Milo(sce)
    milo <- buildGraph(milo, k = k, d = ncol(G_use), reduced.dim = "GRAPH")
    milo <- makeNhoods(milo, prop = prop, k = k, refined = TRUE, reduced_dims = "GRAPH")
    
    vcat(sprintf("Created %d neighborhoods", ncol(nhoods(milo))))
    
    # ---- Get metadata ----
    # milo_md <- as.data.frame(colData(milo)) # 이 라인은 이제 필요 없습니다.

    # ---- Compute counts ----
    vcat("Computing neighborhood-by-sample counts (using official countCells API)...")

    # !! 공식 API인 countCells로 교체 !!
    # 1. countCells는 milo 객체를 수정하고,
    # 2. 'sample' 인자에는 컬럼 이름 (patient_var, 즉 "hos_no")을 전달합니다.
    milo <- countCells(milo, 
                      meta.data = as.data.frame(colData(milo)), 
                      sample = patient_var)

    # 3. 계산된 'counts'를 milo 객체에서 가져옵니다.
    counts <- nhoodCounts(milo)
        
    H <- nrow(counts)
    S <- ncol(counts)
    vcat(sprintf("Count matrix: %d neighborhoods x %d samples", H, S))
    
    # ---- Create design matrix using patient_info ----
    sample_ids <- colnames(counts)
    
    # -------------------------------------------------------------------
    # !! 'match' 오류 해결 !!
    # sample_ids (from counts < sce < sobj)와 
    # patient_info[[patient_var]] (from sobj)가 
    # 모두 동일한 character 소스에서 유래했으므로 매칭이 성공합니다.
    # -------------------------------------------------------------------
    sample_idx <- match(sample_ids, patient_info[[patient_var]])
    
    if (any(is.na(sample_idx))) {
        missing <- sample_ids[is.na(sample_idx)]
        stop(sprintf("Could not find patient info for samples: %s", 
                     paste(missing, collapse=", ")))
    }
    
    sample_design <- patient_info[sample_idx, ]
    rownames(sample_design) <- sample_ids
    
    # Convert to factors
    sample_design[[target_var]] <- factor(sample_design[[target_var]])
    sample_design[[batch_var]] <- factor(sample_design[[batch_var]])
    
    vcat(sprintf("\nSample design summary:"))
    vcat(sprintf("  Target (%s) levels: %s", target_var,
                 paste(levels(sample_design[[target_var]]), collapse=", ")))
    vcat(sprintf("  Batch (%s) levels: %s", batch_var,
                 paste(levels(sample_design[[batch_var]]), collapse=", ")))
    
    # Check design
    design_table <- table(sample_design[[target_var]], sample_design[[batch_var]])
    vcat("\nDesign structure (rows=target, cols=batch):")
    print(design_table)
    
    # ---- Filter neighborhoods ----
    nh_sizes <- rowSums(counts)
    keep_nh <- nh_sizes >= min_cells
    if (sum(!keep_nh) > 0) {
        vcat(sprintf("Filtering %d neighborhoods with < %d cells", sum(!keep_nh), min_cells))
    }

    # ---- Majority cluster ----

    # !! 여기(Majority cluster 섹션 상단)에 milo_md 정의를 다시 추가합니다 !!
    milo_md <- as.data.frame(colData(milo))

    nh_mat <- miloR::nhoods(milo)
    cell_clusters <- milo_md[[cluster_var]]  # <-- 이제 'milo_md'를 찾을 수 있습니다.
    majority_cluster <- vapply(seq_len(ncol(nh_mat)), function(i) {
        idx <- as.logical(nh_mat[, i])
        labs <- cell_clusters[idx]
        if (length(labs) == 0) return(NA_character_)
        tbl <- table(labs)
        names(tbl)[which.max(tbl)]
    }, character(1))
    
    # ---- EdgeR analysis ----
    vcat("\nRunning edgeR analysis...")
    
    # -------------------------------------------------------------------
    # !! 'nlevels' 오류 해결 !!
    # sample_design이 patient_info로부터 올바르게 생성되었으므로,
    # nlevels(sample_design[[target_var]])가 2 이상이 됩니다.
    # -------------------------------------------------------------------
    if (nlevels(sample_design[[target_var]]) < 2) {
        stop(sprintf("Target variable '%s' must have at least 2 levels, found: %s",
                     target_var, 
                     paste(levels(sample_design[[target_var]]), collapse=", ")))
    }
    
    # Simple design with target only (since batch is nested)
    formula_str <- paste("~ 0 +", target_var)
    vcat(sprintf("Using formula: %s", formula_str))
    
    mm <- model.matrix(as.formula(formula_str), data = sample_design)
    vcat(sprintf("Design matrix: %d x %d", nrow(mm), ncol(mm)))
    vcat(sprintf("Column names: %s", paste(colnames(mm), collapse=", ")))
    
    # EdgeR workflow
    dge <- DGEList(counts = counts[keep_nh, , drop = FALSE])
    dge <- calcNormFactors(dge, method = "TMM")
    
    vcat("Estimating dispersions...")
    dge <- estimateDisp(dge, design = mm, robust = TRUE)
    vcat(sprintf("Common dispersion: %.4f", dge$common.dispersion))
    
    vcat("Fitting GLM...")
    fit <- glmQLFit(dge, design = mm, robust = TRUE)
    
    # Test contrast
    if (ncol(mm) == 2) {
        # Two groups - simple contrast
        contrast_name <- paste(colnames(mm)[2], "vs", colnames(mm)[1])
        vcat(sprintf("Testing: %s", contrast_name))
        
        contrast <- numeric(ncol(mm))
        contrast[1] <- -1
        contrast[2] <- 1
        qlf <- glmQLFTest(fit, contrast = contrast)
    } else {
        # Multiple groups
        vcat("Testing all groups")
        qlf <- glmQLFTest(fit, coef = 2:ncol(mm))
    }
    
    # Results
    nh_idx <- which(keep_nh)
    da.res <- tibble(
        nhood = nh_idx,
        logFC = qlf$table$logFC,
        logCPM = qlf$table$logCPM,
        F = qlf$table$F,
        pval = qlf$table$PValue
    )
    
    vcat("Computing spatial FDR...")
    # !! k=k 와 reduced.dim="GRAPH" 를 모두 추가 !!
    da.res <- graphSpatialFDR(milo, da.res, 
                              pvalues = "pval", 
                              indices = "nhood", 
                              k = k,
                              reduced.dim = "GRAPH")
    
    da.res$major_cluster <- majority_cluster[da.res$nhood]
    da.res <- da.res %>%
        mutate(is_sig = spatialFDR < fdr_thresh & abs(logFC) > logfc_thresh) %>%
        arrange(spatialFDR, desc(abs(logFC)))
    
    # Summary
    cluster_summary <- if (sum(da.res$is_sig, na.rm = TRUE) > 0) {
        da.res %>%
            filter(is_sig) %>%
            group_by(major_cluster) %>%
            summarise(
                n_sig_nhoods = n(),
                mean_logFC = mean(logFC),
                median_logFC = median(logFC),
                min_spatialFDR = min(spatialFDR),
                .groups = "drop"
            ) %>%
            arrange(desc(n_sig_nhoods))
    } else {
        tibble()
    }
    
    vcat(sprintf("\nFound %d significant neighborhoods", sum(da.res$is_sig, na.rm = TRUE)))
    if (nrow(cluster_summary) > 0) {
        vcat(sprintf("Affecting %d cell types", nrow(cluster_summary)))
    }
    
    return(list(
        mode = "edgeR_fixed_effects",
        formula = formula_str,
        milo = milo,
        da_nhood_with_cluster = da.res,
        da_by_cluster = cluster_summary,
        design_matrix = mm,
        sample_design = sample_design
    ))
}

milo_opus_final2 <- function(
    sobj,
    target_var,
    batch_var, 
    patient_var,
    cluster_var,
    graph_reduction,              
    layout_reduction = NULL,      
    k = 30,
    d = NULL,                     
    prop = 0.1,
    fixed_effects  = NULL,
    random_effects = NULL,
    glmm_family    = c("nb","poisson"),
    p_adj_method   = "BH",
    fdr_thresh     = 0.1,
    logfc_thresh   = 0.5,
    seed           = 1,
    use_future     = FALSE,
    workers        = max(1, parallel::detectCores()-1),
    verbose        = TRUE,
    min_cells      = 3
){
    suppressPackageStartupMessages({
        library(Seurat)
        library(SingleCellExperiment)
        library(miloR)
        library(Matrix)
        library(dplyr)
        library(tibble)
        library(edgeR)
    })
    set.seed(seed)
    
    vcat <- function(...) if(verbose) cat(..., "\n")
    
    # ---- Checks ----
    req_cols <- c(target_var, batch_var, patient_var, cluster_var)
    stopifnot("Required columns missing" = all(req_cols %in% colnames(sobj@meta.data)))
    stopifnot("Graph reduction not found" = graph_reduction %in% Reductions(sobj))
    
    if (is.null(fixed_effects)) fixed_effects <- c(target_var, batch_var)
    glmm_family <- match.arg(glmm_family)
    use_glmm <- length(random_effects) > 0
    
    vcat("Starting Milo DA analysis...")
    vcat("Mode:", ifelse(use_glmm, "GLMM with random effects", "edgeR fixed effects"))
    
    # ---- 1. 타입 통일 ----
    vcat("Coercing metadata columns to character in source object...")
    sobj@meta.data[[patient_var]] <- as.character(sobj@meta.data[[patient_var]])
    sobj@meta.data[[target_var]]  <- as.character(sobj@meta.data[[target_var]])
    sobj@meta.data[[batch_var]]   <- as.character(sobj@meta.data[[batch_var]])
    sobj@meta.data[[cluster_var]] <- as.character(sobj@meta.data[[cluster_var]])
    
    # ---- 2. 샘플 정보 생성 ----
    vcat("Extracting sample-level information...")
    patient_info <- sobj@meta.data[, c(patient_var, target_var, batch_var)] %>%
        unique() %>%
        as.data.frame()
    
    vcat(sprintf("Found %d unique patients", nrow(patient_info)))
    vcat(sprintf("Target values in data: %s", 
                 paste(sort(unique(patient_info[[target_var]])), collapse=", ")))
    vcat(sprintf("Batch values in data: %s", 
                 paste(sort(unique(patient_info[[batch_var]])), collapse=", ")))
    
    # ---- Embeddings ----
    G_emb <- Embeddings(sobj, reduction = graph_reduction)
    if (is.null(d)) d <- ncol(G_emb)
    d <- min(d, ncol(G_emb))
    G_use <- G_emb[, seq_len(d), drop = FALSE]
    L_use <- NULL
    if (!is.null(layout_reduction)) {
        L_emb <- Embeddings(sobj, reduction = layout_reduction)
        L_use <- L_emb
    }
    vcat(sprintf("Using %d dimensions from %s reduction", d, graph_reduction))
    
    # ---- Create SCE ----
    sce <- as.SingleCellExperiment(sobj, assay = DefaultAssay(sobj))
    reducedDim(sce, "GRAPH") <- as.matrix(G_use)
    if (!is.null(L_use)) reducedDim(sce, "LAYOUT") <- as.matrix(L_use)
    
    # ---- Build Milo ----
    vcat("Building kNN graph and neighborhoods...")
    milo <- Milo(sce)
    milo <- buildGraph(milo, k = k, d = ncol(G_use), reduced.dim = "GRAPH")
    milo <- makeNhoods(milo, prop = prop, k = k, refined = TRUE, reduced_dims = "GRAPH")
    
    vcat(sprintf("Created %d neighborhoods", ncol(nhoods(milo))))
    
    # ---- Get metadata ----
    milo_md <- as.data.frame(colData(milo))
    
    # ---- 3. Counts 계산 (수정된 .fast_count_cells 사용) ----
    vcat("Computing neighborhood-by-sample counts (using fixed .fast_count_cells)...")
    counts <- .fast_count_cells(milo, samples = milo_md[[patient_var]])
    
    H <- nrow(counts)
    S <- ncol(counts)
    vcat(sprintf("Count matrix: %d neighborhoods x %d samples", H, S))
    
    # ---- 4. Design Matrix 생성 (매칭 오류 해결됨) ----
    sample_ids <- colnames(counts)
    sample_idx <- match(sample_ids, patient_info[[patient_var]])
    
    if (any(is.na(sample_idx))) {
        missing <- sample_ids[is.na(sample_idx)]
        stop(sprintf("Could not find patient info for samples: %s", 
                     paste(missing, collapse=", ")))
    }
    
    sample_design <- patient_info[sample_idx, ]
    rownames(sample_design) <- sample_ids
    
    sample_design[[target_var]] <- factor(sample_design[[target_var]])
    sample_design[[batch_var]] <- factor(sample_design[[batch_var]])
    
    vcat(sprintf("\nSample design summary:"))
    vcat(sprintf("  Target (%s) levels: %s", target_var,
                 paste(levels(sample_design[[target_var]]), collapse=", ")))
    vcat(sprintf("  Batch (%s) levels: %s", batch_var,
                 paste(levels(sample_design[[batch_var]]), collapse=", ")))
    
    design_table <- table(sample_design[[target_var]], sample_design[[batch_var]])
    vcat("\nDesign structure (rows=target, cols=batch):")
    print(design_table)
    
    # ---- Filter neighborhoods ----
    nh_sizes <- rowSums(counts)
    keep_nh <- nh_sizes >= min_cells
    if (sum(!keep_nh) > 0) {
        vcat(sprintf("Filtering %d neighborhoods with < %d cells", sum(!keep_nh), min_cells))
    }
    
    # ---- 5. Majority cluster (milo_md 오류 수정됨) ----
    nh_mat <- miloR::nhoods(milo)
    cell_clusters <- milo_md[[cluster_var]]
    majority_cluster <- vapply(seq_len(ncol(nh_mat)), function(i) {
        idx <- as.logical(nh_mat[, i])
        labs <- cell_clusters[idx]
        if (length(labs) == 0) return(NA_character_)
        tbl <- table(labs)
        names(tbl)[which.max(tbl)]
    }, character(1))
    
    # ---- EdgeR analysis ----
    vcat("\nRunning edgeR analysis...")
    
    if (nlevels(sample_design[[target_var]]) < 2) {
        stop(sprintf("Target variable '%s' must have at least 2 levels, found: %s",
                     target_var, 
                     paste(levels(sample_design[[target_var]]), collapse=", ")))
    }
    
    formula_str <- paste("~ 0 +", target_var)
    vcat(sprintf("Using formula: %s", formula_str))
    
    mm <- model.matrix(as.formula(formula_str), data = sample_design)
    vcat(sprintf("Design matrix: %d x %d", nrow(mm), ncol(mm)))
    vcat(sprintf("Column names: %s", paste(colnames(mm), collapse=", ")))
    
    dge <- DGEList(counts = counts[keep_nh, , drop = FALSE])
    dge <- calcNormFactors(dge, method = "TMM")
    
    vcat("Estimating dispersions...")
    dge <- estimateDisp(dge, design = mm, robust = TRUE)
    vcat(sprintf("Common dispersion: %.4f", dge$common.dispersion))
    
    vcat("Fitting GLM...")
    fit <- glmQLFit(dge, design = mm, robust = TRUE)
    
    if (ncol(mm) == 2) {
        contrast_name <- paste(colnames(mm)[2], "vs", colnames(mm)[1])
        vcat(sprintf("Testing: %s", contrast_name))
        contrast <- numeric(ncol(mm)); contrast[1] <- -1; contrast[2] <- 1
        qlf <- glmQLFTest(fit, contrast = contrast)
    } else {
        vcat("Testing all groups")
        qlf <- glmQLFTest(fit, coef = 2:ncol(mm))
    }
    
    nh_idx <- which(keep_nh)
    da.res <- tibble(
        nhood = nh_idx,
        logFC = qlf$table$logFC,
        logCPM = qlf$table$logCPM,
        F = qlf$table$F,
        pval = qlf$table$PValue
    )
    
    # -------------------------------------------------------------------
    # !! 7. [신규 수정] NA/Inf 오류 수정 !!
    # -------------------------------------------------------------------
    vcat(sprintf("GLM results: %d rows", nrow(da.res)))
    da.res <- da.res %>%
        filter(!is.na(pval) & is.finite(logFC))
    
    if(nrow(da.res) == 0) {
        stop("No valid neighborhoods remained after filtering NA/Inf values. Check GLM results.")
    }
    vcat(sprintf("Filtered GLM results (NA/Inf removed): %d rows remain", nrow(da.res)))

    # ---- 6. Spatial FDR (k=k 및 reduced.dim 오류 수정됨) ----
    vcat("Computing spatial FDR...")
    da.res <- graphSpatialFDR(milo, da.res, 
                              pvalues = "pval", 
                              indices = "nhood", 
                              k = k,
                              reduced.dim = "GRAPH") # <-- 모든 수정 사항 적용
    
    da.res$major_cluster <- majority_cluster[da.res$nhood]
    da.res <- da.res %>%
        mutate(is_sig = spatialFDR < fdr_thresh & abs(logFC) > logfc_thresh) %>%
        arrange(spatialFDR, desc(abs(logFC)))
    
    cluster_summary <- if (sum(da.res$is_sig, na.rm = TRUE) > 0) {
        da.res %>%
            filter(is_sig) %>%
            group_by(major_cluster) %>%
            summarise(
                n_sig_nhoods = n(),
                mean_logFC = mean(logFC),
                median_logFC = median(logFC),
                min_spatialFDR = min(spatialFDR),
                .groups = "drop"
            ) %>%
            arrange(desc(n_sig_nhoods))
    } else {
        tibble()
    }
    
    vcat(sprintf("\nFound %d significant neighborhoods", sum(da.res$is_sig, na.rm = TRUE)))
    if (nrow(cluster_summary) > 0) {
        vcat(sprintf("Affecting %d cell types", nrow(cluster_summary)))
    }
    
    return(list(
        mode = "edgeR_fixed_effects",
        formula = formula_str,
        milo = milo,
        da_nhood_with_cluster = da.res,
        da_by_cluster = cluster_summary,
        design_matrix = mm,
        sample_design = sample_design
    ))
}

#' --------------------------------------------------
#' 1. 속도 개선을 위한 보조 함수 (컬럼명 오류 수정됨)
#' --------------------------------------------------
.fast_count_cells <- function(milo, samples) {
  stopifnot(inherits(milo, "Milo"))
  X <- miloR::nhoods(milo)           # cells x H (dgCMatrix)
  stopifnot(nrow(X) == ncol(milo))
  
  if (!is.factor(samples)) samples <- factor(samples)
  
  S <- Matrix::sparse.model.matrix(~ 0 + samples)   # cells x S
  
  # H x S 매트릭스 계산
  counts_matrix <- Matrix::t(X) %*% S |> as.matrix() 
  
  # !! 컬럼명 오류 수정 !!
  if (ncol(counts_matrix) == length(levels(samples))) {
      colnames(counts_matrix) <- levels(samples)
  } else {
      warning("Count matrix columns do not match sample levels. Using default names.")
  }

  return(counts_matrix)
}


#' --------------------------------------------------
#' 2. 최종 Milo Opus 함수 (모든 수정 사항 통합)
#' --------------------------------------------------
milo_opus_final3 <- function(
    sobj,
    target_var,
    batch_var, 
    patient_var,
    cluster_var,
    graph_reduction,              
    layout_reduction = NULL,      
    k = 30,
    d = NULL,                     
    prop = 0.1,
    fixed_effects  = NULL,
    random_effects = NULL,
    glmm_family    = c("nb","poisson"),
    p_adj_method   = "BH",
    fdr_thresh     = 0.1,
    logfc_thresh   = 0.5,
    seed           = 1,
    use_future     = FALSE,
    workers        = max(1, parallel::detectCores()-1),
    verbose        = TRUE,
    min_cells      = 3
){
    suppressPackageStartupMessages({
        library(Seurat)
        library(SingleCellExperiment)
        library(miloR)
        library(Matrix)
        library(dplyr)
        library(tibble)
        library(edgeR)
    })
    set.seed(seed)
    
    vcat <- function(...) if(verbose) cat(..., "\n")
    
    # ---- Checks ----
    req_cols <- c(target_var, batch_var, patient_var, cluster_var)
    stopifnot("Required columns missing" = all(req_cols %in% colnames(sobj@meta.data)))
    stopifnot("Graph reduction not found" = graph_reduction %in% Reductions(sobj))
    
    if (is.null(fixed_effects)) fixed_effects <- c(target_var, batch_var)
    glmm_family <- match.arg(glmm_family)
    use_glmm <- length(random_effects) > 0
    
    vcat("Starting Milo DA analysis...")
    vcat("Mode:", ifelse(use_glmm, "GLMM with random effects", "edgeR fixed effects"))
    
    # ---- 1. 타입 통일 ----
    vcat("Coercing metadata columns to character in source object...")
    sobj@meta.data[[patient_var]] <- as.character(sobj@meta.data[[patient_var]])
    sobj@meta.data[[target_var]]  <- as.character(sobj@meta.data[[target_var]])
    sobj@meta.data[[batch_var]]   <- as.character(sobj@meta.data[[batch_var]])
    sobj@meta.data[[cluster_var]] <- as.character(sobj@meta.data[[cluster_var]])
    
    # ---- 2. 샘플 정보 생성 ----
    vcat("Extracting sample-level information...")
    patient_info <- sobj@meta.data[, c(patient_var, target_var, batch_var)] %>%
        unique() %>%
        as.data.frame()
    
    vcat(sprintf("Found %d unique patients", nrow(patient_info)))
    vcat(sprintf("Target values in data: %s", 
                 paste(sort(unique(patient_info[[target_var]])), collapse=", ")))
    vcat(sprintf("Batch values in data: %s", 
                 paste(sort(unique(patient_info[[batch_var]])), collapse=", ")))
    
    # ---- Embeddings ----
    G_emb <- Embeddings(sobj, reduction = graph_reduction)
    if (is.null(d)) d <- ncol(G_emb)
    d <- min(d, ncol(G_emb))
    G_use <- G_emb[, seq_len(d), drop = FALSE]
    L_use <- NULL
    if (!is.null(layout_reduction)) {
        L_emb <- Embeddings(sobj, reduction = layout_reduction)
        L_use <- L_use
    }
    vcat(sprintf("Using %d dimensions from %s reduction", d, graph_reduction))
    
    # ---- Create SCE ----
    sce <- as.SingleCellExperiment(sobj, assay = DefaultAssay(sobj))
    reducedDim(sce, "GRAPH") <- as.matrix(G_use)
    if (!is.null(L_use)) reducedDim(sce, "LAYOUT") <- as.matrix(L_use)
    
    # ---- Build Milo ----
    vcat("Building kNN graph and neighborhoods...")
    milo <- Milo(sce)
    milo <- buildGraph(milo, k = k, d = ncol(G_use), reduced.dim = "GRAPH")
    milo <- makeNhoods(milo, prop = prop, k = k, refined = TRUE, reduced_dims = "GRAPH")
    
    # ----------------------------------------------------
    # !! [핵심 수정] 누락된 이웃 그래프 생성 단계 추가 !!
    # ----------------------------------------------------
    vcat("Building neighborhood graph...")
    milo <- buildNhoodGraph(milo)
    # ----------------------------------------------------

    vcat(sprintf("Created %d neighborhoods", ncol(nhoods(milo))))
    
    # [DEBUG] 이웃 그래프가 제대로 생성되었는지 확인
    nh_graph_dim <- dim(nhoodGraph(milo))
    vcat(sprintf("[DEBUG] Neighborhood graph dimensions: %d x %d", nh_graph_dim[1], nh_graph_dim[2]))
    if (nh_graph_dim[1] != ncol(nhoods(milo))) {
        stop("Neighborhood graph dimension mismatch!")
    }
    # ---- Get metadata ----
    milo_md <- as.data.frame(colData(milo))
    
    # ---- 3. Counts 계산 (수정된 .fast_count_cells 사용) ----
    vcat("Computing neighborhood-by-sample counts (using fixed .fast_count_cells)...")
    counts <- .fast_count_cells(milo, samples = milo_md[[patient_var]])
    
    H <- nrow(counts)
    S <- ncol(counts)
    vcat(sprintf("Count matrix: %d neighborhoods x %d samples", H, S))
    
    # ---- 4. Design Matrix 생성 (매칭 오류 해결됨) ----
    sample_ids <- colnames(counts)
    sample_idx <- match(sample_ids, patient_info[[patient_var]])
    
    if (any(is.na(sample_idx))) {
        missing <- sample_ids[is.na(sample_idx)]
        stop(sprintf("Could not find patient info for samples: %s", 
                     paste(missing, collapse=", ")))
    }
    
    sample_design <- patient_info[sample_idx, ]
    rownames(sample_design) <- sample_ids
    
    sample_design[[target_var]] <- factor(sample_design[[target_var]])
    sample_design[[batch_var]] <- factor(sample_design[[batch_var]])
    
    vcat(sprintf("\nSample design summary:"))
    vcat(sprintf("  Target (%s) levels: %s", target_var,
                 paste(levels(sample_design[[target_var]]), collapse=", ")))
    vcat(sprintf("  Batch (%s) levels: %s", batch_var,
                 paste(levels(sample_design[[batch_var]]), collapse=", ")))
    
    design_table <- table(sample_design[[target_var]], sample_design[[batch_var]])
    vcat("\nDesign structure (rows=target, cols=batch):")
    print(design_table)
    
    # ---- Filter neighborhoods ----
    nh_sizes <- rowSums(counts)
    keep_nh <- nh_sizes >= min_cells
    if (sum(!keep_nh) > 0) {
        vcat(sprintf("Filtering %d neighborhoods with < %d cells", sum(!keep_nh), min_cells))
    }
    
    # ---- 5. Majority cluster ----
    nh_mat <- miloR::nhoods(milo)
    cell_clusters <- milo_md[[cluster_var]]
    majority_cluster <- vapply(seq_len(ncol(nh_mat)), function(i) {
        idx <- as.logical(nh_mat[, i])
        labs <- cell_clusters[idx]
        if (length(labs) == 0) return(NA_character_)
        tbl <- table(labs)
        names(tbl)[which.max(tbl)]
    }, character(1))
    
    # ---- EdgeR analysis ----
    vcat("\nRunning edgeR analysis...")
    
    if (nlevels(sample_design[[target_var]]) < 2) {
        stop(sprintf("Target variable '%s' must have at least 2 levels, found: %s",
                     target_var, 
                     paste(levels(sample_design[[target_var]]), collapse=", ")))
    }
    
    formula_str <- paste("~ 0 +", target_var)
    vcat(sprintf("Using formula: %s", formula_str))
    
    mm <- model.matrix(as.formula(formula_str), data = sample_design)
    vcat(sprintf("Design matrix: %d x %d", nrow(mm), ncol(mm)))
    vcat(sprintf("Column names: %s", paste(colnames(mm), collapse=", ")))
    
    dge <- DGEList(counts = counts[keep_nh, , drop = FALSE])
    dge <- calcNormFactors(dge, method = "TMM")
    
    vcat("Estimating dispersions...")
    dge <- estimateDisp(dge, design = mm, robust = TRUE)
    vcat(sprintf("Common dispersion: %.4f", dge$common.dispersion))
    
    vcat("Fitting GLM...")
    fit <- glmQLFit(dge, design = mm, robust = TRUE)
    
    if (ncol(mm) == 2) {
        contrast_name <- paste(colnames(mm)[2], "vs", colnames(mm)[1])
        vcat(sprintf("Testing: %s", contrast_name))
        contrast <- numeric(ncol(mm)); contrast[1] <- -1; contrast[2] <- 1
        qlf <- glmQLFTest(fit, contrast = contrast)
    } else {
        vcat("Testing all groups")
        qlf <- glmQLFTest(fit, coef = 2:ncol(mm))
    }
    
    # 1. 필터링된 이웃(M개)에 대한 결과 테이블 생성
    nh_idx <- which(keep_nh)
    da.res_sub <- tibble(
      nhood = nh_idx,
      logFC = qlf$table$logFC,
      logCPM = qlf$table$logCPM,
      F = qlf$table$F,
      pval = qlf$table$PValue
    )
    
    vcat(sprintf("GLM results produced %d rows", nrow(da.res_sub)))

    # 2. (선택적) GLM 결과의 NA/Inf 값 필터링 (필요시)
    #    이전 디버깅에서 모든 행이 통과했으므로, 
    #    da.res_sub를 직접 필터링할 수 있습니다.
    da.res_sub <- da.res_sub %>%
       filter(
           is.finite(logFC) &
           is.finite(logCPM) &
           is.finite(F) &
           !is.na(pval)
       )
    
    if(nrow(da.res_sub) == 0) {
       stop("No valid neighborhoods remained after filtering NA/Inf/NaN values.")
    }
    vcat(sprintf("Filtered GLM results (non-finite removed): %d rows remain", nrow(da.res_sub)))


    # 3. [핵심 수정] 모든 이웃(H개)을 포함하는 전체 테이블 생성
    #    H = ncol(nhoods(milo))
    H_total_nhoods <- ncol(nhoods(milo))
    da.res_full <- tibble(
        nhood = seq_len(H_total_nhoods)
    )
    
    # 4. GLM 결과를 전체 테이블에 left_join
    #    필터링된(keep_nh=FALSE) 이웃과 GLM에서 NA가 된 이웃은
    #    logFC, pval 등이 모두 NA가 됩니다.
    da.res <- left_join(da.res_full, da.res_sub, by = "nhood")
    
vcat(sprintf("Re-expanded DA results to %d total neighborhoods", nrow(da.res)))

    # -------------------------------------------------------------------
    # !! [체계적인 디버그 및 수정 블록] !!
    # -------------------------------------------------------------------
    
    # [DEBUG] 1. 현재 'nhood' 컬럼의 상태 확인
    vcat(sprintf("[DEBUG] Class of 'nhood' column before fix: %s", class(da.res$nhood)))
    if (!is.numeric(da.res$nhood)) {
        vcat("[DEBUG] !! WARNING !! 'nhood' is NOT numeric. This is the likely cause of the error.")
    }
    
    # [FIX] 'nhood' 컬럼을 강제로 integer로 변환하여 인덱싱 오류 방지
    # 'NAs introduced by coercion' 및 'indices out of range' 오류를 해결합니다.
    da.res$nhood <- as.integer(da.res$nhood)
    
    # [DEBUG] 2. 수정한 'nhood' 컬럼의 상태 재확인
    vcat(sprintf("[DEBUG] Class of 'nhood' column after fix: %s", class(da.res$nhood)))
    
    # [DEBUG] 3. 인덱스 유효성 최종 검증
    n_milo_nhoods <- ncol(nhoods(milo))
    vcat(sprintf("[DEBUG] Max index in da.res: %d", max(da.res$nhood, na.rm = TRUE)))
    vcat(sprintf("[DEBUG] Total neighborhoods in milo: %d", n_milo_nhoods))
    
    if (any(is.na(da.res$nhood))) {
         vcat("[DEBUG] !! CRITICAL ERROR !! 'nhood' column contains NAs after coercion.")
         stop("Failed to create valid numeric indices for 'nhood'.")
    }
    if (max(da.res$nhood) > n_milo_nhoods) {
         vcat("[DEBUG] !! CRITICAL ERROR !! Max index > total neighborhoods.")
         stop("Index mismatch between da.res and milo object.")
    }
    vcat("[DEBUG] Index validation passed.")
    # -------------------------------------------------------------------

    # ---- 6. Spatial FDR (수정된 da.res 사용) ----
    vcat("Computing spatial FDR...")
    
    da.res <- graphSpatialFDR(milo, da.res, 
                              pvalues = "pval", 
                              indices = "nhood",  # 인덱스 컬럼 지정
                              k = k,
                              reduced.dim = "GRAPH")
    
    da.res$major_cluster <- majority_cluster[da.res$nhood]
    da.res <- da.res %>%
        mutate(is_sig = spatialFDR < fdr_thresh & abs(logFC) > logfc_thresh) %>%
        arrange(spatialFDR, desc(abs(logFC)))
    
    cluster_summary <- if (sum(da.res$is_sig, na.rm = TRUE) > 0) {
        da.res %>%
            filter(is_sig) %>%
            group_by(major_cluster) %>%
            summarise(
                n_sig_nhoods = n(),
                mean_logFC = mean(logFC),
                median_logFC = median(logFC),
                min_spatialFDR = min(spatialFDR),
                .groups = "drop"
            ) %>%
            arrange(desc(n_sig_nhoods))
    } else {
        tibble()
    }
    
    vcat(sprintf("\nFound %d significant neighborhoods", sum(da.res$is_sig, na.rm = TRUE)))
    if (nrow(cluster_summary) > 0) {
        vcat(sprintf("Affecting %d cell types", nrow(cluster_summary)))
    }
    
    return(list(
        mode = "edgeR_fixed_effects",
        formula = formula_str,
        milo = milo,
        da_nhood_with_cluster = da.res,
        da_by_cluster = cluster_summary,
        design_matrix = mm,
        sample_design = sample_design
    ))
}




#' MiloR DA 분석 전체 워크플로우 래퍼
#'
#' @param sobj Seurat 객체
#' @param output_dir 중간 RDS 파일을 저장할 폴더
#' @param prefix 저장할 파일 이름의 접두사 (예: "ProjectA")
#' @param patient_var 환자/샘플 ID 컬럼명
#' @param cluster_var 세포 유형(클러스터) 컬럼명
#' @param target_var DA 분석 대상 변수 (예: "Condition")
#' @param batch_var 배치 변수 (모델에 포함)
#' @param graph_reduction Nhood 생성용 리덕션
#' @param layout_reduction 시각화용 리덕션 (예: "umap")
#' @param k k-NN 파라미터
#' @param d 리덕션 차원
#' @param prop Nhood 생성 파라미터
#' @param force_rerun_step (선택) 1, 2, 3, 4. 특정 단계부터 강제로 다시 실행합니다.
#'
#' @return DA 결과 플롯 리스트
#'
#' @export 
MILO <- function(
    sobj,
    output_dir,
    prefix,
    patient_var,
    cluster_var,
    target_var,
    batch_var,
    graph_reduction,
    layout_reduction = "umap",
    k = 30,
    d = 30,
    prop = 0.1,
    force_rerun_step = NULL
) {
    
    # --- 0. 라이브러리 로드 ---
    suppressPackageStartupMessages({
        library(Seurat)
        library(SingleCellExperiment)
        library(miloR)
        library(Matrix)
        library(dplyr)
        library(tibble)
        library(patchwork)
        library(ggplot2)
    })
    set.seed(1)
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    # --- 파일 경로 정의 ---
    file_step1 <- file.path(output_dir, paste0(prefix, "_01_nhoods.qs"))
    file_step2 <- file.path(output_dir, paste0(prefix, "_02_distances.qs"))
    file_step3_milo <- file.path(output_dir, paste0(prefix, "_03_tested_milo.qs"))
    file_step3_res  <- file.path(output_dir, paste0(prefix, "_03_tested_results.qs"))

    # --- 1단계: Nhood 생성 ---
    if (!is.null(force_rerun_step) && force_rerun_step <= 1 || !file.exists(file_step1)) {
        message("1단계: Milo 객체 생성 시작...")
        sobj@meta.data[[patient_var]] <- as.character(sobj@meta.data[[patient_var]])
        sobj@meta.data[[target_var]]  <- as.character(sobj@meta.data[[target_var]])
        sobj@meta.data[[batch_var]]   <- as.character(sobj@meta.data[[batch_var]])
        sobj@meta.data[[cluster_var]] <- as.character(sobj@meta.data[[cluster_var]])

        G_emb <- Embeddings(sobj, reduction = graph_reduction)
        G_use <- G_emb[, seq_len(d), drop = FALSE]
        sce <- as.SingleCellExperiment(sobj, assay = DefaultAssay(sobj))
        reducedDim(sce, "GRAPH") <- as.matrix(G_use)
        
        milo <- Milo(sce)
        milo <- buildGraph(milo, k = k, d = d, reduced.dim = "GRAPH")
        milo <- makeNhoods(milo, prop = prop, k = k, refined = TRUE, reduced_dims = "GRAPH")
        
        qs::qsave(milo, file_step1)
        message("1단계 완료 및 저장.")
    } else {
        message("1단계 파일 로드 중...")
        milo <- qs::qread(file_step1)
    }

    # --- 2단계: Nhood 거리 계산 (병목) ---
    if (!is.null(force_rerun_step) && force_rerun_step <= 2 || !file.exists(file_step2)) {
        message("2단계: Nhood 거리 계산 시작...")
        milo <- calcNhoodDistance(milo, d = d)
        qs::qsave(milo, file_step2)
        message("2단계 완료 및 저장.")
    } else {
        message("2단계 파일 로드 중...")
        milo <- qs::qread(file_step2)
    }

    # --- 3단계: DA 테스트 ---
    if (!is.null(force_rerun_step) && force_rerun_step <= 3 || 
        !file.exists(file_step3_milo) || !file.exists(file_step3_res)) {
        
        message("3단계: 카운트 집계 및 DA 테스트 시작...")
        
        # 3.1 카운트
        milo <- countCells(milo, 
                           samples = patient_var, 
                           meta.data = as.data.frame(colData(milo)))
        
        # 3.2 디자인
        patient_info <- as.data.frame(colData(milo))[, c(patient_var, target_var, batch_var)] %>%
            distinct()
        rownames(patient_info) <- patient_info[[patient_var]]
        sample_design <- patient_info[colnames(nhoodCounts(milo)), , drop = FALSE]
        sample_design[[target_var]] <- factor(sample_design[[target_var]])
        sample_design[[batch_var]]  <- factor(sample_design[[batch_var]])
        
        # 3.3 DA 테스트
        da_results <- testNhoods(milo, 
                                 design = as.formula(paste("~", batch_var, "+", target_var)),
                                 design.df = sample_design)
        
        qs::qsave(milo, file_step3_milo)
        qs::qsave(da_results, file_step3_res)
        message("3단계 완료 및 저장.")
    } else {
        message("3단계 파일 로드 중...")
        milo       <- qs::qread(file_step3_milo)
        da_results <- qs::qread(file_step3_res)
    }

    # --- 4단계: 시각화 ---
    message("4단계: 시각화...")
    
    if (is.null(nhoodGraph(milo))) {
        milo <- buildNhoodGraph(milo)
    }
    if (!("UMAP" %in% reducedDimNames(milo)) && (layout_reduction %in% Reductions(sobj))) {
        reducedDim(milo, "UMAP") <- Embeddings(sobj, reduction = layout_reduction)
    }

    nh_mat <- miloR::nhoods(milo)
    cell_clusters <- colData(milo)[[cluster_var]]
    majority_cluster <- vapply(seq_len(ncol(nh_mat)), function(i) {
        idx <- as.logical(nh_mat[, i])
        labs <- cell_clusters[idx]
        if (length(labs) == 0) return(NA_character_)
        tbl <- table(labs)
        names(tbl)[which.max(tbl)]
    }, character(1))
    
    da_results$major_cluster <- majority_cluster[da_results$nhood]

    p_da_graph <- plotNhoodGraphDA(milo, da_results, alpha = 0.1) + 
        ggtitle(paste(prefix, "DA Results (logFC)"))
        
    p_beeswarm <- plotDAbeeswarm(da_results, group.by = "major_cluster") +
        ggtitle(paste(prefix, "DA results by Cell Type"))

    message("워크플로우 완료.")
    
    # 최종 결과 반환
    return(list(
        milo_obj = milo,
        da_results = da_results,
        plot_da_graph = p_da_graph,
        plot_beeswarm = p_beeswarm
    ))
}
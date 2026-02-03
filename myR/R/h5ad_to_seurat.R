#' Load h5ad file and convert to Seurat object, then save as .qs file
#'
#' @param h5ad_path Path to the h5ad file
#' @param condaenv Conda environment with anndata installed (default: "/home/jaecheon/miniconda3/envs/scenvi")
#' @param save_path Directory to save the .qs file (default: current directory)
#' @param save_name Name for the saved file (without extension). If NULL, uses basename of h5ad_path
#' @param use_layer If h5ad has layers, which layer to use for counts (e.g., "counts", "RNA_counts"). If NULL, uses X matrix
#' @param project_name Project name for Seurat object (default: "h5ad_import")
#' @import reticulate Seurat qs Matrix
#' @export
load_h5ad_to_seurat_qs <- function(h5ad_path,
                                   condaenv = "/home/jaecheon/miniconda3/envs/scvi-tools",
                                   save_path = NULL,
                                   save_name = NULL,
                                   use_layer = NULL,
                                   project_name = "h5ad_import") {
  # Python과 anndata 설정
  reticulate::use_condaenv(condaenv)
  anndata <- reticulate::import("anndata")

  # Matrix 패키지 로드 (S4 메서드 디스패치 보장)
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Matrix package is required")
  }
  library(Matrix)
  # Seurat 패키지 로드
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package is required")
  }
  library(Seurat)

  # h5ad 파일 로드
  message("Loading h5ad file: ", h5ad_path)
  adata <- anndata$read_h5ad(h5ad_path)

  # Python builtins 가져오기 (list 변환용)
  builtins <- reticulate::import("builtins")

  # Helper function to safely extract matrix
  # handles cases where reticulate already converted it to S4, or if it's still a python object
  extract_matrix_to_r <- function(mat) {
    # If it is a python object, convert it
    if (inherits(mat, "python.builtin.object")) {
      if (reticulate::py_is_null_xptr(mat)) {
        return(NULL)
      }

      # Try to convert to CSC if available and it's a sparse matrix
      if (reticulate::py_has_attr(mat, "tocsc")) {
        mat <- mat$tocsc()
      }
      return(reticulate::py_to_r(mat))
    }
    # If it's already an R object (e.g. S4 matrix), return it as is
    return(mat)
  }

  # Counts 매트릭스 추출
  # AnnData는 cells x genes 형식이므로 transpose 필요
  # 먼저 행/열 이름 가져오기 (fallback 포함)
  var_names <- tryCatch(
    {
      reticulate::py_to_r(adata$var_names$tolist())
    },
    error = function(e) {
      tryCatch(
        {
          reticulate::py_to_r(adata$var$index$tolist())
        },
        error = function(e2) {
          message("Warning: Could not retrieve var_names")
          NULL
        }
      )
    }
  )

  obs_names <- tryCatch(
    {
      reticulate::py_to_r(adata$obs_names$tolist())
    },
    error = function(e) {
      tryCatch(
        {
          reticulate::py_to_r(adata$obs$index$tolist())
        },
        error = function(e2) {
          message("Warning: Could not retrieve obs_names")
          NULL
        }
      )
    }
  )

  layer_keys <- NULL
  if (!is.null(adata$layers)) {
    # Python의 list()로 KeysView를 변환 후 R 문자벡터로 가져옴
    layer_keys <- tryCatch(
      reticulate::py_to_r(builtins$list(adata$layers$keys())),
      error = function(e) {
        message("Error getting layer keys: ", e$message)
        NULL
      }
    )
  }

  # 매트릭스 추출 및 transpose
  if (!is.null(use_layer) && !is.null(layer_keys) && use_layer %in% layer_keys) {
    message("Using layer: ", use_layer)
    counts_matrix_raw <- extract_matrix_to_r(adata$layers[[use_layer]])
  } else {
    # X 매트릭스 사용 (또는 layers에 "counts"가 있으면 사용)
    if (!is.null(layer_keys) && "counts" %in% layer_keys) {
      message("Using 'counts' layer from h5ad")
      counts_matrix_raw <- extract_matrix_to_r(adata$layers[["counts"]])
    } else {
      message("Using X matrix from h5ad")
      counts_matrix_raw <- extract_matrix_to_r(adata$X)
    }
  }

  # AnnData는 cells x genes 형식이므로 transpose 필요 (Seurat는 genes x cells)
  message("Matrix dimensions before transpose: ", nrow(counts_matrix_raw), " x ", ncol(counts_matrix_raw))
  message("Class of counts_matrix_raw: ", paste(class(counts_matrix_raw), collapse = ", "))

  # sparse matrix 여부에 따라 적절한 transpose 사용
  if (inherits(counts_matrix_raw, "sparseMatrix") || inherits(counts_matrix_raw, "dgCMatrix") || inherits(counts_matrix_raw, "dgRMatrix")) {
    counts_matrix <- Matrix::t(counts_matrix_raw)
  } else {
    counts_matrix <- t(counts_matrix_raw)
  }

  # Seurat는 CsparseMatrix(dgCMatrix)를 선호하므로 dgRMatrix인 경우 변환
  if (inherits(counts_matrix, "dgRMatrix")) {
    counts_matrix <- as(counts_matrix, "CsparseMatrix")
  }
  message("Matrix dimensions after transpose: ", nrow(counts_matrix), " x ", ncol(counts_matrix))

  # 행/열 이름 설정
  if (!is.null(var_names) && length(var_names) == nrow(counts_matrix)) {
    rownames(counts_matrix) <- var_names
  } else if (nrow(counts_matrix) > 0) {
    rownames(counts_matrix) <- paste0("Gene_", 1:nrow(counts_matrix))
    message("Warning: Using default gene names. Loaded names length: ", length(var_names), ", Matrix rows: ", nrow(counts_matrix))
  }

  if (!is.null(obs_names) && length(obs_names) == ncol(counts_matrix)) {
    colnames(counts_matrix) <- obs_names
  } else if (ncol(counts_matrix) > 0) {
    colnames(counts_matrix) <- paste0("Cell_", 1:ncol(counts_matrix))
    message("Warning: Using default cell names. Loaded names length: ", length(obs_names), ", Matrix cols: ", ncol(counts_matrix))
  }

  # 메타데이터 추출 - 모든 cell-level metadata 가져오기
  if (!is.null(adata$obs)) {
    metadata <- reticulate::py_to_r(adata$obs)
    # row.names가 cell names와 일치하는지 확인
    if (!is.null(obs_names) && nrow(metadata) == length(obs_names)) {
      rownames(metadata) <- obs_names
    }
    message("Loaded ", ncol(metadata), " metadata columns")
  } else {
    metadata <- data.frame(row.names = colnames(counts_matrix))
    message("Warning: No metadata found in h5ad file")
  }

  # Seurat 객체 생성
  message("Creating Seurat object...")
  seurat_obj <- CreateSeuratObject(
    counts = counts_matrix,
    project = project_name,
    meta.data = metadata
  )

  # Layers에서 다른 데이터 추가 (예: normalized data)
  if (!is.null(layer_keys)) {
    message("Available layers: ", paste(layer_keys, collapse = ", "))

    # "data" layer가 있으면 normalized data로 추가
    if ("data" %in% layer_keys) {
      data_matrix <- extract_matrix_to_r(adata$layers[["data"]])

      # sparse matrix 여부에 따라 적절한 transpose 사용
      if (inherits(data_matrix, "sparseMatrix") || inherits(data_matrix, "dgCMatrix") || inherits(data_matrix, "dgRMatrix")) {
        data_matrix <- Matrix::t(data_matrix)
      } else {
        data_matrix <- t(data_matrix)
      }

      rownames(data_matrix) <- rownames(counts_matrix)
      colnames(data_matrix) <- colnames(counts_matrix)

      # Seurat는 CsparseMatrix(dgCMatrix)를 선호
      if (inherits(data_matrix, "dgRMatrix")) {
        data_matrix <- as(data_matrix, "CsparseMatrix")
      }
      seurat_obj[["RNA"]]@layers[["data"]] <- data_matrix
      message("Added 'data' layer to Seurat object")
    }
  }

  # 차원 축소 결과 추가 (obsm) - 모든 obsm을 reductions로 변환
  if (!is.null(adata$obsm)) {
    # Python의 list()로 KeysView를 변환 후 R 문자벡터로 가져옴
    obsm_keys <- tryCatch(
      reticulate::py_to_r(builtins$list(adata$obsm$keys())),
      error = function(e) {
        message("Error getting obsm keys: ", e$message)
        NULL
      }
    )

    if (!is.null(obsm_keys)) {
      message("Available obsm keys: ", paste(obsm_keys, collapse = ", "))

      # obsm 이름을 Seurat reduction 이름으로 매핑
      obsm_to_reduction <- list(
        "X_pca" = list(name = "pca", key = "PC_"),
        "X_umap" = list(name = "umap", key = "UMAP_"),
        "X_tsne" = list(name = "tsne", key = "tSNE_"),
        "X_scVI" = list(name = "integrated.scvi", key = "SCVI_"),
        "spatial" = list(name = "spatial", key = "SPATIAL_")
      )

      for (obsm_key in obsm_keys) {
        # 매핑이 있으면 사용, 없으면 obsm_key에서 "X_" 제거하여 사용
        if (obsm_key %in% names(obsm_to_reduction)) {
          reduction_info <- obsm_to_reduction[[obsm_key]]
          reduction_name <- reduction_info$name
          reduction_key <- reduction_info$key
        } else {
          # X_ prefix 제거
          reduction_name <- sub("^X_", "", obsm_key)
          reduction_key <- paste0(toupper(reduction_name), "_")
        }

        tryCatch(
          {
            embeddings <- reticulate::py_to_r(adata$obsm[[obsm_key]])

            # embeddings가 matrix인지 확인
            if (!is.matrix(embeddings)) {
              embeddings <- as.matrix(embeddings)
            }

            # 차원 확인
            if (nrow(embeddings) != ncol(seurat_obj)) {
              message("Warning: ", obsm_key, " embeddings dimension mismatch. Skipping.")
              next
            }

            rownames(embeddings) <- colnames(seurat_obj)
            colnames(embeddings) <- paste0(reduction_key, 1:ncol(embeddings))

            seurat_obj[[reduction_name]] <- CreateDimReducObject(
              embeddings = embeddings,
              key = reduction_key,
              assay = DefaultAssay(seurat_obj)
            )
            message("Added reduction: ", reduction_name, " (from ", obsm_key, ")")
          },
          error = function(e) {
            message("Warning: Failed to add reduction from ", obsm_key, ": ", e$message)
          }
        )
      }
    }
  }

  # 저장 경로 및 파일명 설정
  if (is.null(save_path)) {
    save_path <- dirname(h5ad_path)
  }
  if (is.null(save_name)) {
    save_name <- tools::file_path_sans_ext(basename(h5ad_path))
  }

  # .qs 파일로 저장
  qs_file <- file.path(save_path, paste0(save_name, ".qs"))
  message("Saving Seurat object to: ", qs_file)
  qs::qsave(seurat_obj, qs_file)

  message("Successfully saved Seurat object to: ", qs_file)
  message("Seurat object dimensions: ", nrow(seurat_obj), " genes x ", ncol(seurat_obj), " cells")

  return(seurat_obj)
}

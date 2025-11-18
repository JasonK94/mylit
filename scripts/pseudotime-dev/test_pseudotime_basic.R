# Pseudotime Analysis 기본 테스트 스크립트
#
# 이 스크립트는 다운샘플 데이터를 사용하여 기본적인 pseudotime 분석을 테스트합니다.

# --- 0. 환경 설정 ---
message("=== Pseudotime Analysis Basic Test ===")

# 라이브러리 로드
if (!requireNamespace("qs", quietly = TRUE)) {
  stop("qs package is required. Install with: install.packages('qs')")
}

# 패키지 로드 (작업 디렉토리에서 실행 시 자동 로드됨)
# devtools::load_all("myR")

# --- 1. 데이터 로드 (전처리된 데이터 사용) ---
message("\n--- Step 1: Loading preprocessed data ---")
data_file <- "/data/user3/sobj/IS_scvi_251107_ds2500_preprocessed.qs"

if (!file.exists(data_file)) {
  # 전처리된 파일이 없으면 원본 파일 사용 (경고)
  data_file <- "/data/user3/sobj/IS_scvi_251107_ds2500.qs"
  warning("Preprocessed file not found. Using original file. ",
          "Please run preprocess_pseudotime_data() first!", call. = FALSE)
}

message("Loading data from: ", data_file)
sobj <- qs::qread(data_file)

if (!inherits(sobj, "Seurat")) {
  stop("Loaded object is not a Seurat object.", call. = FALSE)
}

message("Data loaded successfully:")
message("  Cells: ", ncol(sobj))
message("  Genes: ", nrow(sobj))
message("  Metadata columns: ", paste(head(colnames(sobj@meta.data), 10), collapse = ", "))

# --- 2. 주요 Feature 확인 ---
message("\n--- Step 2: Checking key features ---")
key_features <- c("DDIT4", "UTY", "S100B", "XIST", "HLA-B", "CCL4", "HLA-C", "TXNIP")
available_features <- intersect(key_features, rownames(sobj))
missing_features <- setdiff(key_features, rownames(sobj))

message("Available key features: ", length(available_features), "/", length(key_features))
if (length(available_features) > 0) {
  message("  Found: ", paste(available_features, collapse = ", "))
}
if (length(missing_features) > 0) {
  warning("Missing features: ", paste(missing_features, collapse = ", "), call. = FALSE)
}

# --- 3. Metadata 확인 ---
message("\n--- Step 3: Checking metadata ---")
important_metadata <- c("nih_change", "sex", "g3", "arrival_gcs_score", 
                       "age", "ini_nih", "nih_end1", 
                       "icu_adm_dt", "ia_start", "ia_end")

available_metadata <- intersect(important_metadata, colnames(sobj@meta.data))
missing_metadata <- setdiff(important_metadata, colnames(sobj@meta.data))

message("Available important metadata: ", length(available_metadata), "/", length(important_metadata))
if (length(available_metadata) > 0) {
  for (md in available_metadata) {
    md_class <- class(sobj@meta.data[[md]])
    md_na <- sum(is.na(sobj@meta.data[[md]]))
    md_total <- length(sobj@meta.data[[md]])
    message("  ", md, ": ", md_class, " (NA: ", md_na, "/", md_total, ")")
    
    # sex와 g3는 factor인지 확인
    if (md == "sex") {
      if (is.factor(sobj@meta.data[[md]])) {
        message("    Levels: ", paste(levels(sobj@meta.data[[md]]), collapse = ", "))
      } else {
        warning("    sex is not a factor! Please run preprocess_pseudotime_data().", call. = FALSE)
      }
    }
    if (md == "g3") {
      if (is.factor(sobj@meta.data[[md]])) {
        message("    Levels: ", paste(levels(sobj@meta.data[[md]]), collapse = ", "))
      } else {
        warning("    g3 is not a factor! Please run preprocess_pseudotime_data().", call. = FALSE)
      }
    }
  }
}
if (length(missing_metadata) > 0) {
  warning("Missing metadata: ", paste(missing_metadata, collapse = ", "), call. = FALSE)
}

# --- 4. Slingshot 테스트 ---
message("\n--- Step 4: Testing Slingshot trajectory inference ---")

# 클러스터 컬럼 확인
cluster_cols <- grep("cluster|seurat_clusters", colnames(sobj@meta.data), 
                     value = TRUE, ignore.case = TRUE)
if (length(cluster_cols) == 0) {
  warning("No cluster column found. Skipping Slingshot test.", call. = FALSE)
} else {
  cluster_col <- cluster_cols[1]
  message("Using cluster column: ", cluster_col)
  
  # Reduction 확인
  available_reductions <- names(sobj@reductions)
  if (length(available_reductions) == 0) {
    warning("No dimensionality reduction found. Skipping Slingshot test.", call. = FALSE)
  } else {
    reduction_name <- if ("umap" %in% tolower(available_reductions)) {
      available_reductions[tolower(available_reductions) == "umap"][1]
    } else {
      available_reductions[1]
    }
    message("Using reduction: ", reduction_name)
    
    # Slingshot 실행
    tryCatch({
      sce_slingshot <- run_slingshot_from_seurat(
        seurat_obj = sobj,
        cluster_col = cluster_col,
        reduced_dim_name = reduction_name,
        counts_assay_name = "RNA"
      )
      
      if (!is.null(sce_slingshot)) {
        message("Slingshot analysis successful!")
        message("  Cells in SCE: ", ncol(sce_slingshot))
        
        # Pseudotime 확인
        if ("slingshot" %in% SingleCellExperiment::reducedDimNames(sce_slingshot)) {
          pst_matrix <- SingleCellExperiment::reducedDim(sce_slingshot, "slingshot")
          message("  Pseudotime matrix: ", nrow(pst_matrix), " cells x ", ncol(pst_matrix), " lineages")
        }
        
        # 결과 저장
        output_file_sling <- "/data/user3/sobj/IS_scvi_251107_ds2500_slingshot.qs"
        qs::qsave(sce_slingshot, output_file_sling)
        message("  Saved to: ", output_file_sling)
      }
    }, error = function(e) {
      warning("Slingshot test failed: ", e$message, call. = FALSE)
    })
  }
}

# --- 5. Monocle3 테스트 ---
message("\n--- Step 5: Testing Monocle3 trajectory inference ---")

tryCatch({
  cds_monocle3 <- run_monocle3_from_seurat(
    seurat_obj = sobj,
    counts_assay_name = "RNA",
    reduction_method = "UMAP",
    verbose = TRUE
  )
  
  if (!is.null(cds_monocle3)) {
    message("Monocle3 analysis successful!")
    message("  Cells in CDS: ", ncol(cds_monocle3))
    
    # Pseudotime 확인
    pt <- tryCatch({
      monocle3::pseudotime(cds_monocle3)
    }, error = function(e) {
      NULL
    })
    
    if (!is.null(pt)) {
      message("  Pseudotime extracted: ", sum(!is.na(pt)), " cells with pseudotime")
    }
    
    # 결과 저장
    output_file_m3 <- "/data/user3/sobj/IS_scvi_251107_ds2500_monocle3.qs"
    qs::qsave(cds_monocle3, output_file_m3)
    message("  Saved to: ", output_file_m3)
  }
}, error = function(e) {
  warning("Monocle3 test failed: ", e$message, call. = FALSE)
})

message("\n=== Test completed ===")


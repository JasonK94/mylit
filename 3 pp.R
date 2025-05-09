# pp.R
# Seurat 데이터 전처리 스크립트
#
# 기능:
# 1. 10x Genomics 데이터 로드
# 2. Demultiplexing 방식 (SNP 또는 HTO)에 따른 처리 분기
#    - SNP: 제공된 demultiplexing 정보 로드, 메타데이터 추가, 필요시 doublet 제거
#    - HTO: HTO 정보 기반으로 메타데이터 추가 (이미 demultiplexing 되었다고 가정)
# 3. Normalization 방식 (SCT 또는 LogNormalize)에 따른 처리 분기
# 4. 기본적인 Seurat 전처리 (Variable Features, PCA, Clustering)
# 5. 처리된 Seurat 객체 및 raw count matrix 저장
#
# 실행 예시:
# Rscript pp.R --data "/path/to/filtered_feature_bc_matrix" \
#              --output "/path/to/output_prefix" \
#              --multi "SNP" \
#              --nmz "SCT" \
#              --GEM "GEM1" \
#              --demulti "/path/to/demux_posterior.csv" \
#              --doublet_removal_by_demulti TRUE
#
# Rscript pp.R --data "/path/to/hto_count_matrix" \
#              --output "/path/to/output_prefix" \
#              --multi "HTO" \
#              --nmz "SCT" \
#              --GEM "GEM9"

# 라이브러리 로드
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(optparse)))

# 커맨드 라인 인자 정의
option_list <- list(
  make_option(c("-d", "--data"), type="character", default=NULL, help="입력 데이터 경로 (filtered_feature_bc_matrix 또는 HTO count matrix 디렉토리)", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="seurat_preprocessed", help="출력 파일 접두사 (RDS 및 raw counts 저장용)", metavar="character"),
  make_option(c("-m", "--multi"), type="character", default="SNP", help="Demultiplexing 방식 ('SNP' 또는 'HTO')", metavar="character"),
  make_option(c("-n", "--nmz"), type="character", default="SCT", help="Normalization 방식 ('SCT' 또는 'LogNormalize')", metavar="character"),
  make_option(c("-g", "--GEM"), type="character", default="GEM_default", help="GEM 또는 샘플 그룹 정보", metavar="character"),
  make_option(c("--demulti"), type="character", default=NULL, help="SNP-based demultiplexing 결과 CSV 파일 경로", metavar="character"),
  make_option(c("--doublet_removal_by_demulti"), type="logical", default=FALSE, help="SNP-based demultiplexing으로 doublet 제거 여부", action = "store_true")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$data)){
  print_help(opt_parser)
  stop("입력 데이터 경로 (--data)를 지정해야 합니다.", call.=FALSE)
}

# 사용자 정의 함수 (제공된 스크립트 기반, 필요시 실제 구현 필요)
# 예시: get_barcode_mapping, is_doublet 등
# 여기서는 SNP demultiplexing 시 필요한 함수라고 가정합니다.
# get_barcode_mapping <- function(demux_data_df) {
#   # 실제 구현: demux_data_df를 처리하여 Barcode, Best_Sample 등을 포함하는 데이터프레임 반환
#   # 예시:
#   # demux_data_df <- demux_data_df %>%
#   #   mutate(Best_Sample = apply(demux_data_df[,-1], 1, function(x) colnames(demux_data_df[,-1])[which.max(x)])) %>%
#   #   rename(Barcode = BARCODE)
#   return(demux_data_df)
# }
# is_doublet <- function(best_sample_vector) {
#   # 실제 구현: Best_Sample 정보를 바탕으로 doublet 여부 판단
#   # 예시:
#   # return(grepl("DBL", best_sample_vector) | best_sample_vector == "doublet") # 예시 조건
#   return(rep(FALSE, length(best_sample_vector))) # 임시로 모두 singlet으로 처리
# }


# 1. 데이터 로드
message("1. 데이터 로드 중: ", opt$data)
if (opt$multi == "HTO" || file.exists(file.path(opt$data, "matrix.mtx.gz")) || file.exists(file.path(opt$data, "matrix.mtx"))) {
    counts <- Read10X(data.dir = opt$data)
    # HTO 데이터의 경우 'Gene Expression' assay를 사용할 수 있음
    if (is.list(counts) && "Gene Expression" %in% names(counts)) {
        raw_counts_for_soupx <- counts$`Gene Expression` # SoupX를 위해 raw counts 저장
        counts <- counts$`Gene Expression`
        message("HTO 데이터의 'Gene Expression' assay 사용.")
    } else if (is.list(counts) && !is.null(counts[[1]])) {
        # 간혹 CellPlex 데이터가 Assay 객체 리스트로 로드될 때
        raw_counts_for_soupx <- counts[[1]]
        counts <- counts[[1]]
        message("첫 번째 assay 사용.")
    } else {
        raw_counts_for_soupx <- counts # SoupX를 위해 raw counts 저장
    }
} else {
    stop("제공된 --data 경로가 유효한 10x 데이터 디렉토리인지 확인하세요.")
}

# Seurat 객체 생성 (cell ID에 GEM 정보 추가하여 고유성 보장)
s_obj <- CreateSeuratObject(counts = counts, project = opt$GEM, add.cell.ids = opt$GEM)
message("Seurat 객체 생성 완료. 초기 세포 수: ", ncol(s_obj))

# 원본 raw counts 저장 (SoupX 입력용)
# CreateSeuratObject는 필터링을 수행할 수 있으므로, Read10X 직후의 raw_counts_for_soupx 사용
output_raw_counts_path <- paste0(opt$output, "_raw_counts.rds")
saveRDS(raw_counts_for_soupx, file = output_raw_counts_path)
message("Raw counts 매트릭스 저장 완료: ", output_raw_counts_path)


# 2. Demultiplexing 및 메타데이터 추가
s_obj$GEM <- opt$GEM

if (opt$multi == "SNP") {
  message("SNP 기반 demultiplexing 처리 중...")
  if (is.null(opt$demulti)) {
    stop("SNP demultiplexing의 경우 --demulti 경로를 지정해야 합니다.", call.=FALSE)
  }
  demux_data <- read.csv(opt$demulti)
  
  # --- 사용자 제공 스크립트의 demultiplexing 로직 참조 ---
  # colnames(demux_data) 첫번째는 바코드, 나머지는 샘플 확률로 가정
  # 아래는 예시이며, 실제 `get_barcode_mapping` 및 `is_doublet` 함수 구현 필요
  # 또는 사용자의 `barcode_mapping_list[[i]]` 생성 로직 적용
  
  # 바코드 형식 일치 확인 (Seurat 객체는 "BARCODE-1_GEM번호" 형식, demux_data는 "BARCODE" 형식일 수 있음)
  # Seurat 객체의 바코드에서 GEM 접미사 및 '-1' 제거하여 매칭 시도
  original_barcodes <- rownames(s_obj@meta.data)
  s_obj@meta.data$original_barcode_for_demux <- gsub(paste0("-1_", opt$GEM), "", original_barcodes)
  s_obj@meta.data$original_barcode_for_demux <- gsub(paste0("_", opt$GEM), "", s_obj@meta.data$original_barcode_for_demux)


  # demux_data의 바코드 컬럼 이름을 'BARCODE'로 가정
  if (!"BARCODE" %in% colnames(demux_data)) {
      if ("Barcode" %in% colnames(demux_data)) {
          colnames(demux_data)[colnames(demux_data) == "Barcode"] <- "BARCODE"
      } else {
          warning("demux_data에 'BARCODE' 또는 'Barcode' 컬럼이 없습니다. 첫 번째 컬럼을 바코드로 사용합니다.")
          colnames(demux_data)[1] <- "BARCODE"
      }
  }
  
  # Best sample 및 doublet 정보 추가 (사용자 정의 함수 필요)
  # demux_data$Best_Sample <- apply(demux_data[, -which(colnames(demux_data)=="BARCODE")], 1, function(x) colnames(demux_data[, -which(colnames(demux_data)=="BARCODE")])[which.max(x)])
  # demux_data$droplet_demulti <- ifelse(is_doublet(demux_data$Best_Sample), "doublet_demulti", "singlet_demulti")
  
  # 임시: Best_Sample을 확률이 가장 높은 컬럼명으로 지정, doublet은 없다고 가정
  sample_cols <- setdiff(colnames(demux_data), "BARCODE")
  if(length(sample_cols) > 0) {
      best_sample_indices <- apply(demux_data[, sample_cols, drop = FALSE], 1, which.max)
      demux_data$Best_Sample <- sample_cols[best_sample_indices]
      demux_data$droplet_demulti <- "singlet_demulti" # 임시
  } else {
      warning("demux_data에 샘플 확률 컬럼이 없습니다.")
      demux_data$Best_Sample <- NA
      demux_data$droplet_demulti <- NA
  }

  # Seurat 객체 메타데이터와 demux_data 조인
  # AddMetaData는 rownames가 일치해야 하므로, demux_data의 바코드를 s_obj의 바코드와 일치시키거나,
  # s_obj 메타데이터에 demux_data 정보를 매핑하여 추가.
  # 여기서는 s_obj의 original_barcode_for_demux 와 demux_data의 BARCODE를 기준으로 merge
  current_meta <- s_obj@meta.data
  current_meta <- dplyr::left_join(current_meta, demux_data, by = c("original_barcode_for_demux" = "BARCODE"))
  rownames(current_meta) <- original_barcodes # join 후 rownames 복원 중요!
  s_obj@meta.data <- current_meta
  
  if (opt$doublet_removal_by_demulti && "droplet_demulti" %in% colnames(s_obj@meta.data)) {
    s_obj <- subset(s_obj, subset = droplet_demulti == "singlet_demulti")
    message("SNP demultiplexing 기반 doublet 제거 후 세포 수: ", ncol(s_obj))
  }

} else if (opt$multi == "HTO") {
  message("HTO 기반 처리 중... (별도 demultiplexing 정보 파일 사용 안 함)")
  # HTO의 경우, demultiplexing은 이미 수행되었다고 가정.
  # 필요한 경우 HTO.global 이나 HTO.Classification 등의 메타데이터가 이미 s_obj에 있을 수 있음 (Cell Ranger `multi` 파이프라인)
  # 또는 사용자가 별도 컬럼으로 HTO 정보를 추가해야 할 수 있음.
  # 여기서는 GEM 정보를 메타데이터에 추가하는 것 외에는 특별한 HTO demultiplexing 처리는 하지 않음.
  # 예: s_obj$hash.ID <- s_obj$HTO_classification (만약 있다면)
}

# 3. Normalization
message("3. Normalization (", opt$nmz, ") 및 기본 전처리 중...")
if (opt$nmz == "SCT") {
  s_obj <- SCTransform(s_obj, verbose = FALSE)
  DefaultAssay(s_obj) <- "SCT"
} else if (opt$nmz == "LogNormalize") {
  DefaultAssay(s_obj) <- "RNA" # LogNormalize는 RNA assay에 대해 수행
  s_obj <- NormalizeData(s_obj, verbose = FALSE)
  s_obj <- FindVariableFeatures(s_obj, verbose = FALSE)
  s_obj <- ScaleData(s_obj, verbose = FALSE) # LogNormalize의 경우 ScaleData 필요
} else {
  stop("지원하지 않는 Normalization 방식입니다. 'SCT' 또는 'LogNormalize'를 사용하세요.", call.=FALSE)
}

# 4. 기본 전처리 (PCA, Clustering)
s_obj <- FindVariableFeatures(s_obj, verbose = FALSE) # SCT 후에도 HVG를 다시 찾거나, SCT의 variable.features.SCT를 사용
s_obj <- RunPCA(s_obj, verbose = FALSE)
s_obj <- FindNeighbors(s_obj, dims = 1:30, verbose = FALSE) # dims는 데이터에 따라 조절
s_obj <- FindClusters(s_obj, verbose = FALSE)

message("전처리 완료. 최종 세포 수: ", ncol(s_obj))

# 5. 결과 저장
output_rds_path <- paste0(opt$output, ".rds")
saveRDS(s_obj, file = output_rds_path)
message("전처리된 Seurat 객체 저장 완료: ", output_rds_path)

message("pp.R 스크립트 실행 완료.")
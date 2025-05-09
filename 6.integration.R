# integration.R
# 여러 Seurat 객체 통합 스크립트
#
# 기능:
# 1. 통합할 Seurat 객체 파일 목록 (xlsx) 로드
# 2. 각 Seurat 객체 로드
# 3. 지정된 방법 (예: rpca)을 사용하여 데이터 통합
# 4. 통합된 데이터에 대한 후처리 (PCA, UMAP, Clustering)
# 5. 최종 통합된 Seurat 객체 저장
#
# 실행 예시:
# Rscript integration.R --data_list "/path/to/list.xlsx" \
#                       --method "rpca" \
#                       --output "/path/to/integrated_seurat"

# 라이브러리 로드
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(readxl)))
suppressMessages(suppressWarnings(library(optparse)))

# 커맨드 라인 인자 정의
option_list <- list(
  make_option(c("-l", "--data_list"), type="character", default=NULL, help="통합할 Seurat 객체 RDS 파일 경로 목록이 포함된 XLSX 파일", metavar="character"),
  make_option(c("-m", "--method"), type="character", default="rpca", help="통합 방법 ('rpca', 'cca', 'harmony', 등)", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="integrated_data", help="출력 통합 Seurat 객체 RDS 파일 접두사", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$data_list)){
  print_help(opt_parser)
  stop("데이터 리스트 파일 (--data_list) 경로를 지정해야 합니다.", call.=FALSE)
}

# 1. 데이터 리스트 로드
message("1. 통합할 데이터 리스트 로드 중: ", opt$data_list)
if (!file.exists(opt$data_list)) {
  stop("제공된 데이터 리스트 파일이 존재하지 않습니다: ", opt$data_list, call.=FALSE)
}
# XLSX 파일의 첫 번째 시트, 첫 번째 컬럼에 파일 경로가 있다고 가정
file_paths_df <- read_excel(opt$data_list, col_names = FALSE)
seurat_object_paths <- file_paths_df[[1]] # 첫 번째 컬럼 사용

if (length(seurat_object_paths) < 2) {
  stop("통합을 위해서는 최소 2개 이상의 데이터 경로가 필요합니다.", call.=FALSE)
}
message(length(seurat_object_paths), "개의 Seurat 객체를 통합합니다.")

# 2. Seurat 객체 로드
message("2. Seurat 객체 로드 중...")
seurat_list <- lapply(seurat_object_paths, function(path) {
  if (!file.exists(path)) {
    warning("파일을 찾을 수 없습니다: ", path, ". 이 파일은 건너<0xEB><0xB5>니다.")
    return(NULL)
  }
  message("  로드 중: ", path)
  s_obj <- readRDS(path)
  # 통합을 위해 모든 객체가 동일한 DefaultAssay (예: SCT)를 사용하도록 하거나,
  # PrepSCTIntegration 등이 이를 처리하도록 함.
  # 사용자 스크립트에서는 SCTransform된 객체들을 사용.
  # 만약 DefaultAssay가 RNA이고 SCT 데이터가 있다면 변경
  if ("SCT" %in% names(s_obj@assays) && DefaultAssay(s_obj) != "SCT"){
      DefaultAssay(s_obj) <- "SCT"
  } else if ("SCT_soupx" %in% names(s_obj@assays) && DefaultAssay(s_obj) != "SCT_soupx") {
      DefaultAssay(s_obj) <- "SCT_soupx"
  }
  # SCT나 유사한 assay가 없으면 경고
  if (!DefaultAssay(s_obj) %in% c("SCT", "SCT_soupx") && !("SCT" %in% names(s_obj@assays) || "SCT_soupx" %in% names(s_obj@assays))) {
      warning(paste0("객체 ", path, "의 DefaultAssay가 SCT 또는 SCT_soupx가 아닙니다. 통합 결과에 영향이 있을 수 있습니다."))
  }
  return(s_obj)
})
seurat_list <- seurat_list[!sapply(seurat_list, is.null)] # NULL인 항목 제거

if (length(seurat_list) < 2) {
  stop("로드된 Seurat 객체가 2개 미만입니다. 통합을 진행할 수 없습니다.", call.=FALSE)
}

# 3. 데이터 통합
message("3. 데이터 통합 중 (방법: ", opt$method, ")...")

if (opt$method == "rpca") {
  # SCTransform된 데이터 통합 (RPCA 사용)
  # 모든 객체가 SCTransform 되었는지 확인 (또는 유사한 assay 이름 사용)
  # 예: DefaultAssay가 "SCT" 또는 "SCT_soupx" 등인지 확인
  
  # PrepSCTIntegration 전에 variable features를 각 객체에서 다시 계산하거나,
  # 이미 계산된 variable features를 사용하도록 SelectIntegrationFeatures를 먼저 실행.
  # 사용자 스크립트처럼 SelectIntegrationFeatures 사용
  features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000, verbose = FALSE)
  seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = features, verbose = FALSE)
  
  # RPCA를 사용하기 위해 각 객체에 대해 RunPCA 실행 (PrepSCTIntegration 후에도 필요할 수 있음)
  # FindIntegrationAnchors에서 reduction="rpca"를 사용하면 내부적으로 PCA를 활용
  # 사용자 스크립트에서는 k.anchor=20 사용
  anchors <- FindIntegrationAnchors(object.list = seurat_list,
                                    normalization.method = "SCT", # 중요: SCT 데이터임을 명시
                                    anchor.features = features,
                                    reduction = "rpca", # RPCA 사용
                                    k.anchor = 20, # 사용자 설정값
                                    verbose = FALSE)
  integrated_data <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = FALSE)
  # 통합 후 DefaultAssay는 'integrated'가 됨.

} else if (opt$method == "cca") {
  # CCA 기반 통합 (SCT 데이터의 경우)
  features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000, verbose = FALSE)
  seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = features, verbose = FALSE)
  anchors <- FindIntegrationAnchors(object.list = seurat_list,
                                    normalization.method = "SCT",
                                    anchor.features = features,
                                    reduction = "cca", # CCA 사용
                                    verbose = FALSE)
  integrated_data <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = FALSE)

} else if (opt$method == "harmony") {
  # Harmony 통합 (일반적으로 LogNormalized 데이터에 사용, SCT 데이터에도 적용 가능)
  # Harmony는 보통 하나의 Seurat 객체로 merge 후 실행
  # 여기서는 개별 객체에 대해 PCA 등을 실행 후, 그 결과를 이용
  # merged_obj <- merge(seurat_list[[1]], y = seurat_list[-1])
  # merged_obj <- NormalizeData(merged_obj) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
  # library(harmony)
  # integrated_data <- RunHarmony(merged_obj, group.by.vars = "orig.ident") # orig.ident 또는 다른 배치 변수
  stop("Harmony 통합은 이 스크립트에서 아직 완전하게 지원되지 않습니다. RPCA 또는 CCA를 사용하세요.")
} else {
  stop("지원하지 않는 통합 방법입니다. 'rpca', 'cca' 등을 사용하세요.", call.=FALSE)
}

message("데이터 통합 완료.")

# 4. 통합 후 처리
message("4. 통합 후 데이터 처리 중...")
# 통합된 데이터에 대해 PCA, UMAP, Clustering 등 수행
# 만약 SCT 기반으로 통합했다면, DefaultAssay(integrated_data)는 'integrated'
# 사용자 스크립트에서는 통합 후 다시 SCTransform을 수행하는 부분이 있는데,
# 이는 여러 SCT 모델을 병합하기 위한 것이나, IntegrateData(normalization.method = "SCT")가 이를 처리.
# 불필요할 수 있으므로 주석 처리하고, 필요시 해제.
# integrated_data <- SCTransform(integrated_data, assay = "RNA", new.assay.name = "SCT_integrated", verbose = FALSE)
# DefaultAssay(integrated_data) <- "SCT_integrated" # 또는 기존 'integrated' assay 사용

DefaultAssay(integrated_data) <- "integrated" # IntegrateData의 결과 assay 사용
integrated_data <- RunPCA(integrated_data, verbose = FALSE)
integrated_data <- RunUMAP(integrated_data, dims = 1:30, verbose = FALSE) # dims는 데이터에 따라 조절
integrated_data <- FindNeighbors(integrated_data, dims = 1:30, verbose = FALSE)
integrated_data <- FindClusters(integrated_data, verbose = FALSE) # graph.name은 필요시 지정 (예: "integrated_snn")

message("통합 후 처리 완료.")

# (선택적) 마커 유전자 탐색 - 사용자 스크립트 부분
# DefaultAssay(integrated_data) <- "SCT" # 또는 "integrated" assay의 "data" slot 사용
# integrated_data <- PrepSCTFindMarkers(integrated_data) # SCT 데이터용
# markers_list <- list()
# if ("seurat_clusters" %in% colnames(integrated_data@meta.data)) {
#   for(cluster_id in levels(integrated_data$seurat_clusters)){
#     message(paste0("Finding markers for cluster: ", cluster_id))
#     markers <- FindMarkers(integrated_data, ident.1 = cluster_id, verbose = FALSE)
#     markers$gene <- rownames(markers)
#     markers_list[[as.character(cluster_id)]] <- markers
#   }
#   saveRDS(markers_list, file = paste0(opt$output, "_markers.rds"))
#   message("마커 유전자 탐색 및 저장 완료.")
# } else {
#   warning("클러스터 정보 ('seurat_clusters')가 없어 마커 유전자 탐색을 건너<0xEB><0xB5>니다.")
# }


# 5. 결과 저장
output_rds_path <- paste0(opt$output, ".rds")
saveRDS(integrated_data, file = output_rds_path)
message("통합된 Seurat 객체 저장 완료: ", output_rds_path)

message("integration.R 스크립트 실행 완료.")
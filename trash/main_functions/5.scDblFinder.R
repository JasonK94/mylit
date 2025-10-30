# scDblFinder.R
# scDblFinder를 이용한 doublet 제거 스크립트
#
# 기능:
# 1. Seurat 객체 로드 (soupx.R의 출력)
# 2. scDblFinder 실행하여 doublet 예측
# 3. 예측된 singlet만 선택하여 Seurat 객체 subset
# 4. 최종 Seurat 객체 저장
#
# 실행 예시:
# Rscript scDblFinder.R --data "/path/to/soupx_corrected_seurat.rds" \
#                       --output "/path/to/scDblFinder_singlets_seurat"

# 라이브러리 로드
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(scDblFinder)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(optparse)))
suppressMessages(suppressWarnings(library(SingleCellExperiment))) # scDblFinder는 SCE 객체를 반환

# 커맨드 라인 인자 정의
option_list <- list(
  make_option(c("-d", "--data"), type="character", default=NULL, help="입력 Seurat 객체 RDS 파일 경로 (soupx.R의 출력)", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="scDblFinder_singlets", help="출력 Seurat 객체 RDS 파일 접두사", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$data)){
  print_help(opt_parser)
  stop("입력 Seurat 객체 (--data) 경로를 지정해야 합니다.", call.=FALSE)
}

# 1. 데이터 로드
message("1. 데이터 로드 중: ", opt$data)
s_obj <- readRDS(opt$data)
message("Seurat 객체 로드 완료. 초기 세포 수: ", ncol(s_obj))

# 2. scDblFinder 실행
message("2. scDblFinder 실행 중...")
# scDblFinder는 SingleCellExperiment 객체 또는 count matrix를 입력으로 받음
# Seurat 객체의 count 데이터와 클러스터 정보 사용
# 사용자의 원본 스크립트에서는 GetAssayData(..., layer="counts")와 Idents()를 사용
# DefaultAssay의 counts를 사용한다고 가정 (SCT_soupx assay)

if (!DefaultAssay(s_obj) %in% names(s_obj@assays)) {
    stop(paste("Default assay", DefaultAssay(s_obj), "가 객체에 존재하지 않습니다."))
}
counts_for_dblfinder <- GetAssayData(s_obj, layer = "counts") # DefaultAssay의 counts 사용

# 클러스터 정보 (scDblFinder는 클러스터 정보를 활용하여 doublet을 더 잘 찾을 수 있음)
if (!"seurat_clusters" %in% colnames(s_obj@meta.data)) {
    warning("Seurat 객체에 'seurat_clusters' 정보가 없습니다. Idents(s_obj)를 사용합니다.")
    clusters_for_dblfinder <- Idents(s_obj)
    if (length(unique(clusters_for_dblfinder)) <= 1) {
        warning("클러스터 정보가 단일하거나 없습니다. scDblFinder 성능에 영향이 있을 수 있습니다.")
    }
} else {
    clusters_for_dblfinder <- s_obj$seurat_clusters
}

# scDblFinder 실행
# scDblFinder는 SingleCellExperiment 객체를 반환하므로, 결과에서 doublet 정보만 추출
sce <- scDblFinder(counts_for_dblfinder, clusters = clusters_for_dblfinder)
doublet_class <- sce$scDblFinder.class # "singlet" 또는 "doublet"

# Seurat 객체에 doublet 정보 추가
s_obj$scDblFinder_class <- doublet_class
s_obj$scDblFinder_score <- sce$scDblFinder.score # 필요시 score도 추가

message("scDblFinder 실행 완료.")
table(s_obj$scDblFinder_class)

# 3. Singlet만 선택
message("3. Singlet만 선택하여 subset 중...")
s_obj_singlets <- subset(s_obj, subset = scDblFinder_class == "singlet")
message("Doublet 제거 후 세포 수: ", ncol(s_obj_singlets))

# (선택적) Singlet만 남은 객체에 대해 재-클러스터링 등 후처리 가능
# s_obj_singlets <- FindNeighbors(s_obj_singlets, dims = 1:30, verbose = FALSE)
# s_obj_singlets <- FindClusters(s_obj_singlets, verbose = FALSE)

# 4. 결과 저장
output_rds_path <- paste0(opt$output, ".rds")
saveRDS(s_obj_singlets, file = output_rds_path)
message("Doublet 제거된 Seurat 객체 저장 완료: ", output_rds_path)

message("scDblFinder.R 스크립트 실행 완료.")
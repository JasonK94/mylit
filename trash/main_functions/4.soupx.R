# soupx.R
# SoupX를 이용한 ambient RNA 제거 스크립트
#
# 기능:
# 1. 전처리된 Seurat 객체 (pp.R의 출력) 및 raw count matrix 로드
# 2. SoupX 파이프라인 실행 (SoupChannel, setClusters, autoEstCont, adjustCounts)
# 3. SoupX로 보정된 count matrix로 새로운 Seurat 객체 생성 또는 기존 객체 업데이트
# 4. 재-전처리 (Normalization, PCA, Clustering)
# 5. 최종 Seurat 객체 저장
#
# 실행 예시:
# Rscript soupx.R --data "/path/to/output_prefix.rds" \
#                 --raw "/path/to/output_prefix_raw_counts.rds" \
#                 --output "/path/to/soupx_corrected_seurat"

# 라이브러리 로드
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(SoupX)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(optparse)))

# 커맨드 라인 인자 정의
option_list <- list(
  make_option(c("-d", "--data"), type="character", default=NULL, help="입력 Seurat 객체 RDS 파일 경로 (pp.R의 출력)", metavar="character"),
  make_option(c("-r", "--raw"), type="character", default=NULL, help="입력 raw count matrix RDS 파일 경로 (pp.R의 출력)", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="soupx_corrected", help="출력 Seurat 객체 RDS 파일 접두사", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$data) || is.null(opt$raw)){
  print_help(opt_parser)
  stop("입력 Seurat 객체 (--data)와 raw count matrix (--raw) 경로를 모두 지정해야 합니다.", call.=FALSE)
}

# 1. 데이터 로드
message("1. 데이터 로드 중...")
s_obj_prev <- readRDS(opt$data)
message("Seurat 객체 로드 완료: ", opt$data)
raw_matrix <- readRDS(opt$raw) # Read10X로 읽은, CreateSeuratObject 전의 matrix
message("Raw count matrix 로드 완료: ", opt$raw)

# SoupX는 필터링되지 않은 전체 droplet 정보 (tod)와 필터링된 세포의 count 정보 (toc)를 필요로 함.
# toc는 현재 s_obj_prev의 counts (RNA assay의 counts 또는 SCT assay의 counts).
# tod는 raw_matrix (모든 바코드 포함).

# toc (table of counts) 준비: s_obj_prev의 count 데이터 사용
# SCTransform을 사용했다면 s_obj_prev@assays$SCT$counts를 사용할 수 있으나,
# SoupX 문서에서는 raw UMI counts (RNA assay)를 권장하는 경우가 많음.
# 사용자의 원본 스크립트에서는 SCT counts를 사용했으므로 그 방식을 따름.
if ("SCT" %in% names(s_obj_prev@assays)) {
    toc <- GetAssayData(s_obj_prev, assay = "SCT", layer = "counts")
    message("SCT assay의 counts를 toc로 사용.")
} else if ("RNA" %in% names(s_obj_prev@assays)) {
    toc <- GetAssayData(s_obj_prev, assay = "RNA", layer = "counts")
    message("RNA assay의 counts를 toc로 사용.")
} else {
    stop("Seurat 객체에 SCT 또는 RNA assay의 count 데이터가 없습니다.")
}

# tod (table of droplets) 준비: pp.R에서 저장한 raw_counts 사용
# tod와 toc의 유전자 목록(rownames) 일치시키기
common_genes <- intersect(rownames(raw_matrix), rownames(toc))
tod <- raw_matrix[common_genes, ]
toc <- toc[common_genes, colnames(toc)] # toc의 cell은 s_obj_prev의 cell과 동일

# 2. SoupX 실행
message("2. SoupX 실행 중...")
# SoupChannel 생성 시 tod는 모든 droplet, toc는 분석 대상 cell의 count.
# tod의 colnames는 모든 바코드, toc의 colnames는 s_obj_prev의 세포 바코드.
# raw_matrix의 바코드와 toc의 바코드가 일치하는지 확인 필요.
# CreateSeuratObject 시 add.cell.ids로 GEM이 추가되었으므로, raw_matrix의 바코드도 그에 맞게 조정되었거나,
# 혹은 SoupX가 알아서 매칭할 수 있도록 바코드 형식을 통일해야 함.
# 여기서는 raw_matrix가 원본 바코드, toc가 GEM이 추가된 바코드라고 가정. 이 경우 문제가 될 수 있음.
# pp.R 에서 저장한 raw_counts_for_soupx 는 CreateSeuratObject 이전의 것이므로, 원래 바코드임.
# toc의 바코드는 s_obj_prev의 바코드임. (예: "AAACCCAAGAGCATC-1_GEM1")

# toc의 바코드에서 GEM 접미사를 제거하여 tod의 바코드와 매칭 시도 (만약 tod가 원본 바코드라면)
# 하지만, SoupX는 toc의 메타데이터 (클러스터 정보)를 사용하므로 toc의 바코드는 s_obj_prev와 일치해야 함.
# tod는 toc에 있는 모든 cell을 포함해야 함.
# 따라서 tod는 raw_matrix 전체를 사용하고, SoupX가 내부적으로 처리하도록 함.
# toc의 세포들이 tod에 존재하는지 확인
cells_in_tod <- colnames(toc)[colnames(toc) %in% colnames(tod)]
if(length(cells_in_tod) != ncol(toc)){
    warning("일부 세포가 raw matrix (tod)에 존재하지 않습니다. SoupX 결과에 영향을 줄 수 있습니다.")
    # toc를 tod에 있는 세포들로만 구성해야 할 수 있음.
    # s_obj_prev도 그에 맞게 subset해야 함. 이는 복잡하므로,
    # pp.R에서 raw_matrix를 저장할 때 CreateSeuratObject 직전의, 필터링 안된 count를 저장하는 것이 중요.
}


sc <- SoupChannel(tod = tod, toc = toc, calcSoupProfile = FALSE) # calcSoupProfile은 autoEstCont에서 수행

# 클러스터 정보 설정
# s_obj_prev의 seurat_clusters 또는 Idents() 사용
if (!"seurat_clusters" %in% colnames(s_obj_prev@meta.data)) {
    # s_obj_prev가 클러스터링되지 않았다면 임시 클러스터링 또는 에러 처리
    # 여기서는 FindClusters를 이미 했다고 가정
    stop("Seurat 객체에 'seurat_clusters' 메타데이터가 없습니다. pp.R에서 클러스터링이 수행되었는지 확인하세요.")
}
clusters <- s_obj_prev@meta.data[colnames(toc), "seurat_clusters"] # toc의 순서에 맞게 클러스터 정보 가져오기
sc <- setClusters(sc, clusters = clusters)
message("SoupChannel 생성 및 클러스터 정보 설정 완료.")

# 오염도 자동 추정 및 count 보정
# 사용자의 파라미터 사용 (tfidfMin, soupQuantile, forceAccept)
sc <- autoEstCont(sc, tfidfMin = 0.05, soupQuantile = 0.8, forceAccept = TRUE, verbose = FALSE)
message("오염도 추정 완료.")
out_matrix <- adjustCounts(sc, roundToInt = TRUE, verbose = FALSE) # 보정된 count matrix (정수형 권장)
message("Ambient RNA 보정 완료.")

# 3. 보정된 count로 Seurat 객체 생성
message("3. 보정된 count로 새로운 Seurat 객체 생성 중...")
# 기존 메타데이터를 유지하면서 count matrix만 교체
s_obj_soupx <- CreateSeuratObject(counts = out_matrix, meta.data = s_obj_prev@meta.data)
# 만약 assay가 여러 개였다면, DefaultAssay를 RNA로 설정하고 RNA assay의 counts를 교체하는 방식도 고려
# s_obj_soupx <- s_obj_prev
# s_obj_soupx <- SetAssayData(s_obj_soupx, assay="RNA", slot="counts", new.data=out_matrix)

# 4. 재-전처리
message("4. SoupX 보정 후 재-전처리 중...")
# pp.R에서 사용한 normalization 방식과 동일하게 적용 (또는 사용자가 선택)
# 여기서는 SCTransform을 기본으로 가정 (사용자 스크립트처럼)
# SCTransform은 raw counts에 적용되므로, out_matrix (보정된 raw counts) 사용
DefaultAssay(s_obj_soupx) <- "RNA" # SCTransform은 RNA assay의 counts를 사용
s_obj_soupx <- SCTransform(s_obj_soupx, assay = "RNA", new.assay.name = "SCT_soupx", verbose = FALSE)
DefaultAssay(s_obj_soupx) <- "SCT_soupx"

s_obj_soupx <- FindVariableFeatures(s_obj_soupx, verbose = FALSE)
s_obj_soupx <- RunPCA(s_obj_soupx, verbose = FALSE)
s_obj_soupx <- FindNeighbors(s_obj_soupx, dims = 1:30, verbose = FALSE) # dims는 데이터에 따라 조절
s_obj_soupx <- FindClusters(s_obj_soupx, verbose = FALSE)

message("재-전처리 완료. 세포 수: ", ncol(s_obj_soupx))

# 5. 결과 저장
output_rds_path <- paste0(opt$output, ".rds")
saveRDS(s_obj_soupx, file = output_rds_path)
message("SoupX로 보정된 Seurat 객체 저장 완료: ", output_rds_path)

message("soupx.R 스크립트 실행 완료.")
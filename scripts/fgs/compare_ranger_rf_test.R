# 환경 변수 설정
Sys.setenv(
    FGS_MAX_CPU_CORES = "4",
    FGS_BLAS_THREADS = "1",
    FGS_DISABLE_PARALLEL = "FALSE"
)

source("/home/user3/data_user3/git_repo/_wt/fgs/scripts/fgs/init_fgs_env.R")
renv::activate("/home/user3/GJC_KDW_250721")

library(qs)
library(Seurat)

# 직접 source하여 최신 코드 사용
source("/home/user3/data_user3/git_repo/_wt/fgs/myR/R/signature.R")

# 데이터 로드
message("Loading dataset...")
is5a <- qs::qread("/data/user3/sobj/IS6_sex_added_251110.qs")

# Naive B-cells 선택 (적당한 크기)
cluster_name <- "Naive B-cells"
sobj <- subset(is5a, anno3.scvi == cluster_name)
if (ncol(sobj) > 300) {
    set.seed(42)
    sobj <- subset(sobj, cells = sample(Cells(sobj), 300))
}
message(sprintf("Subset created for %s: %d cells", cluster_name, ncol(sobj)))

# FGS 실행 (RF vs Ranger)
message("Running FGS with random_forest and random_forest_ranger...")
fgsa <- find_gene_signature_v5.4(
    sobj,
    target_var = "g3",
    control_vars = "hos_no",
    n_features = 50, # 적은 수로 테스트
    method = c("random_forest", "random_forest_ranger"),
    min_cells = 10
)

# 결과 비교
sig_rf <- fgsa$random_forest
sig_ranger <- fgsa$random_forest_ranger

if (is.null(sig_rf)) message("random_forest signature is NULL")
if (is.null(sig_ranger)) message("random_forest_ranger signature is NULL")

if (!is.null(sig_rf) && !is.null(sig_ranger)) {
    genes_rf <- names(sig_rf)
    genes_ranger <- names(sig_ranger)

    common <- base::intersect(genes_rf, genes_ranger)
    jaccard <- length(common) / length(base::union(genes_rf, genes_ranger))

    message(sprintf("\n=== Comparison: RF vs Ranger ==="))
    message(sprintf("RF genes: %d", length(genes_rf)))
    message(sprintf("Ranger genes: %d", length(genes_ranger)))
    message(sprintf("Common genes: %d", length(common)))
    message(sprintf("Jaccard Index: %.4f", jaccard))

    message("\nTop 10 RF genes:")
    print(head(sort(sig_rf, decreasing = TRUE), 10))

    message("\nTop 10 Ranger genes:")
    print(head(sort(sig_ranger, decreasing = TRUE), 10))

    if (jaccard > 0.5) {
        message("\n✓ Result: High similarity between RF and Ranger.")
    } else {
        message("\n⚠ Result: Low similarity. Check if parameters match.")
    }
}

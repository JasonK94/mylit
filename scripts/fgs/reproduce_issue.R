# 환경 변수 설정 (init_fgs_env.R 이전에)
Sys.setenv(
    FGS_MAX_CPU_CORES = "16", # 최대 코어
    FGS_BLAS_THREADS = "4", # BLAS 4스레드
    FGS_DISABLE_PARALLEL = "FALSE" # 병렬 처리 활성화
)

# FGS 환경 초기화
if (file.exists("/home/user3/data_user3/git_repo/_wt/fgs/scripts/fgs/init_fgs_env.R")) {
    source("/home/user3/data_user3/git_repo/_wt/fgs/scripts/fgs/init_fgs_env.R")
}

renv::activate("/home/user3/GJC_KDW_250721")

# 필요한 패키지 로드
library(qs)
library(pROC)

# FGS 함수 로드 (devtools 권장)
# devtools::load_all("/home/user3/data_user3/git_repo/_wt/fgs/myR", quiet = TRUE)
source("/home/user3/data_user3/git_repo/_wt/fgs/myR/R/signature.R")
# Source other potentially needed files if functions are missing
# source('/home/user3/data_user3/git_repo/_wt/fgs/myR/R/utils_data.R')

# 데이터 로드
data_path <- "/data/user3/sobj/IS6_sex_added_251110_ds2500.qs" # Using downsampled version
if (!file.exists(data_path)) {
    stop(paste("Data file not found:", data_path))
}
is5a <- qs::qread(data_path)

# 에러 로그 파일 경로
error_log_file <- "fgs_errors.log"

# Test on full dataset (not subsetting by cluster)
print("Using full dataset (all clusters combined)")
print(paste("Data dimensions:", paste(dim(is5a), collapse = "x")))
print("Target variable (g3) distribution:")
print(table(is5a$g3, useNA = "always"))

# FGS execution - testing ranger and nmf_loadings (the ones we fixed!)
fgs_methods <- c("random_forest_ranger", "nmf_loadings")
print(paste("Running FGS with methods:", paste(fgs_methods, collapse = ", ")))

fgs_list <- list()
for (m in fgs_methods) {
    print(paste("\n--- Running Method:", m, "---"))
    tryCatch(
        {
            res_list <- find_gene_signature_v5_impl(
                data = is5a,
                target_var = "g3",
                control_vars = NULL,
                method = m,
                n_features = 50,
                min_cells = 10,
                min_pct = 0.01,
                method_impl = "v5.4"
            )
            # find_gene_signature_v5_impl returns a list of results keyed by method name
            res <- res_list[[m]]
            fgs_list[[m]] <- res
            print(paste("✓", m, "SUCCEEDED"))
            print(paste("  Genes found:", length(res$genes)))
            print("  Top 5 genes:")
            print(head(res$genes, 5))
        },
        error = function(e) {
            print(paste("✗", m, "FAILED:", e$message))
            traceback()
        }
    )
}

message("\n=== FGS Test Results ===")
print(list(
    methods_attempted = fgs_methods,
    methods_succeeded = names(fgs_list),
    ranger_ok = "random_forest_ranger" %in% names(fgs_list),
    nmf_ok = "nmf_loadings" %in% names(fgs_list)
))

fgsa <- fgs_list
message("FGS completed successfully.")

# TML7 execution
message("Running TML7...")
tmla <- TML7(
    l1_signatures = fgsa,
    holdout_data = is5a, # Use full dataset
    target_var = "g3",
    cv_group_var = NULL, # Removed to simplify
    l2_methods = c("ranger", "earth", "glm")
)

message("TML7 completed successfully.")

# Test CMGI
message("Testing CMGI...")

# Test GLM (should work)
cmgi_glm <- try(compute_meta_gene_importance(tmla, target_model = "glm"), silent = FALSE)
if (!inherits(cmgi_glm, "try-error")) {
    message("✓ CMGI for GLM succeeded")
    print(head(cmgi_glm$gene_summary))
} else {
    message("✗ CMGI for GLM failed")
}

# Test ranger (this is what we're debugging)
cmgi_ranger <- try(compute_meta_gene_importance(tmla, target_model = "ranger"), silent = FALSE)
if (!inherits(cmgi_ranger, "try-error")) {
    message("✓ CMGI for ranger succeeded")
    print(head(cmgi_ranger$gene_summary))
} else {
    message("✗ CMGI for ranger failed")
}

# Test earth (may fail gracefully as documented)
cmgi_earth <- try(compute_meta_gene_importance(tmla, target_model = "earth"), silent = FALSE)
if (!inherits(cmgi_earth, "try-error") && !is.null(cmgi_earth)) {
    message("✓ CMGI for earth succeeded")
    print(head(cmgi_earth$gene_summary))
} else {
    message("✗ CMGI for earth failed or returned NULL (expected)")
}

message("\n=== All tests completed ===")
print(list(
    fgs_methods_run = names(fgsa),
    l2_models_trained = names(tmla$trained_models),
    cmgi_glm_ok = !inherits(cmgi_glm, "try-error"),
    cmgi_ranger_ok = !inherits(cmgi_ranger, "try-error"),
    cmgi_earth_ok = !inherits(cmgi_earth, "try-error") && !is.null(cmgi_earth)
))

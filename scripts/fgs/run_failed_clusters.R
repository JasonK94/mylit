# 환경 변수 설정
Sys.setenv(
    FGS_MAX_CPU_CORES = "16",
    FGS_BLAS_THREADS = "4",
    FGS_DISABLE_PARALLEL = "FALSE"
)

# FGS 환경 초기화
source("/home/user3/data_user3/git_repo/_wt/fgs/scripts/fgs/init_fgs_env.R")

renv::activate("/home/user3/GJC_KDW_250721")

library(qs)
library(pROC)
library(Seurat) # Seurat 객체 처리를 위해 필요

# FGS 함수 로드
devtools::load_all("/home/user3/data_user3/git_repo/_wt/fgs/myR", quiet = TRUE)

# 데이터 로드
message("Loading full dataset...")
is5a <- qs::qread("/data/user3/sobj/IS6_sex_added_251110.qs")

# 에러 로그 파일 경로
error_log_file <- "/data/user3/sobj/fgs/fgs_errors.log"

# 실패한 것들 (이미 생성된 것들 제외)
fgs_dirs <- dir("/data/user3/sobj/fgs/fgsas", full.names = TRUE)
# 정규 표현식 (파일명 패턴에 맞게 조정)
# 예: fgsa_Activated Monocytes _ Macrophages_25-11-21-10-00.qs
regex_pattern <- "^.*fgsa_(.*)_\\d{2}-\\d{2}-\\d{2}-\\d{2}-\\d{2}\\.qs$"

if (length(fgs_dirs) > 0) {
    extracted_names <- sapply(fgs_dirs, function(p) {
        # 파일명만 추출
        fname <- basename(p)
        # 패턴 매칭
        if (grepl(regex_pattern, fname)) {
            sub(pattern = regex_pattern, replacement = "\\1", x = fname)
        } else {
            NA
        }
    })
    extracted_names <- extracted_names[!is.na(extracted_names)]
    # URL decoding might be needed if names were encoded, but assuming simple match for now
    # Also handle potential space replacements if any
} else {
    extracted_names <- character(0)
}

all_clusters <- unique(is5a$anno3.scvi)
fgs_keys <- all_clusters[!all_clusters %in% extracted_names]

message(sprintf("Total clusters: %d", length(all_clusters)))
message(sprintf("Already processed: %d", length(unique(extracted_names))))
message(sprintf("Clusters to process: %d", length(fgs_keys)))

if (length(fgs_keys) == 0) {
    message("No failed/missing clusters to process.")
    quit(save = "no")
}

# 각 클러스터별로 FGS + TML7 실행
results <- lapply(fgs_keys, function(x) {
    tryCatch(
        {
            message(sprintf("\n=== Processing cluster: %s ===\n", x))

            # 클러스터 서브셋
            sobj <- subset(is5a, anno3.scvi == x)

            # 데이터 체크
            if (ncol(sobj) < 10) {
                warning(sprintf("Cluster %s has too few cells (%d). Skipping.", x, ncol(sobj)))
                return(list(cluster = x, success = FALSE, error = "Too few cells"))
            }

            # FGS 실행
            message(sprintf("Running FGS for cluster %s...", x))
            # min_cells=3으로 설정하여 작은 클러스터도 허용
            fgsa <- find_gene_signature_v5.4(
                sobj,
                target_var = "g3",
                control_vars = "hos_no",
                n_features = 200,
                min_cells = 3
            )

            # FGS 저장
            # 파일명에 특수문자가 있을 수 있으므로 안전하게 처리하되, 전체 경로에는 make.names 쓰지 않음
            safe_x <- gsub("/", "_", x) # 슬래시만 언더바 처리
            timestamp <- format(Sys.time(), "%y-%m-%d-%H-%M")

            fgs_dir <- "/data/user3/sobj/fgs/fgsas"
            if (!dir.exists(fgs_dir)) dir.create(fgs_dir, recursive = TRUE)
            fgs_file <- file.path(fgs_dir, paste0("fgsa_", safe_x, "_", timestamp, ".qs"))

            qs::qsave(fgsa, fgs_file)
            message(sprintf("FGS saved: %s", fgs_file))

            # TML7 실행
            message(sprintf("Running TML7 for cluster %s...", x))
            tmla <- TML7(
                l1_signatures = fgsa,
                holdout_data = sobj,
                target_var = "g3",
                cv_group_var = "hos_no"
            )

            # TML7 저장
            tml_dir <- "/data/user3/sobj/fgs/tmlas"
            if (!dir.exists(tml_dir)) dir.create(tml_dir, recursive = TRUE)
            tml_file <- file.path(tml_dir, paste0("tmla_", safe_x, "_", timestamp, ".qs"))

            qs::qsave(tmla, tml_file)
            message(sprintf("TML7 saved: %s", tml_file))

            return(list(cluster = x, fgs_file = fgs_file, tml_file = tml_file, success = TRUE))
        },
        error = function(e) {
            error_msg <- sprintf("ERROR processing cluster %s: %s", x, conditionMessage(e))
            message(error_msg)

            # 에러 로그 파일에 기록
            log_entry <- sprintf(
                "[%s] Cluster: %s | Error: %s\n",
                format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                x,
                conditionMessage(e)
            )

            # 파일에 추가
            cat(log_entry, file = error_log_file, append = TRUE)

            return(list(cluster = x, success = FALSE, error = conditionMessage(e)))
        }
    )
})

# 최종 요약 출력
cat("\n=== Processing Summary ===\n")
successful <- sum(sapply(results, function(r) isTRUE(r$success)))
failed <- length(results) - successful
cat(sprintf("Successful: %d, Failed: %d\n", successful, failed))

if (failed > 0) {
    cat("\nFailed clusters:\n")
    for (r in results) {
        if (!isTRUE(r$success)) {
            cat(sprintf("  - %s: %s\n", r$cluster, r$error))
        }
    }
    cat(sprintf("\nSee error log: %s\n", error_log_file))
}

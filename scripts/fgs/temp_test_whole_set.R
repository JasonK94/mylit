source("/home/user3/data_user3/git_repo/_wt/fgs/myR/R/signature.R")
source("/home/user3/data_user3/git_repo/mylit/myR/R/fgs_core.R")
Sys.setenv(
    FGS_MAX_CPU_CORES = "16", # 최대 코어
    FGS_BLAS_THREADS = "4", # BLAS 4스레드
    FGS_DISABLE_PARALLEL = "FALSE" # 병렬 처리 활성화
)
source("/home/user3/data_user3/git_repo/_wt/fgs/scripts/fgs/init_fgs_env.R")

message("Reading data...")
is5s <- qs::qread("/data/user3/sobj/is2_IS_3_clustered.qs")
message("Running find_gene_signature_v5_impl...")
fgss <- find_gene_signature_v5_impl(is5s, target_var = "g3", control_vars = "hos_no", n_features = 200)
output_fgss <- paste0("/data/user3/sobj/fgs/fgss_", format(Sys.time(), "%y-%m-%d-%H-%M"), ".qs")
qs::qsave(fgss, output_fgss)
message("Saved fgss to ", output_fgss)

# fgss=qs::qread("/data/user3/sobj/fgs/fgss_25-11-18-08-56")
message("Running TML7...")
tmls <- TML7(fgss, is5s, "g3")
output_tmls <- paste0("/data/user3/sobj/fgs/tmls", format(Sys.time(), "%y-%m-%d-%H-%M"), ".qs")
qs::qsave(tmls, output_tmls)
message("Saved tmls to ", output_tmls)

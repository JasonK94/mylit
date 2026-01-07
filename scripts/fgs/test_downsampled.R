source("/home/user3/data_user3/git_repo/_wt/fgs/myR/R/signature.R")
source("/home/user3/data_user3/git_repo/_wt/fgs/myR/R/fgs_core.R")
Sys.setenv(
    FGS_MAX_CPU_CORES = "16", # 최대 코어
    FGS_BLAS_THREADS = "4", # BLAS 4스레드
    FGS_DISABLE_PARALLEL = "FALSE" # 병렬 처리 활성화
)
source("/home/user3/data_user3/git_repo/_wt/fgs/scripts/fgs/init_fgs_env.R")

message("Reading downsampled data (0.1x)...")
# Using the downsampled data path from guide.md
is5s_down <- qs::qread("/data/user3/sobj/IS6_sex_added_0.1x_251110.qs")
message("Data read. Dimensions: ", paste(dim(is5s_down), collapse = " x "))

message("Running find_gene_signature on downsampled data...")
fgss <- find_gene_signature(is5s_down, target_var = "g3", control_vars = "hos_no", n_features = 200)

output_fgss <- paste0("/data/user3/sobj/fgs/fgss_downsampled_", format(Sys.time(), "%y-%m-%d-%H-%M"), ".qs")
qs::qsave(fgss, output_fgss)
message("Saved fgss to ", output_fgss)

message("Running TML7 on downsampled data...")
tmls <- TML7(fgss, is5s_down, "g3")
output_tmls <- paste0("/data/user3/sobj/fgs/tmls_downsampled_", format(Sys.time(), "%y-%m-%d-%H-%M"), ".qs")
qs::qsave(tmls, output_tmls)
message("Saved tmls to ", output_tmls)

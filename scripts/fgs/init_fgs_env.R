#!/usr/bin/env Rscript
# FGS-specific environment initialization
# Uses start.R but with parallel processing disabled via environment variables
# This prevents cascade parallelization when child processes load start.R

# CPU configuration: can be overridden by environment variables
# Set FGS_MAX_CPU_CORES to control maximum CPU cores (default: 16)
# Set FGS_BLAS_THREADS to control BLAS/LAPACK threads (default: 1, for sequential)
max_cpu_cores <- as.integer(Sys.getenv("FGS_MAX_CPU_CORES", "16"))
blas_threads <- as.integer(Sys.getenv("FGS_BLAS_THREADS", "1"))
disable_parallel <- isTRUE(as.logical(Sys.getenv("FGS_DISABLE_PARALLEL", "TRUE")))

# Set environment variables BEFORE sourcing start.R
# These will be inherited by child processes and prevent them from spawning workers
Sys.setenv(
  # Disable parallel processing in start.R (can be overridden)
  MYLIT_DISABLE_PARALLEL = if (disable_parallel) "TRUE" else "FALSE",
  # Set memory limit (GB) - can be overridden
  MYLIT_FUTURE_MEMORY_GB = Sys.getenv("FGS_MEMORY_GB", "200"),
  # BLAS/LAPACK threads (1 = sequential, higher = parallel)
  OMP_NUM_THREADS = as.character(blas_threads),
  OPENBLAS_NUM_THREADS = as.character(blas_threads),
  MKL_NUM_THREADS = as.character(blas_threads),
  VECLIB_MAXIMUM_THREADS = as.character(blas_threads),
  NUMEXPR_NUM_THREADS = as.character(blas_threads),
  # Future configuration
  R_FUTURE_AVAILABLE_CORES = if (disable_parallel) "1" else as.character(max_cpu_cores),
  R_FUTURE_PLAN = if (disable_parallel) "sequential" else "multisession",
  # Prevent start.R from loading all packages
  KDW_START_LOAD_ALL_PACKAGES = "FALSE",
  KDW_START_AUTOLOAD_QS = "FALSE",
  # Store FGS settings for child processes
  FGS_MAX_CPU_CORES = as.character(max_cpu_cores),
  FGS_BLAS_THREADS = as.character(blas_threads),
  FGS_DISABLE_PARALLEL = if (disable_parallel) "TRUE" else "FALSE"
)
options(mc.cores = if (disable_parallel) 1L else max_cpu_cores)

# Change to directory with start.R
original_wd <- getwd()
if (file.exists("/home/user3/GJC_KDW_250721/start.R")) {
  setwd("/home/user3/GJC_KDW_250721")
  message("✓ Changed to directory: /home/user3/GJC_KDW_250721")
}

# Source start.R (it will respect MYLIT_DISABLE_PARALLEL=TRUE)
tryCatch({
  source("start.R")
  message("✓ start.R sourced (parallel processing disabled)")
}, error = function(e) {
  warning("Failed to source start.R: ", conditionMessage(e))
  # Fallback: just activate renv
  if (file.exists("renv/activate.R")) {
    source("renv/activate.R")
    message("✓ renv activated (fallback)")
  }
})

# Ensure sequential plan is set (double-check)
if (requireNamespace("future", quietly = TRUE)) {
  tryCatch({
    future::plan(future::sequential)
  }, error = function(e) {
    # Ignore if already set
  })
}

# Disable all parallel backends
if (requireNamespace("foreach", quietly = TRUE)) {
  tryCatch({
    foreach::registerDoSEQ()
  }, error = function(e) {
    # Ignore
  })
}

if (requireNamespace("doParallel", quietly = TRUE)) {
  tryCatch({
    doParallel::stopImplicitCluster()
  }, error = function(e) {
    # Ignore if no cluster
  })
}

message("\n=== FGS Environment Ready ===")
message(sprintf("  - start.R: sourced with MYLIT_DISABLE_PARALLEL=%s", if (disable_parallel) "TRUE" else "FALSE"))
message(sprintf("  - Memory limit: %sGB", Sys.getenv("MYLIT_FUTURE_MEMORY_GB", "200")))
if (disable_parallel) {
  message("  - Parallel workers: 0 (sequential only)")
} else {
  message(sprintf("  - Max CPU cores: %d", max_cpu_cores))
  message(sprintf("  - BLAS/LAPACK threads: %d", blas_threads))
}
message("  - Child processes will inherit these settings")
message("\nTo customize CPU usage, set environment variables before sourcing:")
message("  - FGS_MAX_CPU_CORES: maximum CPU cores (default: 16)")
message("  - FGS_BLAS_THREADS: BLAS/LAPACK threads (default: 1)")
message("  - FGS_DISABLE_PARALLEL: disable parallel (default: TRUE)")
message("  - FGS_MEMORY_GB: memory limit in GB (default: 200)")


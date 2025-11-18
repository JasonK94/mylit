#!/usr/bin/env Rscript
# FGS-specific environment initialization
# Uses start.R but with parallel processing disabled via environment variables
# This prevents cascade parallelization when child processes load start.R

# Set environment variables BEFORE sourcing start.R
# These will be inherited by child processes and prevent them from spawning workers
Sys.setenv(
  # Disable parallel processing in start.R
  MYLIT_DISABLE_PARALLEL = "TRUE",
  # Set memory limit (GB)
  MYLIT_FUTURE_MEMORY_GB = "200",
  # Force single-threaded BLAS/LAPACK
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1",
  # Prevent future from auto-detecting cores
  R_FUTURE_AVAILABLE_CORES = "1",
  R_FUTURE_PLAN = "sequential",
  # Prevent start.R from loading all packages
  KDW_START_LOAD_ALL_PACKAGES = "FALSE",
  KDW_START_AUTOLOAD_QS = "FALSE"
)
options(mc.cores = 1L)

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
message("  - start.R: sourced with MYLIT_DISABLE_PARALLEL=TRUE")
message("  - Memory limit: 200GB")
message("  - Parallel workers: 0 (sequential only)")
message("  - Child processes will inherit these settings")


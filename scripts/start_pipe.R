# Lightweight start.R for pipeline only
# Loads only essential packages for reading & integration pipeline

if (isTRUE(getOption("pipe.start.initialised", FALSE))) {
  message("✨ start_pipe.R already ran; skipping re-initialisation.")
  invisible(TRUE)
} else {
  options(pipe.start.initialised = TRUE)
  
  message("✨ Initialising pipeline session (start_pipe.R)")
  
  # Change to renv project directory and load
  renv_project <- "/home/user3/GJC_KDW_250721"
  if (file.exists(renv_project)) {
    old_wd <- getwd()
    setwd(renv_project)
    
    if (file.exists(".Rprofile")) {
      tryCatch({
        source(".Rprofile")
      }, error = function(e) {
        # Silently fail
      })
    }
    
    if (file.exists("renv") && requireNamespace("renv", quietly = TRUE)) {
      tryCatch({
        renv::load()
      }, error = function(e) {
        warning("Could not load renv: ", e$message)
      })
    }
    
    setwd(old_wd)
  }
  
  # ---- Essential packages for pipeline ----
  essential_packages <- c(
    "optparse",      # Command line argument parsing
    "Seurat",        # Core Seurat functions
    "SeuratWrappers", # Integration methods (RPCA, scVI)
    "SoupX",         # Ambient RNA correction
    "scDblFinder",  # Doublet detection
    "SingleCellExperiment", # For scDblFinder
    "dplyr",         # Data manipulation
    "stringr",       # String operations
    "qs",            # Fast serialization
    "future",        # Parallel processing
    "reticulate",    # Python integration (for scVI)
    "readr",         # CSV reading
    "ggplot2",       # Plotting
    "devtools"       # For loading myR
  )
  
  # Helper to load packages safely
  load_package <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      warning(sprintf("Package '%s' is not installed; skipping.", pkg), call. = FALSE)
      return(FALSE)
    }
    suppressPackageStartupMessages(
      library(pkg, character.only = TRUE)
    )
    TRUE
  }
  
  message("• Loading essential packages …")
  loaded <- vapply(essential_packages, load_package, logical(1))
  if (any(!loaded)) {
    message("  - Missing packages: ",
            paste(essential_packages[!loaded], collapse = ", "))
  }
  
  # ---- Load myR package ----
  myr_candidates <- c(
    "/home/user3/data_user3/git_repo/mylit/myR",
    "/data/user3/git_repo/_wt/pipe/../../myR"
  )
  myr_path <- NULL
  for (candidate in myr_candidates) {
    if (dir.exists(candidate)) {
      myr_path <- candidate
      break
    }
  }
  
  if (!is.null(myr_path) && requireNamespace("devtools", quietly = TRUE)) {
    message("• Loading local package 'myR' from ", myr_path, " …")
    tryCatch(
      devtools::load_all(myr_path, quiet = TRUE),
      error = function(e) {
        warning("Failed to load myR: ", conditionMessage(e), call. = FALSE)
        # Fallback: source essential files
        if (file.exists(file.path(myr_path, "R/utils_demulti.R"))) {
          source(file.path(myr_path, "R/utils_demulti.R"))
        }
      }
    )
  } else {
    warning("Could not locate myR package directory; skipped load_all.", call. = FALSE)
  }
  
  # ---- Configure future for parallel processing ----
  if (requireNamespace("future", quietly = TRUE)) {
    memory_limit <- as.numeric(Sys.getenv("PIPE_FUTURE_MEMORY_GB", "200"))
    options(future.globals.maxSize = memory_limit * 1024^3)
    message("• Configured future with memory limit: ", memory_limit, " GB")
  }
  
  message("✅ start_pipe.R initialisation complete.")
  invisible(TRUE)
}


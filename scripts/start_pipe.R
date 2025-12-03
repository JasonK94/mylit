# Lightweight start.R for pipeline only
# Loads step-specific packages based on STEP environment variable or step argument

# Get step number from environment variable or argument
step_num <- as.integer(Sys.getenv("PIPE_STEP", NA))
if (is.na(step_num)) {
  # Try to get from command line args
  cmd_args <- commandArgs(trailingOnly = TRUE)
  step_idx <- which(cmd_args == "--output_step" | cmd_args == "-o")
  if (length(step_idx) > 0 && step_idx < length(cmd_args)) {
    step_num <- as.integer(cmd_args[step_idx + 1])
  }
}

if (isTRUE(getOption("pipe.start.initialised", FALSE))) {
  message("✨ start_pipe.R already ran; skipping re-initialisation.")
  invisible(TRUE)
} else {
  options(pipe.start.initialised = TRUE)
  
  message("✨ Initialising pipeline session (start_pipe.R)")
  if (!is.na(step_num)) {
    message(sprintf("  Step %d detected", step_num))
  }
  
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
  
  # ---- Step-specific package lists ----
  # Common packages needed by all steps
  common_packages <- c(
    "optparse",      # Command line argument parsing
    "qs",            # Fast serialization
    "readr"          # CSV reading
  )
  
  # Step-specific packages
  step_packages <- list(
    step1 = c(
      "Seurat",        # Core Seurat functions
      "dplyr",         # Data manipulation
      "stringr"        # String operations
    ),
    step2 = c(
      "Seurat",        # Core Seurat functions
      "dplyr"          # Data manipulation
    ),
    step3 = c(
      "Seurat",        # Core Seurat functions
      "SoupX",         # Ambient RNA correction
      "ggplot2"        # Plotting
    ),
    step4 = c(
      "Seurat",        # Core Seurat functions
      "dplyr"          # Data manipulation
    ),
    step5 = c(
      "Seurat",        # Core Seurat functions
      "scDblFinder",  # Doublet detection
      "SingleCellExperiment", # For scDblFinder
      "ggplot2"        # Plotting
    ),
    step6 = c(
      "Seurat",        # Core Seurat functions
      "future"         # Parallel processing
    ),
    step6_scvi = c(
      "Seurat",        # Core Seurat functions
      "SeuratWrappers", # Integration methods (scVI)
      "future",        # Parallel processing
      "reticulate"     # Python integration (for scVI)
    ),
    step6_rpca = c(
      "Seurat",        # Core Seurat functions
      "future"         # Parallel processing
    )
  )
  
  # Determine which packages to load
  packages_to_load <- common_packages
  
  if (!is.na(step_num)) {
    step_key <- paste0("step", step_num)
    
    # Add step-specific packages
    if (step_key %in% names(step_packages)) {
      packages_to_load <- c(packages_to_load, step_packages[[step_key]])
    }
    
    # For step 6, check method if available
    if (step_num == 6) {
      cmd_args <- commandArgs(trailingOnly = TRUE)
      method_idx <- which(cmd_args == "--method" | cmd_args == "-m")
      if (length(method_idx) > 0 && method_idx < length(cmd_args)) {
        method <- tolower(cmd_args[method_idx + 1])
        method_key <- paste0("step6_", method)
        if (method_key %in% names(step_packages)) {
          packages_to_load <- c(packages_to_load, step_packages[[method_key]])
        }
      }
    }
  } else {
    # If step not detected, load all packages (fallback)
    message("  Step not detected; loading all packages (fallback mode)")
    all_step_packages <- unique(unlist(step_packages))
    packages_to_load <- c(packages_to_load, all_step_packages)
  }
  
  # Remove duplicates
  packages_to_load <- unique(packages_to_load)
  
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
  
  message(sprintf("• Loading %d packages for step %s …", 
                 length(packages_to_load), 
                 ifelse(is.na(step_num), "?", as.character(step_num))))
  loaded <- vapply(packages_to_load, load_package, logical(1))
  if (any(!loaded)) {
    message("  - Missing packages: ",
            paste(packages_to_load[!loaded], collapse = ", "))
  }
  
  # ---- Load myR package ----
  # Only load myR for steps that need it (Step 1, 6)
  load_myr <- is.na(step_num) || step_num %in% c(1, 6)
  
  if (load_myr) {
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
  }
  
  # ---- Configure future for parallel processing ----
  # Only for steps that need it (Step 6)
  if (step_num == 6 && requireNamespace("future", quietly = TRUE)) {
    memory_limit <- as.numeric(Sys.getenv("PIPE_FUTURE_MEMORY_GB", "200"))
    options(future.globals.maxSize = memory_limit * 1024^3)
    message("• Configured future with memory limit: ", memory_limit, " GB")
  }
  
  message("✅ start_pipe.R initialisation complete.")
  invisible(TRUE)
}

# ==============================================================================
# FGS/TML Main Entry Point (Legacy Support)
# ==============================================================================
# This file is maintained for backward compatibility.
# The actual implementations have been moved to:
# - fgs_core.R
# - tml_core.R
# - tml_utils.R
# - seurat_utils.R
# ==============================================================================

# Note: When loaded as a package (devtools::load_all or library(myR)),
# all files in R/ are automatically loaded, so this file doesn't need to do anything.

# However, if this file is sourced directly via source(), we need to manually
# source the split files to maintain compatibility.

# Determine directory
current_dir <- NULL
if (sys.nframe() > 0) {
  try(
    {
      if (!is.null(sys.frame(1)$ofile)) {
        current_dir <- dirname(sys.frame(1)$ofile)
      }
    },
    silent = TRUE
  )
}

if (is.null(current_dir)) {
  # Fallback: check if we are in project root
  if (file.exists("myR/R/fgs_core.R")) {
    current_dir <- "myR/R"
  } else if (file.exists("fgs_core.R")) {
    current_dir <- "."
  }
}

if (!is.null(current_dir)) {
  # List of files to source
  files_to_source <- c(
    "seurat_utils.R",
    "fgs_core.R",
    "tml_core.R",
    "tml_utils.R"
  )

  for (f in files_to_source) {
    f_path <- file.path(current_dir, f)
    if (file.exists(f_path)) {
      source(f_path, local = FALSE)
    } else {
      # Try without path if in same dir
      if (file.exists(f)) {
        source(f, local = FALSE)
      } else {
        # Only warn if we are sure about the directory but file is missing
        if (current_dir != ".") {
          warning(sprintf("Could not find %s to source from signature.R (dir: %s)", f, current_dir))
        }
      }
    }
  }
}

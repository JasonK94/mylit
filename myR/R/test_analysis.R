# ------------------------------------------------------------------------------
# Deprecated: Differential expression helpers relocated to myR/R/analysis.R
# ------------------------------------------------------------------------------
# This shim keeps backwards compatibility by sourcing the consolidated
# implementation from `analysis.R`. Please update scripts to reference
# `myR/R/analysis.R` directly.

.test_analysis_source <- function() {
  candidates <- character(0)
  this_file <- tryCatch(normalizePath(sys.frames()[[1]]$ofile),
                        error = function(...) NA_character_)
  if (!is.na(this_file)) {
    candidates <- c(candidates, file.path(dirname(this_file), "analysis.R"))
  }
  system_path <- system.file("R", "analysis.R", package = "myR")
  candidates <- c(
    candidates,
    "analysis.R",
    file.path("R", "analysis.R"),
    file.path("myR", "R", "analysis.R"),
    file.path(getwd(), "myR", "R", "analysis.R"),
    system_path
  )
  analysis_path <- NULL
  for (path in unique(candidates)) {
    if (!is.null(path) && file.exists(path)) {
      analysis_path <- normalizePath(path)
      break
    }
  }
  if (is.null(analysis_path)) {
    stop("Unable to locate analysis.R while sourcing test_analysis.R. ",
         "Please source myR/R/analysis.R directly.")
  }
  message("⚠️  test_analysis.R is deprecated. Sourcing ", analysis_path)
  source(analysis_path, local = parent.frame())
}
.test_analysis_source()
rm(.test_analysis_source)

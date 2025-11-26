local({
  root_dir <- normalizePath(".", winslash = "/", mustWork = FALSE)

  activate <- file.path(root_dir, "renv", "activate.R")
  if (file.exists(activate)) {
    try(source(activate), silent = TRUE)
  }

  start_file <- file.path(root_dir, "st", "start.R")
  if (interactive() && file.exists(start_file)) {
    tryCatch(
      source(start_file, local = FALSE),
      error = function(e) warning("Failed to source st/start.R: ", conditionMessage(e), call. = FALSE)
    )
  }
})


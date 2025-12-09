# Helper script to load renv environment - use lightweight start_pipe.R instead
source(file.path(dirname(normalizePath(commandArgs(trailingOnly = FALSE)[grep("--file", commandArgs(trailingOnly = FALSE))])), "start_pipe.R"))


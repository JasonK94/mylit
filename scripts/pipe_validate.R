#!/usr/bin/env Rscript
# Validation script for pipeline configuration
# Usage: Rscript pipe_validate.R --config <config_path> [--run_id <run_id>]

# Load lightweight pipeline environment
# Find script directory
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- cmd_args[grepl("--file=", cmd_args)]
if (length(file_arg) > 0) {
  script_path <- sub("--file=", "", file_arg)
  script_dir <- dirname(normalizePath(script_path))
} else {
  script_dir <- getwd()
}
source(file.path(script_dir, "start_pipe.R"))

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
})

# Parse command line arguments
option_list <- list(
  make_option(c("--config", "-c"), type = "character", default = NULL,
              help = "Path to config.csv file", metavar = "character"),
  make_option(c("--run_id", "-r"), type = "character", default = "run1",
              help = "Run ID", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$config)) {
  print_help(opt_parser)
  stop("--config argument is required", call. = FALSE)
}

# Source utility functions
source(file.path(dirname(opt$config), "../myR/R/pipe_utils.R"))

# Setup logging
log_list <- setup_logging(opt$run_id)
log_message(sprintf("Starting validation for config: %s", opt$config), log_list)

# Load config
config_list <- load_config(opt$config)
config <- config_list$config
defaults <- config_list$defaults

errors <- character(0)
warnings <- character(0)

# Validation 1: Check required columns in config.csv
required_cols <- c("no", "name", "sample_name", "gem_name", "multiplex_method", 
                   "dir_input_filtered_barcode_matrix")
missing_cols <- setdiff(required_cols, colnames(config))
if (length(missing_cols) > 0) {
  errors <- c(errors, sprintf("사전 검토 단계에서, config.csv에 필수 컬럼이 없습니다: %s", 
                                paste(missing_cols, collapse = ", ")))
}

# Validation 2: Check file existence
log_message("Checking file existence...", log_list)
for (i in seq_len(nrow(config))) {
  row <- config[i, ]
  sample_name <- row$sample_name
  
  # Check filtered barcode matrix
  if (!is.na(row$dir_input_filtered_barcode_matrix) && row$dir_input_filtered_barcode_matrix != "") {
    if (!dir.exists(row$dir_input_filtered_barcode_matrix)) {
      errors <- c(errors, sprintf("사전 검토 단계에서, 파일이 실제 존재하지 않습니다: sample_name='%s', dir_input_filtered_barcode_matrix='%s'",
                                  sample_name, row$dir_input_filtered_barcode_matrix))
    }
  }
  
  # Check raw barcode matrix (optional but warn if missing)
  if (!is.na(row$dir_input_raw_barcode_matrix) && row$dir_input_raw_barcode_matrix != "") {
    if (!dir.exists(row$dir_input_raw_barcode_matrix)) {
      warnings <- c(warnings, sprintf("사전 검토 단계에서, raw barcode matrix가 없습니다 (경고): sample_name='%s', dir_input_raw_barcode_matrix='%s'",
                                       sample_name, row$dir_input_raw_barcode_matrix))
    }
  }
  
  # Check demultiplex output (required for SNP, optional for HTO)
  if (row$multiplex_method == "SNP") {
    if (is.na(row$dir_demultiplex_output) || row$dir_demultiplex_output == "") {
      errors <- c(errors, sprintf("사전 검토 단계에서, SNP 샘플인데 dir_demultiplex_output이 없습니다: sample_name='%s'", sample_name))
    } else if (!file.exists(row$dir_demultiplex_output)) {
      errors <- c(errors, sprintf("사전 검토 단계에서, 파일이 실제 존재하지 않습니다: sample_name='%s', dir_demultiplex_output='%s'",
                                  sample_name, row$dir_demultiplex_output))
    }
  }
  
  # Check metadata file
  if (!is.na(row$dir_meta_data) && row$dir_meta_data != "") {
    if (!file.exists(row$dir_meta_data)) {
      errors <- c(errors, sprintf("사전 검토 단계에서, 파일이 실제 존재하지 않습니다: sample_name='%s', dir_meta_data='%s'",
                                  sample_name, row$dir_meta_data))
    } else {
      # Check if metadata file can be read and has required columns
      tryCatch({
        if (grepl("\\.csv$", row$dir_meta_data, ignore.case = TRUE)) {
          meta <- read.csv(row$dir_meta_data, nrows = 1, stringsAsFactors = FALSE)
        } else if (grepl("\\.xlsx$", row$dir_meta_data, ignore.case = TRUE)) {
          if (!requireNamespace("readxl", quietly = TRUE)) {
            warnings <- c(warnings, sprintf("readxl 패키지가 없어 xlsx 파일 검증을 건너뜁니다: %s", row$dir_meta_data))
          } else {
            meta <- readxl::read_excel(row$dir_meta_data, n_max = 1)
          }
        }
        
        if (exists("meta")) {
          if (!is.na(row$sample_col_meta_data) && !row$sample_col_meta_data %in% colnames(meta)) {
            errors <- c(errors, sprintf("사전 검토 단계에서, '%s' 변수가 '%s' 오브젝트에 존재하지 않습니다: sample_name='%s', column='%s'",
                                        row$sample_col_meta_data, row$dir_meta_data, sample_name, row$sample_col_meta_data))
          }
          if (!is.na(row$gem_col_meta_data) && !row$gem_col_meta_data %in% colnames(meta)) {
            errors <- c(errors, sprintf("사전 검토 단계에서, '%s' 변수가 '%s' 오브젝트에 존재하지 않습니다: sample_name='%s', column='%s'",
                                        row$gem_col_meta_data, row$dir_meta_data, sample_name, row$gem_col_meta_data))
          }
        }
      }, error = function(e) {
        warnings <- c(warnings, sprintf("메타데이터 파일을 읽을 수 없습니다 (경고): %s, 오류: %s", row$dir_meta_data, e$message))
      })
    }
  }
  
  # Check BAM file (optional, warn if missing)
  if (!is.na(row$dir_input_bam) && row$dir_input_bam != "") {
    if (!file.exists(row$dir_input_bam)) {
      warnings <- c(warnings, sprintf("사전 검토 단계에서, BAM 파일이 없습니다 (경고, 선택사항): sample_name='%s', dir_input_bam='%s'",
                                      sample_name, row$dir_input_bam))
    }
  }
}

# Validation 3: Check GEM name consistency
log_message("Checking GEM name consistency...", log_list)
gem_mapping <- list()
for (i in seq_len(nrow(config))) {
  row <- config[i, ]
  gem_name <- row$gem_name
  filtered_dir <- row$dir_input_filtered_barcode_matrix
  
  if (!gem_name %in% names(gem_mapping)) {
    gem_mapping[[gem_name]] <- filtered_dir
  } else {
    # Same GEM name but different directory - check if this is expected
    if (gem_mapping[[gem_name]] != filtered_dir) {
      # This might be OK if samples are pre-separated (1:1 mapping)
      # But we should warn
      warnings <- c(warnings, sprintf("같은 gem_name('%s')인데 dir_input_filtered_barcode_matrix가 다릅니다: '%s' vs '%s'. 샘플이 미리 분리된 경우 정상일 수 있습니다.",
                                     gem_name, gem_mapping[[gem_name]], filtered_dir))
    }
  }
}

# Validation 4: Check demultiplex method compatibility
log_message("Checking demultiplex method compatibility...", log_list)
for (i in seq_len(nrow(config))) {
  row <- config[i, ]
  if (row$multiplex_method == "SNP") {
    if (is.na(row$demultiplex_method) || row$demultiplex_method == "") {
      errors <- c(errors, sprintf("SNP 샘플인데 demultiplex_method가 지정되지 않았습니다: sample_name='%s'", row$sample_name))
    } else if (row$demultiplex_method != "demuxalot") {
      warnings <- c(warnings, sprintf("SNP 샘플인데 demultiplex_method가 'demuxalot'이 아닙니다: sample_name='%s', method='%s'",
                                     row$sample_name, row$demultiplex_method))
    }
  } else if (row$multiplex_method %in% c("HTO", "CMO")) {
    if (is.na(row$demultiplex_method) || row$demultiplex_method == "") {
      warnings <- c(warnings, sprintf("HTO/CMO 샘플인데 demultiplex_method가 지정되지 않았습니다. 기본값 'HTODemux'를 사용합니다: sample_name='%s'",
                                      row$sample_name))
    }
  }
}

# Print results
if (length(errors) > 0) {
  log_message("=== VALIDATION ERRORS ===", log_list, level = "ERROR")
  for (err in errors) {
    log_message(err, log_list, level = "ERROR")
  }
  close_logging(log_list)
  stop("Validation failed with errors. Please fix the issues above.")
}

if (length(warnings) > 0) {
  log_message("=== VALIDATION WARNINGS ===", log_list, level = "WARNING")
  for (warn in warnings) {
    log_message(warn, log_list, level = "WARNING")
  }
}

log_message("Validation completed successfully!", log_list)
close_logging(log_list)

cat("Validation passed!\n")
if (length(warnings) > 0) {
  cat(sprintf("Note: %d warning(s) were issued. Please review them.\n", length(warnings)))
}


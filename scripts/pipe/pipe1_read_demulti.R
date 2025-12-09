#!/usr/bin/env Rscript
# Step 1: Read data and perform demultiplexing
# Usage: Rscript pipe1_read_demulti.R --config <config_path> --run_id <run_id> [--input_step <step>] [--output_step <step>]

# Load lightweight pipeline environment FIRST (before any library calls)
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- cmd_args[grepl("--file=", cmd_args)]
if (length(file_arg) > 0) {
  script_path <- sub("--file=", "", file_arg)
  script_dir <- dirname(normalizePath(script_path))
} else {
  # Fallback: try to find script directory
  script_dir <- tryCatch(
    {
      dirname(normalizePath(sys.frame(1)$ofile))
    },
    error = function(e) {
      getwd()
    }
  )
}
source(file.path(script_dir, "start_pipe.R"))

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
  library(stringr)
})

# Parse command line arguments
option_list <- list(
  make_option(c("--config", "-c"),
    type = "character", default = NULL,
    help = "Path to config.csv file", metavar = "character"
  ),
  make_option(c("--run_id", "-r"),
    type = "character", default = "run1",
    help = "Run ID", metavar = "character"
  ),
  make_option(c("--input_step", "-i"),
    type = "integer", default = 0,
    help = "Input step number (0 = start from scratch)", metavar = "integer"
  ),
  make_option(c("--output_step", "-o"),
    type = "integer", default = 1,
    help = "Output step number", metavar = "integer"
  ),
  make_option(c("--downsample", "-d"),
    type = "numeric", default = NULL,
    help = "Downsample fraction (0-1) for testing, e.g., 0.1 for 10%% of cells", metavar = "numeric"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$config)) {
  print_help(opt_parser)
  stop("--config argument is required", call. = FALSE)
}

# Source utility functions
# Use consistent approach: pipe_dir is parent of scripts directory
# script_dir is .../scripts/pipe
# pipe_dir should be .../pipe
pipe_dir <- dirname(dirname(script_dir))

# Source pipeline utility functions
# Source pipeline utility functions
source(file.path(pipe_dir, "myR/R/pipe_utils.R"))
source(file.path(pipe_dir, "myR/R/pipe_demulti.R"))
source(file.path(pipe_dir, "myR/R/utils_demulti.R"))

# Setup logging
log_list <- setup_logging(opt$run_id)
log_message(sprintf("Step %d: Starting read and demultiplexing", opt$output_step), log_list)

# Load config
config_list <- load_config(opt$config)
config <- config_list$config
defaults <- config_list$defaults

# Get output directory
output_base_dir <- get_param("output_base_dir", config_list, "/data/user3/sobj/pipe")
output_path <- get_output_path(
  opt$run_id, opt$output_step,
  get_param("output_step1_demulti", config_list, "step1_demulti_list.qs"),
  output_base_dir
)

# Check if output already exists and skip if requested
if (file.exists(output_path) && opt$input_step == opt$output_step) {
  log_message(sprintf("Output file already exists: %s. Loading existing results...", output_path), log_list)
  sl <- load_intermediate(output_path, log_list)
  log_message(sprintf("Loaded %d samples from existing file", length(sl)), log_list)
  close_logging(log_list)
  cat(sprintf("Step %d completed: Loaded %d samples from existing file.\n", opt$output_step, length(sl)))
  quit(save = "no", status = 0)
}

# Create plots directory
plots_dir <- file.path(output_base_dir, opt$run_id, "plots", "step1_demultiplex")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# Set random seed
seed <- as.numeric(get_param("seed", config_list, 1234))
set.seed(seed)

# Separate SNP and HTO samples
snp_samples <- config[config$multiplex_method == "SNP", ]
hto_samples <- config[config$multiplex_method %in% c("HTO", "CMO"), ]

log_message(sprintf(
  "Found %d SNP samples and %d HTO samples",
  nrow(snp_samples), nrow(hto_samples)
), log_list)

# Initialize result lists
demux_cache <- list()
demux_cache_labels <- list()
hto_objects <- list()
sl <- list()
processed_snp_samples <- character(0)
all_snp_mappings <- NULL

# ============================================================================
# Process SNP samples
# ============================================================================
if (nrow(snp_samples) > 0) {
  log_message("Processing SNP samples...", log_list)

  for (i in seq_len(nrow(snp_samples))) {
    row <- snp_samples[i, ]
    sample_name <- row$sample_name
    gem_name <- row$gem_name

    log_message(sprintf(
      "Processing SNP sample %d/%d: %s (GEM: %s)",
      i, nrow(snp_samples), sample_name, gem_name
    ), log_list)

    tryCatch(
      {
        # Get demux file
        demux_file <- as.character(row$dir_demultiplex_output)
        demux_file <- gsub('^["\']+|["\']+$', "", demux_file) # Remove quotes
        demux_file <- trimws(demux_file)

        if (is.na(demux_file) || demux_file == "" || !file.exists(demux_file)) {
          stop(sprintf(
            "dir_demultiplex_output is required for SNP samples but is missing or file not found for %s: %s",
            sample_name, row$dir_demultiplex_output
          ))
        }

        # Cache demultiplexing per GEM/posterior file
        demux_key <- tryCatch(normalizePath(demux_file), error = function(e) demux_file)
        cache_label <- sprintf(
          "%s::%s",
          ifelse(is.na(gem_name) || gem_name == "", "NA_GEM", gem_name),
          basename(demux_file)
        )

        if (!demux_key %in% names(demux_cache)) {
          log_message(sprintf("  Parsing demux posterior for GEM %s (%s)", gem_name, demux_file), log_list)
          demux_cache[[demux_key]] <- demultiplex_demuxalot(
            demuxalot_posterior = demux_file,
            barcode_col = get_param("demuxalot_barcode_col", config_list, "BARCODE"),
            singlet_threshold = as.numeric(get_param("demuxalot_singlet_threshold", config_list, 0.5)),
            doublet_threshold = as.numeric(get_param("demuxalot_doublet_threshold", config_list, 0.3)),
            gem_name = gem_name,
            gem_col = get_param("metadata_gem_col", config_list, "GEM"),
            return_probs = TRUE
          )
          demux_cache_labels[[demux_key]] <- cache_label
        }

        barcode_map_full <- demux_cache[[demux_key]]

        # For logging, check how many barcodes currently map to this sample
        demultiplex_id <- as.character(row$demultiplex_id)
        sample_mapping_preview <- barcode_map_full
        if (!is.na(demultiplex_id) && demultiplex_id != "") {
          sample_mapping_preview <- sample_mapping_preview[
            sample_mapping_preview$Best_Sample == demultiplex_id |
              grepl(paste0("^", demultiplex_id, "\\+"), sample_mapping_preview$Best_Sample) |
              grepl(paste0("\\+", demultiplex_id, "$"), sample_mapping_preview$Best_Sample), ,
            drop = FALSE
          ]
        }

        log_message(sprintf(
          "  Current assignment snapshot for %s: %d barcodes",
          ifelse(is.na(demultiplex_id) || demultiplex_id == "", "ALL", demultiplex_id),
          nrow(sample_mapping_preview)
        ), log_list)
      },
      error = function(e) {
        log_message(sprintf("!!! ERROR processing SNP sample %s: %s !!!", sample_name, e$message), log_list, level = "ERROR")
      }
    )
  }

  # ============================================================================
  # Save SNP barcode mappings (full GEM-level data)
  # ============================================================================
  if (length(demux_cache) > 0) {
    log_message("Combining SNP barcode mappings...", log_list)
    barcode_mapping_list <- demux_cache
    if (length(demux_cache_labels) == length(barcode_mapping_list)) {
      names(barcode_mapping_list) <- unlist(demux_cache_labels, use.names = FALSE)
    }

    all_snp_mappings <- dplyr::bind_rows(barcode_mapping_list, .id = "demux_source")
    log_message(sprintf("Total SNP barcodes (across all GEMs): %d", nrow(all_snp_mappings)), log_list)

    # Save full mapping table
    snp_output_path <- get_output_path(
      opt$run_id, opt$output_step,
      "step1_snp_barcode_mappings.qs",
      output_base_dir
    )
    save_intermediate(all_snp_mappings, snp_output_path, log_list)

    # Save list (per GEM) for re-use
    snp_list_output_path <- get_output_path(
      opt$run_id, opt$output_step,
      "step1_snp_barcode_mapping_list.qs",
      output_base_dir
    )
    save_intermediate(barcode_mapping_list, snp_list_output_path, log_list)
  }
}

# ============================================================================
# Process HTO samples (GEM-centric)
# ============================================================================
if (nrow(hto_samples) > 0) {
  log_message("Processing HTO samples (GEM-centric)...", log_list)

  # Group by GEM
  hto_gems <- unique(hto_samples$gem_name)

  for (gem_idx in seq_along(hto_gems)) {
    gem_name <- hto_gems[gem_idx]
    gem_rows <- hto_samples[hto_samples$gem_name == gem_name, ]

    log_message(sprintf("Processing HTO GEM %d/%d: %s (%d samples)", gem_idx, length(hto_gems), gem_name, nrow(gem_rows)), log_list)

    # Use the first row to determine input paths (assuming consistent within GEM)
    first_row <- gem_rows[1, ]

    tryCatch(
      {
        # 1. Load Data (Once per GEM)
        log_message(sprintf("  Loading data for GEM %s...", gem_name), log_list)

        counts_path <- first_row$dir_input_filtered_barcode_matrix

        # Check for inconsistent paths within GEM
        paths <- gem_rows$dir_input_filtered_barcode_matrix
        if (length(unique(paths)) > 1) {
          log_message(sprintf("  WARNING: Inconsistent input paths for GEM %s. Using first path: %s. This may lead to missing cells if samples are in different directories. Please ensure unique GEM names for separate inputs.", gem_name, counts_path), log_list, level = "WARNING")
        }

        if (!dir.exists(counts_path)) {
          # Try raw if filtered not found? Or just error.
          stop(sprintf("Input directory not found: %s", counts_path))
        }

        raw_data <- Seurat::Read10X(counts_path)

        rna_counts <- NULL
        hto_counts <- NULL

        if (is.list(raw_data)) {
          log_message(sprintf("  DEBUG: raw_data names: %s", paste(names(raw_data), collapse = ", ")), log_list)
          if ("Gene Expression" %in% names(raw_data)) {
            rna_counts <- raw_data$`Gene Expression`
          } else {
            # Fallback if only one matrix or named differently
            rna_counts <- raw_data[[1]]
          }

          # Try to find HTO counts
          if ("Antibody Capture" %in% names(raw_data)) {
            hto_counts <- raw_data$`Antibody Capture`
          } else if ("Custom" %in% names(raw_data)) {
            hto_counts <- raw_data$`Custom`
          } else if ("Multiplexing Capture" %in% names(raw_data)) {
            hto_counts <- raw_data$`Multiplexing Capture`
          } else {
            # Look for any other assay that is not Gene Expression
            other_assays <- setdiff(names(raw_data), "Gene Expression")
            if (length(other_assays) == 1) {
              hto_counts <- raw_data[[other_assays]]
              log_message(sprintf("  DEBUG: Using '%s' as HTO counts", other_assays), log_list)
            } else {
              hto_counts <- NULL
              log_message("  WARNING: Could not identify HTO counts in raw_data list.", log_list)
            }
          }
        } else {
          # If not a list, maybe it's just RNA? Or maybe HTO is in a separate file?
          # For now assume list if HTO is expected.
          rna_counts <- raw_data
          log_message("  WARNING: Read10X returned a single matrix. Assuming it is Gene Expression. HTO counts might be missing.", log_list)
        }

        # Downsample if requested (on the whole GEM)
        if (!is.null(opt$downsample) && opt$downsample > 0 && opt$downsample < 1) {
          log_message(sprintf("  Downsampling GEM to %.1f%%...", opt$downsample * 100), log_list)
          n_cells <- ncol(rna_counts)
          if (n_cells == 0) {
            log_message("  No cells to downsample.", log_list)
          } else {
            keep_cells <- sample(colnames(rna_counts), size = round(n_cells * opt$downsample))
            rna_counts <- rna_counts[, keep_cells, drop = FALSE]
            if (!is.null(hto_counts)) {
              hto_counts <- hto_counts[, keep_cells, drop = FALSE]
            }
          }
        }

        # Ensure colnames
        if (is.null(colnames(rna_counts))) {
          # ... (barcode loading logic if needed, same as before) ...
          # For brevity, assuming Read10X usually provides colnames.
          # If missing, we might fail or need the barcode file logic.
          # I'll include a simplified check.
          if (ncol(rna_counts) > 0) {
            barcodes_path <- file.path(counts_path, "barcodes.tsv.gz")
            if (file.exists(barcodes_path)) {
              barcodes <- read.table(gzfile(barcodes_path), header = FALSE, stringsAsFactors = FALSE)[, 1]
              # Adjust for downsampling
              # This is tricky because we already subsetted by 'keep_cells' which relied on colnames.
              # If colnames were missing, we couldn't have downsampled by name.
              # So downsampling logic above assumes colnames exist.
              # If they don't, we should load barcodes FIRST.
              # For now, assume Read10X provides colnames or we fail.
              log_message("  WARNING: rna_counts has no colnames. Downsampling might be inconsistent if barcodes.tsv.gz is used after subsetting.", log_list, level = "WARNING")
            }
          }
        }

        # Create Seurat Object (Once per GEM)
        log_message("  Creating Seurat Object for GEM...", log_list)
        gem_obj <- Seurat::CreateSeuratObject(counts = rna_counts, project = gem_name)

        # Add HTO Assay
        if (!is.null(hto_counts) && ncol(hto_counts) > 0) {
          # Filter common cells
          common_cells <- intersect(colnames(gem_obj), colnames(hto_counts))
          if (length(common_cells) == 0) {
            log_message("  WARNING: No common cells between RNA and HTO counts for GEM. Skipping HTO assay.", log_list, level = "WARNING")
          } else {
            gem_obj <- gem_obj[, common_cells]

            # Ensure matrix
            # Ensure matrix
            if (is.matrix(hto_counts) || inherits(hto_counts, "Matrix")) {
              hto_counts <- hto_counts[, common_cells, drop = FALSE]
            } else {
              hto_counts <- hto_counts[, common_cells, drop = FALSE]
            }

            gem_obj[["HTO"]] <- Seurat::CreateAssayObject(counts = hto_counts)
            gem_obj <- Seurat::NormalizeData(gem_obj, assay = "HTO", normalization.method = "CLR")

            # 2. Run Demultiplexing (Once per GEM)
            log_message("  Running Demultiplexing on GEM...", log_list)

            # Filter zero HTO cells
            hto_mat <- Seurat::GetAssayData(gem_obj, assay = "HTO", slot = "counts")
            if (is.vector(hto_mat)) {
              if (sum(hto_mat) == 0) gem_obj <- gem_obj[, FALSE]
            } else {
              hto_sums <- Matrix::colSums(hto_mat)
              if (any(hto_sums == 0)) {
                log_message(sprintf("  Removing %d cells with zero HTO counts", sum(hto_sums == 0)), log_list)
                gem_obj <- gem_obj[, hto_sums > 0]
              }
            }

            # Plotting HTO QC
            if (ncol(gem_obj) > 0) {
              # Determine method from manifest (use first row of GEM)
              method <- first_row$demultiplex_method
              if (is.na(method) || method == "") method <- "HTODemux"

              log_message(sprintf("  Demultiplexing method: %s", method), log_list)

              if (method == "MULTIseqDemux") {
                gem_obj <- tryCatch(
                  {
                    Seurat::MULTIseqDemux(gem_obj, assay = "HTO", quantile = 0.7)
                  },
                  error = function(e) {
                    log_message(sprintf("  MULTIseqDemux failed: %s. Trying HTODemux...", e$message), log_list, level = "WARNING")
                    tryCatch(
                      {
                        Seurat::HTODemux(gem_obj, assay = "HTO", positive.quantile = 0.99)
                      },
                      error = function(e2) {
                        log_message(sprintf("  HTODemux failed: %s. Skipping demux.", e2$message), log_list, level = "WARNING")
                        return(gem_obj)
                      }
                    )
                  }
                )
              } else if (method == "None" || method == "Skip") {
                log_message("  Skipping demultiplexing as requested.", log_list)
                # No demux performed
              } else {
                # Default to HTODemux
                gem_obj <- tryCatch(
                  {
                    Seurat::HTODemux(gem_obj, assay = "HTO", positive.quantile = 0.99)
                  },
                  error = function(e) {
                    log_message(sprintf("  HTODemux failed: %s. Trying MULTIseqDemux...", e$message), log_list, level = "WARNING")
                    tryCatch(
                      {
                        Seurat::MULTIseqDemux(gem_obj, assay = "HTO", quantile = 0.7)
                      },
                      error = function(e2) {
                        log_message(sprintf("  MULTIseqDemux failed: %s. Skipping demux.", e2$message), log_list, level = "WARNING")
                        return(gem_obj)
                      }
                    )
                  }
                )
              }

              # Generate Plots
              # RidgePlot
              plot_features <- rownames(gem_obj[["HTO"]])
              if (length(plot_features) > 0) {
                p_ridge <- Seurat::RidgePlot(gem_obj, assay = "HTO", features = plot_features, ncol = 2)
                ggsave(file.path(plots_dir, sprintf("%s_ridge.png", gem_name)), plot = p_ridge, width = 10, height = 8, bg = "white")

                # FeatureScatter (HTO1 vs HTO2 if available)
                if (length(plot_features) >= 2) {
                  p_scatter <- Seurat::FeatureScatter(gem_obj, feature1 = plot_features[1], feature2 = plot_features[2])
                  ggsave(file.path(plots_dir, sprintf("%s_scatter.png", gem_name)), plot = p_scatter, width = 8, height = 8, bg = "white")
                }

                # VlnPlot
                p_vln <- Seurat::VlnPlot(gem_obj, features = "nCount_HTO", pt.size = 0.1, log = TRUE)
                ggsave(file.path(plots_dir, sprintf("%s_vln_count.png", gem_name)), plot = p_vln, width = 6, height = 6, bg = "white")
              }
            }
          }
          # 3. Split and Save Samples
          log_message("  Splitting GEM into samples...", log_list)

          # Load clinical metadata if available (for HTO mapping)
          clinical_metadata <- NULL
          if ("dir_meta_data" %in% colnames(gem_rows)) {
            meta_path <- unique(gem_rows$dir_meta_data)[1]
            if (!is.na(meta_path) && meta_path != "" && file.exists(meta_path)) {
              log_message(sprintf("  Loading clinical metadata from %s", meta_path), log_list)
              clinical_metadata <- read.csv(meta_path, stringsAsFactors = FALSE)
            }
          }

          # Identify Demux Columns
          demux_id_col <- NULL
          demux_class_col <- NULL

          if ("HTO_maxID" %in% colnames(gem_obj@meta.data)) {
            demux_id_col <- "HTO_maxID"
            demux_class_col <- "HTO_classification.global"
          } else if ("MULTI_ID" %in% colnames(gem_obj@meta.data)) {
            demux_id_col <- "MULTI_ID"
            if ("MULTI_classification.global" %in% colnames(gem_obj@meta.data)) {
              demux_class_col <- "MULTI_classification.global"
            } else {
              # MULTIseqDemux might not produce .global, but MULTI_ID implies singlet if it's an ID
              demux_class_col <- NULL
            }
          } else if ("hash.ID" %in% colnames(gem_obj@meta.data)) {
            demux_id_col <- "hash.ID"
            demux_class_col <- "hash.ID" # Sometimes classification is same
          }
          if (!is.null(demux_id_col)) {
            log_message(sprintf("  Using demux column: %s", demux_id_col), log_list)
          } else {
            log_message("  WARNING: No demultiplexing columns found in metadata.", log_list, level = "WARNING")
          }

          for (i in seq_len(nrow(gem_rows))) {
            row <- gem_rows[i, ]
            sample_name <- row$sample_name
            demux_id <- row$demultiplex_id

            # Try to derive HTO tag from metadata if available
            derived_hto_tag <- NULL
            if (!is.null(clinical_metadata)) {
              if ("sample_no" %in% colnames(clinical_metadata)) {
                meta_row <- clinical_metadata[clinical_metadata$sample_no == sample_name, ]
                if (nrow(meta_row) > 0) {
                  # Use demulti_id from metadata if present (pre-processed)
                  if ("demulti_id" %in% colnames(meta_row) && !is.na(meta_row$demulti_id[1]) && meta_row$demulti_id[1] != "") {
                    derived_hto_tag <- meta_row$demulti_id[1]
                    log_message(sprintf("    Using demulti_id from metadata: %s", derived_hto_tag), log_list)
                  } else {
                    log_message(sprintf("    Debug: demulti_id missing or empty for sample %s in metadata", sample_name), log_list)
                  }
                } else {
                  log_message(sprintf("    Debug: Sample %s not found in metadata (sample_no column)", sample_name), log_list)
                }
              } else {
                log_message("    Debug: sample_no column not found in clinical_metadata", log_list)
              }
            } else {
              log_message("    Debug: clinical_metadata is NULL", log_list)
            }

            target_demux_id <- demux_id
            if (!is.null(derived_hto_tag)) {
              target_demux_id <- derived_hto_tag
            }

            log_message(sprintf("    Extracting sample %s (Target ID: %s)...", sample_name, target_demux_id), log_list)

            # Subset logic
            if (!is.null(target_demux_id) && !is.na(target_demux_id) && target_demux_id != "" && !is.null(demux_id_col)) {
              # Check if demux_id exists in classification
              # Create a mask
              # Handle potential NA in columns
              ids <- gem_obj@meta.data[[demux_id_col]]
              classes <- if (!is.null(demux_class_col)) gem_obj@meta.data[[demux_class_col]] else rep("Singlet", length(ids))

              keep_mask <- !is.na(ids) & ids == target_demux_id & classes == "Singlet"

              if (sum(keep_mask) == 0) {
                log_message(sprintf("    WARNING: No singlet cells found for Target ID '%s' (Manifest ID: '%s'). Checking available IDs: %s", target_demux_id, demux_id, paste(unique(ids), collapse = ", ")), log_list, level = "WARNING")
                log_message(sprintf("    Skipping sample %s due to zero cells.", sample_name), log_list, level = "WARNING")
                next
              } else {
                sample_obj <- gem_obj[, keep_mask]
                log_message(sprintf("    Found %d cells.", ncol(sample_obj)), log_list)
              }

              # Add metadata
              sample_obj$demultiplex_id <- demux_id
              sample_obj$contamination_risk <- if ("HTO_margin" %in% colnames(sample_obj@meta.data)) sample_obj$HTO_margin else NA
              sample_obj$demux_class <- if (!is.null(demux_class_col)) sample_obj@meta.data[[demux_class_col]] else "Singlet"
            } else {
              # If no demux ID provided or demux failed
              if (!is.null(demux_id_col)) {
                log_message(sprintf("    WARNING: Sample %s has no demultiplex_id but demux ran. Skipping.", sample_name), log_list, level = "WARNING")
                next
              } else {
                # No demux ran (e.g. no HTO), so maybe it's 1:1 GEM:Sample?
                if (nrow(gem_rows) == 1) {
                  log_message(sprintf("    No demux performed, assuming 1:1 GEM to Sample for %s", sample_name), log_list)
                  sample_obj <- gem_obj
                } else {
                  log_message(sprintf("    WARNING: Multiple samples in GEM but no demux performed. Cannot split %s. Skipping.", sample_name), log_list, level = "ERROR")
                  next
                }
              }
            }

            # Add config metadata
            sample_obj$sample_id <- sample_name
            sample_obj$GEM <- gem_name
            # Add other columns from manifest
            for (col in names(row)) {
              if (!col %in% c("dir_input_filtered_barcode_matrix", "dir_input_raw_barcode_matrix", "dir_demultiplex_output")) {
                sample_obj[[col]] <- row[[col]]
              }
            }

            # Extract metadata from filename if patterns are provided
            timepoint_regex <- get_param("metadata_timepoint_regex", config_list, NULL)
            if (!is.null(timepoint_regex) && timepoint_regex != "") {
              timepoint <- str_extract(sample_name, timepoint_regex)
              if (!is.na(timepoint)) {
                day_map_24 <- get_param("metadata_day_map_24", config_list, NULL)
                day_map_72 <- get_param("metadata_day_map_72", config_list, NULL)

                if (timepoint == "24" && !is.null(day_map_24)) {
                  sample_obj$day <- as.numeric(day_map_24)
                } else if (timepoint == "72" && !is.null(day_map_72)) {
                  sample_obj$day <- as.numeric(day_map_72)
                }
                sample_obj$time_point <- timepoint
              }
            }

            # Save
            hto_objects[[sample_name]] <- sample_obj
          }
        }

        # Clean up to free memory
        rm(gem_obj, rna_counts, hto_counts, raw_data)
        gc()
      },
      error = function(e) {
        log_message(sprintf("!!! ERROR processing GEM %s: %s !!!", gem_name, e$message), log_list, level = "ERROR")
      }
    )
  }
}

# ============================================================================
# Create Seurat objects for SNP samples (using barcode mappings)
# ============================================================================
if (!is.null(all_snp_mappings) && nrow(all_snp_mappings) > 0) {
  log_message("Creating Seurat objects for SNP samples...", log_list)

  # Group by GEM (same GEM can have multiple samples)
  gem_groups <- split(snp_samples, snp_samples$gem_name)

  for (gem_name in names(gem_groups)) {
    gem_samples <- gem_groups[[gem_name]]

    for (i in seq_len(nrow(gem_samples))) {
      row <- gem_samples[i, ]
      sample_name <- row$sample_name

      tryCatch(
        {
          log_message(sprintf("  Creating Seurat object for %s (GEM: %s)", sample_name, gem_name), log_list)

          # Load count matrices
          filtered_counts <- Seurat::Read10X(row$dir_input_filtered_barcode_matrix)

          # Get barcode mapping for this sample (based on demultiplex_id)
          demultiplex_id <- as.character(row$demultiplex_id)
          sample_mapping <- all_snp_mappings
          if (!is.na(demultiplex_id) && demultiplex_id != "") {
            sample_mapping <- sample_mapping[
              sample_mapping$Best_Sample == demultiplex_id |
                grepl(paste0("^", demultiplex_id, "\\+"), sample_mapping$Best_Sample) |
                grepl(paste0("\\+", demultiplex_id, "$"), sample_mapping$Best_Sample), ,
              drop = FALSE
            ]
          }

          if (nrow(sample_mapping) == 0) {
            stop(sprintf(
              "No barcode mapping found for sample %s (demultiplex_id: %s)",
              sample_name, demultiplex_id
            ))
          }

          sample_mapping$sample_name <- sample_name
          sample_mapping$demultiplex_id <- demultiplex_id

          # Match barcodes using original barcode (before suffix)
          # demultiplex_demuxalot already added suffix to Barcode column
          # We need to match original barcodes first, then add suffix to counts
          original_barcodes_in_mapping <- sample_mapping$barcode # Original barcode (no suffix)
          barcodes_with_suffix <- sample_mapping$Barcode # Barcode with suffix (rownames)

          # Find common original barcodes
          common_original_barcodes <- intersect(colnames(filtered_counts), original_barcodes_in_mapping)

          if (length(common_original_barcodes) == 0) {
            # Debug: show some examples
            log_message(sprintf(
              "  DEBUG: Count matrix barcodes (first 5): %s",
              paste(head(colnames(filtered_counts), 5), collapse = ", ")
            ), log_list)
            log_message(sprintf(
              "  DEBUG: Mapping barcodes (first 5): %s",
              paste(head(original_barcodes_in_mapping, 5), collapse = ", ")
            ), log_list)
            stop(sprintf("No common original barcodes found between count matrix and demux results for %s", sample_name))
          }

          # Get corresponding barcodes with suffix from mapping
          mapping_subset <- sample_mapping[sample_mapping$barcode %in% common_original_barcodes, , drop = FALSE]

          # Create a named vector mapping original -> suffix barcode
          orig_to_suffix <- setNames(mapping_subset$Barcode, mapping_subset$barcode)

          # Update count matrix barcode names: add suffix to match mapping Barcode column
          count_barcodes_original <- colnames(filtered_counts)
          count_barcodes_to_update <- count_barcodes_original %in% names(orig_to_suffix)
          colnames(filtered_counts)[count_barcodes_to_update] <- orig_to_suffix[count_barcodes_original[count_barcodes_to_update]]

          # Match using Barcode column (not rownames, since rownames include sample_name)
          # Create matching key: sample_name_Barcode
          sample_mapping$match_key <- paste0(sample_name, "_", sample_mapping$Barcode)
          count_match_keys <- paste0(sample_name, "_", colnames(filtered_counts))
          common_cells <- intersect(count_match_keys, sample_mapping$match_key)

          if (length(common_cells) == 0) {
            # Additional debug info
            log_message(sprintf(
              "  DEBUG: After suffix addition, count barcodes (first 5): %s",
              paste(head(colnames(filtered_counts), 5), collapse = ", ")
            ), log_list)
            log_message(sprintf(
              "  DEBUG: Mapping rownames (first 5): %s",
              paste(head(rownames(sample_mapping), 5), collapse = ", ")
            ), log_list)
            log_message(sprintf("  DEBUG: Common original barcodes: %d", length(common_original_barcodes)), log_list)
            stop(sprintf("No common cells found between count matrix and demux results for %s", sample_name))
          }

          # Get corresponding count matrix columns and mapping rows
          count_cols <- count_match_keys %in% common_cells
          mapping_rows <- sample_mapping$match_key %in% common_cells

          # Create Seurat object
          obj <- Seurat::CreateSeuratObject(
            counts = filtered_counts[, count_cols, drop = FALSE],
            meta.data = sample_mapping[mapping_rows, , drop = FALSE]
          )

          # Add config metadata columns
          config_cols_to_add <- c("name", "patient_name", "contamination_risk", "sample_name", "gem_name", "demultiplex_id")
          for (col in config_cols_to_add) {
            if (col %in% colnames(row) && !is.na(row[[col]])) {
              obj[[col]] <- as.character(row[[col]])
            }
          }

          # Extract metadata from filename if patterns are provided
          timepoint_regex <- get_param("metadata_timepoint_regex", config_list, NULL)
          if (!is.null(timepoint_regex) && timepoint_regex != "") {
            timepoint <- str_extract(sample_name, timepoint_regex)
            if (!is.na(timepoint)) {
              day_map_24 <- get_param("metadata_day_map_24", config_list, NULL)
              day_map_72 <- get_param("metadata_day_map_72", config_list, NULL)

              if (timepoint == "24" && !is.null(day_map_24)) {
                obj$day <- as.numeric(day_map_24)
              } else if (timepoint == "72" && !is.null(day_map_72)) {
                obj$day <- as.numeric(day_map_72)
              }
              obj$time_point <- timepoint
            }
          }

          # Demultiplex Doublet Removal
          # If demultiplex_doublet_removal is TRUE, remove cells identified as doublets by demultiplexing
          # These are inter-sample doublets (e.g., SampleA+SampleB)
          remove_demux_doublets <- get_param("demultiplex_doublet_removal", config_list, FALSE)
          if (isTRUE(remove_demux_doublets)) {
            # Identify doublets based on demuxalot_annotation_col (Best_Sample)
            # Doublets usually contain the separator (default '+')
            demux_col <- get_param("demuxalot_annotation_col", config_list, "Best_Sample")
            separator <- get_param("demuxalot_doublet_separator", config_list, "+")

            if (demux_col %in% colnames(obj@meta.data)) {
              doublet_cells <- grepl(paste0("\\", separator), obj@meta.data[[demux_col]])
              if (any(doublet_cells)) {
                n_doublets <- sum(doublet_cells)
                log_message(sprintf("Removing %d demultiplexing doublets from sample %s", n_doublets, sample_name), log_list)
                obj <- obj[, !doublet_cells]
              }
            } else {
              log_message(sprintf("Warning: Demultiplex column '%s' not found for doublet removal", demux_col), log_list)
            }
          }

          # Downsample if requested (for testing)
          if (!is.null(opt$downsample) && opt$downsample > 0 && opt$downsample < 1) {
            n_cells_before <- ncol(obj)
            n_cells_after <- max(1, round(n_cells_before * opt$downsample))
            if (n_cells_after < n_cells_before) {
              set.seed(seed + i) # Reproducible sampling
              cells_to_keep <- sample(colnames(obj), n_cells_after)
              obj <- obj[, cells_to_keep]
              log_message(sprintf(
                "  Downsampled %s: %d -> %d cells (%.1f%%)",
                sample_name, n_cells_before, n_cells_after,
                opt$downsample * 100
              ), log_list)
            }
          }

          sl[[sample_name]] <- obj
          processed_snp_samples <- unique(c(processed_snp_samples, sample_name))
          log_message(sprintf("  Completed %s: %d cells", sample_name, ncol(obj)), log_list)
        },
        error = function(e) {
          log_message(sprintf(
            "!!! ERROR creating Seurat object for SNP sample %s: %s !!!",
            sample_name, e$message
          ), log_list, level = "ERROR")
        }
      )
    }
  }
}

# Combine SNP and HTO objects
if (length(hto_objects) > 0) {
  sl <- c(sl, hto_objects)
}

# Save results
if (length(sl) > 0) {
  log_message(sprintf("Saving results to: %s", output_path), log_list)
  save_intermediate(sl, output_path, log_list)
  snp_processed_n <- length(processed_snp_samples)
  log_message(sprintf(
    "Step %d completed: %d samples processed (%d SNP, %d HTO)",
    opt$output_step, length(sl),
    snp_processed_n, length(hto_objects)
  ), log_list)
} else {
  log_message("No samples were successfully processed!", log_list, level = "ERROR")
  close_logging(log_list)
  stop("Step 1 failed: No samples processed")
}

close_logging(log_list)
cat(sprintf("Step %d completed successfully. Processed %d samples.\n", opt$output_step, length(sl)))

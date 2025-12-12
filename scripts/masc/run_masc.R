#!/usr/bin/env Rscript
# MASC CLI runner (optparse)
#
# - Seurat .qs 입력을 받아 MASC(클러스터 abundance 차이) 분석 수행
# - 변수명/옵션은 모두 CLI 플래그로 입력
# - BMI 계산(ht, wt -> bmi)은 하드코딩(있으면 생성)
#
# Usage example:
#   Rscript scripts/masc/run_masc.R \
#     -i /data/user3/sobj/is2_IS_3_clustered.qs \
#     -o /data/user3/sobj/masc/stroke_complex_cli \
#     --cluster_var anno3 \
#     --contrast_var g3 \
#     --random_effects hos_no \
#     --fixed_effects GEM,SET,age,sex,bmi,hx_smok,hx_alcohol \
#     --prefix masc_anno3_complex

suppressPackageStartupMessages({
  library(optparse)
  library(qs)
  library(Seurat)
})

.split_csv <- function(x) {
  if (is.null(x) || is.na(x) || x == "") return(NULL)
  parts <- unlist(strsplit(x, ",", fixed = TRUE))
  parts <- trimws(parts)
  parts <- parts[nzchar(parts)]
  if (!length(parts)) return(NULL)
  parts
}

.calc_bmi <- function(ht, wt) {
  ht_m <- ht / 100
  wt / (ht_m^2)
}

option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "Input Seurat .qs path"),
  make_option(c("-o", "--output_dir"), type = "character", help = "Output directory"),
  make_option(c("--prefix"), type = "character", default = "masc", help = "Output prefix"),
  make_option(c("--cluster_var"), type = "character", help = "Cluster column name in meta.data"),
  make_option(c("--contrast_var"), type = "character", help = "Contrast variable (must have >=2 levels)"),
  make_option(c("--random_effects"), type = "character", default = "", help = "Random effects (comma-separated). At least one is recommended; e.g. hos_no"),
  make_option(c("--fixed_effects"), type = "character", default = "", help = "Fixed effects (comma-separated); e.g. GEM,SET,age,sex,bmi,hx_smok,hx_alcohol"),
  make_option(c("--min_cluster_cells"), type = "integer", default = 10L, help = "Minimum total cells per cluster to keep"),
  make_option(c("--force_run"), action = "store_true", default = FALSE, help = "Force recompute (ignore cached qs)"),
  make_option(c("--no_plot"), action = "store_true", default = FALSE, help = "Disable plotting"),
  make_option(c("--save_models"), action = "store_true", default = FALSE, help = "Save fitted model objects (qs)"),
  make_option(c("--no_fdr"), action = "store_true", default = FALSE, help = "Disable FDR adjustment"),
  make_option(c("--seed"), type = "integer", default = 123, help = "Seed for any sampling steps"),
  make_option(c("--max_cells"), type = "integer", default = 0L, help = "Optional global downsample (0=off). If >0 and < nCells, randomly sample that many cells."),
  make_option(c("--masc_r_path"), type = "character",
              default = "/home/user3/data_user3/git_repo/mylit/Git_Repo/_wt/masc/myR/R/masc.R",
              help = "Path to masc.R (defaults to current worktree)"),
  make_option(c("--start_r_path"), type = "character",
              default = "/home/user3/GJC_KDW_250721/start.R",
              help = "Optional start.R path (loads env/renv). If missing, ignored.")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

if (is.null(opt$input) || is.null(opt$output_dir) || is.null(opt$cluster_var) || is.null(opt$contrast_var)) {
  print_help(parser)
  stop("Missing required args: --input, --output_dir, --cluster_var, --contrast_var", call. = FALSE)
}

if (!file.exists(opt$input)) stop("Input not found: ", opt$input, call. = FALSE)
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

if (file.exists(opt$start_r_path)) {
  suppressWarnings(try(source(opt$start_r_path), silent = TRUE))
}
if (!file.exists(opt$masc_r_path)) stop("masc.R not found: ", opt$masc_r_path, call. = FALSE)
source(opt$masc_r_path)

set.seed(opt$seed)

cat(sprintf("Loading Seurat object: %s\n", opt$input))
sobj <- qs::qread(opt$input)

if (!inherits(sobj, "Seurat")) stop("Input is not a Seurat object.", call. = FALSE)

meta <- sobj@meta.data

# BMI (hardcoded)
if (!"bmi" %in% colnames(meta) && all(c("ht", "wt") %in% colnames(meta))) {
  meta$bmi <- .calc_bmi(meta$ht, meta$wt)
  sobj@meta.data <- meta
  cat("Added bmi from ht/wt\n")
}

# Ensure hos_no and other ID-like columns are not numeric if chosen
random_effects <- .split_csv(opt$random_effects)
fixed_effects <- .split_csv(opt$fixed_effects)

if (!is.null(random_effects)) {
  for (v in random_effects) {
    if (v %in% colnames(sobj@meta.data)) {
      sobj@meta.data[[v]] <- as.character(sobj@meta.data[[v]])
    }
  }
}

# Ensure contrast is factor
if (!opt$contrast_var %in% colnames(sobj@meta.data)) stop("contrast_var missing: ", opt$contrast_var, call. = FALSE)
sobj@meta.data[[opt$contrast_var]] <- as.factor(sobj@meta.data[[opt$contrast_var]])

# Basic global downsample if requested
if (!is.null(opt$max_cells) && opt$max_cells > 0 && opt$max_cells < ncol(sobj)) {
  keep <- sample(colnames(sobj), opt$max_cells)
  sobj <- subset(sobj, cells = keep)
  cat(sprintf("Downsampled to %d cells\n", opt$max_cells))
}

# Filter clusters by total cells
if (!opt$cluster_var %in% colnames(sobj@meta.data)) stop("cluster_var missing: ", opt$cluster_var, call. = FALSE)
cluster_counts <- table(sobj@meta.data[[opt$cluster_var]])
keep_clusters <- names(cluster_counts[cluster_counts >= opt$min_cluster_cells])
if (length(keep_clusters) < 2) stop("Not enough clusters after filtering.", call. = FALSE)
sobj <- subset(sobj, cells = colnames(sobj)[sobj@meta.data[[opt$cluster_var]] %in% keep_clusters])

# Keep only available covariates
if (!is.null(fixed_effects)) fixed_effects <- fixed_effects[fixed_effects %in% colnames(sobj@meta.data)]
if (!is.null(random_effects)) random_effects <- random_effects[random_effects %in% colnames(sobj@meta.data)]

cat(sprintf("cluster_var=%s, contrast_var=%s\n", opt$cluster_var, opt$contrast_var))
cat(sprintf("random_effects=%s\n", paste(random_effects, collapse = ",")))
cat(sprintf("fixed_effects=%s\n", paste(fixed_effects, collapse = ",")))

res <- run_masc_pipeline(
  seurat_obj = sobj,
  cluster_var = opt$cluster_var,
  contrast_var = opt$contrast_var,
  random_effects = random_effects,
  fixed_effects = fixed_effects,
  save = TRUE,
  output_dir = opt$output_dir,
  prefix = opt$prefix,
  force_run = opt$force_run,
  plotting = !opt$no_plot,
  save_models = opt$save_models,
  adjust_pvalue = !opt$no_fdr,
  verbose = TRUE
)

cat("\n=== MASC done ===\n")
print(res$masc_results)



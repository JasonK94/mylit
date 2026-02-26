# Real Data Test Script for DEG Consensus (Stoke PBMC)
# Dataset: Downsampled IS6 (is5s)

suppressPackageStartupMessages({
    library(Seurat)
    library(qs)
    library(dplyr)
    if (dir.exists("myR")) {
        devtools::load_all("myR")
    }
})

# Explicitly source if not found
if (!exists("run_deg_consensus")) {
    message("run_deg_consensus not found in package, sourcing manually...")
    source("myR/R/deg_consensus/run_deg_consensus.R")
    source("myR/R/deg_consensus/deg_methods_base.R")
    source("myR/R/deg_consensus/deg_methods_edger.R")
    source("myR/R/deg_consensus/deg_methods_deseq2.R")
    source("myR/R/deg_consensus/deg_methods_limma.R")
    source("myR/R/deg_consensus/deg_consensus_analysis.R")
}

# 1. Load Data
qs_path <- "/data/user3/sobj/IS6_sex_added_0.1x_251110.qs"
message("Loading data from: ", qs_path)
sobj <- qs::qread(qs_path)

# 2. Preprocessing
message("Preprocessing metadata...")
meta <- sobj@meta.data

# Calculate BMI: wt (kg) / (ht (cm) / 100)^2
if (all(c("wt", "ht") %in% colnames(meta))) {
    meta$bmi <- meta$wt / (meta$ht / 100)^2
    message("BMI calculated. Mean BMI: ", round(mean(meta$bmi, na.rm = TRUE), 2))
} else {
    warning("wt or ht columns missing. BMI not calculated.")
}

# Ensure factors and formatting
meta$g3 <- factor(meta$g3)
meta$sex <- factor(meta$sex)
meta$GEM <- factor(meta$GEM)
meta$set <- factor(meta$set)

# Update Seurat object
sobj@meta.data <- meta

# 3. Define Analysis Parameters
contrast_str <- "2 - 1"
group_col <- "g3"
cluster_col <- "anno3.scvi"
sample_col <- "hos_no"
# Covariates to adjust
# Decision: GEM (8 levels) consumes too many degrees of freedom (7 DF) given N=23.
# We use 'set' (2 levels, 1 DF) which is nested in GEM but sufficient for testing.
covariates <- c("sex", "age", "bmi", "set")

methods_to_test <- c(
    "edgeR-LRT",
    "limma-voom",
    "nebula" # Temporarily disabled to speed up debugging of basic pipeline
)

# 4. Run Consensus Analysis
message("\n=== Running DEG Consensus Analysis ===")
results <- run_deg_consensus(
    sobj = sobj,
    contrast = contrast_str,
    methods = methods_to_test,
    cluster_id = cluster_col,
    sample_id = sample_col,
    group_id = group_col,
    covar_effects = covariates,
    n_cores = 4,
    verbose = TRUE
)

# 5. Save Results
out_file <- "stats/deg_consensus_test_real_is5s.qs"
dir.create("stats", showWarnings = FALSE)
qs::qsave(results, out_file)
message("Results saved to: ", out_file)

# 6. Basic Validation
message("\n=== Validation ===")
print(names(results$results))

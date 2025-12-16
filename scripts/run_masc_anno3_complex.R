#!/usr/bin/env Rscript
# MASC Analysis for Stroke Data (Complex Model, anno3)
# Tests cluster abundance differences by g3 (target variable) with full covariates

suppressPackageStartupMessages({
    library(Seurat)
    library(qs)
    library(dplyr)
})

# Source MASC functions
masc_r_path <- "/home/user3/data_user3/git_repo/mylit/Git_Repo/_wt/masc/myR/R/masc.R"
if (!file.exists(masc_r_path)) stop("MASC.R not found")
source(masc_r_path)

# Configuration
DATA_PATH <- "/data/user3/sobj/is2_IS_3_clustered.qs"
OUTPUT_DIR <- "/data/user3/sobj/masc/stroke_complex"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("=== MASC Analysis: anno3 with Complex Model ===\n\n")

# Helper function to calculate BMI
calculate_bmi <- function(ht, wt) {
    ht_m <- ht / 100
    bmi <- wt / (ht_m^2)
    bmi
}

# Load Data
cat(sprintf("Loading data from: %s\n", DATA_PATH))
seurat_obj <- qs::qread(DATA_PATH)

# --- Data Preprocessing ---

# 1. Cluster Variable: anno3
cluster_var <- "anno3"
if (!cluster_var %in% colnames(seurat_obj@meta.data)) {
    stop(sprintf("Cluster variable '%s' not found.\n", cluster_var))
}

# 2. Target Variable: g3
if (!"g3" %in% colnames(seurat_obj@meta.data)) stop("Target variable 'g3' not found")
seurat_obj@meta.data$g3 <- as.factor(seurat_obj@meta.data$g3)

# 3. Random Effect: hos_no
seurat_obj@meta.data$hos_no <- as.character(seurat_obj@meta.data$hos_no)

# 4. Covariates (Fixed Effects)
# Calculate BMI (kept, but BMI is omitted from the model by default due to missingness/type artifacts)
if ("ht" %in% colnames(seurat_obj@meta.data) && "wt" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj@meta.data$bmi <- calculate_bmi(seurat_obj@meta.data$ht, seurat_obj@meta.data$wt)
    cat("Calculated BMI from ht and wt\n")
}

# Clean hx_alcohol
if ("hx_alcohol" %in% colnames(seurat_obj@meta.data)) {
    cat("Cleaning hx_alcohol variable...\n")
    alcohol_raw <- seurat_obj@meta.data$hx_alcohol
    alcohol_clean <- ifelse(is.na(alcohol_raw), "Missing",
                   ifelse(grepl("None|^NA$", alcohol_raw, ignore.case = TRUE), "None",
                   ifelse(grepl("Light|Infrequent", alcohol_raw, ignore.case = TRUE), "Light",
                   ifelse(grepl("Moderate", alcohol_raw, ignore.case = TRUE), "Moderate",
                   ifelse(grepl("Heavy|Excessive", alcohol_raw, ignore.case = TRUE), "Heavy",
                   as.character(alcohol_raw))))))
    seurat_obj@meta.data$hx_alcohol <- as.factor(alcohol_clean)
}

# Clean other Yes/No variables
yes_no_vars <- c("hx_smok")
for (var in yes_no_vars) {
    if (var %in% colnames(seurat_obj@meta.data)) {
        values <- seurat_obj@meta.data[[var]]
        values_clean <- ifelse(is.na(values), "Missing",
                      ifelse(grepl("^[Yy]", values), "Yes",
                      ifelse(grepl("^[Nn]", values), "No", values)))
        seurat_obj@meta.data[[var]] <- as.factor(values_clean)
    }
}

# Model choice note (nested structure):
# - hos_no is fully contained within GEM, and GEM is fully contained within SET.
# - Avoid including both GEM and SET together (collinearity/rank deficiency).
# - BMI often contains NULL/character artifacts; omit by default.
random_effects <- "hos_no"

# Reduced/stable model (recommended)
fixed_effects <- c("age", "sex")
if ("GEM" %in% colnames(seurat_obj@meta.data)) {
    fixed_effects <- c(fixed_effects, "GEM")
} else if ("SET" %in% colnames(seurat_obj@meta.data)) {
    fixed_effects <- c(fixed_effects, "SET")
}

# Filter to only available columns (defensive)
fixed_effects <- fixed_effects[fixed_effects %in% colnames(seurat_obj@meta.data)]

cat(sprintf("Fixed effects: %s\n", paste(fixed_effects, collapse = ", ")))
cat(sprintf("Random effects: %s\n", paste(random_effects, collapse = ", ")))

# Filter clusters
cluster_counts <- table(seurat_obj@meta.data[[cluster_var]], seurat_obj@meta.data$g3)
min_cells <- 10
valid_clusters <- rownames(cluster_counts)[rowSums(cluster_counts) >= min_cells]
cat(sprintf("Using %d clusters (>= %d cells)\n", length(valid_clusters), min_cells))

# Subset
seurat_obj_subset <- subset(seurat_obj, cells = colnames(seurat_obj)[seurat_obj@meta.data[[cluster_var]] %in% valid_clusters])

# Run Pipeline
cat("\nRunning MASC pipeline...\n")
results <- run_masc_pipeline(
    seurat_obj = seurat_obj_subset,
    cluster_var = cluster_var,
    contrast_var = "g3",
    random_effects = random_effects,
    fixed_effects = fixed_effects,
    save = TRUE,
    output_dir = OUTPUT_DIR,
    prefix = "masc_anno3_complex",
    force_run = TRUE,
    plotting = TRUE,
    save_models = FALSE,
    adjust_pvalue = TRUE,
    verbose = TRUE
)

cat("\n=== MASC Results (Complex Model) ===\n")
print(results$masc_results)

# Save summary
summary_path <- file.path(OUTPUT_DIR, "masc_anno3_complex_summary.txt")
sink(summary_path)
cat(sprintf("Model: %s ~ g3 + %s + (1|%s)\n\n", cluster_var, paste(fixed_effects, collapse=" + "), random_effects))
print(results$masc_results)
sink()
cat(sprintf("\nResults saved to: %s\n", summary_path))


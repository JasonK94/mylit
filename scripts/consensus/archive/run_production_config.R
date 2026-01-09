# Production Analysis Script for Stroke PBMC (IS5)
# Levels: Global, Broad (anno3big), Fine (anno3 within anno3big)

# 0. Setup Environment
renv_lib <- "/home/user3/GJC_KDW_250721/renv/library/R-4.3/x86_64-pc-linux-gnu"
if (dir.exists(renv_lib)) .libPaths(c(renv_lib, .libPaths()))

suppressPackageStartupMessages({
    library(Seurat)
    library(qs)
    library(dplyr)
    if (dir.exists("myR")) {
        devtools::load_all("myR")
    }
})

# Explicitly source if not found (Safety)
if (!exists("run_deg_consensus")) {
    message("Sourcing myR manually...")
    source("myR/R/deg_consensus/run_deg_consensus.R")
    source("myR/R/deg_consensus/deg_methods_base.R")
    source("myR/R/deg_consensus/deg_methods_edger.R")
    source("myR/R/deg_consensus/deg_methods_deseq2.R")
    source("myR/R/deg_consensus/deg_methods_limma.R")
    source("myR/R/deg_consensus/deg_consensus_analysis.R")
}

# 1. Load Data
qs_path <- "/data/user3/sobj/IS6_sex_added_251110_bu.qs"
message("Loading data: ", qs_path)
sobj_full <- qs::qread(qs_path)

# 2. Preprocessing
meta <- sobj_full@meta.data
if (all(c("wt", "ht") %in% colnames(meta))) {
    meta$bmi <- meta$wt / (meta$ht / 100)^2
}
meta$g3 <- factor(meta$g3)
meta$sex <- factor(meta$sex)
meta$set <- factor(meta$set)
meta$global_group <- "All_Cells"
sobj_full@meta.data <- meta

# Output Base Dir
base_out_dir <- "/data/user3/sobj/consensus"
dir.create(base_out_dir, recursive = TRUE, showWarnings = FALSE)

# Parameters
contrast_str <- "2 - 1"
group_col <- "g3"
sample_col <- "hos_no"
covariates <- c("sex", "age", "bmi", "set")
methods_to_test <- c("edgeR-LRT", "limma-voom", "nebula")
cores <- 6 # Adjust based on availability

# ==============================================================================
# Helper Function
# ==============================================================================
run_analysis <- function(obj, label, cluster_col, out_subdir) {
    out_dir <- file.path(base_out_dir, out_subdir)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    out_file <- file.path(out_dir, paste0("results_", label, ".qs"))

    message(sprintf("\n>>> Starting Analysis: %s (Cluster: %s)", label, cluster_col))

    if (file.exists(out_file)) {
        message("File exists, skipping: ", out_file)
        return(NULL)
    }

    tryCatch(
        {
            res <- run_deg_consensus(
                sobj = obj,
                contrast = contrast_str,
                methods = methods_to_test,
                cluster_id = cluster_col,
                sample_id = sample_col,
                group_id = group_col,
                covar_effects = covariates,
                n_cores = cores,
                verbose = TRUE
            )
            qs::qsave(res, out_file)
            message("Saved: ", out_file)

            # Optional: Generate Consensus Table CSV immediately
            if (!is.null(res$consensus_table)) {
                write.csv(res$consensus_table, file.path(out_dir, paste0("consensus_", label, ".csv")))
            }
        },
        error = function(e) {
            message("Error in ", label, ": ", e$message)
        }
    )
}

# ==============================================================================
# L1. Global Analysis
# ==============================================================================
run_analysis(sobj_full, "Global", "global_group", "Level1_Global")

# ==============================================================================
# L2. Broad Analysis (anno3big)
# ==============================================================================
# "Tc, Bc, DC, Mono, NKc, PLT, others"
run_analysis(sobj_full, "Broad", "anno3big", "Level2_Broad")

# ==============================================================================
# L3. Fine Analysis (anno3) - Split by anno3big
# ==============================================================================
message("\n=== Starting Level 3: Fine Analysis (Split by anno3big) ===")
# Get broad categories
broad_types <- unique(na.omit(as.character(sobj_full$anno3big)))
message("Broad Types found: ", paste(broad_types, collapse = ", "))

for (bt in broad_types) {
    label <- gsub("[^A-Za-z0-9]", "_", bt) # Sanitize

    # Check if result already exists to avoid subsetting overhead
    target_file <- file.path(base_out_dir, paste0("Level3_Fine_", label), paste0("results_", label, ".qs"))
    if (file.exists(target_file)) {
        message("Skipping existing: ", label)
        next
    }

    message("Subsetting for: ", bt)
    # Subset
    cells <- colnames(sobj_full)[sobj_full$anno3big == bt & !is.na(sobj_full$anno3big)]
    if (length(cells) < 100) {
        message("Too few cells, skipping: ", bt)
        next
    }

    sobj_sub <- subset(sobj_full, cells = cells)

    # Run
    # Use 'anno3' as cluster_id
    # Clean up unused clusters
    sobj_sub$anno3 <- droplevels(factor(sobj_sub$anno3))

    run_analysis(sobj_sub, label, "anno3", paste0("Level3_Fine_", label))

    rm(sobj_sub)
    gc()
}

message("\n=== All Analyses Completed ===")

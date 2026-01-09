# Production Analysis Script for Stroke PBMC (IS5) - Refined
# Levels: Global, Broad (anno3big), Fine (anno3 within anno3big)

# 0. Setup Environment
renv_lib <- "/home/user3/GJC_KDW_250721/renv/library/R-4.3/x86_64-pc-linux-gnu"
if (dir.exists(renv_lib)) .libPaths(c(renv_lib, .libPaths()))

suppressPackageStartupMessages({
    library(Seurat)
    library(qs)
    library(dplyr)
    library(doFuture)
    library(future)
    if (dir.exists("myR")) {
        devtools::load_all("myR")
    }
})

# Setup Parallel
registerDoFuture()
plan(multisession, workers = 6)

# Explicitly source if not found
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

# Set anno3big manually based on standard PBMC mapping if missing
if (!"anno3big" %in% colnames(meta)) {
    message("Maps 'anno3.scvi' to 'anno3big' using string matching...")
    # Heuristic mapping for PBMC
    # T Cells: CD4, CD8, Treg, T-cell
    # B Cells: B-cell, Plasma
    # Monocytes: Mono, CD14, CD16
    # NK: NK
    # DC: DC, pDC, cDC
    # Platelet: Platelet, Megakaryocyte
    # Other: Proliferating, Erythrocytes, etc.

    a3 <- as.character(meta$anno3.scvi)
    a3b <- rep("Others", length(a3))

    a3b[grep("T-cell|Treg|CD4|CD8|MAIT|gdT", a3, ignore.case = TRUE)] <- "Tc"
    a3b[grep("B-cell|Plasma|B naive|B memory", a3, ignore.case = TRUE)] <- "Bc"
    a3b[grep("Monocyte|Mono|Macrophage", a3, ignore.case = TRUE)] <- "Mono"
    a3b[grep("NK", a3, ignore.case = TRUE)] <- "NKc"
    a3b[grep("DC|Dendritic", a3, ignore.case = TRUE)] <- "DC"
    a3b[grep("Platelet|Megakaryocyte", a3, ignore.case = TRUE)] <- "PLT"
    a3b[grep("Erythrocyte|RBC", a3, ignore.case = TRUE)] <- "RBC"

    meta$anno3big <- factor(a3b)
    message("Derived anno3big levels: ", paste(levels(meta$anno3big), collapse = ", "))
}
sobj_full@meta.data <- meta

# Output Base Dir
base_out_dir <- "/data/user3/sobj/consensus"
dir.create(base_out_dir, recursive = TRUE, showWarnings = FALSE)

# Parameters
contrast_str <- "2 - 1"
group_col <- "g3"
sample_col <- "hos_no"
# Removed BMI as per user request
covariates <- c("sex", "age", "set")
cores <- 6

# ==============================================================================
# Helper Function
# ==============================================================================
run_analysis <- function(obj, label, cluster_col, out_subdir, methods) {
    out_dir <- file.path(base_out_dir, out_subdir)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    out_file <- file.path(out_dir, paste0("results_", label, ".qs"))

    message(sprintf("\n>>> Starting Analysis: %s (Cluster: %s) methods: %s", label, cluster_col, paste(methods, collapse = ",")))

    if (FALSE && file.exists(out_file)) {
        # Check file size to avoid empty skip
        info <- file.info(out_file)
        if (info$size > 10000) {
            message("File exists and seems valid (>10KB), skipping: ", out_file)
            return(NULL)
        } else {
            message("File exists but too small (<10KB), overwriting: ", out_file)
        }
    }

    tryCatch(
        {
            res <- run_deg_consensus(
                sobj = obj,
                contrast = contrast_str,
                methods = methods,
                cluster_id = cluster_col,
                sample_id = sample_col,
                group_id = group_col,
                covar_effects = covariates,
                n_cores = cores,
                verbose = TRUE
            )

            # --- Aggregation Steps ---
            message("Aggregating consensus results...")
            std <- lapply(res$methods_run, function(m) {
                standardize_deg_results(res$results[[m]], m)
            })
            names(std) <- res$methods_run
            res$standardized_results <- std

            # Matrices (Using logic from pipeline: pvalue mode with 0.05, FDR threshold 0.1 for building)
            # Adjust thresholds as needed
            mat <- build_deg_matrices(std, fdr_threshold = 0.1, significance_mode = "pvalue", pvalue_threshold = 0.05)
            res$deg_matrices <- mat

            # Consensus
            agree <- compute_agreement_scores(mat$significance)
            cons_scores <- compute_consensus_scores(mat, agree)
            res$consensus_scores <- cons_scores

            # Consensus List (Final Filtering)
            cons_list <- generate_consensus_deg_list(cons_scores, fdr_threshold = 0.05)
            res$consensus_deg_list <- cons_list
            res$consensus_table <- cons_list

            qs::qsave(res, out_file)
            message("Saved: ", out_file)

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
# Full pseudobulk suite
pb_methods <- c("edgeR-LRT", "edgeR-QLF", "limma-voom", "limma-trend", "DESeq2-Wald", "DESeq2-LRT", "dream")
run_analysis(sobj_full, "Global", "global_group", "Level1_Global",
    methods = pb_methods
)

# ==============================================================================
# L2. Broad Analysis (anno3big)
# ==============================================================================
run_analysis(sobj_full, "Broad", "anno3big", "Level2_Broad",
    methods = pb_methods
)

# ==============================================================================
# L3. Fine Analysis (anno3) - Split by anno3big
# ==============================================================================
message("\n=== Starting Level 3: Fine Analysis (Split by anno3big) ===")
# NEBULA included in separate process, here we run full PB suite
methods_fine <- pb_methods

broad_types <- unique(na.omit(as.character(sobj_full$anno3big)))
message("Broad Types found: ", paste(broad_types, collapse = ", "))

for (bt in broad_types) {
    label <- gsub("[^A-Za-z0-9]", "_", bt) # Sanitize

    # Check if result already exists
    target_file <- file.path(base_out_dir, paste0("Level3_Fine_", label), paste0("results_", label, ".qs"))
    if (file.exists(target_file)) {
        info <- file.info(target_file)
        if (info$size > 10000) {
            message("Skipping existing (valid): ", label)
            next
        }
    }

    message("Subsetting for: ", bt)
    cells <- colnames(sobj_full)[sobj_full$anno3big == bt & !is.na(sobj_full$anno3big)]
    if (length(cells) < 50) {
        message("Too few cells, skipping: ", bt)
        next
    }

    sobj_sub <- subset(sobj_full, cells = cells)
    sobj_sub$anno3 <- droplevels(factor(sobj_sub$anno3.scvi)) # Use anno3.scvi as fine cluster

    run_analysis(sobj_sub, label, "anno3", paste0("Level3_Fine_", label), methods_fine)

    rm(sobj_sub)
    gc()
}

message("\n=== All Analyses Completed ===")

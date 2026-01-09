# Production Analysis - NEBULA Only
# Split due to memory requirements and long runtime

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

if (!exists("run_deg_consensus")) {
    source("myR/R/deg_consensus/run_deg_consensus.R")
    source("myR/R/deg_consensus/deg_methods_base.R")
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
# Impute BMI for NAs (Median Imputation)
if (any(is.na(meta$bmi))) {
    med_bmi <- median(meta$bmi, na.rm = TRUE)
    meta$bmi[is.na(meta$bmi)] <- med_bmi
    message("Imputed BMI for NAs with median: ", round(med_bmi, 2))
}

meta$g3 <- factor(meta$g3)
meta$sex <- factor(meta$sex)
meta$set <- factor(meta$set)
meta$global_group <- "All_Cells"

# Set anno3big manually
if (!"anno3big" %in% colnames(meta)) {
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
}
sobj_full@meta.data <- meta

# Output Base Dir
base_out_dir <- "/data/user3/sobj/consensus"
dir.create(base_out_dir, recursive = TRUE, showWarnings = FALSE)

# Parameters
contrast_str <- "2 - 1"
group_col <- "g3"
sample_col <- "hos_no"
covariates <- c("sex", "age", "set")
cores <- 6

# Helper
run_nebula_safe <- function(obj, label, cluster_col, out_subdir) {
    out_dir <- file.path(base_out_dir, out_subdir)
    out_file <- file.path(out_dir, paste0("results_", label, "_nebula.qs"))

    if (file.exists(out_file)) {
        message("Skipping existing: ", out_file)
        return(NULL)
    }

    message(sprintf("\n>>> Starting NEBULA: %s (Cluster: %s)", label, cluster_col))

    tryCatch(
        {
            # Only run nebula
            res <- run_deg_consensus(
                sobj = obj,
                contrast = contrast_str,
                methods = c("nebula"),
                cluster_id = cluster_col,
                sample_id = sample_col,
                group_id = group_col,
                covar_effects = covariates,
                n_cores = cores,
                verbose = TRUE
            )
            qs::qsave(res, out_file)
        },
        error = function(e) {
            message("Error in ", label, ": ", e$message)
        }
    )
}

# Levels
run_nebula_safe(sobj_full, "Global", "global_group", "Level1_Global")
run_nebula_safe(sobj_full, "Broad", "anno3big", "Level2_Broad")

broad_types <- unique(na.omit(as.character(sobj_full$anno3big)))
for (bt in broad_types) {
    label <- gsub("[^A-Za-z0-9]", "_", bt)
    cells <- colnames(sobj_full)[sobj_full$anno3big == bt & !is.na(sobj_full$anno3big)]
    if (length(cells) < 50) next

    sobj_sub <- subset(sobj_full, cells = cells)
    sobj_sub$anno3 <- droplevels(factor(sobj_sub$anno3.scvi))

    run_nebula_safe(sobj_sub, label, "anno3", paste0("Level3_Fine_", label))
    rm(sobj_sub)
    gc()
}

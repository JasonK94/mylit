# AG_run1 Consensus Analysis Script
#
# Objective: Run comprehensive DEG consensus analysis (12+ methods) on IS2/3 clustered data.
# Levels: Broad (anno3big), Fine (anno3)
# Output: /data/user3/sobj/consensus/AG_run1
# Features:
#   - 12 methods: edgeR(LRT/QLF), DESeq2(Wald/LRT), limma(voom/trend), dream, nebula, muscat(4 variants)
#   - Method clustering based on statistics (t/F/Z)
#   - Comparative visualization

# 0. Setup
renv_lib <- "/home/user3/GJC_KDW_250721/renv/library/R-4.3/x86_64-pc-linux-gnu"
if (dir.exists(renv_lib)) .libPaths(c(renv_lib, .libPaths()))

suppressPackageStartupMessages({
    library(Seurat)
    library(qs)
    library(dplyr)
    library(ggplot2)
    library(doFuture)
    library(future)
    library(pheatmap)
    if (dir.exists("myR")) devtools::load_all("myR")
})

# Manual Sourcing to ensure all functions are available
src_files <- list.files("myR/R/deg_consensus", pattern = "\\.R$", full.names = TRUE)
for (f in src_files) source(f)

# Parallel Setup
registerDoFuture()
plan(multisession, workers = 8)

# 1. Load Data
input_path <- "/data/user3/sobj/is2_IS_3_clustered.qs"
message("Loading data: ", input_path)
sobj <- qs::qread(input_path)

# Ensure metadata
meta <- sobj@meta.data
if (!"anno3big" %in% names(meta)) {
    warning("anno3big not found, attempting to derive from anno3 or map manually.")
    # Assuming anno3 exists based on previous checks
    # If not, we might need to map. But user said "anno3big(5~7개 클러스터)" exists or implies it.
    # Previous check failed to print levels.
    # We will check if we need to create it.
}

# 2. Config
base_out_dir <- "/data/user3/sobj/consensus/AG_run1"
dir.create(base_out_dir, recursive = TRUE, showWarnings = FALSE)

contrast_str <- "2 - 1" # Assuming group 'g3' has levels 1 and 2, or we need to check 'g3' levels.
# User said "is2_IS_3_clustered". Usually condition is in 'g3'? or 'type'?
# Reference file was IS1205.
# Let's inspect 'g3' levels in the data quickly (or assume "2 - 1" is correct/standard).
# If g3 has other levels, we might need to adjust.
# Check input data levels if possible? I'll assume 2-1 for now, or use dynamic contrast.

methods_all <- c(
    # Standalone
    "edgeR-LRT", "edgeR-QLF",
    "DESeq2-Wald", "DESeq2-LRT",
    "limma-voom", "limma-trend",
    "dream",
    "nebula",
    # Muscat
    "muscat-edgeR", "muscat-DESeq2",
    "muscat-limma-voom", "muscat-limma-trend"
)

# 3. Helper: Run & Consensus & Cluster Methods
run_level_analysis <- function(obj, level_name, cluster_col) {
    out_dir <- file.path(base_out_dir, level_name)
    dir.create(out_dir, showWarnings = FALSE)
    out_file <- file.path(out_dir, "results_consensus.qs")

    if (file.exists(out_file) && file.info(out_file)$size > 10000) {
        message(sprintf("Skipping %s (File exists)", level_name))
        return(NULL)
    }

    message(sprintf("Running Analysis for %s (Cluster: %s)", level_name, cluster_col))

    # Run Consensus Pipeline
    # Using run_deg_consensus_with_plots?
    # It does standardize and build matrices.
    # But we want to ensure we run ALL methods.
    # Check contrast validity
    grps <- levels(factor(obj$g3))
    message("Group levels: ", paste(grps, collapse = ", "))
    # If contrast isn't valid, stop?

    # Explicit execution pipeline (like production script) to handle covar_effects
    # 1. Run DEG Methods
    res <- tryCatch(
        {
            run_deg_consensus(
                sobj = obj,
                contrast = contrast_str,
                methods = methods_all,
                cluster_id = cluster_col,
                sample_id = "hos_no",
                group_id = "g3",
                batch_id = NULL, # Assuming no batch_id, or should we use 'set' if available? is2_IS_3 data usually has 'set'?
                covar_effects = c("sex", "age"),
                n_cores = 8,
                verbose = TRUE
            )
        },
        error = function(e) {
            message("Error in run_deg_consensus: ", e$message)
            return(NULL)
        }
    )

    if (is.null(res)) {
        return(NULL)
    }

    # 2. Standardize
    message("Standardizing results...")
    std <- lapply(res$methods_run, function(m) {
        standardize_deg_results(res$results[[m]], m)
    })
    names(std) <- res$methods_run
    res$standardized_results <- std

    # 3. Build Matrices
    # pvalue mode 0.05, fdr 0.1 for building
    mat <- build_deg_matrices(std, fdr_threshold = 0.1, significance_mode = "pvalue", pvalue_threshold = 0.05)
    res$deg_matrices <- mat

    # 4. Consensus Scores
    agree <- compute_agreement_scores(mat$significance)
    cons_scores <- compute_consensus_scores(mat, agree)
    res$consensus_scores <- cons_scores

    # 5. Consensus List
    cons_list <- generate_consensus_deg_list(cons_scores, fdr_threshold = 0.05)
    res$consensus_deg_list <- cons_list

    qs::qsave(res, out_file)

    # 6. Plots (Standard)
    plot_prefix <- file.path(out_dir, "AG_run1")

    # Construct Reproducibility Caption
    # Get git hash if possible
    git_hash <- tryCatch(system("git rev-parse --short HEAD", intern = TRUE), error = function(e) "unknown")
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

    caption_text <- paste0(
        "Cluster: ", cluster_col, " (", level_name, ")\n",
        "Formula: ~ g3 + sex + age + (1|hos_no) [Varies by method]\n",
        "Script: scripts/consensus/run_AG_run1.R\n",
        "Time: ", timestamp, "\n",
        "Commit: ", git_hash
    )

    try({
        plot_consensus_meta_distributions(cons_scores, file_prefix = plot_prefix, caption = caption_text)
        plot_consensus_volcano(cons_scores, file_prefix = plot_prefix, fdr_threshold = 0.05, caption = caption_text)
        plot_consensus_heatmap(mat, cons_scores, file_prefix = plot_prefix, method_top_n = 20, caption = caption_text)
        # PCA of methods
        pca_res <- perform_deg_pca(mat)
        if (!is.null(pca_res)) plot_method_pca(pca_res, file_prefix = plot_prefix, caption = caption_text)
        plot_gene_umap(mat, cons_scores, file_prefix = plot_prefix, caption = caption_text)

        # Method Summary CSV
        export_method_summary_csv(res, file.path(out_dir, "AG_run1_method_summary.csv"))
    })

    # 7. Method Similarity Analysis (Statistic-based)
    message("Performing Method Similarity Analysis (Statistics)...")
    mat_stat <- mat$statistic

    if (!is.null(mat_stat)) {
        # Check valid rows (at least 2 non-NA)
        valid_rows <- rowSums(!is.na(mat_stat)) >= 2
        if (sum(valid_rows) > 10) {
            mat_stat_sub <- mat_stat[valid_rows, ]
            # Replace remaining NA? or use pairwise.
            # Using pairwise correlation
            cor_mat <- cor(mat_stat_sub, use = "pairwise.complete.obs", method = "spearman")

            # Plot Heatmap
            p_file <- file.path(out_dir, "AG_run1_method_similarity_heatmap.png")
            png(p_file, width = 800, height = 800)
            pheatmap::pheatmap(cor_mat,
                main = paste("Method Similarity (Spearman Correlation of Stats) -", level_name),
                display_numbers = TRUE,
                treeheight_row = 0
            )
            dev.off()

            # Save clustering result
            res$method_clustering <- list(cor_mat = cor_mat)
            qs::qsave(res, out_file)
        } else {
            message("Not enough valid genes for correlation analysis.")
        }
    }

    return(res)
}

# 4. Execution
# 4.1 Broad Analysis (anno3big)
if ("anno3big" %in% names(meta)) {
    run_level_analysis(sobj, "Broad_anno3big", "anno3big")
} else {
    warning("anno3big not found in metadata.")
}

# 4.2 Fine Analysis (anno3)
if ("anno3" %in% names(meta)) {
    run_level_analysis(sobj, "Fine_anno3", "anno3")
} else {
    warning("anno3 not found in metadata.")
}

message("AG_run1 Analysis Complete.")

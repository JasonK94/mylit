# Test CMGI Dev (Normalization Methods)

source("/home/user3/data_user3/git_repo/_wt/fgs/scripts/fgs/init_fgs_env.R")
renv::activate("/home/user3/GJC_KDW_250721")
library(qs)

# Load dev function
source("/home/user3/data_user3/git_repo/_wt/fgs/myR/R/signature_dev.R")

# Load a TML file (use one from previous runs or failed clusters run if available)
# Let's look for a file
tml_files <- dir("/data/user3/sobj/fgs/tmlas", pattern = "tmla_.*\\.qs", full.names = TRUE)
if (length(tml_files) == 0) {
    # Try run2 folder
    tml_files <- dir("/data/user3/sobj/fgs/tmlas/run2", pattern = "tmla_.*\\.qs", full.names = TRUE)
}

if (length(tml_files) == 0) {
    stop("No TML files found for testing.")
}

tml_file <- tml_files[1]
message(sprintf("Testing with TML file: %s", basename(tml_file)))
tmla <- qread(tml_file)

# Test methods
methods <- c("max_abs", "min_max", "softmax", "rank", "z_score")

results <- list()

for (m in methods) {
    message(sprintf("\nTesting normalization: %s", m))
    res <- try(compute_meta_gene_importance_v2(tmla, normalization_method = m), silent = TRUE)

    if (inherits(res, "try-error")) {
        message(sprintf("Failed: %s", as.character(res)))
    } else if (is.null(res)) {
        message("Result is NULL")
    } else {
        message("Success!")
        top_genes <- head(res$gene_summary$gene, 5)
        message(sprintf("Top 5 genes: %s", paste(top_genes, collapse = ", ")))
        results[[m]] <- res
    }
}

# Compare top genes across methods
if (length(results) > 1) {
    message("\n=== Comparison of Top 10 Genes ===")
    top10_df <- data.frame(Rank = 1:10)
    for (m in names(results)) {
        top10_df[[m]] <- head(results[[m]]$gene_summary$gene, 10)
    }
    print(top10_df)
}

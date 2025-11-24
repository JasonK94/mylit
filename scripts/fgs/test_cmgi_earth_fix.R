#!/usr/bin/env Rscript
# Test CMGI Earth Fix

# 1. Activate Environment
cat("Activating renv...\n")
renv::activate("/home/user3/GJC_KDW_250721")

# 2. Check Earth Package
cat("Checking earth package...\n")
if (!requireNamespace("earth", quietly = TRUE)) {
    cat("Package 'earth' is NOT available. Attempting to install...\n")
    # Try installing from CRAN if missing (though renv should handle it if in lockfile)
    # For now, just report it.
} else {
    cat("Package 'earth' is available.\n")
}

library(qs)
source("/home/user3/data_user3/git_repo/_wt/fgs/myR/R/signature.R")

# 3. Load TML with Earth model
# Using a file known to have earth model from previous steps
tml_file <- "/data/user3/sobj/fgs/tmlas/run2/tmla_CD4+ T-cells_25-11-20-10-40.qs"

if (!file.exists(tml_file)) {
    stop("TML file not found: ", tml_file)
}

cat("Loading TML file:", basename(tml_file), "\n")
tmla <- qread(tml_file)

# 4. Test CMGI for Earth
if ("earth" %in% names(tmla$trained_models)) {
    cat("Testing CMGI for 'earth' model...\n")

    res <- tryCatch(
        {
            compute_meta_gene_importance(tmla, target_model = "earth")
        },
        error = function(e) {
            cat("Error in CMGI:", conditionMessage(e), "\n")
            return(NULL)
        }
    )

    if (!is.null(res)) {
        cat("✓ CMGI for Earth SUCCEEDED!\n")
        cat("  Top 3 genes:\n")
        print(head(res$gene_summary, 3))
    } else {
        cat("✗ CMGI for Earth FAILED or returned NULL.\n")
    }
} else {
    cat("Earth model not found in this TML file.\n")
}

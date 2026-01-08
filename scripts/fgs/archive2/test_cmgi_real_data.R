#!/usr/bin/env Rscript
# Test CMGI with real TML data

library(qs)
source("/home/user3/data_user3/git_repo/_wt/fgs/myR/R/signature.R")

# List of TML files to test
tml_files <- list.files("/data/user3/sobj/fgs/tmlas/run2",
    pattern = "tmla.*\\.qs$",
    full.names = TRUE
)

cat("Found", length(tml_files), "TML files\n\n")

# Test a few different ones
test_files <- tml_files[c(1, 5, 9)] # Sample 3 files

results <- list()

for (tml_file in test_files) {
    cat("\n========================================\n")
    cat("Testing:", basename(tml_file), "\n")
    cat("========================================\n")

    tryCatch(
        {
            # Load TML
            tmla <- qread(tml_file)

            cat("Loaded TML successfully\n")
            cat("Trained models:", paste(names(tmla$trained_models), collapse = ", "), "\n")
            cat("Best model:", tmla$best_model_name, "\n\n")

            # Test CMGI for each trained model
            for (model_name in names(tmla$trained_models)) {
                cat("--- Testing CMGI for model:", model_name, "---\n")

                cmgi_result <- try(
                    compute_meta_gene_importance(tmla, target_model = model_name),
                    silent = FALSE
                )

                if (!inherits(cmgi_result, "try-error") && !is.null(cmgi_result)) {
                    cat("✓ CMGI SUCCESS for", model_name, "\n")
                    cat("  Signature importance count:", length(cmgi_result$signature_importance), "\n")
                    cat("  Gene importance rows:", nrow(cmgi_result$gene_importance), "\n")
                    cat("  Top 3 genes:\n")
                    print(head(cmgi_result$gene_summary, 3))

                    results[[paste(basename(tml_file), model_name, sep = "_")]] <- "SUCCESS"
                } else {
                    cat("✗ CMGI FAILED for", model_name, "\n")
                    if (inherits(cmgi_result, "try-error")) {
                        cat("  Error:", as.character(cmgi_result), "\n")
                    }
                    results[[paste(basename(tml_file), model_name, sep = "_")]] <- "FAILED"
                }
                cat("\n")
            }
        },
        error = function(e) {
            cat("ERROR loading/processing TML:", conditionMessage(e), "\n")
            results[[basename(tml_file)]] <- paste("ERROR:", conditionMessage(e))
        }
    )
}

cat("\n\n========================================\n")
cat("SUMMARY\n")
cat("========================================\n")
print(results)

# Count successes
success_count <- sum(grepl("SUCCESS", results))
total_count <- length(results)
cat(
    "\nSuccess rate:", success_count, "/", total_count,
    sprintf("(%.1f%%)\n", 100 * success_count / total_count)
)

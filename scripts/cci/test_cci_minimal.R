# Minimal test to check if functions load correctly
# This can be run with: Rscript scripts/cci/test_cci_minimal.R

cat("=== CCI Tool Minimal Test ===\n")

# Check if required packages are available
cat("Checking packages...\n")
pkgs <- c("Seurat", "dplyr", "qs")
for (pkg in pkgs) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat("  ✓", pkg, "\n")
  } else {
    cat("  ✗", pkg, "(NOT INSTALLED)\n")
  }
}

# Try to source functions
cat("\nLoading functions...\n")
func_files <- c(
  "/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/prepare_cci_data.R",
  "/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/utils_cci.R",
  "/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/save_cci_results.R",
  "/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/run_cci_analysis.R"
)

for (f in func_files) {
  if (file.exists(f)) {
    tryCatch({
      source(f)
      cat("  ✓", basename(f), "\n")
    }, error = function(e) {
      cat("  ✗", basename(f), "- Error:", e$message, "\n")
    })
  } else {
    cat("  ✗", basename(f), "- File not found\n")
  }
}

# Check if run_nichenet_analysis is available
cat("\nChecking run_nichenet_analysis...\n")
if (file.exists("/home/user3/data_user3/git_repo/mylit/myR/R/CCI.R")) {
  tryCatch({
    source("/home/user3/data_user3/git_repo/mylit/myR/R/CCI.R")
    if (exists("run_nichenet_analysis")) {
      cat("  ✓ run_nichenet_analysis loaded\n")
    } else {
      cat("  ✗ run_nichenet_analysis not found after sourcing\n")
    }
  }, error = function(e) {
    cat("  ✗ Error loading CCI.R:", e$message, "\n")
  })
} else {
  cat("  ✗ CCI.R file not found\n")
}

# Check function availability
cat("\nChecking function availability...\n")
funcs <- c("validate_cci_inputs", "extract_receiver_degs", "identify_sender_clusters", 
           "run_cci_analysis", "save_cci_final")
for (func in funcs) {
  if (exists(func)) {
    cat("  ✓", func, "\n")
  } else {
    cat("  ✗", func, "(NOT FOUND)\n")
  }
}

cat("\n=== Test Complete ===\n")


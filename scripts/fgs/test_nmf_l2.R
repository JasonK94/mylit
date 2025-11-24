# scripts/fgs/test_nmf_l2.R

# --- Environment Setup ---
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
renv::activate("/home/user3/GJC_KDW_250721")

# Load necessary libraries
suppressPackageStartupMessages({
    library(Seurat)
    library(caret)
    library(dplyr)
    library(myR)
})

# Source latest signature.R if needed (for dev)
source("/home/user3/data_user3/git_repo/_wt/fgs/myR/R/signature.R")
source("/home/user3/data_user3/git_repo/_wt/fgs/myR/R/utils_tml.R")

# --- Data Setup ---
# Create a small dummy Seurat object or load existing one
data_path <- "/data/user3/sobj/is5a_0926.qs"

if (file.exists(data_path)) {
    message("Loading real data subset...")
    sobj_full <- qs::qread(data_path)
    # Subset to small size: 200 cells
    set.seed(42)
    cells <- sample(Cells(sobj_full), 200)
    sobj <- subset(sobj_full, cells = cells)

    # Ensure target variable exists and has 2 levels
    if (!"g3" %in% colnames(sobj@meta.data)) {
        sobj$g3 <- sample(c("GroupA", "GroupB"), ncol(sobj), replace = TRUE)
    }
} else {
    message("Creating dummy data...")
    # Create dummy data
    counts <- matrix(rpois(200 * 500, lambda = 10), nrow = 500, ncol = 200)
    rownames(counts) <- paste0("Gene", 1:500)
    colnames(counts) <- paste0("Cell", 1:200)
    sobj <- CreateSeuratObject(counts = counts)
    sobj$g3 <- sample(c("GroupA", "GroupB"), 200, replace = TRUE)
    sobj$hos_no <- sample(c("H1", "H2"), 200, replace = TRUE)
    sobj <- NormalizeData(sobj)
    sobj <- FindVariableFeatures(sobj)
    sobj <- ScaleData(sobj)
}

message("Data ready. Cells: ", ncol(sobj))

# --- Test 1: FGS with NMF ---
message("\n=== Test 1: FGS with NMF and Ranger ===")
fgs_methods <- c("random_forest_ranger", "nmf_loadings", "lasso")

tryCatch(
    {
        fgsa <- find_gene_signature_v5.4(
            sobj,
            target_var = "g3",
            control_vars = "hos_no",
            n_features = 50,
            method = fgs_methods,
            min_cells = 10
        )
        message("FGS completed successfully.")
        print(names(fgsa))
    },
    error = function(e) {
        message("FGS Failed: ", e$message)
        stop(e)
    }
)

# --- Test 2: TML7 with All L2 Methods ---
message("\n=== Test 2: TML7 with All L2 Methods ===")
# Define L2 methods to test (including those that failed before)
l2_methods <- c("glm", "ranger", "svmRadial", "mlp", "earth", "xgbTree")

tryCatch(
    {
        tmla <- TML7(
            l1_signatures = fgsa,
            holdout_data = sobj, # Using same data for test simplicity
            target_var = "g3",
            l2_methods = l2_methods,
            cv_folds = 3 # Small folds for speed
        )
        message("TML7 completed successfully.")

        # Check trained models
        print("Trained models:")
        print(names(tmla$trained_models))

        # Check CV folds info
        if (!is.null(tmla$cv_folds)) {
            message("CV folds info present.")
            print(str(tmla$cv_folds))
        } else {
            warning("CV folds info MISSING.")
        }
    },
    error = function(e) {
        message("TML7 Failed: ", e$message)
        stop(e)
    }
)

message("\n=== All Tests Completed ===")

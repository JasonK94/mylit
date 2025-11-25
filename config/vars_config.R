# ==============================================================================
# Dataset Configuration for FGS/TML Testing
# ==============================================================================
# Based on guide.md specifications
# Usage: source("scripts/fgs/vars_config.R"); config <- get_dataset_config("stroke")

get_dataset_config <- function(dataset_name = "stroke", level = "full") {
    if (dataset_name == "stroke") {
        # --- Stroke Dataset (IS6) ---
        # Target: g3 ("1", "2", "NA") - NA is string!
        # Patient: hos_no (8-digit number)
        # Cluster: anno3.scvi

        base_path <- "/data/user3/sobj"

        data_path <- switch(level,
            "full" = file.path(base_path, "IS6_sex_added_251110.qs"),
            "downsampled" = file.path(base_path, "IS6_sex_added_0.1x_251110.qs"), # 0.1x
            "ds2500" = file.path(base_path, "IS_scvi_251107_ds2500.qs"), # Legacy test set
            stop(sprintf("Unknown level '%s' for stroke dataset", level))
        )

        config <- list(
            name = "stroke",
            path = data_path,

            # Variables
            patient_var = "hos_no",
            target_var = "g3",
            cluster_var = "anno3.scvi",

            # Covariates (Fixed effects)
            covariates = c("GEM", "set"),

            # Target classes to use (exclude NA string)
            target_classes = c("1", "2"),
            target_na_value = "NA" # String "NA"
        )
    } else if (dataset_name == "ibd") {
        # --- IBD Dataset ---
        # Placeholder based on guide.md

        config <- list(
            name = "ibd",
            sobj_path_geomx = "/data/user3/sobj/data_seurat_251104.qs",
            sobj_path_sc = "/data/user3/sobj/IBD_biologics_251125.qs",

            # Variables (Placeholder - update when known)
            patient_var = "emrid",
            target_var = "response",
            cluster_var = "ck",
            covariates = NULL,
            target_classes = c("R", "NR"),
            target_na_value = "NA" # String "NA"
        )
    } else {
        stop(sprintf("Unknown dataset name: %s", dataset_name))
    }

    return(config)
}


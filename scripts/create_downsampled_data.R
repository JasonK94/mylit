# Script to create downsampled data for debugging
library(Seurat)
library(qs)
library(dplyr)

# Config
input_file <- "/data/user3/sobj/pipe/run_ag1/step5/step5_doublet_list.qs"
output_dir <- "/data/user3/sobj/pipe/run_ag_ds/step5"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_file <- file.path(output_dir, "step5_doublet_list.qs")

print(paste("Reading input:", input_file))
if (file.exists(input_file)) {
    obj_list <- qread(input_file)
    print(paste("Loaded list with", length(obj_list), "objects"))

    # Downsample each object
    downsampled_list <- lapply(names(obj_list), function(name) {
        obj <- obj_list[[name]]
        if (inherits(obj, "Seurat")) {
            n_cells <- ncol(obj)
            if (n_cells > 500) {
                # print(paste("Downsampling", name, "from", n_cells, "to 500"))
                obj <- subset(obj, cells = sample(Cells(obj), 500))
            } else {
                # print(paste("Keeping", name, "as is (", n_cells, "cells)"))
            }
        }
        return(obj)
    })
    names(downsampled_list) <- names(obj_list)

    print("Saving downsampled list...")
    qsave(downsampled_list, output_file)
    print(paste("Saved to:", output_file))
} else {
    stop("Input file not found!")
}

print("Downsampling complete.")

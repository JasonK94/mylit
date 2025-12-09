# Script to update manifest_stroke.csv with new columns
library(dplyr)
library(readr)

manifest_path <- "config/manifest_stroke.csv"
manifest <- read_csv(manifest_path, show_col_types = FALSE)

# Add new columns as requested
# gem_name, set_name, batch_name, project_name
# sample_col, gem_col, set_col, batch_col, project_col

manifest <- manifest %>%
    mutate(
        # Mapping logic for names
        set_name = SET_DETAIL, # Use SET_DETAIL as set_name (SET1_1, etc.)
        batch_name = SET, # Use SET as batch_name (SET1, SET2, etc.)
        project_name = "Stroke", # Assuming project name

        # Meta columns (pointers)
        sample_col = "sample_name",
        gem_col = "gem_name", # Or "GEM" if referring to Seurat column? User example says "GEM"
        set_col = "set_name", # Or "SET"?
        batch_col = "batch_name", # Or "BATCH"?
        project_col = "project_name"
    ) %>%
    # Reorder columns to match example roughly
    select(
        no, name, sample_name, patient_name, gem_name, set_name, batch_name, project_name,
        demultiplex_id, sample_col, gem_col, set_col, batch_col, project_col,
        everything()
    )

# Write back
write_csv(manifest, manifest_path)
print(paste("Updated", manifest_path))

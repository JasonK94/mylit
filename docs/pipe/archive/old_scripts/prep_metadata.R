# Read the original metadata file
meta_file <- "config/meta_data_prep_251204_original.csv"
if (!file.exists(meta_file)) {
    stop("Metadata file not found: ", meta_file)
}

meta_data <- read.csv(meta_file, stringsAsFactors = FALSE)

# Update demulti_id for HTO samples
# Logic: if hto_no is present (not NA), demulti_id = HTO + (250 + hto_no)
# Assuming hto_no is numeric or convertible to numeric

if ("hto_no" %in% colnames(meta_data)) {
    # Convert hto_no to numeric, suppressing warnings for non-numeric values (which become NA)
    hto_nos <- suppressWarnings(as.numeric(meta_data$hto_no))

    # Identify rows where hto_no is valid
    valid_hto <- !is.na(hto_nos)

    # Update demulti_id for these rows
    # Note: The user said "$demulti_id에 매칭시켜서". This implies overwriting the existing demulti_id column.
    # We should be careful not to overwrite SNP demulti_ids if they happen to have hto_no (unlikely based on previous view).

    new_ids <- paste0("HTO", 250 + hto_nos[valid_hto])

    # Check if we are overwriting anything important
    # For HTO samples, demulti_id in original file seemed to be sample ID-like (e.g., 5_SAH_24).
    # Overwriting it with HTO tag is what the user requested for the pipeline to work.

    meta_data$demulti_id[valid_hto] <- new_ids

    cat(sprintf("Updated %d rows with HTO tags.\n", sum(valid_hto)))
    print(head(meta_data[valid_hto, c("sample_no", "hto_no", "demulti_id")]))
} else {
    warning("Column 'hto_no' not found in metadata.")
}

# Save the new metadata file
output_file <- "config/meta_data_prep_251204.csv"
write.csv(meta_data, output_file, row.names = FALSE, quote = TRUE)
cat("Saved updated metadata to:", output_file, "\n")

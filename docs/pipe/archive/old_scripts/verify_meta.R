meta <- read.csv("config/meta_data_prep_251204.csv", stringsAsFactors = FALSE)
row <- meta[meta$sample_no == "1_3", ]
print(row)
if (nrow(row) > 0) {
    cat("demulti_id:", row$demulti_id, "\n")
    cat("hto_no:", row$hto_no, "\n")
} else {
    cat("Sample 1_3 not found.\n")
}

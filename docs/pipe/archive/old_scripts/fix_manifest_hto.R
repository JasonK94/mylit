manifest <- read.csv("config/manifest_hto_fixed.csv", stringsAsFactors = FALSE)
# Update all rows since this manifest only contains the HTO samples
manifest$gem_name <- paste0(manifest$gem_name, "_", manifest$sample_name)
manifest$demultiplex_method <- "None"
write.csv(manifest, "config/manifest_hto_fixed.csv", row.names = FALSE, quote = TRUE)
